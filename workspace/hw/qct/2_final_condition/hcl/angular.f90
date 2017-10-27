module angular
  implicit none
  
contains

  subroutine rotn(iflag,n,Q,v,m,AM,EROT,womega,abc)
    implicit none
    integer,intent(in) :: iflag,n
    real(kind=8),intent(in) :: v(n,3),Q(n,3),m(n)
    real(kind=8),intent(inout) :: AM(3)
    real(kind=8),intent(out) :: EROT
    real(kind=8),intent(out) :: womega(3)
    ! ::::::::::::::::::::::::::::
    integer :: i,j,k,lwork,info
    real(kind=8) :: moi(3,3),moi_inv(3,3),AM4,P(n,3)
    real(kind=8) :: xl(3,n),x0(n,3)
    real(kind=8) :: tmas,det
    real(kind=8) :: abc(3)!To store const ABC for polyatomic mol
    real(kind=8) :: eig(3),work(1:40*size(eig)),evect(3,3)
    real(kind=8) :: sr,dum,cbmax,bmax,detcb1,detcb2,detcb3,detcb,convert
    real(kind=8) :: atoau,wavenm,pi,EROTC(2),w_abc(3),rconst(3)
    real(kind=8) :: WX,WY,WZ
    character(len=100) :: formatstring
    
abc = 0
!   write(*,*) 'iflag =',iflag
!   write(*,*) 'n =',n
!   write(*,*) 'Q =', Q
!   write(*,*) 'v =', v
!   write(*,*) 'm =', m

    atoau=1.66053886d-27/9.1093826d-31
    wavenm=2.1947463067d+05 !! multiply by this to convert Hartree to cm-1
    lwork=40*size(eig)
    pi = dacos(-1D0)
    rconst=0.0d0
    
    if(iflag.eq.0) then
       AM=0.0d0
       EROT=0.0d0
       if(n.eq.1) return
       do i=1,n
          P(i,1:3) = m(i)*v(i,1:3)
          AM(1) =  AM(1) +  (Q(i,2)*P(i,3)-Q(i,3)*P(i,2))
          AM(2) =  AM(2) +  (Q(i,3)*P(i,1)-Q(i,1)*P(i,3))
          AM(3) =  AM(3) +  (Q(i,1)*P(i,2)-Q(i,2)*P(i,1))
       end do
       AM4 = DSQRT(AM(1)**2+AM(2)**2+AM(3)**2)
       if(n.eq.2) then
          moi(1,1)=0.0d0
          do i=1,n
            moi(1,1)=moi(1,1)+(Q(i,1)**2+Q(i,2)**2+Q(i,3)**2)*m(i)
          end do
       EROT=AM4**2/moi(1,1)/2.0d0

       return
       end if !! n.eq.2
    end if  !! iflag.eq.0
    
    moi=0.0d0
    moi_inv=0.0d0
    !! construct moment of inertia
    do i=1,n
     tmas=m(i)
     moi(1,1) = moi(1,1) + tmas*(Q(i,2)**2+Q(i,3)**2)
     moi(2,2) = moi(2,2) + tmas*(Q(i,1)**2+Q(i,3)**2)
     moi(3,3) = moi(3,3) + tmas*(Q(i,2)**2+Q(i,1)**2)
     moi(1,2) = moi(1,2) - tmas*(Q(i,1)*Q(i,2))
     moi(2,3) = moi(2,3) - tmas*(Q(i,2)*Q(i,3))
     moi(1,3) = moi(1,3) - tmas*(Q(i,1)*Q(i,3))
    end do
     moi(2,1)=moi(1,2)
     moi(3,2)=moi(2,3)
     moi(3,1)=moi(1,3)
    
    !! save moment of inertia
    evect=moi
    !! diagonalize moment of inertia matrix
    call dsyev('v','u',3,evect,3,eig,work,lwork,info)
    if(info.ne.0) stop 'Moment of intertia diagonalization failed'
    
    !! calculate determinant
          det = ( moi(1,1)*moi(2,2)*moi(3,3) + &
                  moi(1,2)*moi(2,3)*moi(3,1) + &
                  moi(1,3)*moi(2,1)*moi(3,2) ) - &
                ( moi(3,1)*moi(2,2)*moi(1,3) + &
                  moi(3,2)*moi(2,3)*moi(1,1) + &
                  moi(3,3)*moi(2,1)*moi(1,2) )
    
    !! calculate inverse of the intertia tensor
         if(abs(det).ge.0.001d0) then
                moi_inv(1,1) =  (moi(2,2)*moi(3,3)-moi(2,3)*moi(3,2))/det
                moi_inv(1,2) =  (moi(1,3)*moi(3,2)-moi(1,2)*moi(3,3))/det
                moi_inv(1,3) =  (moi(1,2)*moi(2,3)-moi(1,3)*moi(2,2))/det
    
                moi_inv(2,1) =  (moi(2,3)*moi(3,1)-moi(2,1)*moi(3,3))/det
                moi_inv(2,2) =  (moi(1,1)*moi(3,3)-moi(1,3)*moi(3,1))/det
                moi_inv(2,3) =  (moi(1,3)*moi(2,1)-moi(1,1)*moi(2,3))/det
    
                moi_inv(3,1) =  (moi(2,1)*moi(3,2)-moi(2,2)*moi(3,1))/det
                moi_inv(3,2) =  (moi(1,2)*moi(3,1)-moi(1,1)*moi(3,2))/det
                moi_inv(3,3) =  (moi(1,1)*moi(2,2)-moi(1,2)*moi(2,1))/det
    
    !! angular velocity
                womega(1:3) = matmul(moi_inv,AM(1:3))
    
                WX=moi_inv(1,1)*AM(1)+moi_inv(1,2)*AM(2)+moi_inv(1,3)*AM(3)
                WY=moi_inv(2,1)*AM(1)+moi_inv(2,2)*AM(2)+moi_inv(2,3)*AM(3)
                WZ=moi_inv(3,1)*AM(1)+moi_inv(3,2)*AM(2)+moi_inv(3,3)*AM(3)
    
                w_abc(1:3)=matmul(womega,evect)
          else
    !! calculate rotational energy 
                womega(1:3)=0.0d0
                if(moi(1,1).ne.0.0) womega(1)=AM(1)/moi(1,1)
                if(moi(2,2).ne.0.0) womega(2)=AM(2)/moi(2,2)
                if(moi(3,3).ne.0.0) womega(3)=AM(3)/moi(3,3)
                moi(1,1)=0.0d0
             do i=1,n
               sr=0.0d0
               do j=1,3
                  sr=sr+Q(i,j)*Q(i,j)
               end do
                  moi(1,1)=moi(1,1)+sr*m(i)
             end do
             AM4=DSQRT(AM(1)**2+AM(2)**2+AM(3)**2)
             EROT=AM4**2/moi(1,1)/2.0d0
            return
          end if
    
     EROT=(womega(1)*AM(1)+womega(2)*AM(2)+womega(3)*AM(3))/2.0d0
    
     EROTC(1)=(eig(1)*w_abc(1)**2+eig(2)*w_abc(2)**2)/2.0
     EROTC(2)=(eig(3)*w_abc(3)**2)/2.0

abc(:) = 0.5/eig(:)

! write(*,*) 'Eigenvec'
!write(*,*) evect
!write(*,*) 'Eigenvalue',eig,abc
!write(*,*) 'Verify by evect*I*transpose(evect) = eigenvlue'
!write(*,*) matmul(transpose(evect),matmul(moi,evect))
    return
  end subroutine rotn

  !*****************************************************!
  ! Eckart frame       JCP, 122, 124103 (2005)          !
  !*****************************************************!
  subroutine eckart_coor(natom,rr,refStruc,mass,vv,tmin)
    double precision,parameter::threshold=0.25d0
    integer,intent(in)::natom
    double precision,dimension(3,natom),intent(inout)::rr,vv
    double precision,dimension(3,natom),intent(in)::refStruc
    double precision,dimension(natom),intent(in)::mass
    double precision,dimension(3,3),intent(out)::tmin
    !::::::::::::::::::::::::::::::::::::::::::::::::
    integer::i,j,k,ierr,i1,i2,i3
    double precision::totMass,temp,dxmin
    double precision,dimension(3,3)::a,a1,a2,u,v,t,vTemp
    double precision,dimension(3,natom)::rrTemp,vvTemp,rrmin,vvmin
    double precision,dimension(3)::d,e
  
    ! total mass of the system
    totMass = 0.d0
    do i = 1, natom
       totMass = totMass + mass(i)
    end do
  
    ! shifting the origin into the nuclear center of mass
    do i = 1, 3
       temp = 0.d0
       do j = 1, natom
          temp = temp + mass(j) * rr(i,j)
       end do
       temp = temp/totMass
       do j = 1, natom
          rr(i,j) = rr(i,j) - temp
       end do
    end do
  
    ! constructing a, a1 and a2 matrices
    do i = 1, 3
       do j = 1, 3
          temp = 0.d0
          do k = 1, natom
             temp = temp + mass(k) * rr(i,k) * refStruc(j,k)
          end do
          a(i,j) = temp
       end do
    end do
  
    do i = 1, 3
       do j = i, 3
          temp = 0.d0
          do k = 1, 3
             temp = temp + a(i,k) * a(j,k)
          end do
          a1(i,j) = temp
          a1(j,i) = temp
          temp = 0.d0
          do k = 1, 3
             temp = temp + a(k,i) * a(k,j)
          end do
          a2(i,j) = temp
          a2(j,i) = temp
       end do
    end do
  
    ! a1 and a2 eigenproblems
    call subtred2(3,3,a1,d,e,u)
    call tql2(3,3,d,e,u,ierr)
    call subtred2(3,3,a2,d,e,v)
    call tql2(3,3,d,e,v,ierr)

    dxmin = 1000000.0d0 
    rrTemp = rr
    vvTemp = vv
    vTemp = v
    do i1 = 1, 2
       do i2 = 1, 2
          do i3 = 1, 2

             do j = 1, 3
                v(j,1) = (-1.d0)**i1*vTemp(j,1)
                v(j,2) = (-1.d0)**i2*vTemp(j,2)
                v(j,3) = (-1.d0)**i3*vTemp(j,3)
             end do
  
             ! constructing t
             do i = 1, 3
                do j = 1, 3
                   temp = 0.d0
                   do k = 1, 3
                      temp = temp + v(i,k) * u(j,k)
                   end do
                   t(i,j) = temp
                end do
             end do
  
             ! rotating the Descartes position vectors of the nuclei
!             rrTemp = rr
             do i = 1, 3
                do j = 1, natom
                   temp = 0.d0
                   do k = 1, 3
                      temp = temp + t(i,k) * rrTemp(k,j)
                   end do
                   rr(i,j) = temp
                end do
             end do
  
!             vvTemp = vv
             do i = 1, 3
                do j = 1, natom
                   temp = 0.d0
                   do k = 1, 3
                      temp = temp + t(i,k) * vvTemp(k,j)
                   end do
                   vv(i,j) = temp
                end do
             end do
  
             ! testing the rotational Eckart conditions
!             temp = 0.d0
!             do j = 1, natom
!                temp = temp + mass(j) * (rr(2,j) * refStruc(3,j) - &
!                       rr(3,j) * refStruc(2,j))
!             end do
!             if (dabs(temp) > threshold) then
!                write(*,*) 'x error', temp
!             end if
!
!             temp = 0.d0
!             do j = 1, natom
!                temp = temp + mass(j) * (rr(3,j) * refStruc(1,j) - &
!                       rr(1,j) * refStruc(3,j))
!             end do
!             if (dabs(temp) > threshold) then
!                write(*,*) 'y error', temp
!             end if
!
!             temp = 0.d0
!             do j = 1, natom
!                temp = temp + mass(j) * (rr(1,j) * refStruc(2,j) - &
!                       rr(2,j)*refStruc(1,j))
!             end do
!             if (dabs(temp)>threshold) then
!                write(*,*) 'z error', temp
!             end if
  
             ! computing rr - ref
             temp = 0.0d0
             do i = 1, 3
                do j = 1, natom
                   temp = temp + (rr(i,j) - refStruc(i,j))**2
                end do
             end do
             if (temp < dxmin) then
                dxmin = temp  
                rrmin = rr
                vvmin = vv
                tmin = t 
             end if
  
          end do
       end do
    end do
   
    rr = rrmin
    vv = vvmin
  
    !determinant of the transformation matrix  
!    det=tmin(1,1)*(tmin(2,2)*tmin(3,3)-tmin(2,3)*tmin(3,2))
!    det=det-tmin(1,2)*(tmin(2,1)*tmin(3,3)-tmin(2,3)*tmin(3,1))
!    det=det+tmin(1,3)*(tmin(2,1)*tmin(3,2)-tmin(2,2)*tmin(3,1)) 
  
    return
  end subroutine eckart_coor

  !=======================================!
  ! diagnolize the interial tensor        !
  ! and set J=0                           !
  !=======================================!
  subroutine diag_it(mass,x,v)
    real,dimension(:),intent(in)::mass
    real,dimension(:),intent(inout)::x
    real,dimension(:),intent(inout)::v
    ! ::::::::::::::::::::
    real,dimension(1:3,1:3)::it
    integer::i,n
    real::tmas
    real,dimension(1:3)::cord,velc,J,omega
    integer::info,lwork
    real,dimension(1:3)::egvu !Eigen Vaule
    real,dimension(1:12)::work !just temporary space
    lwork=12
    
    n=size(mass,1)
    
    it=0    
    do i=1,n
       tmas = mass(i)
       cord = x(3*i-2:3*i)
       it(1,1) = it(1,1) + tmas*( cord(2)**2 + cord(3)**2 )
       it(2,2) = it(2,2) + tmas*( cord(1)**2 + cord(3)**2 )
       it(3,3) = it(3,3) + tmas*( cord(2)**2 + cord(1)**2 )
       it(1,2) = it(1,2) - tmas*( cord(1)*cord(2) )
       it(2,3) = it(2,3) - tmas*( cord(2)*cord(3) )
       it(1,3) = it(1,3) - tmas*( cord(1)*cord(3) )
    end do

    it(2,1)=it(1,2)
    it(3,2)=it(2,3)
    it(3,1)=it(1,3)

    call dsyev('v','u',3,it,3,egvu,work,lwork,info)
   
    do i=1,n
       cord = x(3*i-2:3*i)
       x(3*i-2:3*i) = matmul(transpose(it),cord)
       velc = v(3*i-2:3*i)
       v(3*i-2:3*i) = matmul(transpose(it),velc)
    end do

    !================!
    ! set J=0        !
    !================!
    J=0; velc=0; omega=0; cord=0
    
    do i=1,n
       tmas = mass(i)
       cord = x(3*i-2:3*i)
       velc = v(3*i-2:3*i)
       
       J(1) = J(1) + tmas*(cord(2)*velc(3)-cord(3)*velc(2))
       J(2) = J(2) + tmas*(cord(3)*velc(1)-cord(1)*velc(3))
       J(3) = J(3) + tmas*(cord(1)*velc(2)-cord(2)*velc(1))
    end do
    
    omega=J/egvu

    do i=1,n
       cord = x(3*i-2:3*i)
       velc = v(3*i-2:3*i)
       v(3*i-2) = velc(1) + cord(2)*omega(3) - cord(3)*omega(2)
       v(3*i-1) = velc(2) + cord(3)*omega(1) - cord(1)*omega(3)
       v(3*i)   = velc(3) + cord(1)*omega(2) - cord(2)*omega(1)
    end do

    return
  end subroutine diag_it

  !=====================================!
  ! Calculate the angular momentum      !
  !=====================================!
  subroutine monitor_J(mass,x,v)
    real,dimension(:),intent(in)::mass,x,v
!    real,dimension(1;3),intent(out)::J
    !:::::::::::::::::::
    real,dimension(1:3)::xx,vv,J
    integer::natm,i
    real::tmas

    natm = size(mass)

    J = 0.0
    do i=1,natm
       tmas = mass(i)
       xx = x(3*i-2:3*i)
       vv = v(3*i-2:3*i)
       
       J(1) = J(1) + tmas*(xx(2)*vv(3)-xx(3)*vv(2))
       J(2) = J(2) + tmas*(xx(3)*vv(1)-xx(1)*vv(3))
       J(3) = J(3) + tmas*(xx(1)*vv(2)-xx(2)*vv(1))
    end do

    write(15,'(3F15.8)') J

  end subroutine

end module angular
