module normal_mode
  use pes_shell
  use angular
  implicit none

contains

  subroutine e_rot(xe,mass,x,v,secmode)
    !=================================================!
    ! given the initial geometry, gives out current rotational energy
    !=================================================!
    real,dimension(:),intent(in)   :: xe,mass
    real,dimension(:),intent(in)   :: secmode 
    ! A vector showing which mode should add the zpe to
    real,dimension(:),intent(inout):: x,v
    ! ::::::::::::::::::::
    real,dimension(1:size(v))::qx,qv,dx ! position and velocity in normal coord.
    real,dimension(size(mass),3)::xx,vv
    integer::dim,natm,i,imax
    real,dimension(3)::js,dj,omega
    real::zpe,pot,pot0,e1,e0,erot,gamma,sdum

    dim=size(secmode)
    natm=size(mass)

    zpe=sum(freq*secmode/fmd_aucm)/2.0
    e0 = zpe
  write(*,*) 'zpe(cm-1)=',zpe*219474.63


    ! calculate the potential of the equilibrium
    pot0 = fho(x)
!    pot0 = f(x)

    ! assign ZPE to each mode
    qx=0.0
    qv=0.0
    call aszpe(freq,qx,qv)
    call norm2card(qv,mass,nvec,v)
    call norm2card(qx,mass,nvec,dx)
    x=xe+dx

    fmd_time = 0.0

    ! 1. calculate suprious angular momentum (js)
    !  1.1 move to center of mass
50  call xcom(mass,x)
    call vcom(mass,v)

    !  1.2 calculate angular momentum vector
    js=0.0
    do i=1,natm
       xx(i,:)=x(3*i-2:3*i)
       vv(i,:)=v(3*i-2:3*i)
    end do

    do i=1,natm
       js(1) = js(1) + mass(i)*(xx(i,2)*vv(i,3)-xx(i,3)*vv(i,2))
       js(2) = js(2) + mass(i)*(xx(i,3)*vv(i,1)-xx(i,1)*vv(i,3))
       js(3) = js(3) + mass(i)*(xx(i,1)*vv(i,2)-xx(i,2)*vv(i,1))
    end do

    ! 2. substract js from the desired angular momentum j0
    !    and calculate the omega needed to add in
    dj = -js
    call rotn(1,natm,xx,vv,mass,dj,erot,omega)

    ! 3. add in dj so that J=js+dj=j0
    do i=1,natm
       v(3*i-2) = vv(i,1) + xx(i,3)*omega(2) - xx(i,2)*omega(3)
       v(3*i-1) = vv(i,2) + xx(i,1)*omega(3) - xx(i,3)*omega(1)
       v(3*i)   = vv(i,3) + xx(i,2)*omega(1) - xx(i,1)*omega(2)
    end do

    do i=1,natm
       xx(i,:)=x(3*i-2:3*i)
       vv(i,:)=v(3*i-2:3*i)
    end do
    call rotn(0,natm,xx,vv,mass,js,erot,omega)

    ! 4. calculate the energy
    pot = fho(x)
!    pot = f(x)
    e1 = pot + calc_kine(mass,v) - pot0 ! e1 is the current internal energy

    ! 5. scale
    sdum=dabs(e1-e0)/e0
    i=0; imax=30
    if (sdum .gt. 0.00001d0) then
       i=i+1
       if (i .gt. imax) stop 'Too many scaling cycles'
       gamma=dsqrt(e0/e1)

       x = xe + (x-xe)*gamma
       v = v * gamma
       go to 50
    end if

    return
  end subroutine nm_sample

  subroutine nm_sample(xe,mass,x,v,secmode)
    !=================================================!
    ! given the initial equilibrium geometry and mass,!
    ! do the Normal Mode Analysis and return the      !
    ! frequency, normal vectors and initial velocities!
    !=================================================!
    real,dimension(:),intent(in)   :: xe,mass
    real,dimension(:),intent(in)   :: secmode 
    ! A vector showing which mode should add the zpe to
    real,dimension(:),intent(inout):: x,v
    ! ::::::::::::::::::::
    real,dimension(1:size(v))::qx,qv,dx ! position and velocity in normal coord.
    real,dimension(size(mass),3)::xx,vv
    integer::dim,natm,i,imax
    real,dimension(3)::js,dj,omega
    real::zpe,pot,pot0,e1,e0,erot,gamma,sdum

    dim=size(secmode)
    natm=size(mass)

    zpe=sum(freq*secmode/fmd_aucm)/2.0
    e0 = zpe
  write(*,*) 'zpe(cm-1)=',zpe*219474.63


    ! calculate the potential of the equilibrium
    pot0 = fho(x)
!    pot0 = f(x)

    ! assign ZPE to each mode
    qx=0.0
    qv=0.0
    call aszpe(freq,qx,qv)
    call norm2card(qv,mass,nvec,v)
    call norm2card(qx,mass,nvec,dx)
    x=xe+dx

    fmd_time = 0.0

    ! 1. calculate suprious angular momentum (js)
    !  1.1 move to center of mass
50  call xcom(mass,x)
    call vcom(mass,v)

    !  1.2 calculate angular momentum vector
    js=0.0
    do i=1,natm
       xx(i,:)=x(3*i-2:3*i)
       vv(i,:)=v(3*i-2:3*i)
    end do

    do i=1,natm
       js(1) = js(1) + mass(i)*(xx(i,2)*vv(i,3)-xx(i,3)*vv(i,2))
       js(2) = js(2) + mass(i)*(xx(i,3)*vv(i,1)-xx(i,1)*vv(i,3))
       js(3) = js(3) + mass(i)*(xx(i,1)*vv(i,2)-xx(i,2)*vv(i,1))
    end do

    ! 2. substract js from the desired angular momentum j0
    !    and calculate the omega needed to add in
    dj = -js
    call rotn(1,natm,xx,vv,mass,dj,erot,omega)

    ! 3. add in dj so that J=js+dj=j0
    do i=1,natm
       v(3*i-2) = vv(i,1) + xx(i,3)*omega(2) - xx(i,2)*omega(3)
       v(3*i-1) = vv(i,2) + xx(i,1)*omega(3) - xx(i,3)*omega(1)
       v(3*i)   = vv(i,3) + xx(i,2)*omega(1) - xx(i,1)*omega(2)
    end do

    do i=1,natm
       xx(i,:)=x(3*i-2:3*i)
       vv(i,:)=v(3*i-2:3*i)
    end do
    call rotn(0,natm,xx,vv,mass,js,erot,omega)

    ! 4. calculate the energy
    pot = fho(x)
!    pot = f(x)
    e1 = pot + calc_kine(mass,v) - pot0 ! e1 is the current internal energy

    ! 5. scale
    sdum=dabs(e1-e0)/e0
    i=0; imax=30
    if (sdum .gt. 0.00001d0) then
       i=i+1
       if (i .gt. imax) stop 'Too many scaling cycles'
       gamma=dsqrt(e0/e1)

       x = xe + (x-xe)*gamma
       v = v * gamma
       go to 50
    end if

    return
  end subroutine nm_sample


  subroutine aszpe(freq,qx,qv)
    !====================================!
    ! ASsign Zero-Point Energy to        !
    ! the ith normal mode                !
    !====================================!
    real,dimension(:),intent(in)   ::freq  ! Frequency of a specific normal mode
    real,dimension(:),intent(inout)::qx,qv ! Velocities in Normal Coordinates
    ! ::::::::::::::::::::
    real,dimension(1:size(freq))::zpe,w
    real::ran,pi
    integer::dim,i
    
    pi=acos(-1.0)

    w=freq/fmd_aucm
    zpe=0.5*abs(w)*secmode   !E=(n+1/2)hv, where n is the quantum number
 
    dim=size(freq)
    do i=1,dim
       if(w(i)==0.0) cycle
       call random_number(ran)
       qv(i)=-sqrt(2*zpe(i))*sin(2*pi*ran)
       qx(i)= sqrt(2*zpe(i))/w(i)*cos(2*pi*ran)
    end do
    
    return
  end subroutine aszpe
  
  subroutine card2norm(x,xe,mass,nvec,bQ)
    !================================================!
    ! convert the usual cardesian coordinates        !
    ! to the normal coordinates                      !
    !================================================!
    real,dimension(:),intent(in)::x,xe,mass
    real,dimension(:,:),intent(in)::nvec
    real,dimension(:),intent(out)::bQ
    ! ::::::::::::::::::::
    real,dimension(1:size(x))::q
    
    call card2mswt(x,xe,mass,q)
    call mswt2norm(q,nvec,bQ)
    
    return
  end subroutine card2norm

  subroutine norm2card(bQ,mass,nvec,x)
    !==============================================!
    ! convert the normal coordinates to the        !
    ! usual cardesian coordinates                  !
    !==============================================!
    real,dimension(:),intent(in)::bQ !Normal Coordinates
    real,dimension(:),intent(in)::mass
    real,dimension(:,:),intent(in)::nvec
    real,dimension(:),intent(out)::x
    ! ::::::::::::::::::::
    real,dimension(1:size(bQ))::q !mass weighted coordinates
    integer::i
    
    call norm2mswt(bQ,nvec,q)
    call mswt2card(q,mass,x)

    return
  end subroutine norm2card

  subroutine mswt2norm(q,L,bQ)
    !====================================================!
    ! convert the mass weighted coordinate to the        !
    ! Normal Coordinate                                  !
    !====================================================!
    real,dimension(:),intent(in)::q
    real,dimension(:,:),intent(in)::L
    real,dimension(:),intent(out)::bQ
    ! ::::::::::::::::::::
    
    ! Q=L'q
    bQ=matmul(transpose(L),q)

    return
  end subroutine mswt2norm
  
  subroutine norm2mswt(bQ,nvec,q)
    !=============================================!
    ! convert the normal coordinate to the        !
    ! mass weighted coordinate                    !
    !=============================================!
    real,dimension(:),intent(in)::bQ
    real,dimension(:,:),intent(in)::nvec ! L
    real,dimension(:),intent(out)::q
    ! ::::::::::::::::::::
    
    ! q=LQ
    q=matmul(nvec,bQ)
    
    return
  end subroutine norm2mswt

  subroutine mswt2card(q,mass,x)
    !=================================================!
    ! convert the mass weighted coordinates to        !
    ! the cardesian coordinates                       !
    !=================================================!
    real,dimension(:),intent(in)::q,mass
    real,dimension(:),intent(out)::x
    ! ::::::::::::::::::::
    integer::i,natm

    natm=size(mass)

    do i=1,natm
       x(3*i-2:3*i)=q(3*i-2:3*i)/sqrt(mass(i))
    end do
    
    return
  end subroutine mswt2card

  subroutine card2mswt(x,xe,mass,q)
    !=============================================!
    ! convert the Cardisian Coordinates to        !
    ! mass weighted coordinates                   !
    !=============================================!
    real,dimension(:),intent(in) ::x,xe,mass
    real,dimension(:),intent(out)::q
    ! ::::::::::::::::::::
    real,dimension(1:size(x))::d
    integer::i,natm

    natm=size(mass)
 
    d=x-xe
    do i=1,natm
       q(3*i-2:3*i)=sqrt(mass(i))*d(3*i-2:3*i)
    end do
    
    return
  end subroutine card2mswt

  subroutine mw_hessian(p,mass,H)
    real,dimension(:),intent(in)::p
    real,dimension(:,:),intent(inout)::H
    real,dimension(:),intent(in)::mass
    !:::::::::::::::::::::::::::::
    real,dimension(3,size(p)/3)::xx
    real::f_ff,f_fb,f_bf,f_bb,fp
    real::rmass
    real,dimension(1:size(p))::tp
    integer::dim,i,j,ii,natm
    real,parameter::det=5.0d-3

    natm=size(p)/3

    tp=p
    fp=f(tp)

    dim=size(p)

    do i=1,dim-1
       do j=i+1,dim
          rmass=sqrt(mass(ceiling(i/3.0)))*sqrt(mass(ceiling(j/3.0)))

          call pt(tp,i,j, 1, 1)
          f_ff=f(tp)  ! f( 1, 1)
          tp=p

          call pt(tp,i,j, 1,-1)
          f_fb=f(tp)  ! f( 1,-1)
          tp=p 

          call pt(tp,i,j,-1, 1)
          f_bf=f(tp)  ! f(-1, 1)
          tp=p

          call pt(tp,i,j,-1,-1)
          f_bb=f(tp)  ! f(-1,-1)
          tp=p

          H(i,j)=0.25*(f_ff-f_fb-f_bf+f_bb)/det/rmass/det
          H(j,i)=H(i,j)
       end do
    end do

    do i=1,dim
       rmass=mass(ceiling(i/3.0))
       call pt(tp,i,i, 0, 1)
       f_ff=f(tp)  ! f( 1)
       tp=p

       call pt(tp,i,i, 0,-1)
       f_bb=f(tp)  ! f(-1)
       tp=p
       
       H(i,i)=(f_ff-2*fp+f_bb)/det/det/rmass
    end do

    return 
  end subroutine mw_hessian

  !==================================================!
  ! move the point p in i,j direction m and n        !
  ! steps(step length is equal to opt_det            !
  !==================================================!
  subroutine pt(p,i,j,m,n)
    real,dimension(:),intent(inout)::p
    integer::i,j,m,n
    ! ::::::::::::::::::::
    real::det = 5.0d-3

    p(i)=p(i)+m*det
    p(j)=p(j)+n*det
    return
  end subroutine pt

  !============================================
  ! diagonalize the hessian matrix and return        
  ! the eigen value and eigenvectors.         
  ! The original Hessian will be destroied
  !============================================
  subroutine diag_hessian(H,w)
    real,dimension(:,:),intent(inout)::H
    real,dimension(:),intent(out)::w
    ! ::::::::::::::::::::
    real,dimension(:),allocatable::work
    integer::dim,lwork,info,i,j
    
    dim=size(H,1)
    lwork=dim*dim*10;
    allocate(work(1:lwork))
    
    call dsyev('v','u',dim,H,dim,w,work,lwork,info) 
    
    do i=1,dim
       w(i)=sign(sqrt(abs(w(i)))*fmd_aucm,w(i))
    end do
    
    return
  end subroutine diag_hessian

  !==========================================!
  !  calculate the zero-th order potential   !
  !==========================================!
  function fho(x)
    real,dimension(:),intent(in)::x ! Cartesian coordinates
    real::fho
    !::::::::::::::::::::
    real,dimension(1:size(x))::bQ ! Normal coordinate
    real,dimension(1:size(x))::xx,vv
    real::pot,t(3,3)
    integer::i

    xx = x
    vv = fmd_v
 
    call eckart_coor(fmd_natm,xx,fmd_xe,fmd_mass,vv,t)
    ! Transform to normal coordinates
    call card2norm(xx,fmd_xe,fmd_mass,nvec,bQ)

    pot = 0.d0
    do i = 7, size(x)
       pot = pot + 0.5*((freq(i)/fmd_aucm)**2)*(bQ(i)**2)
    end do

    fho = pot

    return
  end function

  !================================================!
  ! calculate the kinetic energy of current        !
  ! configuration                                  !
  !================================================!
  function calc_kine(mass,v) result(te)
    real,dimension(:)::mass
    real,dimension(:)::v
    real::te
    ! ::::::::::::::::::::
    integer::i,n,dim
    real,dimension(1:size(v,1))::sqv
    
    dim=size(v)
    
    sqv=v*v
    
    te=0
    do i=1,dim
       te=te+0.5*mass(ceiling(i/3.0))*sqv(i)
    end do
    
    return
  end function calc_kine

  !=============================!
  ! set the velocity to         !
  ! center of mass frame        !
  !=============================!
  subroutine vcom(mass,v)
    real,dimension(:),intent(in)::mass
    real,dimension(:),intent(inout)::v
    ! ::::::::::::::::::::
    integer::i,n,dim
    real,dimension(1:3)::com_velc

    dim=size(v,1)
    n=dim/3
    
    com_velc=0
    
    do i=1,n
       com_velc=com_velc+mass(i)*v(3*i-2:3*i)
    end do
    
    com_velc=com_velc/sum(mass)

    do i=1,n
       v(3*i-2:3*i)=v(3*i-2:3*i)-com_velc
    end do
    
    return
  end subroutine vcom
  
  !=============================!
  ! set the molecule to         !
  ! center of mass frame        !
  !=============================!
  subroutine xcom(mass,x)
    real,dimension(:),intent(in)::mass
    real,dimension(:),intent(inout)::x
    ! ::::::::::::::::::::
    integer::i,n,dim
    real,dimension(1:3)::com_cord !com coordinate

    dim=size(x,1)
    n=dim/3
    
    com_cord=0                  !init the com coordinates:set to 0
    
    do i=1,n
       com_cord=com_cord+mass(i)*x(3*i-2:3*i)
    end do
    
    com_cord=com_cord/sum(mass)

    do i=1,n
       x(3*i-2:3*i)=x(3*i-2:3*i)-com_cord
    end do

    return
  end subroutine xcom

  subroutine mode_pote(x,xe,mass,pote)
    !==========================================!
    ! calculate the potential energy of        !
    ! a specific normal mode                   !
    !==========================================!
    real,dimension(:),intent(in)::x,xe,mass
    real,dimension(:),intent(out)::pote
    ! ::::::::::::::::::::
    real,dimension(1:size(x))::Q,tQ,tx,dx
    real::pot1,pot0
    integer::dim,i

    dim=size(x)

    call card2norm(x,xe,mass,nvec,Q)

    do i=1,dim
       tQ=0.0
       tQ(i)=Q(i)
       call norm2card(tQ,mass,nvec,dx)
       tx=xe+dx

       pot1 = f(tx)*switch(fmd_time) + (1-switch(fmd_time))*fho(tx)
       pot0 = f(xe)*switch(fmd_time) + (1-switch(fmd_time))*fho(xe)
       pote(i) = pot1 - pot0
    end do
    
    return
  end subroutine mode_pote

  subroutine mode_kine(vx,mass,kine)
    !=========================================!
    ! calculate the kinetic energy of         !
    ! every normal mode                       !
    !=========================================!
    real,dimension(:),intent(in)::vx,mass
    real,dimension(:),intent(out)::kine
    ! ::::::::::::::::::::
    real,dimension(1:size(vx))::vQ,zero

    zero=0.0
    call card2norm(vx,zero,mass,nvec,vQ)
    kine=0.5*vQ*vQ
    
    return
  end subroutine mode_kine

  subroutine mode_enrg(vx,x,xe,mass,enrg)
    ! ==========================================!
    ! calculate the energy for each normal mode !
    ! ==========================================!
    real,dimension(:),intent(in)::vx,x,xe,mass
    real,dimension(:),intent(out)::enrg
    ! ::::::::::::::::::::
    real,dimension(1:size(enrg))::pote,kine
    real,dimension(1:size(x))::xx,vv
    real,dimension(3,3)::t

    xx = x
    vv = vx
    call eckart_coor(fmd_natm,xx,fmd_xe,fmd_mass,vv,t)
    call mode_pote(xx,xe,mass,pote)
    call mode_kine(vv,mass,kine)
    enrg=pote+kine
    
    return
  end subroutine mode_enrg

  subroutine monitor_vib(x,xe,v,mass)
    !==================================!
    ! monitor the energy in each mode  !
    !==================================!
    real,dimension(:),intent(in)::x,xe,v,mass
    ! ::::::::::::::::::::
    real,dimension(1:size(x))::enrg

    call mode_enrg(v,x,xe,mass,enrg)
    write(16,*) fmd_time,enrg(7:)*fmd_aucm

    return
  end subroutine monitor_vib

  subroutine monitor_coord()
    !======================================!
    ! Monitor the normal coodinates before !
    ! and after the Eckart transformation  !
    !======================================!
    !:::::::::::::::::::::
    real,dimension(1:size(fmd_x))::xx
    real,dimension(1:size(fmd_v))::vv
    real,dimension(3,3)::t
    real,dimension(1:size(fmd_x))::bQ1,bQ2

    xx = fmd_x
    vv = fmd_v

    call card2norm(xx,fmd_xe,fmd_mass,nvec,bQ1)
    write(17,'(A,15F7.1)') "Before:",bQ1
    
    call eckart_coor(fmd_natm,xx,fmd_xe,fmd_mass,vv,t)
    call card2norm(xx,fmd_xe,fmd_mass,nvec,bQ2)
    write(17,'(A,15F7.1)') "After :",bQ2

  end subroutine

end module normal_mode
