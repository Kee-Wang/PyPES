module final_condition
!This is mod that contains tools for calcualte finaly condition. Initial
!desigend for HCl-WWW system
use constants
use angular
implicit none
!=====

contains
real function j221(abc)
real::abc(3),a,b,c
A = abc(1)
B = abc(2)
C = abc(3)
j221 = 4*A + B + C
end function

real function j211(abc)
real::abc(3),a,b,c
A = abc(1)
B = abc(2)
C = abc(3)
j211 = A + 4*B + C
end function

real function j220(abc)
real::abc(3),a,b,c
A = abc(1)
B = abc(2)
C = abc(3)
j220 = 2*A + 2*B +2*C + 2*sqrt((B-C)**2 + (A-C)*(A-B)) 
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
subroutine fin_cond(mass,x,v,speed, Ekine, Erot,js,abc)
!All units are in a.u.
real,dimension(:),intent(in)::mass!Mass
real,dimension(:),intent(inout)::x!Cart of frag in 1d. When out, in com frame
real,dimension(:),intent(inout)::v!Velocity of frag in 1d. When out, in com frame.
real,intent(out)::speed!Speed of frag
real,intent(out)::Erot!Rotational energy of frag
real,dimension(3),intent(out)::js!angular momentum
!:::::::::::::::;
real,dimension(3)::com_velc !Center of mass speed
integer::i,natm
real,dimension(size(mass),3)::xx,vv
real,dimension(3)::omega,abc
real::Ekine
!Find com velocity of calculate speed.
natm = size(mass)
Erot=0
com_velc=0

    do i=1,natm
       com_velc=com_velc+mass(i)*v(3*i-2:3*i)
    end do
com_velc=com_velc/sum(mass)
speed =  sqrt((dot_product(com_velc,com_velc)))
Ekine = calc_kine(mass, v)

    call xcom(mass,x)
    call vcom(mass,v)

    do i=1,natm
       xx(i,:)=x(3*i-2:3*i)
       vv(i,:)=v(3*i-2:3*i)
    end do

    call rotn(0,natm,xx,vv,mass,js,erot,omega,abc)




end subroutine


function vec_cross(v1,v2) result(j)
real,dimension(3)::v1,v2,j

j(1) =    v1(2)*v2(3) - v1(3)*v2(2)
j(2) = -( v1(1)*v2(3) - v1(3)*v2(1) ) 
j(3) =    v1(1)*v2(2) - v1(2)*v2(1)

end function


real function f_kine(m,v)
real::m
real,dimension(3)::v
f_kine = 0.5 * m * ( v(1)**2 + v(2)**2 + v(3)**2)
 end function



subroutine com_cart(xx,mass,com_x)
real,dimension(:,:),intent(inout)::xx
real,dimension(:),intent(in)::mass
    ! ::::::::::::::::::::
integer::i,n
real,dimension(3),intent(out)::com_x !com coordinate,(x,y,z)

    com_x=0                  !init the com coordinates:set to 0
   do i=1,3 !(for x,y,z)
    com_x(i) = dot_product(mass,xx(i,:))
  end do
    com_x=com_x/sum(mass)
end subroutine


subroutine com_velo(velo,mass,com_v)
real,dimension(:,:),intent(inout)::velo
real,dimension(:),intent(in)::mass
    ! ::::::::::::::::::::
integer::i,n
real,dimension(3),intent(out)::com_v !com coordinate,(x,y,z)

    com_v=0                  !init the com coordinates:set to 0
   do i=1,3 !(for x,y,z)
    com_v(i) = dot_product(mass,velo(i,:))
  end do
    com_v=com_v/sum(mass)
end subroutine



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
!write(*,*) 'v',v

    return
  end subroutine vcom
end module
