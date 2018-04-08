program main
! For water
! and ZPE will be checked.

use final_condition
  use pes_shell
  use constants
  use angular
  implicit none
real,parameter::Be_HCl=10.59341
real,parameter::aums=2187695.176776!1au = aums m/s
real,parameter::zpe_w=4637.098762
real,parameter::zpe_hww=12220.164040!DMC
real,parameter::E_w_ref =8462.8673097354421770
real,parameter::E_hww_ref = 3537.7663972273066975
real,dimension(:,:),allocatable::xx,gradd,velo
character(len=2),dimension(:),allocatable::sym
character(len=2),dimension(3)::sym_w
character(len=32)::filename
integer::k,i,natm,ierr,j,vio
real::pot,diss, E_trans,E_T, E_V, ab_j, M_w, E_rot, speed, Evib
real::Evib_hww, speed_hww
real,dimension(:),allocatable::mass
real,dimension(:),allocatable::x
real,dimension(9)::x_w
real,dimension(9)::v_w
real,dimension(24)::x_hww
real,dimension(24)::v_hww
real,dimension(3,3)::xxw
integer::j221c,j321c

!For H2O fragment
real,dimension(3)::mass_w !Redefine mass for w
real,dimension(8)::mass_hww !Redefine mass for w
real,dimension(3)::abc,abc_hww

real,dimension(3,3)::w_r, w_v !Coordiante for w and velocity for w
real,dimension(3,3)::w_r_prime, w_v_prime !Coordiante and velocity in com

real,dimension(8,3)::hww_r, hww_v !Coordiante for w and velocity for w
real,dimension(8,3)::hww_r_prime, hww_v_prime !Coordiante and velocity in com

real,dimension(3)::w_com_r, w_com_v!Vectos to define com coord and velocity
real,dimension(3)::j_w, omega !Angular momentum of HCl

real,dimension(3)::hww_com_r, hww_com_v!Vectos to define com coord and velocity
real,dimension(3)::j_hww, omega_hww !Angular momentum of HWW

real,dimension(3,3)::xxx,vv
real,dimension(3,8)::xxhww
real::Erot
real::Erot_hww
integer:: iflag
real::J1(3),J1_hww(3)

vio = 0
j221c=0
j321c=0

  call getarg(1,filename)
  open(21,status='old',file=filename)
  i=index(filename,'.',.true.)
  filename=filename(1:i-1)
  open(22,status='unknown',file=trim(filename)//".data")
  open(23,status='unknown',file=trim(filename)//".abc")
  open(24,status='unknown',file=trim(filename)//".j221")
  open(25,status='unknown',file=trim(filename)//".j321")

  open(26,status='unknown',file=trim(filename)//".zpe_data")
  open(27,status='unknown',file=trim(filename)//".zpe_abc")
  open(28,status='unknown',file=trim(filename)//".zpe_j221")
  open(29,status='unknown',file=trim(filename)//".zpe_j321")
  write(22,*) '#Speed(m/s)      #Erot(cm-1)'
  write(26,*) '#Speed(m/s)      #Erot(cm-1)'

  write(23,*) '#a       #b      #c'
  write(27,*) '#a       #b      #c'

  read(21,*) natm
  read(21,*) 

  allocate(xx(3,natm))
  allocate(gradd(3,natm))
  allocate(velo(3,natm))
  allocate(sym(natm))
  allocate(mass(natm))
  allocate(x(natm*3))

  do i=1,natm
     read(21,*) sym(i),xx(:,i),gradd(:,i),velo(:,i)!(i*3-2:i*3)
       select case(sym(i))
       case("C")
          mass(i)=c_mass
       case("H")
          mass(i)=h_mass
       case("O")
          mass(i)=o_mass
       case("Cl")
          mass(i)=cl_mass
       end select
    end do
  rewind(21)
mass = mass*emass




!Initialize HCl fragment
call hclwat_init(1,2) !Init hww
call wat_init(1) !Init w

mass_w(1:2) = mass(1:2)
mass_w(3) = mass(7)

mass_hww(1:4) = mass(3:6)
mass_hww(5:8) = mass(8:11)


  k = 0!record total number of configs
  do
     read(21,*,iostat=ierr) natm
     if (ierr.ne.0) exit
     k = k + 1
     read(21,*)
     do i=1,natm
        read(21,*) sym(i),xx(:,i),gradd(:,i),velo(:,i)!(i*3-2:i*3)
     end do
     xx = xx/auang !xx in Bohr

!Assign cooridnate and velocity to w, as 1d
do j=1,2
x_w(3*j-2:3*j) = xx(:,j)
v_w(3*j-2:3*j) = velo(:,j)
end do
j=3
x_w(3*j-2:3*j) = xx(:,7)
v_w(3*j-2:3*j) = velo(:,7)

do j=1,4

x_hww(3*j-2:3*j) = xx(:,j+2)
v_hww(3*j-2:3*j) = velo(:,j+2)

x_hww(3*(j+4)-2:3*(j+4)) = xx(:,j+7)
v_hww(3*(j+4)-2:3*(j+4)) = velo(:,j+7)

end do

xxw(:,1:2) = xx(:,1:2)
xxw(:,3) = xx(:,7)

xxhww(:,1:4) = xx(:,3:6)
xxhww(:,5:8) = xx(:,8:11)
!Calculate final condition of HCl


call fin_cond(mass_w,x_w,v_w,speed,Erot,j_w,abc)
write(*,*) 'j subroutine:\n', j_w
J1 = 0
do i=1,8
J1 =J1+  vec_cross(x_w(3*i-2:3*i),mass_w(i)*v_w(3*i-2:3*i))
end do

write(*,*) 'j cross product',J1
!stop

!To test that indeed these 3 levels have rigorous energy order
if (j331(abc)>j321(abc) .and. j331(abc)>j322(abc) .and. j321(abc)>j322(abc)) &
then
else
write(*,*) '3,2,1: ',j331(abc)*aucm, j321(abc)*aucm, j322(abc)*aucm
cycle
end if


!calculate vibrational energy, so that to determin if zpe violated
nwat = 1
nhcl = 0
Evib = calc_kine(mass_w,v_w)- Erot 
Evib = Evib*aucm + (f(xxw)*aucm - E_w_ref )

nwat = 2
nhcl = 1
call fin_cond(mass_hww,x_hww,v_hww,speed_hww,Erot_hww,j_hww,abc_hww)
Evib_hww = calc_kine(mass_hww,v_hww)- Erot_hww
Evib_hww = Evib_hww*aucm + (f(xxhww)*aucm - E_hww_ref )

iflag = 1
if (Evib < zpe_w .or. Evib_hww < zpe_hww) then
vio = vio + 1
iflag = 0
end if
if (iflag .eq. 1) then
write(26,*) speed*aums, Erot*aucm
write(27,*) abc(1)*aucm, abc(2)*aucm, abc(3)*aucm

! For J211
if (((Erot - (j221(abc)+j211(abc))/2)>0) .and. (Erot-(j221(abc)+j220(abc))/2<0)) then
j221c = j221c + 1
write(28,*) speed*aums
end if

! For J321
if (((Erot - (j321(abc)+j322(abc))/2)>0) .and. (Erot-(j321(abc)+j331(abc))/2<0)) then
j321c = j321c + 1
write(29,*) speed*aums
end if

end if

write(22,*) speed*aums, Erot*aucm
write(23,*) abc(1)*aucm, abc(2)*aucm, abc(3)*aucm

! For J211
if (((Erot - (j221(abc)+j211(abc))/2)>0) .and. (Erot-(j221(abc)+j220(abc))/2<0)) then
j221c = j221c + 1
write(24,*) speed*aums
end if

! For J321
if (((Erot - (j321(abc)+j322(abc))/2)>0) .and. (Erot-(j321(abc)+j331(abc))/2<0)) then
j321c = j321c + 1
write(25,*) speed*aums
end if



  enddo

write(*,*) 'Total number of configs/ZPE violation: ', k, ' / ', vio
write(*,*) 'Total number of configs/j221: ', k, ' / ', j221c
write(*,*) 'Total number of configs/j321: ', k, ' / ', j321c

end program main
