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
!ZPE calcualted use our own PES
!HO ZPE using John's calcualtion
real,parameter::zpe_w= 4713.0629051 !HO=4713.0629051 !DMC = 4637.098762
real,parameter::zpe_hww= 12528.93462455! HO = !DMC = 12220.164040
real,parameter::zpe_hwww = 17683.110522 !DMC = 17683.110522
real,parameter::E_w_ref =8462.8673097354421770
real,parameter::E_hww_ref = 3537.7663972273066975
real,parameter::E_given =  21233.1103893624
real:: E_tot
real,dimension(:,:),allocatable::xx,gradd,velo
character(len=2),dimension(:),allocatable::sym
character(len=2),dimension(8)::sym_hww
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
integer::j221c,j321c, D0count(4), iter
real:: D0, red_w, red_hww, kin_hwww, x_hwww(33), v_hwww(33), com_velc(3)
real:: Ekine, Ekine_hww

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

real,dimension(3,3)::xxx,vv,evect
real,dimension(3,8)::xxhww
real::Erot
real::Erot_hww
integer:: iflag(4), jmoi(3), jall,jsum(4),jcount(4),jac(3)
real::J1(3),J1_hww(3),vsum(4), Jall_real_sum(4)
real::Jall_real,jmoi_real(3)
  call getarg(1,filename)
  open(21,status='old',file=filename)
  i=index(filename,'.',.true.)
  filename=filename(1:i-1)
  open(22,status='unknown',file=trim(filename)//".hww")!No restriction

  read(21,*) natm
  read(21,*) 

  allocate(xx(3,natm))
  allocate(gradd(3,natm))
  allocate(velo(3,natm))
  allocate(sym(natm))
  allocate(x(natm*3))
rewind(21)
sym_hww = (/'H','H','H','H','O','O','H','Cl'/)
write(*,*) sym_hww
  call clust_break(sym_hww)
  call pes_init()

  k = 0!record total number of configs
do
     read(21,*,iostat=ierr) natm
     if (ierr.ne.0) exit
     k = k + 1
     read(21,*)
     do i=1,natm
        read(21,*) sym(i),xx(:,i),gradd(:,i),velo(:,i)
     end do
     xx = xx/auang !xx in Bohr

!For potential, 2D coordiante
xxw(:,1:2) = xx(:,1:2)
xxw(:,3) = xx(:,7)
xxhww(:,1:4) = xx(:,3:6)
xxhww(:,5:8) = xx(:,8:11)

write(22,*) f(xxhww)*aucm - E_hww_ref
!write(*,*) f(xxhww)*aucm
!if (k .eq. 3) then
!stop
!end if

end do


end program main
