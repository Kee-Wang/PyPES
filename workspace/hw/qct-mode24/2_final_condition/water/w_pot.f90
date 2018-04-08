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
Jall_real_sum = 0
D0count = 0
vio = 0
j221c=0
j321c=0

  call getarg(1,filename)
  open(21,status='old',file=filename)
  i=index(filename,'.',.true.)
  filename=filename(1:i-1)
  open(22,status='unknown',file=trim(filename)//".wat")!No restriction

end do

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


  call clust_break(sym(1:3))
  call pes_init()



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
        read(21,*) sym(i),xx(:,i),gradd(:,i),velo(:,i)
     end do
     xx = xx/auang !xx in Bohr

!For potential, 2D coordiante
xxw(:,1:2) = xx(:,1:2)
xxw(:,3) = xx(:,7)
xxhww(:,1:4) = xx(:,3:6)
xxhww(:,5:8) = xx(:,8:11)

!For whole HWWW, Flaten to 1D
do j=1,11
x_hwww(3*j-2:3*j) = xx(:,j)
v_hwww(3*j-2:3*j) = velo(:,j)
end do

!Kinetic energy of HWWW, verify this hould be 0
com_velc=0
    do i=1,natm
       com_velc=com_velc+mass(i)*v_hwww(3*i-2:3*i)
    end do
com_velc=com_velc/sum(mass)
kin_hwww = calc_kine(mass,com_velc)

call xcom(mass,x_hwww)
call vcom(mass,v_hwww)

!Assign cooridnate and velocity to w, as 1d
!For W
do j=1,2
x_w(3*j-2:3*j) = x_hwww(3*j-2:3*j)
v_w(3*j-2:3*j) = v_hwww(3*j-2:3*j)
end do
j=3
x_w(3*j-2:3*j) = x_hwww(3*7-2:3*7)
v_w(3*j-2:3*j) = v_hwww(3*7-2:3*7)
!For HWW
do j=1,4
x_hww(3*j-2:3*j) = x_hwww(3*(j+2)-2:3*(j+2))
v_hww(3*j-2:3*j) = v_hwww(3*(j+2)-2:3*(j+2))
x_hww(3*(j+4)-2:3*(j+4)) =  x_hwww(3*(j+7)-2:3*(j+7))
v_hww(3*(j+4)-2:3*(j+4)) =  v_hwww(3*(j+7)-2:3*(j+7))
end do

!To this point, x_w, x_hww, v_w, v_www are all in COM frame.

!Calculate rotational energy
!For water
call fin_cond(mass_w,x_w,v_w,speed,Ekine,Erot,j_w,abc,evect)
call fin_cond(mass_hww,x_hww,v_hww,speed_hww,Ekine_hww,Erot_hww,j_hww,abc_hww,evect)

!calculate vibrational energy, so that to determin if zpe violated
!nwat = 1
!nhcl = 0
!call hclwat_init(nhcl,nwat) !Init hww
!call hcl_init(nhcl)
!call wat_init(nwat) !Init w
!Evib = calc_kine(mass_w,v_w)- Erot 
!Evib = Evib*aucm + (f(xxw,1)*aucm - E_w_ref ) !Evib for water
write(*,*) f(xxw,1)*aucm, '****'
!nwat = 2
!nhcl = 1
!call hclwat_init(nhcl,nwat) !Init hww
!call hcl_init(nhcl)
!call wat_init(nwat) !Init w
!Evib_hww = calc_kine(mass_hww,v_hww)- Erot_hww
!Evib_hww = Evib_hww*aucm + (f(xxhww,0)*aucm - E_hww_ref ) !Evib for hww

!write(*,*) f(xxhww,0)*aucm 
if (k .eq. 3) then
stop
end if


!Classification
iflag = 0
iflag(1) = 1
if (Evib > zpe_w) then !Only water
iflag(2) = 1
end if
if (Evib + Evib_hww > zpe_hww + zpe_w) then !Soft ZPE
iflag(3) = 1
end if
if (Evib > zpe_w .and. Evib_hww > zpe_hww) then !Hard ZPE
iflag(4) = 1
end if

!Rotational energy
Erot = Erot * aucm
Erot_hww = Erot_hww * aucm

!Excessive vibrational energy
red_hww = Evib_hww - zpe_hww
red_w = Evib - zpe_w


!D0 according to current result
D0 = (E_given - zpe_hwww) - (Erot + Erot_hww) - (red_w + red_hww)&
-(Ekine + Ekine_hww)


do i = 1,4 
if (iflag(i) .eq. 1)  then!not violate zpe of hcl
Jsum(i) = Jsum(i) + jall
Jcount(i) = Jcount(i) + 1
vsum(i) = vsum(i) + speed*aums
Jall_real_sum(i) = Jall_real_sum(i) + Jall_real
D0count(i) = D0count(i) + D0

!write(21+i,*) Jall_real, Erot
write(21+i,'(6(F15.2))')  Erot, Erot_hww, &
Evib - zpe_w, Evib_hww-zpe_hww, D0, speed*aums!jall
!write(21+i,'(I2,8(F12.2))') int(jall), Erot, Erot_hww, &
!Evib, Evib_hww, speed*aums, kin_hwww, E_tot, D0


!Check J321 state for soft and no constraint
!if (i .eq. 4 .or. i .eq. 1) then

!For j=321
!if ( (abs(Jall-3) .lt. 1e-5) .and. (abs(Jmoi(1)-2) .lt. 1e-5) .and. (abs(Jmoi(3)-1) .lt.1e-5)) then
!write(25+i,*) speed*aums, Erot
!write(25+i,*) Jall_real, Erot
!end if

!For j=505

end if
end do


  enddo

write(*,*) 'Jcount,no/hard on hcl/soft/hard',Jcount
write(*,'(A20,4F8.2)') 'Rejection Rate:',1-real(Jcount)/real(Jcount(1))
write(*,'(A20,4F8.2)') 'Average J: ',Jsum(:)/real(Jcount(:))
write(*,'(A20,4F8.2)') 'Average Speed (m/s): ',vsum(:)/real(Jcount(:))
write(*,'(A20,4F8.2)') 'Average D0 (cm-1): ',D0count(:)/real(Jcount(:))


end program main
