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
real,parameter::E_w_ref =8462.8673097354421770 !Ref to min
real,parameter::E_hww_ref = 3537.7663972273066975 !Ref to min
real,parameter::E_given =  21233.1103893624 !Total energy
real,dimension(:,:),allocatable::xx,gradd,velo
character(len=2),dimension(:),allocatable::sym
character(len=2),dimension(3)::sym_w
character(len=32)::filename
integer::k,i,natm,ierr,j,vio
real::pot,diss, E_trans,E_T, E_V, M_w, E_rot, speed, Evib
real::Evib_hww, speed_hww
real,dimension(:),allocatable::mass
real,dimension(:),allocatable::x
real,dimension(9)::x_w
real,dimension(9)::v_w
real,dimension(24)::x_hww
real,dimension(24)::v_hww
real,dimension(3,3)::xxw
integer::j221c,j321c, D0count(4), iter, ab_j
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
real::pot_w, pot_hww, pot_hwww

  call getarg(1,filename)
open(19, status='old',file='watpot/result_HWW_W.wat') !Read water potential
open(18, status='old',file='hwwpot/result_HWW_W.hww') !Read HWW potential
  open(21,status='old',file=filename)
  i=index(filename,'.',.true.)
  filename=filename(1:i-1)

!These files are for recording
  open(34,status='unknown',file=trim(filename)//"_s4.geom")!Hard restriction
  open(22,status='unknown',file=trim(filename)//"_s1.txt")!No restriction
  open(23,status='unknown',file=trim(filename)//"_s2.txt")!water restriction
  open(24,status='unknown',file=trim(filename)//"_s3.txt")!Soft restriction
  open(25,status='unknown',file=trim(filename)//"_s4.txt")!Hard restriction

  open(26,status='unknown',file=trim(filename)//"_s1_j221.txt")!J=221
  open(27,status='unknown',file=trim(filename)//"_s2_j221.txt")
  open(28,status='unknown',file=trim(filename)//"_s3_j221.txt")
  open(29,status='unknown',file=trim(filename)//"_s4_j221.txt")
  open(30,status='unknown',file=trim(filename)//"_s1_j321.txt")!J=321
  open(31,status='unknown',file=trim(filename)//"_s2_j321.txt")
  open(32,status='unknown',file=trim(filename)//"_s3_j321.txt")
  open(33,status='unknown',file=trim(filename)//"_s4_j321.txt")

write(22,'(3A15)')  '#No_Constraint','#No_Constraint','#No_Constraint'
write(23,'(3A15)')  '#Hard_water','#Hard_water','#Hard_water'
write(24,'(3A15)')  '#Soft_ZPE', '#Soft_ZPE', '#Soft_ZPE'
write(25,'(3A15)')  '#Hard_ZPE','#Hard_ZPE','#Hard_ZPE'

write(26,'(3A15)')  '#No_Constraint','#No_Constraint','#No_Constraint'
write(27,'(3A15)')  '#Hard_water','#Hard_water','#Hard_water'
write(28,'(3A15)')  '#Soft_ZPE', '#Soft_ZPE', '#Soft_ZPE'
write(29,'(3A15)')  '#Hard_ZPE','#Hard_ZPE','#Hard_ZPE'

write(30,'(3A15)')  '#No_Constraint','#No_Constraint','#No_Constraint'
write(31,'(3A15)')  '#Hard_water','#Hard_water','#Hard_water'
write(32,'(3A15)')  '#Soft_ZPE', '#Soft_ZPE', '#Soft_ZPE'
write(33,'(3A15)')  '#Hard_ZPE','#Hard_ZPE','#Hard_ZPE'

do j=1,12
!write(21+j,*) '#Water ZPE (DMC): 4713.0629051, HWW ZPE (DMC): 12528.93462455'
!write(21+j,*) '#E_given = DMC ZPE + 3550 =21233.1103893624.'
!write(21+j,*) '#Dvib = Evib - ZPE'
!write(21+j,*) "#D0 = E_given - Erot(W) - Erot(HWW) - Evib(W) + Evib(HWW) - E_COM,&
! E_COM=0"
!write(21+j,*)

write(21+j,'(3A15)')  '#Erot(W)cm-1','#J',  '#Speed(m/s)'
!write(21+j,'(7A15,A2)')  '#Erot(W)','#Erot(HWW)',&
!'#Dvib(W)', '#Dvib(HWW)','D0', '#Speed'!,'#J'
end do


!Process begins here
!initialize
Jall_real_sum = 0
D0count = 0

  read(21,*) natm
  read(21,*) 

  allocate(xx(3,natm))
  allocate(gradd(3,natm))
  allocate(velo(3,natm))
  allocate(sym(natm))
  allocate(mass(natm))
  allocate(x(natm*3))

  do i=1,natm
     read(21,*) sym(i),xx(:,i),gradd(:,i),velo(:,i)
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
mass_w(1:2) = mass(1:2)
mass_w(3) = mass(7)
mass_hww(1:4) = mass(3:6)
mass_hww(5:8) = mass(8:11)

  call clust_break(sym)
  call pes_init()



!Loop throught all configurations

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

!For whole HWWW, Flaten to 1D
do j=1,11
x_hwww(3*j-2:3*j) = xx(:,j)
v_hwww(3*j-2:3*j) = velo(:,j)
end do
pot_hwww = f(xx) * aucm
!D0 = f(xx) * aucm + calc_kine(mass,v_hwww)*aucm
!write(*,*) D0
!Total energy using full potential!Checked this is conserved

!Kinetic energy of HWWW, verify this hould be 0
com_velc=0
    do i=1,natm
       com_velc=com_velc+mass(i)*v_hwww(3*i-2:3*i)
    end do
com_velc=com_velc/sum(mass)
!kin_hwww = calc_kine(mass,com_velc)
kin_hwww = calc_kine(mass,v_hwww)*aucm
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

!Calculate rotational energy Erot and translational energy Ekine
!For water
call fin_cond(mass_w,x_w,v_w,speed,Ekine,Erot,j_w,abc,evect)
ab_j = nint(sqrt(0.25 + sum(j_w**2))-0.5)
call fin_cond(mass_hww,x_hww,v_hww,speed_hww,Ekine_hww,Erot_hww,j_hww,abc_hww,evect)

!calculate vibrational energy, so that to determin if zpe violated
read(19,*) pot_w
Evib = calc_kine(mass_w,v_w)- Erot 
Evib = Evib*aucm + pot_w 

read(18,*) pot_hww
Evib_hww = calc_kine(mass_hww,v_hww)- Erot_hww
Evib_hww = Evib_hww*aucm + pot_hww 




!Classification for different ZPE condition
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
!Print ZPE configs
write(34,*) 11
write(34,*) k
do i=1,natm
write(34,*) sym(i),xx(:,i)*auang
end do


end if

!Rotational energy
Erot = Erot * aucm
Erot_hww = Erot_hww * aucm

!Excessive vibrational energy, red=redundent
red_hww = Evib_hww - zpe_hww
red_w = Evib - zpe_w


!D0 according to current result
!De
!D0 = (E_given) - (Erot + Erot_hww) - (Evib + Evib_hww)&
!-(Ekine + Ekine_hww)

D0 =  pot_hwww - pot_w - pot_hww
!D0 = E_given - ( kin_hwww + pot_w + pot_hww)


do i = 1,4 !Record according to different ZPE conditon

if (iflag(i) .eq. 1)  then
!For J221
if ((abs(ab_j-2.0) <=1d-5) .and. (115.669905 < Erot) .and. (Erot < 135.532765) )then
write(25+i,'(F15.2,I15,F15.2)')  Erot,ab_j, speed*aums
end if
!For J321
if ((abs(ab_j-3.0) <=1d-5) .and. (209.22888 < Erot) .and. (Erot < 248.78747) )then
write(29+i,'(F15.2,I15,F15.2)')  Erot,ab_j, speed*aums
end if


Jsum(i) = Jsum(i) + jall
Jcount(i) = Jcount(i) + 1
vsum(i) = vsum(i) + speed*aums
Jall_real_sum(i) = Jall_real_sum(i) + Jall_real
D0count(i) = D0count(i) + D0

write(21+i,'(F15.2,I15,F15.2)')  Erot,ab_j, speed*aums
!write(21+i,'(6(F15.2))')  Erot, Erot_hww, &
!Evib - zpe_w, Evib_hww-zpe_hww, D0, speed*aums



end if
end do


  enddo

write(*,*) 'Jcount,no/hard on hcl/soft/hard',Jcount
write(*,'(A20,4F8.2)') 'Rejection Rate:',1-real(Jcount)/real(Jcount(1))
write(*,'(A20,4F8.2)') 'Average J: ',Jsum(:)/real(Jcount(:))
write(*,'(A20,4F8.2)') 'Average Speed (m/s): ',vsum(:)/real(Jcount(:))
write(*,'(A20,4F8.2)') 'Average D0 (cm-1): ',D0count(:)/real(Jcount(:))


end program main
