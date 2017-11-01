program main
! This program is for HWWW -> H + WWW reaction. |J| for HCl will be calculated
! and ZPE will be checked.

use final_condition
  use pes_shell
  use constants
  use angular
  implicit none
real,parameter::Be_HCl=10.59341 !cm-1
real,parameter::aums=2187695.176776!1au = aums m/s
real,parameter::zpe_hcl=1483.743446 !From DMC
real,parameter::zpe_www=15602.536472 !From DMC
real,parameter::zpe_hwww = 17683.110522 !DMC = 17683.110522

real,parameter::E_hcl_ref =8463.2867115594708594 
real,parameter::E_www_ref =2796.6585222145845364 
real,parameter::E_given =  21233.1103893624
real,dimension(:,:),allocatable::xx,gradd,velo
character(len=2),dimension(:),allocatable::sym
character(len=32)::filename
integer::k,i,natm,ierr,j,vio
real::pot,diss, E_trans,E_T, E_V,  M_hcl, E_rot, speed, Evib
real::speed_www, Erot_www,Evib_www
real,dimension(:),allocatable::mass
integer::ab_j
real,dimension(:),allocatable::x
real::j_temp
real,dimension(27)::x_www
real,dimension(27)::v_www
real,dimension(6)::x_hcl
real,dimension(6)::v_hcl
real::x_hwww(33), v_hwww(33), kin_hwww, E_tot, D0
real::abc(3), abc_www(3), D0count(4), red_w, red_h

!For HCl fragment
real,dimension(2)::mass_hcl !Redefine mass for hcl
real,dimension(9)::mass_www !Redefine mass for www

real,dimension(3,9)::www_r, www_v !Coordiante for hcl and velocity for hcl
real,dimension(3,9)::www_r_prime, www_v_prime !Coordiante and velocity in com

real,dimension(3,2)::hcl_r, hcl_v !Coordiante for hcl and velocity for hcl
real,dimension(3,2)::hcl_r_prime, hcl_v_prime !Coordiante and velocity in com

real,dimension(3)::www_com_r, www_com_v!Vectos to define com coord and velocity
real,dimension(3)::j_www, omega_www !Angular momentum of www

real,dimension(3)::hcl_com_r, hcl_com_v!Vectos to define com coord and velocity
real,dimension(3)::j_hcl, omega, j_hcl_calc !Angular momentum of HCl
real::j_calc
real,dimension(2,3)::xxx,vv
real::Erot,J1,J1_sum,J2,J2_sum,J3,J3_sum,com_velc(3), kine_h, kine_www
integer::count1,count2,count3
integer::Jcount(4),iflag(4),Jrot(4),Jsum(4),vsum(4)
!real::speed(4)
Jcount = 0
Jsum = 0
Jrot = 0
vsum=0
D0count = 0


  call getarg(1,filename)
  open(21,status='old',file=filename)
  i=index(filename,'.',.true.)
  filename=filename(1:i-1)
  open(34,status='unknown',file=trim(filename)//"_s4.geom")!Hard restriction


  open(22,status='unknown',file=trim(filename)//"_s1.txt")!No restriction
  open(23,status='unknown',file=trim(filename)//"_s2.txt")!HCl restriction
  open(24,status='unknown',file=trim(filename)//"_s3.txt")!Soft restriction
  open(25,status='unknown',file=trim(filename)//"_s4.txt")!Hard restriction

  open(26,status='unknown',file=trim(filename)//"_s1_j4.txt")!J=4
  open(27,status='unknown',file=trim(filename)//"_s2_j4.txt")
  open(28,status='unknown',file=trim(filename)//"_s3_j4.txt")
  open(29,status='unknown',file=trim(filename)//"_s4_j4.txt")
  open(30,status='unknown',file=trim(filename)//"_s1_j6.txt")!J=6
  open(31,status='unknown',file=trim(filename)//"_s2_j6.txt")
  open(32,status='unknown',file=trim(filename)//"_s3_j6.txt")
  open(33,status='unknown',file=trim(filename)//"_s4_j6.txt")

!write(22,*) '# No ZPE restriction, Energy (cm-1), speed (m/s)'

!write(23,*) '# HCl ZPE restriction,  Energy (cm-1), speed (m/s)'
!write(24,*) '# Soft ZPE restriction, Energy (cm-1), speed (m/s)'
!write(25,*) '# Hard ZPE restriction,  Energy (cm-1), speed (m/s)'
write(22,'(3A15)')  '#No_Constraint','#No_Constraint','#No_Constraint'
write(23,'(3A15)')  '#Hard_on_HCl', '#Hard_on_HCl','#Hard_on_HCl'
write(24,'(3A15)')  '#Soft_ZPE','#Soft_ZPE','#Soft_ZPE'
write(25,'(3A15)')  '#Hard_ZPE','#Hard_ZPE','#Hard_ZPE'

write(26,'(3A15)')  '#No_Constraint','#No_Constraint','#No_Constraint'
write(27,'(3A15)')  '#Hard_on_HCl', '#Hard_on_HCl','#Hard_on_HCl'
write(28,'(3A15)')  '#Soft_ZPE','#Soft_ZPE','#Soft_ZPE'
write(29,'(3A15)')  '#Hard_ZPE','#Hard_ZPE','#Hard_ZPE'

write(30,'(3A15)')  '#No_Constraint','#No_Constraint','#No_Constraint'
write(31,'(3A15)')  '#Hard_on_HCl', '#Hard_on_HCl','#Hard_on_HCl'
write(32,'(3A15)')  '#Soft_ZPE','#Soft_ZPE','#Soft_ZPE'
write(33,'(3A15)')  '#Hard_ZPE','#Hard_ZPE','#Hard_ZPE'

do j=1,12
!write(21+j,*) '#HCl ZPE (DMC): 1483.743446, WWW ZPE (DMC): 15602.536472'
!write(21+j,*) '#E_given = DMC ZPE + 3550 =21233.1103893624.'
!write(21+j,*) '#Dvib = Evib - ZPE'
!write(21+j,*) "#D0 = E_given - Erot(W) - Erot(HWW) - Evib(W) + Evib(HWW) -E_COM,&
! E_COM=0"
!write(21+j,*)
!For self readiing
!write(21+j,'(8A15,A2)')  '#Erot(HCl)','#Erot(WWW)',&
!'#Dvib(HCl)', '#Dvib(WWW)','D0', '#Speed','#J'
!For export
write(21+j,'(3A15)')  '#Erot(HCl)cm-1','#J','#Speed(m/s)'
!write(21+j,'(3A15)')  '#Erot(HCl) (cm-1)','#J','#Speed(m/s)'


end do

  do i=1,4
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

call hcl_init(1)
call wat_init(3)

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

!For whole
do j=1,11
x_hwww(3*j-2:3*j) = xx(:,j)
v_hwww(3*j-2:3*j) = velo(:,j)
end do





! Kinetic energy of HWWW, this should be 0
com_velc=0

    do i=1,natm
       com_velc=com_velc+mass(i)*v_hwww(3*i-2:3*i)
    end do
com_velc=com_velc/sum(mass)
kin_hwww = calc_kine(mass,com_velc)

call xcom(mass,x_hwww)
call vcom(mass,v_hwww)



!Assign cooridnate and velocity to hcl, as 1d
do j=1,2
x_hcl(3*j-2:3*j) = x_hwww(3*(9+j)-2:3*(9+j))
v_hcl(3*j-2:3*j) = v_hwww(3*(9+j)-2:3*(9+j))
end do

!For water
do j=1,9
x_www(3*j-2:3*j) = x_hwww(3*j-2:3*j)
v_www(3*j-2:3*j) = v_hwww(3*j-2:3*j) 
end do



!Calculate final condition of HCl

call fin_cond(mass(10:11),x_hcl,v_hcl,speed,kine_h, Erot,j_hcl,abc)
call fin_cond(mass(1:9),x_www,v_www,speed_www, kine_www, Erot_www,j_www,abc_www)

!j_hcl_calc = &
!vec_cross(x_hcl(1:3),v_hcl(1:3)*mass(10))+vec_cross(x_hcl(4:6),v_hcl(4:6)*mass(11))
!j_calc = sqrt(sum(j_hcl_calc**2))
!write(*,*) '1', Erot*aucm
!write(*,*) '2', (j_calc+1)*j_calc*Be_HCl
!Be_HCl

! Calculate the rotation constant |J|

!ab_j = nint(sqrt(0.25 + ab_j)-0.5)
ab_j = nint(sqrt(sum(j_hcl**2))-0.5) !Another way to get ab_j


!j_calc = sqrt(0.25 + ab_j)-0.5
!write(*,*) 'compare', ab_j
Erot = ab_j*(ab_j+1)*Be_HCl 
!write(*,*) Erot*aucm,j_calc*(j_calc+1)*Be_HCl


!calculate vibrational energy, so that to determin if zpe violat
nwat = 0
nhcl = 1!Calculate 1 hcl only
Evib =calc_kine(mass(10:11),v_hcl) - Erot/aucm
Evib = Evib*aucm + (f(xx(:,10:11))*aucm - E_hcl_ref)
!For www
nwat = 3!Calculate water trimer
nhcl = 0
Evib_www =calc_kine(mass(1:9),v_www) - Erot_www
Evib_www = Evib_www*aucm + (f(xx(:,1:9))*aucm - E_www_ref)

!if (k .eq. 3) then
!stop
!end if


!Claasification
iflag = 0
iflag(1) = 1
if (Evib > zpe_hcl) then !Only HCl
iflag(2) = 1
end if
if (Evib + Evib_www > zpe_www + zpe_hcl) then !Soft ZPE
iflag(3) = 1
end if
if (Evib > zpe_hcl .and. Evib_www > zpe_www) then !Hard ZPE
iflag(4) = 1


end if

!Rotational energy
!Erot = Erot * aucm
Erot_www = Erot_www * aucm

!Excessive vibrational energy
red_w = Evib_www - zpe_www
red_h = Evib - zpe_hcl


!ikin_hwww = kin_hwww * aucm
!E_tot = Erot+Erot_www+Evib+Evib_www+kin_hwww
D0 = E_given-zpe_hwww   - red_h - red_w-Erot-Erot_www&
-kine_h-kine_www
write(*,*) D0
! Write files for different J condition
do i = 1,4
if (iflag(i) .eq. 1)  then!not violate zpe of hcl
Jsum(i) = Jsum(i) + ab_j!
Jcount(i) = Jcount(i) + 1
vsum(i) = vsum(i) + speed*aums
D0count(i) = D0count(i) + D0

!Record
! Record overal speed distrubtion and j distribution
!write(21+j,'(3A15)')  '#Erot(HCl) (cm-1)','#J','#Speed(m/s)'
!For export
write(21+i,'(F15.2,I15,F15.2)')  Erot,ab_j, speed*aums

!For selfuse
!write(21+i,'(6(F15.2),I15)')  Erot, Erot_www, &
!Evib - zpe_hcl, Evib_www-zpe_www, D0, speed*aums, ab_j

!For J=4 configs
if (abs(ab_j-4.0) <=1d-5) then
write(25+i,'(F15.2,I15,F15.2)')  Erot,ab_j, speed*aums
if (iflag(4) .eq. 1 .and. i .eq. 4) then
!Find one of the good ZPE result
!write(*,*) k, ab_j
write(34,*) 11
write(34,*) k
do j=1,natm
write(34,*) sym(j),xx(:,j)*auang
end do
end if
!write(25+i,'(6(F15.2),I15)')  Erot, Erot_www, &
!Evib - zpe_hcl, Evib_www-zpe_www, D0, speed*aums, ab_j
end if

!For J=6 configs
if (abs(ab_j-6.0) <=1d-5) then
write(29+i,'(F15.2,I15,F15.2)')  Erot,ab_j, speed*aums
!write(29+i,'(6(F15.2),I15)')  Erot, Erot_www, &
!Evib - zpe_hcl, Evib_www-zpe_www, D0, speed*aums, ab_j
end if

end if
end do
end do

write(*,*) 'Jcount,no/hard on hcl/soft/hard',Jcount
write(*,'(A20,4F8.2)') 'Rejection Rate:',1-real(Jcount)/real(Jcount(1))
write(*,'(A20,4F8.2)') 'Average J: ',Jsum(:)/real(Jcount(:))
write(*,'(A20,4F8.2)') 'Average Speed (m/s): ',vsum(:)/real(Jcount(:))
write(*,'(A20,4F8.2)') 'Average D0 (cm-1): ',D0count(:)/real(Jcount(:))


end program main
