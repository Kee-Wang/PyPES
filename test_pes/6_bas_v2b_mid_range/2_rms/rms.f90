program main
use pes_shell
  implicit none

!  integer,parameter::wp=selected_real_kind(12,300)
  real::auang=0.5291772083
  real::aucm=219474.63
  real::xx(3,6),pot0,diff
  character::symb(6)
  character(len=32)::filename
  integer::i,j,natm,ierr,n
  real::pot,a,b
  real::rms
real::switch1=0,switch2=0
real,allocatable::ab(:),pott(:),arms(:)

  call getarg(1,filename)

  open(21,status='old',file=trim(filename))
  i=index(filename,'.',.true.)
  filename=filename(1:i-1)
  open(22,status='unknown',file=trim(filename)//".eng")

  call pes_init()

  rms=0.0
  n = 0
  do
     read(21,*,iostat=ierr) natm
     if (ierr < 0) exit
     n = n + 1
     read(21,*) pot0
     do i=1,natm
        read(21,*) symb(i),xx(:,i)
     end do
     xx=xx/auang

     pot=f(xx) ! Calculate the potential
     diff = abs(pot-pot0)*aucm
     write(22,'(2F15.8,F13.2)') pot0,pot,diff

     rms=rms+(pot0-pot)**2

  end do

  rms=sqrt(rms/real(n))*219474.63
  write(*,*) 'RMS=',rms,'cm-1'
close (unit=22)


allocate(ab(N))
allocate(pott(N))
allocate(arms(N))


open(11,status='old',file=trim(filename)//".eng")

do i=1,N
read(11,*) ab(i),pott(i),arms(i)
end do

do i=1,N
 do j=i,N
   if ( ab(i) - ab(j) > 0) then
   switch1=ab(i)
   switch2=arms(i)
   ab(i)=ab(j)
   arms(i)=arms(j)
   ab(j)=switch1
   arms(j)=switch2
end if
end do
end do

a=0
b=0
open(12,status='unknown',file='rms_vs_energy.out')
write(12,'(2A25)') '#ab initio E(cm-1)','cummulative RMS(cm-1)'

do i=1,N
a=a+arms(i)**2

if (b>100) then
write(12,'(2F25.10)') ab(i)*219474.63,sqrt(a/real(i))
b = 0
end if
b=b+1
end do




end program main
