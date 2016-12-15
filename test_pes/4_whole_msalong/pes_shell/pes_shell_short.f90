MODULE pes_shell_short
use pes, wp=>pes_wp

private
public :: pes_init_short, f_short
save
integer, parameter :: &
  nki(0:2)=(/3,2,1/), nk=6, iord(0:nk-1)=(/1,2,3,4,5,6/)
real (kind=wp) :: x1_cf, y1_cf, z1_cf
type (cx_t) :: &
  x3y2z1_pc = cx_null
real (kind=wp), allocatable ::   &
  x3y2z1_cf(:)

CONTAINS
  !=========================!
  ! pes_init()              !
  !=========================!
  SUBROUTINE pes_init_short()
  integer :: iun, nb
  logical :: b0
  character (len=255) :: chd
  character (len=255) :: dirname

  dirname = '../pes_shell/coef_short/'
  call pes_getiun (iun)
  b0 = dirname(len_trim(dirname):len_trim(dirname)).eq.'/'
  if (b0) then
   chd = dirname
  else
   chd = trim(dirname)//'/'
  endif
  
  write (*,*) 'Principal data directory: ', chd(1:len_trim(chd))
  write (*,*) ' reading pcf-x3y2z1.dat'
  open (iun, status='old', file=trim(chd)//'pcf-x3y2z1.dat')
  read (iun,*) x3y2z1_pc
  read (iun,*) nb
  if (nb.ne.pes_x3y2z1_nb(x3y2z1_pc%dg)) then
   stop 'pes_init_short: x3y2z1 dimension error'
  endif
  allocate (x3y2z1_cf(0:nb-1))
  if (1.le.nb) then
   read (iun,*) x3y2z1_cf
  endif
  close (iun)

  return
  END SUBROUTINE pes_init_short

  !=================================!
  ! function to calculate potential !
  !=================================!
  FUNCTION f_short(xn)
  implicit none
  real (kind=wp), dimension(:,:), intent (in) :: xn
  real (kind=wp) :: f_short
  !:::::::::::::::::::::::
  real (kind=wp) :: r(0:nk-1,0:nk-1)
  
  call pes_dists(xn, r)
  f_short = cx_f321(nki, r, x3y2z1_pc, x3y2z1_cf) 
  
  return
  END FUNCTION f_short

END MODULE pes_shell_short
