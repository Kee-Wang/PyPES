module pes_shell
  use pes_shell_short
  use pes_shell_long

contains
  subroutine pes_init()
    call pes_init_short()
    call pes_init_long()  !a0=7.0
  return
  end subroutine pes_init

  function f(x)
    integer,parameter::wp=selected_real_kind(12,300)
    real(kind=wp),dimension(3,6)::x
    real(kind=wp)::f,a,b,rCO,s,p,auang,aucm
    auang = 0.5292
    aucm = 219474.63
    
!    f = f_short(x)
    a = 5.0/auang
    b = 6.0/auang
   !write(*,*) 'a=',a
    rCO = sqrt(abs((x(1,3)-x(1,6))**2+(x(2,3)-x(2,6))**2+(x(3,3)-x(3,6))**2)) 
    if (rCO <= a) then
      f = f_short(x)
        write(*,*) 'dis>a, f_short='!,f*aucm
    else if (rCO < b) then
      p = (rCO-a)/(b-a)
      s = 10*p**3 - 15*p**4 + 6*p**5
      f= (1-s)*f_short(x) + s*f_long(x)
        write(*,*) 'a<dis<b, s,l,swith'!,f_short(x)*aucm,f_long(x)*aucm,f*aucm
    else
      f = f_long(x)
 
        write(*,*) 'dis>b, f_long='!,f*aucm
    end if
    return

  end function f

end module
