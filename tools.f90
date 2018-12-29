module tools
  use global
  implicit none
contains
    ! --------------------------------------------------------------------------------------
  !
  !              simpson()
  !
  ! --------------------------------------------------------------------------------------
  function  simpson(m,f)
    implicit none
    type (t_mesh)::m
    double precision::simpson
    double precision :: f(:)
    integer::i
    simpson=0.0
    do i=1,m%N-2,2
       simpson=simpson+f(i)*f(i)+4*f(i+1)*f(i+1)+f(i+2)*f(i+2)
    end do
    simpson=m%dv*simpson/3.0
  end function simpson
  ! --------------------------------------------------------------------------------------
  !
  !              trapz()
  !
  ! --------------------------------------------------------------------------------------
  function  trapz(m,f)
    implicit none
    type (t_mesh)::m
    double precision::trapz
    double precision :: f(:)
    integer::i
     trapz=0.0
    do i=1,m%N-1
        trapz=trapz+f(i)*f(i)+f(i+1)*f(i+1)
    end do
    trapz=0.5*m%dv*trapz
  end function trapz
  ! --------------------------------------------------------------------------------------
  !
  !              norm()
  !
  ! --------------------------------------------------------------------------------------
  subroutine  norm(m,evec) 
    implicit none
    double precision :: evec(:),normloc
    double precision, external :: ddot
    type(t_mesh)::m
    if(m%dim.eq.3) then
       normloc=1.0/sqrt(m%dv*ddot(m%N,evec(:),1,evec(:),1))
    else if(m%dim.eq.2) then
       normloc=1.0/sqrt(m%dv*ddot(m%N,evec(:),1,evec(:),1))
    else    if(m%dim.eq.1) then
       !       normloc=1.0/sqrt(trapz(m,evec))
       !       print *,normloc,sqrt(trapz(m,evec))
       normloc=1.0/sqrt(simpson(m,evec))
    else
       print *,' STOP in norm(): dimension=',m%dim,' not yet implemented!'
       stop
    end if
    call dscal(m%N,normloc,evec(:),1)
  end subroutine norm

end module tools
