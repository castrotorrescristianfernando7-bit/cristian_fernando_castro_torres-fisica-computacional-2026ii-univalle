! ============================================================
! Programa:  ejemplo_newton_falla
! Propósito: efecto del valor inicial en el método de Newton
!            para f(x) = tan(pi x) - x - 6.
! Salida:    tabla en pantalla
! ============================================================
program ejemplo_newton_falla
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: pi = 4*atan(1.0_dp)
  real(dp) :: x(3)                          ! una estimación por columna
  integer :: n

  x = [0.48_dp, 0.40_dp, 0.01_dp]           ! tres valores iniciales
  print '(a)', ' n    x0 = 0.48     x0 = 0.40     x0 = 0.01'
  do n = 1, 11
     x = x - f(x)/df(x)                     ! una iteración de las tres columnas
     print '(i2,3f14.7)', n, x
  end do

contains

  elemental real(dp) function f(x)
    real(dp), intent(in) :: x
    f = tan(pi*x) - x - 6.0_dp
  end function f

  ! derivada analítica: pi/cos(pi x)**2 - 1
  elemental real(dp) function df(x)
    real(dp), intent(in) :: x
    df = pi/cos(pi*x)**2 - 1.0_dp
  end function df

end program ejemplo_newton_falla
