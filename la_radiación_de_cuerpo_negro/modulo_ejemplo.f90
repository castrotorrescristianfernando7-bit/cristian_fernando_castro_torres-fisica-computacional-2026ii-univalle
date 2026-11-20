! ============================================================
! Módulo:    modulo_ejemplo
! Propósito: función de ejemplo f(x) = x**6 - x - 1 y su derivada.
! ============================================================
module modulo_ejemplo
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: raiz_exacta = 1.134724138401519_dp   ! solo para medir errores
contains

  ! f: función cuyo cero se busca
  real(dp) function f(x)
    real(dp), intent(in) :: x
    f = x**6 - x - 1.0_dp
  end function f

  ! df: derivada analítica de f, f'(x) = 6 x**5 - 1
  real(dp) function df(x)
    real(dp), intent(in) :: x
    df = 6.0_dp*x**5 - 1.0_dp
  end function df

end module modulo_ejemplo
