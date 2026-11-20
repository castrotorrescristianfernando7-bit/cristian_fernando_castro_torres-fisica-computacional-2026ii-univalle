! ============================================================
! Programa:  ejercicio_secante   (Ejercicio de calentamiento 1)
! Propósito: raíz de f con el método de la secante.
! Salida:    tabla en pantalla; errores_secante.dat
! ============================================================
program ejercicio_secante
  use modulo_ejemplo, only: dp, f, raiz_exacta
  implicit none
  real(dp), parameter :: tol = 1.0e-10_dp
  integer, parameter :: max_iter = 15
  real(dp) :: x0, x1, x_nuevo
  integer :: n, u

  x0 = 1.0_dp
  x1 = 2.0_dp
  open(newunit=u, file='errores_secante.dat', status='replace', action='write')
  print '(a)', ' n      x_n             error'
  do n = 2, max_iter
     if (abs(f(x1) - f(x0)) < tiny(1.0_dp)) stop 'secante: denominador nulo'
     x_nuevo = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0))   ! fórmula de la sección 4.1
     x0 = x1                                           ! se descarta la estimación más vieja
     x1 = x_nuevo
     print '(i2,f17.12,es12.3)', n, x_nuevo, abs(x_nuevo - raiz_exacta)
     write(u, '(i4,es14.6)') n, abs(x_nuevo - raiz_exacta)
     if (abs(x1 - x0) < tol) exit                      ! dos estimaciones consecutivas coinciden
  end do
  close(u)
end program ejercicio_secante
