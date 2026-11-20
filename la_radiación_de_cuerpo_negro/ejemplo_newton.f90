! ============================================================
! Programa:  ejemplo_newton
! Propósito: hallar la raíz de f partiendo de x0 = 1.5 con el
!            método de Newton-Raphson.
! Salida:    tabla en pantalla; errores_newton.dat
! ============================================================
program ejemplo_newton
  use modulo_ejemplo, only: dp, f, df, raiz_exacta
  implicit none
  real(dp), parameter :: tol = 1.0e-10_dp
  integer, parameter :: max_iter = 50
  real(dp) :: x, paso, dfx
  integer :: n, u

  x = 1.5_dp                                ! estimación inicial x0
  open(newunit=u, file='errores_newton.dat', status='replace', action='write')
  write(u, '(i4,es14.6)') 0, abs(x - raiz_exacta)
  print '(a)', ' n       x_n              f(x_n)     x_n - x_(n-1)   alfa - x_(n-1)'
  print '(i2,f17.12,es13.3)', 0, x, f(x)
  do n = 1, max_iter
     dfx = df(x)
     if (abs(dfx) < tiny(1.0_dp)) stop 'Newton: derivada nula'
     paso = -f(x)/dfx                       ! corrección de Newton
     x = x + paso
     print '(i2,f17.12,es13.3,2es15.3)', n, x, f(x), paso, raiz_exacta - (x - paso)
     write(u, '(i4,es14.6)') n, abs(x - raiz_exacta)
     if (abs(paso) < tol) exit              ! dos iteraciones ya coinciden
  end do
  close(u)
  if (n > max_iter) print '(a)', 'Newton: no convergió en max_iter iteraciones'
  print '(a,es10.2)', 'Residuo final |f(x)|: ', abs(f(x))    ! revisión del resultado
end program ejemplo_newton
