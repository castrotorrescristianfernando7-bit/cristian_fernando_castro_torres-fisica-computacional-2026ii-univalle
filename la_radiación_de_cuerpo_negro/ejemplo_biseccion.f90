! ============================================================
! Programa:  ejemplo_biseccion
! Propósito: hallar la raíz de f en [1, 2] por bisección.
! Salida:    tabla en pantalla; errores_biseccion.dat
! ============================================================
program ejemplo_biseccion
  use modulo_ejemplo, only: dp, f, raiz_exacta
  implicit none
  real(dp), parameter :: tol = 1.0e-3_dp
  integer, parameter :: max_iter = 100      ! protección contra ciclos sin fin
  real(dp) :: a, b, c
  integer :: n, u

  a = 1.0_dp
  b = 2.0_dp
  ! --- Paso 1: el intervalo debe contener un cambio de signo
  if (f(a)*f(b) >= 0.0_dp) stop 'bisección: f(a) y f(b) no tienen signos distintos'

  open(newunit=u, file='errores_biseccion.dat', status='replace', action='write')
  print '(a)', '  n        a        b        c       b-c       f(c)'
  do n = 1, max_iter
     c = 0.5_dp*(a + b)                     ! paso 2: punto medio
     print '(i3,3f9.4,f10.5,f11.4)', n, a, b, c, b - c, f(c)
     write(u, '(i4,es14.6)') n, abs(c - raiz_exacta)
     if (b - c <= tol) exit                 ! paso 3: criterio de parada
     if (f(a)*f(c) <= 0.0_dp) then          ! paso 4: mitad con cambio de signo
        b = c
     else
        a = c
     end if
  end do
  close(u)
  print '(/,a,f8.4,a,i0,a)', 'Raíz aproximada: ', c, ' en ', n, ' iteraciones'
  print '(a,es9.2)', 'Error real:      ', abs(c - raiz_exacta)
end program ejemplo_biseccion
