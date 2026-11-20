! ============================================================
! Programa:  ejercicio_barrido   (Ejercicio de calentamiento 2)
! Propósito: localizar y refinar todas las raíces reales de f.
! Salida:    en pantalla
! ============================================================
program ejercicio_barrido
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: n = 81              ! puntos de la malla
  real(dp), parameter :: tol = 1.0e-3_dp
  real(dp) :: x(n), y(n), a, b, c
  integer :: i

  integer :: it

  ! --- Paso 1: malla de 81 puntos en [-4, 4] (paso 0.1) y valores de f
  x = [(-4.0_dp + 8.0_dp*real(i - 1, dp)/real(n - 1, dp), i = 1, n)]
  y = f(x)                                   ! f es elemental: admite arreglos

  ! --- Paso 2: cada par vecino con signos distintos es un intervalo de arranque
  print '(a)', 'intervalo de arranque              raíz (bisección)'
  do i = 1, n - 1
     if (y(i)*y(i+1) < 0.0_dp) then
        a = x(i);  b = x(i+1)
        do it = 1, 100                       ! máximo de iteraciones
           c = 0.5_dp*(a + b)
           if (b - c <= tol) exit
           if (f(a)*f(c) <= 0.0_dp) then
              b = c
           else
              a = c
           end if
        end do
        print '(a,f6.2,a,f6.2,a,f10.4)', '[', x(i), ',', x(i+1), ']        ', c
     end if
  end do

contains

  ! f: polinomio de grado 3 (elemental: admite un arreglo como argumento)
  elemental real(dp) function f(x)
    real(dp), intent(in) :: x
    f = x**3 + 2.0_dp*x**2 - 3.0_dp*x - 1.0_dp
  end function f

end program ejercicio_barrido
