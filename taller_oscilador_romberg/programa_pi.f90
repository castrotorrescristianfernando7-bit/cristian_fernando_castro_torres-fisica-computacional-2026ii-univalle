! ============================================================
! Programa:  programa_pi
! Propósito: calcula pi integrando 4/(1+x^2) en [0,1] con el mismo módulo
! Salida:    resultado_pi.dat
! ============================================================
program programa_pi
  use modulo_romberg, only: dp, romberg
  implicit none

  real(dp), parameter :: pi_exacto = 4.0_dp * atan(1.0_dp)
  real(dp), allocatable :: S(:,:)
  real(dp) :: a, b, tol, pi_calc, err_abs
  integer :: kmax, k_final, n_eval

  a = 0.0_dp
  b = 1.0_dp
  tol = 1.0e-5_dp
  kmax = 12

  allocate(S(0:kmax, 0:kmax))
  S = 0.0_dp

  ! Llamar a romberg en modo 2 (integrando pi)
  call romberg(a, b, tol, kmax, 2, S, k_final, n_eval)
  pi_calc = S(k_final, k_final)
  err_abs = abs(pi_calc - pi_exacto)

  ! Guardar resultado_pi.dat
  open(unit=30, file='resultado_pi.dat', status='replace', action='write')
  write(30, '(i5, i7, 2es24.15)') k_final, n_eval, pi_calc, err_abs
  close(30)

  deallocate(S)
end program programa_pi
