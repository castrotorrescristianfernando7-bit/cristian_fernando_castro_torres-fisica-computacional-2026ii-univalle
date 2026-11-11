! ============================================================
! Programa:  programa_energia
! Propósito: calcula E_dis via Romberg y la compara con E(0) - E(T)
! Salida:    tabla_energia.dat y resultado_energia.dat
! ============================================================
program programa_energia
  use modulo_romberg, only: dp, romberg, gamma_param, omega_param
  implicit none

  real(dp), allocatable :: S(:,:)
  real(dp) :: a, b, tol, e_dis, e0, eT, e_ref, dif_abs, w02
  integer :: kmax, k_final, n_eval, k, j

  a = 0.0_dp
  b = 3.0_dp
  tol = 1.0e-6_dp
  kmax = 12

  allocate(S(0:kmax, 0:kmax))
  S = 0.0_dp

  ! Calcular Romberg (modo 1 = energía)
  call romberg(a, b, tol, kmax, 1, S, k_final, n_eval)
  e_dis = S(k_final, k_final)

  ! Cálculo analítico E(0) - E(T)
  w02 = omega_param**2 + gamma_param**2
  e0 = 0.5_dp * (-gamma_param)**2 + 0.5_dp * w02 * (1.0_dp**2)
  eT = 0.5_dp * (-exp(-gamma_param*b)*(gamma_param*cos(omega_param*b) + omega_param*sin(omega_param*b)))**2 &
       + 0.5_dp * w02 * (exp(-gamma_param*b)*cos(omega_param*b))**2

  e_ref = e0 - eT
  dif_abs = abs(e_dis - e_ref)

  ! Guardar tabla_energia.dat
  open(unit=10, file='tabla_energia.dat', status='replace', action='write')
  do k = 0, k_final
     write(10, '(13es24.15)') (S(k, j), j = 0, k)
  end do
  close(10)

  ! Guardar resultado_energia.dat
  open(unit=20, file='resultado_energia.dat', status='replace', action='write')
  write(20, '(i5, i7, 3es24.15)') k_final, n_eval, e_dis, e_ref, dif_abs
  close(20)

  deallocate(S)
end program programa_energia
