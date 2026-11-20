! ============================================================
! Programa:  programa_sigma
! Propósito: constante de Stefan-Boltzmann con Romberg, con y sin
!            cambio de variable x = t/(1 - t).
! Salida:    resultado_sigma.dat, tabla en pantalla
! ============================================================
program programa_sigma
  use modulo_romberg
  implicit none
  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  real(dp), parameter :: h_planck = 6.62607015e-34_dp
  real(dp), parameter :: c_luz    = 299792458.0_dp
  real(dp), parameter :: k_boltz  = 1.380649e-23_dp
  real(dp), parameter :: tol = 1.0e-8_dp
  real(dp), parameter :: i_exacta = pi**4/15.0_dp
  real(dp) :: factor, integral, sigma
  integer  :: u, fila, k_final, n_eval

  factor = 2.0_dp*pi*k_boltz**4/(h_planck**3*c_luz**2)
  open(newunit=u, file='resultado_sigma.dat', status='replace', action='write')
  print '(a)', ' caso                      k  evaluaciones        integral          sigma'
  do fila = 1, 2
     if (fila == 1) then
        call romberg(1, 0.0_dp, 1.0_dp, tol, integral, k_final, n_eval)     ! con cambio de variable
     else
        call romberg(2, 0.0_dp, 100.0_dp, tol, integral, k_final, n_eval)   ! sin cambio de variable
     end if
     sigma = factor*integral
     write(u, '(2i6,2es24.16)') k_final, n_eval, integral, sigma
     print '(a,i4,i8,es20.12,es16.8)', merge(' con cambio de variable ', ' sin cambio de variable ', fila == 1), &
           k_final, n_eval, integral, sigma
  end do
  close(u)
  print '(/,a,es20.12)', 'pi**4/15 = ', i_exacta
  print '(a,es16.8)',    'sigma exacta = ', factor*i_exacta
end program programa_sigma
