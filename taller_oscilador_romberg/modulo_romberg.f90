! ============================================================
! Módulo:    modulo_romberg
! Propósito: integración de Romberg general con refinamiento
!            eficiente del trapecio e integrandos del Taller 7.
! ============================================================
module modulo_romberg
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: gamma_param = 0.5_dp
  real(dp), parameter :: omega_param = 8.0_dp * atan(1.0_dp) ! 2*pi

contains

  ! Integrando 1: 2*gamma * v(t)^2
  real(dp) function integrando_energia(t)
    real(dp), intent(in) :: t
    real(dp) :: v
    v = -exp(-gamma_param*t) * (gamma_param*cos(omega_param*t) + omega_param*sin(omega_param*t))
    integrando_energia = 2.0_dp * gamma_param * v**2
  end function integrando_energia

  ! Integrando 2: 4 / (1 + x^2)
  real(dp) function integrando_pi(x)
    real(dp), intent(in) :: x
    integrando_pi = 4.0_dp / (1.0_dp + x**2)
  end function integrando_pi

  ! Evaluación del integrando seleccionado por modo
  real(dp) function evaluar_f(x, modo)
    real(dp), intent(in) :: x
    integer, intent(in) :: modo
    if (modo == 1) then
       evaluar_f = integrando_energia(x)
    else if (modo == 2) then
       evaluar_f = integrando_pi(x)
    else
       evaluar_f = 0.0_dp
    end if
  end function evaluar_f

  ! Subrutina principal de Romberg
  subroutine romberg(a, b, tol, kmax, modo, S, k_final, n_eval)
    real(dp), intent(in)  :: a, b, tol
    integer, intent(in)   :: kmax, modo
    real(dp), intent(out) :: S(0:, 0:)
    integer, intent(out)  :: k_final, n_eval

    real(dp) :: h, suma, x
    integer :: k, j, i, n_nuevos

    h = b - a
    S(0,0) = 0.5_dp * h * (evaluar_f(a, modo) + evaluar_f(b, modo))
    n_eval = 2
    k_final = 0

    do k = 1, kmax
       h = h / 2.0_dp
       n_nuevos = 2**(k - 1)
       suma = 0.0_dp
       do i = 1, n_nuevos
          x = a + real(2*i - 1, dp) * h
          suma = suma + evaluar_f(x, modo)
       end do
       n_eval = n_eval + n_nuevos

       S(k,0) = 0.5_dp * S(k-1,0) + h * suma

       do j = 1, k
          S(k,j) = S(k,j-1) + (S(k,j-1) - S(k-1,j-1)) / (4.0_dp**j - 1.0_dp)
       end do

       k_final = k
       if (abs(S(k,k) - S(k-1,k-1)) < tol) exit
    end do
  end subroutine romberg

end module modulo_romberg
