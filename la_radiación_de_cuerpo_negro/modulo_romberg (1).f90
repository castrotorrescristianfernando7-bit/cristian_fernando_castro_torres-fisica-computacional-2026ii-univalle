! ============================================================
! Módulo:    modulo_romberg
! Propósito: integración de Romberg con el criterio de la diagonal
!            |R(k,k) - R(k-1,k-1)| < tol, e integrandos de Planck.
! Nivel k:   k = 0 es el trapecio inicial; el nivel k usa 2**k + 1 nodos.
! Integrandos (sel):
!   1: g(t) = x**3/(e**x - 1) / (1 - t)**2,  x = t/(1 - t)   [0,1]
!   2: x**3/(e**x - 1)                                        [0,inf)
! ============================================================
module modulo_romberg
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: k_max = 20             ! límite de niveles
  real(dp), parameter :: x_corte = 700.0_dp    ! e**x desborda cerca de 709
contains

  ! planck: x**3/(e**x - 1) con límite 0 en x = 0 y cero para x grande
  real(dp) function planck(x)
    real(dp), intent(in) :: x
    if (x <= 0.0_dp .or. x > x_corte) then
       planck = 0.0_dp
    else
       planck = x**3/(exp(x) - 1.0_dp)
    end if
  end function planck

  ! integrando: selecciona la función a integrar
  real(dp) function integrando(sel, s)
    integer,  intent(in) :: sel
    real(dp), intent(in) :: s
    real(dp) :: x
    select case (sel)
    case (1)                                  ! s es t en [0,1]
       if (s <= 0.0_dp .or. s >= 1.0_dp) then
          integrando = 0.0_dp                 ! límites en t = 0 y t = 1
       else
          x = s/(1.0_dp - s)
          integrando = planck(x)/(1.0_dp - s)**2
       end if
    case (2)                                  ! s es x
       integrando = planck(s)
    case default
       stop 'integrando: selector inválido'
    end select
  end function integrando

  ! romberg: integral de integrando(sel, .) en [a, b]
  subroutine romberg(sel, a, b, tol, integral, k_final, n_eval)
    integer,  intent(in)  :: sel
    real(dp), intent(in)  :: a, b, tol
    real(dp), intent(out) :: integral
    integer,  intent(out) :: k_final, n_eval
    real(dp) :: R(0:k_max, 0:k_max), h, suma, factor
    integer  :: k, j, i

    h = b - a
    R(0,0) = 0.5_dp*h*(integrando(sel, a) + integrando(sel, b))   ! trapecio inicial
    integral = R(0,0);  k_final = 0
    do k = 1, k_max
       h = 0.5_dp*h
       suma = 0.0_dp
       do i = 1, 2**(k-1)                      ! nodos nuevos a + (2i-1) h
          suma = suma + integrando(sel, a + real(2*i - 1, dp)*h)
       end do
       R(k,0) = 0.5_dp*R(k-1,0) + h*suma       ! trapecio refinado
       factor = 1.0_dp
       do j = 1, k                             ! extrapolación de Richardson
          factor = 4.0_dp*factor
          R(k,j) = (factor*R(k,j-1) - R(k-1,j-1))/(factor - 1.0_dp)
       end do
       integral = R(k,k);  k_final = k
       if (abs(R(k,k) - R(k-1,k-1)) < tol) exit
    end do
    if (k_final == k_max .and. abs(R(k_max,k_max) - R(k_max-1,k_max-1)) >= tol) &
         print '(a)', 'romberg: no convergió en k_max niveles'
    n_eval = 2**k_final + 1
  end subroutine romberg

end module modulo_romberg
