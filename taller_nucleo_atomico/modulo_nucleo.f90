! Taller 6 -- Módulo de la fórmula semiempírica
module modulo_nucleo
  implicit none
  
  ! Constantes del modelo
  integer, parameter :: dp = kind(1.0d0)
  real(dp), parameter :: a_V = 15.8_dp
  real(dp), parameter :: a_S = 18.3_dp
  real(dp), parameter :: a_C = 0.71_dp
  real(dp), parameter :: a_A = 92.7_dp
  real(dp), parameter :: M_n = 939.57_dp
  real(dp), parameter :: M_H = 938.78_dp

contains

  ! Término de Volumen
  real(dp) function e_vol(Z, A)
    integer, intent(in) :: Z, A
    e_vol = -a_V
  end function e_vol

  ! Término de Superficie
  real(dp) function e_sup(Z, A)
    integer, intent(in) :: Z, A
    e_sup = a_S * (real(A, dp)**(-1.0_dp/3.0_dp))
  end function e_sup

  ! Término de Coulomb
  real(dp) function e_coul(Z, A)
    integer, intent(in) :: Z, A
    e_coul = a_C * (real(Z, dp)**2) * (real(A, dp)**(-4.0_dp/3.0_dp))
  end function e_coul

  ! Término de Asimetría
  real(dp) function e_asim(Z, A)
    integer, intent(in) :: Z, A
    e_asim = a_A * (real(Z, dp) - real(A, dp)/2.0_dp)**2 / (real(A, dp)**2)
  end function e_asim

  ! Energía de enlace por nucleón total (E_B)
  real(dp) function E_B(Z, A)
    integer, intent(in) :: Z, A
    E_B = e_vol(Z, A) + e_sup(Z, A) + e_coul(Z, A) + e_asim(Z, A)
  end function E_B

  ! Masa atómica aproximada (M)
  real(dp) function masa(Z, A)
    integer, intent(in) :: Z, A
    masa = real(Z, dp)*M_H + real(A - Z, dp)*M_n + real(A, dp)*E_B(Z, A)
  end function masa

  ! Búsqueda del isótopo más estable para un A dado
  subroutine Z_estable(A, Z_min, EB_min)
    integer, intent(in) :: A
    integer, intent(out) :: Z_min
    real(dp), intent(out) :: EB_min
    
    integer :: Z
    real(dp) :: eb_actual

    ! Inicialización usando Z=1 como referencia base
    EB_min = E_B(1, A)
    Z_min = 1
    
    ! Búsqueda del menor E_B (mayor estabilidad)
    do Z = 2, A
       eb_actual = E_B(Z, A)
       if (eb_actual < EB_min) then
          EB_min = eb_actual
          Z_min = Z
       end if
    end do
  end subroutine Z_estable

end module modulo_nucleo
