! ============================================================
! Módulo:    integracion
! Propósito: reglas compuestas del trapecio y de Simpson sobre
!            valores muestreados en puntos equiespaciados.
!            Se compila como librería de Python con F2PY.
! ============================================================
module integracion
  implicit none
contains

  ! trapecio: regla compuesta del trapecio
  ! Entradas: n, número de puntos; y(n), valores de la función; h, paso
  ! Salida:   s, integral aproximada
  subroutine trapecio(n, y, h, s)
    integer, intent(in)  :: n
    real(8), intent(in)  :: y(n), h
    real(8), intent(out) :: s
    if (n < 2) then
       s = 0.0d0
    else
       s = h*(0.5d0*(y(1) + y(n)) + sum(y(2:n-1)))
    end if
  end subroutine trapecio

  ! simpson: regla compuesta de Simpson (n impar, n >= 3)
  ! Entradas: n, número de puntos (impar); y(n), valores; h, paso
  ! Salida:   s, integral aproximada. Si n es par o menor que 3
  !           devuelve NaN (no se usa stop: mataría el intérprete)
  subroutine simpson(n, y, h, s)
    integer, intent(in)  :: n
    real(8), intent(in)  :: y(n), h
    real(8), intent(out) :: s
    real(8) :: cero
    if (n < 3 .or. mod(n, 2) == 0) then
       cero = 0.0d0
       s = cero/cero                              ! NaN: entrada inválida
    else
       s = (h/3.0d0)*(y(1) + y(n) + 4.0d0*sum(y(2:n-1:2)) + 2.0d0*sum(y(3:n-2:2)))
    end if
  end subroutine simpson

end module integracion
