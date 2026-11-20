! ============================================================
! Módulo:    modulo_ceros
! Propósito: F(x) = x - n(1 - e^-x), su derivada y tres métodos
!            para hallar ceros: bisección, Newton y secante.
! Estado:    0 = convergió; 1 = no convergió en max_iter;
!            2 = sin cambio de signo (bisección); 3 = denominador nulo.
! ============================================================
module modulo_ceros
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: max_iter = 100      ! protección contra ciclos sin fin
contains

  ! F: ecuación de Wien, x = n(1 - e^-x), escrita como F(x) = 0
  real(dp) function F(x, n)
    real(dp), intent(in) :: x
    integer,  intent(in) :: n
    F = x - real(n, dp)*(1.0_dp - exp(-x))
  end function F

  ! dF: derivada analítica, F'(x) = 1 - n e^-x
  real(dp) function dF(x, n)
    real(dp), intent(in) :: x
    integer,  intent(in) :: n
    dF = 1.0_dp - real(n, dp)*exp(-x)
  end function dF

  ! Bisección en [a, b]; hist(k) guarda el punto medio c_k
  subroutine biseccion(n, a0, b0, tol, raiz, iter, estado, hist)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: a0, b0, tol
    real(dp), intent(out) :: raiz, hist(max_iter)
    integer,  intent(out) :: iter, estado
    real(dp) :: a, b, c
    a = a0;  b = b0;  raiz = a0;  iter = 0;  estado = 1;  hist = 0.0_dp
    if (F(a, n)*F(b, n) >= 0.0_dp) then      ! paso 1: cambio de signo
       estado = 2
       print '(a)', 'bisección: F(a) y F(b) no tienen signos distintos'
       return
    end if
    do iter = 1, max_iter
       c = 0.5_dp*(a + b)                    ! paso 2: punto medio
       hist(iter) = c
       if (b - c <= tol) then                ! paso 3: criterio de parada
          raiz = c;  estado = 0;  return
       end if
       if (F(a, n)*F(c, n) <= 0.0_dp) then   ! paso 4: mitad con cambio de signo
          b = c
       else
          a = c
       end if
    end do
    iter = max_iter
    raiz = c
    print '(a)', 'bisección: no convergió en max_iter iteraciones'
  end subroutine biseccion

  ! Newton-Raphson desde x0; hist(k) guarda x_k
  subroutine newton(n, x0, tol, raiz, iter, estado, hist)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: x0, tol
    real(dp), intent(out) :: raiz, hist(max_iter)
    integer,  intent(out) :: iter, estado
    real(dp) :: x, paso, dfx
    x = x0;  raiz = x0;  iter = 0;  estado = 1;  hist = 0.0_dp
    do iter = 1, max_iter
       dfx = dF(x, n)
       if (abs(dfx) < tiny(1.0_dp)) then
          estado = 3
          print '(a)', 'Newton: derivada nula'
          return
       end if
       paso = -F(x, n)/dfx                   ! corrección de Newton
       x = x + paso
       hist(iter) = x
       if (abs(paso) < tol) then
          raiz = x;  estado = 0;  return
       end if
    end do
    iter = max_iter
    raiz = x
    print '(a)', 'Newton: no convergió en max_iter iteraciones'
  end subroutine newton

  ! Secante con dos valores iniciales; hist(k) guarda x_(k+1)
  subroutine secante(n, x0, x1, tol, raiz, iter, estado, hist)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: x0, x1, tol
    real(dp), intent(out) :: raiz, hist(max_iter)
    integer,  intent(out) :: iter, estado
    real(dp) :: xa, xb, xn, den
    xa = x0;  xb = x1;  raiz = x1;  iter = 0;  estado = 1;  hist = 0.0_dp
    do iter = 1, max_iter
       den = F(xb, n) - F(xa, n)
       if (abs(den) < tiny(1.0_dp)) then
          estado = 3
          print '(a)', 'secante: denominador nulo'
          return
       end if
       xn = xb - F(xb, n)*(xb - xa)/den
       hist(iter) = xn
       xa = xb;  xb = xn
       if (abs(xb - xa) < tol) then
          raiz = xb;  estado = 0;  return
       end if
    end do
    iter = max_iter
    raiz = xb
    print '(a)', 'secante: no convergió en max_iter iteraciones'
  end subroutine secante

end module modulo_ceros
