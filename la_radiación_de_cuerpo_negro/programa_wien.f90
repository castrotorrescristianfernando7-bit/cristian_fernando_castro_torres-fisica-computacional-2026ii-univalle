! ============================================================
! Programa:  programa_wien
! Propósito: raíces x5 y x3 de x = n(1 - e^-x) con tres métodos
!            y constante de Wien b = hc/(kB x5).
! Salida:    resultado_wien.dat, iteraciones_wien.dat, constante_wien.dat
! ============================================================
program programa_wien
  use modulo_ceros
  implicit none
  real(dp), parameter :: tol = 1.0e-10_dp
  real(dp), parameter :: h_planck = 6.62607015e-34_dp
  real(dp), parameter :: c_luz    = 299792458.0_dp
  real(dp), parameter :: k_boltz  = 1.380649e-23_dp
  real(dp) :: raiz, hist(max_iter), x5, x3, b
  integer  :: u_res, u_it, u_cte, iter, estado, k, m, i, n

  open(newunit=u_res, file='resultado_wien.dat',   status='replace', action='write')
  open(newunit=u_it,  file='iteraciones_wien.dat', status='replace', action='write')

  do i = 1, 2
     n = merge(5, 3, i == 1)                 ! primero n = 5 y luego n = 3
     do m = 1, 3                             ! 1 bisección, 2 Newton, 3 secante
        select case (m)
        case (1); call biseccion(n, 1.0_dp, 10.0_dp, tol, raiz, iter, estado, hist)
        case (2); call newton(n, 10.0_dp, tol, raiz, iter, estado, hist)
        case (3); call secante(n, 10.0_dp, 9.0_dp, tol, raiz, iter, estado, hist)
        end select
        if (estado /= 0) stop 'programa_wien: un método no convergió'
        write(u_res, '(2i6,es24.16,i6)') n, m, raiz, iter
        if (n == 5) then
           if (m == 2) x5 = raiz             ! x5 por Newton
           do k = 1, iter
              write(u_it, '(2i6,es24.16)') m, k, hist(k)
           end do
        else if (m == 2) then
           x3 = raiz
        end if
        print '(a,i1,a,i1,a,f18.12,a,i0)', 'n=', n, ' método=', m, '  raíz=', raiz, '  iteraciones=', iter
     end do
  end do

  ! Newton con n = 5 desde x0 = 0.5: converge a la raíz trivial x = 0
  call newton(5, 0.5_dp, tol, raiz, iter, estado, hist)
  if (estado /= 0) stop 'programa_wien: Newton desde 0.5 no convergió'
  write(u_res, '(2i6,es24.16,i6)') 5, 4, raiz, iter
  print '(a,es12.3,a,i0)', 'Newton desde x0=0.5: raíz=', raiz, '  iteraciones=', iter
  close(u_res);  close(u_it)

  b = h_planck*c_luz/(k_boltz*x5)            ! constante de Wien [m K]
  open(newunit=u_cte, file='constante_wien.dat', status='replace', action='write')
  write(u_cte, '(3es24.16)') x5, x3, b
  close(u_cte)
  print '(/,a,f14.10,a,f14.10)', 'x5 = ', x5, '   x3 = ', x3
  print '(a,es16.9,a)', 'b  = ', b, ' m K'
end program programa_wien
