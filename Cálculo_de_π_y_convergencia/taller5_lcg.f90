! Taller 5 -- pi por acierto y error con el LCG (57, 1, 256) en lugar de random_number.
! Salida: lcg.dat (N, pi_est).
program taller5_lcg
  implicit none
  integer, parameter :: dp = kind(1.0d0)
  integer, parameter :: a_lcg = 57, c_lcg = 1, M_lcg = 256
  integer :: k, r, u_dat
  integer(8) :: i, n_total, n_aciertos
  real(dp) :: x, y, pi_est

  open(newunit=u_dat,file='lcg.dat',status='replace')
  do k = 2, 7                                   ! N = 10^2, ..., 10^7
     n_total = 10_8**k
     r = 10                                     ! semilla del LCG, la misma para cada N
     pi_est = -1.0_dp                           ! valor provisional

     ! TODO C1: obtenga pi_est con n_total puntos cuyas coordenadas salen del LCG (sección 5).
     n_aciertos = 0
     do i = 1, n_total
        ! --- Se avanza el estado para la coordenada X ---
        r = mod(a_lcg * r + c_lcg, M_lcg)             
        ! Normalizamos de [0, M) a [0, 1) y luego desplazamos a [-1, 1)
        x = 2.0_dp * real(r, dp) / real(M_lcg, dp) - 1.0_dp
        
        ! --- Se avanza el estado de nuevo para la coordenada Y ---
        r = mod(a_lcg * r + c_lcg, M_lcg)             
        y = 2.0_dp * real(r, dp) / real(M_lcg, dp) - 1.0_dp
        
        ! Evaluamos si cae dentro del círculo de radio 1
        if (x*x + y*y < 1.0_dp) n_aciertos = n_aciertos + 1
     end do
     ! Estimación final de pi para este valor de N
     pi_est = 4.0_dp * real(n_aciertos, dp) / real(n_total, dp)

     write(u_dat,'(i10,es25.16e3)') n_total,pi_est
     print '(i10,f12.6)', n_total,pi_est
  end do
  close(u_dat)
end program taller5_lcg
