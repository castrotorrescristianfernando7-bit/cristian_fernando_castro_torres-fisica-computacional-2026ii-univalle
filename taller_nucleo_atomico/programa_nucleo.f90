! Taller 6 -- Programa principal
program programa_nucleo
  use modulo_nucleo
  implicit none
  
  integer :: u_nuc, u_val, u_max, i, A, Z
  integer :: Z_valle, A_global, Z_global
  real(dp) :: EB_valle, EB_global
  
  ! Datos para la primera pregunta: 244Pu, 70Ge, 179Hf, 38Ar, 56Fe, 40Ca, 235U, 238U
  integer, dimension(8) :: arr_Z = [94, 32, 72, 18, 26, 20, 92, 92]
  integer, dimension(8) :: arr_A = [244, 70, 179, 38, 56, 40, 235, 238]

  ! 1. Generar nucleos.dat
  open(newunit=u_nuc, file='nucleos.dat', status='replace')
  do i = 1, 8
     Z = arr_Z(i)
     A = arr_A(i)
     ! 8 columnas: Z, A, E_vol, E_sup, E_coul, E_asim, E_B, M
     write(u_nuc, '(2I5, 5F15.8, F15.4)') Z, A, e_vol(Z,A), e_sup(Z,A), e_coul(Z,A), e_asim(Z,A), E_B(Z,A), masa(Z,A)
  end do
  close(u_nuc)

  ! 2 y 3. Generar valle.dat y encontrar el núcleo globalmente más ligado
  open(newunit=u_val, file='valle.dat', status='replace')
  EB_global = 1.0e10_dp ! Valor inicial arbitrariamente alto para asegurar reemplazo
  
  do A = 4, 260
     call Z_estable(A, Z_valle, EB_valle)
     write(u_val, '(2I5, F15.8)') A, Z_valle, EB_valle
     
     ! Filtrar el mínimo global a lo largo de toda la iteración
     if (EB_valle < EB_global) then
        EB_global = EB_valle
        A_global = A
        Z_global = Z_valle
     end if
  end do
  close(u_val)

  ! Escribir maximo.dat con el resultado global
  open(newunit=u_max, file='maximo.dat', status='replace')
  write(u_max, '(2I5, F15.8)') A_global, Z_global, EB_global
  close(u_max)

end program programa_nucleo
