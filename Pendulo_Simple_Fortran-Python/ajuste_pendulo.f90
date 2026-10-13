program ajuste_pendulo
  use iso_fortran_env, only: real64, iostat_end ! Para precisión y detectar el fin del archivo
  implicit none                                 ! Obliga a declarar todo (buena práctica)

  ! Constante pi (usando arcoseno para mayor exactitud)
  real(real64), parameter :: PI = acos(-1.0_real64)

  ! Variables para manejar los archivos y el estado de lectura
  integer :: arch_in, unidad_salida, error_lectura
  integer :: id_fila, oscilaciones, n

  ! Variables temporales para leer cada línea
  real(real64) :: l_cm, ang_deg, t_total
  real(real64) :: val_x, val_y, periodo

  ! Sumatorias para hacer el ajuste por mínimos cuadrados
  real(real64) :: sum_x, sum_y, sum_xx, sum_xy
  real(real64) :: denominador, a, b, g
  real(real64) :: y_estimada, suma_residuos, suma_total, r2

  ! Inicializamos contadores y sumatorias en cero para empezar limpios
  n = 0
  sum_x  = 0.0_real64
  sum_y  = 0.0_real64
  sum_xx = 0.0_real64
  sum_xy = 0.0_real64

  ! 1. Abrimos el archivo de datos. Si falla, detenemos todo de una.
  open(newunit=arch_in, file='pendulo_limpio.dat', status='old', &
       action="read", iostat=error_lectura)

  if (error_lectura /= 0) error stop '¡Pilas! No pude abrir el archivo pendulo_limpio.dat'

  ! 2 y 3. Leemos el archivo línea por línea con un ciclo infinito hasta que se acabe
  do
      ! Leemos las 5 columnas
      read(arch_in, *, iostat=error_lectura) id_fila, l_cm, ang_deg, oscilaciones, t_total

      if (error_lectura == iostat_end) exit  ! Terminamos de leer todo bien
      if (error_lectura /= 0) error stop '¡Uy! Hay un problema con el formato en alguna fila.'

      ! Hacemos las conversiones que pide el profe
      val_x = l_cm / 100.0_real64                    ! Pasamos la longitud a metros
      periodo = t_total / real(oscilaciones, real64) ! Calculamos el periodo T
      val_y = periodo**2                             ! Elevamos T al cuadrado

      ! 4 y 5. Vamos sumando los valores de esta fila a nuestros acumuladores
      n = n + 1
      sum_x  = sum_x  + val_x
      sum_y  = sum_y  + val_y
      sum_xx = sum_xx + val_x*val_x
      sum_xy = sum_xy + val_x*val_y
  end do

  ! 11. Validamos que haya suficientes datos (mínimo 2 para una recta)
  if (n < 2) error stop 'Necesito por lo menos 2 mediciones para hacer el ajuste lineal.'

  ! 6. Calculamos la parte de abajo de la fórmula (denominador)
  denominador = real(n, real64)*sum_xx - sum_x*sum_x

  ! Revisamos que no sea cero para que no explote la división
  if (abs(denominador) <= tiny(denominador)) error stop 'El denominador dio cero, no se puede calcular a.'

  ! Calculamos pendiente e intercepto
  a = (real(n, real64)*sum_xy - sum_x*sum_y) / denominador
  b = (sum_y - a*sum_x) / real(n, real64)

  ! 7. Calculamos la gravedad despejando de la fórmula del péndulo
  g = (4.0_real64 * PI**2) / a

  ! 8. Toca volver a leer para sacar el coeficiente de determinación R^2
  suma_residuos = 0.0_real64
  suma_total = 0.0_real64

  ! Nos devolvemos a la primera línea del archivo
  rewind(arch_in)

  do
      read(arch_in, *, iostat=error_lectura) id_fila, l_cm, ang_deg, oscilaciones, t_total
      if (error_lectura == iostat_end) exit

      ! Repetimos el cálculo de x y y para esta fila
      val_x = l_cm / 100.0_real64
      periodo = t_total / real(oscilaciones, real64)
      val_y = periodo**2

      ! Predicción del modelo y diferencias al cuadrado
      y_estimada = a*val_x + b
      suma_residuos = suma_residuos + (val_y - y_estimada)**2
      suma_total = suma_total + (val_y - sum_y / real(n, real64))**2
  end do

  ! Ya no necesitamos más el archivo, lo cerramos
  close(arch_in)

  ! Calculamos R^2 protegiendo por si la variación total es nula
  if (suma_total > 0.0_real64) then
      r2 = 1.0_real64 - (suma_residuos / suma_total)
  else
      r2 = 0.0_real64
  end if

  ! 9. Imprimimos el resumen de resultados en pantalla
  print *, '      RESULTADOS DEL EXPERIMENTO DEL PÉNDULO       '
  print *, ''
  print '(A, I5)',      ' Cantidad de datos (N)  = ', n
  print '(A, F12.6)', ' Pendiente (a) [s^2/m]  = ', a
  print '(A, F12.6)', ' Intercepto (b) [s^2]   = ', b
  print '(A, F12.6)', ' Confianza (R^2)        = ', r2
  print '(A, F12.6)', ' Gravedad (g) [m/s^2]   = ', g

  ! 10. Código obligatorio del profe para el archivo PASS/FAIL
  open(newunit=unidad_salida, file='resultados_ajuste.dat', &
       status='replace', action='write')
  write(unidad_salida, *) n, a, b, r2, g
  close(unidad_salida)

end program ajuste_pendulo
