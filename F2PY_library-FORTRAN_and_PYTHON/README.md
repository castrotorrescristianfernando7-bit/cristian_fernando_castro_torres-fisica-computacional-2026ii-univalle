# Taller 10 — Librería de integración en Fortran usada desde Python (F2PY)

**Curso:** Física Computacional (106018C) — Universidad del Valle  
**Estudiante:** Cristian Fernando Castro Torres (Código: 202518873)  
**Fecha de entrega:** 27 de noviembre de 2026  

## Descripción
Implementación de un módulo de integración numérica en **Fortran 2008** (`integracion`) compilado como librería compartida (`integ`) mediante **F2PY** para su uso desde **Python**. El módulo provee las reglas compuestas del trapecio ($\mathcal{O}(h^2)$) y de Simpson ($\mathcal{O}(h^4)$), las cuales se emplean para verificar los órdenes de convergencia y calcular la constante de Stefan-Boltzmann ($\sigma$) a partir de la integral de radiación de cuerpo negro de Planck.

## Contenido del repositorio (`taller_f2py/`)
* `modulo_integracion.f90`: Código fuente en Fortran con las subrutinas `trapecio` y `simpson` en doble precisión (`real(8)`), sin sentencias `stop`.
* `F2PY_library_students-5.ipynb`: Cuaderno de Jupyter con la compilación vía F2PY, la ejecución de las actividades 1 a 4, la celda de verificación en **3 de 3 PASS** y las respuestas de interpretación de la sección D.
* `README.md`: Documentación general del taller y guía de compilación.

## Compilación y ejecución
Desde la carpeta `taller_f2py/`, compilar el módulo de Fortran para generar la extensión de Python:

```bash
python -m numpy.f2py -c modulo_integracion.f90 -m integ
