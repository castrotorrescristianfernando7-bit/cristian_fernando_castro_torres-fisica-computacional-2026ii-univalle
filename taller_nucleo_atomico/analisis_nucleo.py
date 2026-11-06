"""Taller 6 -- Módulo de análisis en Python."""
import numpy as np
import matplotlib.pyplot as plt

def leer_valle(archivo):
    """
    Lee la tabla valle.dat generada por Fortran.
    """
    datos = np.loadtxt(archivo)
    return datos[:, 0], datos[:, 1], datos[:, 2]

def encontrar_maximo(A, Z, E_B):
    """
    Encuentra el núcleo más ligado (menor valor de E_B) usando un ciclo explícito.
    """
    EB_min = 1.0e10
    A_max = 0
    Z_max = 0
    
    for i in range(len(E_B)):
        if E_B[i] < EB_min:
            EB_min = E_B[i]
            A_max = A[i]
            Z_max = Z[i]
            
    return A_max, Z_max, EB_min

def escribir_maximo(archivo, A_max, Z_max, EB_min):
    """
    Escribe el núcleo más ligado en un archivo.
    """
    with open(archivo, 'w') as f:
        f.write(f"{int(A_max):5d}{int(Z_max):5d}{EB_min:15.8f}\n")

def graficar_curva(A, E_B, experimentales, archivo_salida):
    """
    Genera y guarda la gráfica de B/A = -E_B contra A, con datos experimentales.
    """
    B_A_modelo = -E_B
    
    plt.figure(figsize=(9, 5))
    plt.plot(A, B_A_modelo, '-', color="#162A54", label="Modelo (Curva de máximo enlace)")
    
    A_exp = [dato[2] for dato in experimentales]
    B_A_exp = [dato[3] for dato in experimentales]
    
    plt.plot(A_exp, B_A_exp, 'o', color="#C9972C", label="Datos experimentales (AME2020)")
    
    plt.xlabel("Número de nucleones (A)")
    plt.ylabel("Energía de enlace por nucleón B/A [MeV]")
    plt.title("Energía de Enlace por Nucleón vs A (Modelo de Gota Líquida)")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(archivo_salida)
    plt.close()
