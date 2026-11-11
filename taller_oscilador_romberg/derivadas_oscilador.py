import numpy as np
import matplotlib.pyplot as plt

G, W = 0.5, 2 * np.pi
def x_ex(t): return np.exp(-G * t) * np.cos(W * t)
def v_ex(t): return -np.exp(-G * t) * (G * np.cos(W * t) + W * np.sin(W * t))
def a_ex(t): return np.exp(-G * t) * ((G**2 - W**2) * np.cos(W * t) + 2 * G * W * np.sin(W * t))

def procesar_paso(dt):
    t = np.arange(0, 3 + dt/2, dt)
    x = x_ex(t)
    
    v_adel = (x[1:] - x[:-1]) / dt
    v_atras = (x[1:] - x[:-1]) / dt
    v_cent = (x[2:] - x[:-2]) / (2 * dt)
    v_5p = (-x[4:] + 8*x[3:-1] - 8*x[1:-3] + x[:-4]) / (12 * dt)
    a_cent = (x[2:] - 2*x[1:-1] + x[:-2]) / dt**2
    
    err_adel = np.max(np.abs(v_adel - v_ex(t[:-1])))
    err_atras = np.max(np.abs(v_atras - v_ex(t[1:])))
    err_cent = np.max(np.abs(v_cent - v_ex(t[1:-1])))
    err_5p = np.max(np.abs(v_5p - v_ex(t[2:-2])))
    err_a_cent = np.max(np.abs(a_cent - a_ex(t[1:-1])))
    
    return t, x, v_cent, a_cent, [dt, err_adel, err_atras, err_cent, err_5p, err_a_cent]

t_01, x_01, v_cent_01, a_cent_01, errs_01 = procesar_paso(0.01)
_, _, _, _, errs_005 = procesar_paso(0.005)

# Corrección: Guardar los datos con notación exponencial de alta precisión
datos_dat = np.column_stack((t_01[1:-1], x_01[1:-1], v_cent_01, a_cent_01))
np.savetxt("derivadas_oscilador.dat", datos_dat, fmt="%.15e")

errores_mat = np.array([errs_01, errs_005])
np.savetxt("errores_derivadas.dat", errores_mat, fmt="%14.6e")

# Gráfica
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
ax1.plot(t_01, v_ex(t_01), 'k-', label=r'$v_{exacta}$', alpha=0.7)
ax1.plot(t_01[1:-1], v_cent_01, 'r--', label=r'$v_{central}$')
ax1.set_ylabel('Velocidad [m/s]')
ax1.legend(); ax1.grid(True, alpha=0.3)

ax2.plot(t_01, a_ex(t_01), 'k-', label=r'$a_{exacta}$', alpha=0.7)
ax2.plot(t_01[1:-1], a_cent_01, 'b--', label=r'$a_{central}$')
ax2.set_xlabel('Tiempo [s]')
ax2.set_ylabel('Aceleración [m/s²]')
ax2.legend(); ax2.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("derivadas_oscilador.png")
