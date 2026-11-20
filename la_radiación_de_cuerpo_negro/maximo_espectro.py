"""Máximo del espectro solar (T = 5772 K): malla, derivada central y spline.

Entrada : constante_wien.dat  (x5, x3, b)
Salida  : maximo_espectro.dat, espectro_solar.png
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

H, C, KB = 6.62607015e-34, 299792458.0, 1.380649e-23
T = 5772.0
AZUL, DORADO, ROJO, VERDE, GRIS = "#162A54", "#C9972C", "#B02A20", "#2E7D32", "#7F7F7F"


def planck(lam_nm):
    """B_lambda en W m^-2 sr^-1 m^-1 para lam en nm."""
    lam = lam_nm*1e-9
    return 2*H*C**2/lam**5/np.expm1(H*C/(lam*KB*T))


def derivada_central(y, h):
    """Derivada central en los puntos interiores (programada, sin np.gradient)."""
    return (y[2:] - y[:-2])/(2*h)


def cero_maximo(s, cercano_a):
    """Ceros del spline s con segunda derivada negativa; el más cercano a 'cercano_a'."""
    r = s.roots(extrapolate=False)
    r = r[np.isfinite(r)]
    r = r[s(r, 1) < 0]                      # solo máximos (descarta mínimos y ceros espurios)
    if len(r) == 0:
        raise ValueError("no se encontró candidato a máximo")
    return r[np.argmin(np.abs(r - cercano_a))]


x5, x3, b = np.loadtxt("constante_wien.dat")
lam_wien = b/T*1e9                          # lambda_max = b/T [nm]
lam_frec = H*C/(KB*T*x3)*1e9                # lambda del máximo en frecuencia [nm]
print(f"b = {b:.9e} m K   lambda_max (Wien) = {lam_wien:.4f} nm   lambda (max. en frecuencia) = {lam_frec:.2f} nm")

filas = []
for dl in (100.0, 50.0, 25.0, 12.5):
    lam = np.arange(100.0, 3000.0 + dl/2, dl)
    B = planck(lam)
    B = B/B.max()                           # escala 0..1: no cambia la posición de los ceros
    la = lam[np.argmax(B)]                                              # (a) mayor valor
    lb = cero_maximo(CubicSpline(lam[1:-1], derivada_central(B, dl)), la)   # (b) spline de la derivada central
    lc = cero_maximo(CubicSpline(lam, B).derivative(), la)              # (c) derivada del spline de B
    e = [abs(la - lam_wien), abs(lb - lam_wien), abs(lc - lam_wien)]
    filas.append([dl, *e])
    print(f"dl = {dl:6.1f} nm   (a) {la:9.3f}   (b) {lb:9.3f}   (c) {lc:9.3f}   "
          f"errores [nm]: {e[0]:8.3f} {e[1]:8.3f} {e[2]:8.4f}")

filas = np.array(filas)
np.savetxt("maximo_espectro.dat", filas, fmt="%.16e")
print("\nrazón entre errores consecutivos (b, c):")
print(np.round(filas[:-1, 2:]/filas[1:, 2:], 2))

# ---- figura: malla de 25 nm
dl = 25.0
lam = np.arange(100.0, 3000.0 + dl/2, dl)
B = planck(lam)
lf = np.linspace(100, 3000, 1500)
fig, ax = plt.subplots(figsize=(8, 4.4))
ax.axvspan(380, 750, color=DORADO, alpha=0.12, label="visible (380-750 nm)")
ax.plot(lf, CubicSpline(lam, B)(lf)*1e-12, color=GRIS, lw=1, alpha=0.7)
ax.plot(lam, B*1e-12, "o", ms=3, color=AZUL, label="$B_\\lambda$, malla de 25 nm")
ax.axvline(lam_wien, color=ROJO, ls="--", label=f"$\\lambda_{{\\max}} = b/T$ = {lam_wien:.1f} nm")
ax.axvline(lam_frec, color=VERDE, ls=":", lw=2, label=f"máximo en frecuencia: {lam_frec:.0f} nm")
ax.set_xlabel("longitud de onda $\\lambda$ [nm]")
ax.set_ylabel("$B_\\lambda$ [$10^{12}$ W m$^{-2}$ sr$^{-1}$ m$^{-1}$]")
ax.set_title("Espectro del Sol como cuerpo negro, $T = 5772$ K")
ax.legend(); fig.tight_layout(); fig.savefig("espectro_solar.png", dpi=150)
print("figura: espectro_solar.png")
