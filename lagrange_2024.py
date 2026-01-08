import sympy as sp

# ==================================================
# SYMBOLDEFINITIONER
# ==================================================

# Generaliserede koordinater og hastigheder
h, theta = sp.symbols('h theta', real=True)
v, omega = sp.symbols('v omega', real=True)   # givet i opgaven

# Parametre
m1, m2 = sp.symbols('m1 m2', positive=True)
l, s = sp.symbols('l s', positive=True)
J2 = sp.symbols('J2', positive=True)
g = sp.symbols('g', positive=True)

# ==================================================
# OPGAVE (2): KINETISK OG POTENTIEL ENERGI
# ==================================================

# --------------------------------------------------
# Legeme 1: Glidende blok (ren translation)
# --------------------------------------------------

T1 = sp.Rational(1, 2) * m1 * v**2
V1 = m1 * g * h

# --------------------------------------------------
# Legeme 2: Link AB
# --------------------------------------------------
# Hastighed af tyngdepunkt C2
vC2_sq = (
    v**2
    + (omega * l / 2)**2
    + v * omega * l * sp.cos(theta)
)

# Kinetisk energi
T2 = (
    sp.Rational(1, 2) * m2 * vC2_sq
    + sp.Rational(1, 2) * J2 * omega**2
)

# Potentiel energi
V2 = m2 * g * (h + s + l/2 * sp.sin(theta))

# Samlet energi
T = sp.simplify(T1 + T2)
V = sp.simplify(V1 + V2)

print("\n--- KINETISK ENERGI ---")
sp.pretty_print(sp.collect(sp.expand(T), [v, omega]))

print("\n--- POTENTIEL ENERGI ---")
sp.pretty_print(V)

# ==================================================
# OPGAVE (3): LAGRANGE-FUNKTION
# ==================================================

Lagrangian = T - V

print("\n--- LAGRANGE-FUNKTION L = T - V ---")
sp.pretty_print(Lagrangian)

# ==================================================
# OPGAVE (4): LAGRANGE-LIGNINGER
# ==================================================
# Generaliserede kræfter:
# Q_h     = F  (ydre kraft på blokken)
# Q_theta = tau (ydre moment i leddet)

F, tau = sp.symbols('F tau', real=True)

# Da v = dh/dt og omega = dtheta/dt er givet,
# differentierer vi mht. v og omega direkte

# Kraft-ligning (translation)
dL_dv = sp.diff(Lagrangian, v)
eq_F = sp.Eq(sp.diff(dL_dv, h), F)

# Moment-ligning (rotation)
dL_domega = sp.diff(Lagrangian, omega)
eq_tau = sp.Eq(sp.diff(dL_domega, theta), tau)

print("\n--- LAGRANGE-LIGNING (kraft F) ---")
sp.pretty_print(eq_F)

print("\n--- LAGRANGE-LIGNING (moment tau) ---")
sp.pretty_print(eq_tau)

# ==================================================
# LATEX (KLAR TIL WORD)
# ==================================================

print("\n--- LATEX: T ---")
print(sp.latex(T))

print("\n--- LATEX: V ---")
print(sp.latex(V))

print("\n--- LATEX: L ---")
print(sp.latex(Lagrangian))
