import sympy as sp

sp.init_printing()  # enables pretty printing

# ============================
# SYMBOLS
# ============================

t = sp.symbols('t')

# Generalized coordinates
x = sp.Function('x')(t)
theta = sp.Function('theta')(t)

# Time derivatives
xd = sp.diff(x, t)
xdd = sp.diff(x, t, 2)

thetad = sp.diff(theta, t)
thetadd = sp.diff(theta, t, 2)

# Parameters
m1, m2, g, d, l2, I2 = sp.symbols('m1 m2 g d l2 I2')
F, tau = sp.symbols('F tau')

# ============================
# POSITIONS OF MASS CENTERS
# ============================

# sliding block C1
xC1 = x
yC1 = d

# pivot A
xA = x
yA = 0

# midpoint of link AB = C2
xC2 = x + (l2/2)*sp.sin(theta)
yC2 = - (l2/2)*sp.cos(theta)

# ============================
# VELOCITIES
# ============================

vxC2 = sp.diff(xC2, t)
vyC2 = sp.diff(yC2, t)

vC2_sq = sp.simplify(vxC2**2 + vyC2**2)
vC1_sq = sp.simplify(xd**2)

# ============================
# ENERGIES
# ============================

# Kinetic energy
T1 = sp.Rational(1,2)*m1*vC1_sq
T2 = sp.Rational(1,2)*m2*vC2_sq + sp.Rational(1,2)*I2*thetad**2
T = sp.simplify(T1 + T2)

# Potential energy
V1 = m1*g*yC1
V2 = m2*g*yC2
V = sp.simplify(V1 + V2)

# Lagrangian
L = sp.simplify(T - V)

# ============================
# LAGRANGE EQUATIONS
# ============================

Qx = F
Qtheta = tau

EOMx = sp.diff(sp.diff(L, xd), t) - sp.diff(L, x) - Qx
EOMtheta = sp.diff(sp.diff(L, thetad), t) - sp.diff(L, theta) - Qtheta

# ============================
# STRONG SIMPLIFICATION FUNCTION
# ============================

def clean(expr):
    return sp.simplify(
        sp.trigsimp(
            sp.factor(
                sp.expand(expr)
            )
        )
    )

# apply simplification
vC2_sq_s = clean(vC2_sq)
T_s = clean(T)
V_s = clean(V)
L_s = clean(L)
EOMx_s = clean(EOMx)
EOMtheta_s = clean(EOMtheta)

# ============================
# OUTPUT
# ============================

print("\nVelocity of C2 squared:")
sp.pprint(vC2_sq_s)

print("\nKinetic Energy T:")
sp.pprint(T_s)

print("\nPotential Energy V:")
sp.pprint(V_s)

print("\nLagrangian L = T - V:")
sp.pprint(L_s)

print("\nLagrange Equation in x (force F):")
sp.pprint(EOMx_s)

print("\nLagrange Equation in theta (torque tau):")
sp.pprint(EOMtheta_s)
