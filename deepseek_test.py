import sympy as sp
import math as m
import numpy as np

# ======================
# KONSTANTER
# ======================
l1 = 0.5
l2 = 0.3
m1 = 3.0
m2 = 2.0
I1 = 0.5
I2 = 0.3
g = 9.8  # Fasit bruger 9.8

# ======================
# 1. OPGAVE 1 - KINEMATIK
# ======================
time = 0.1

# Fra dit script (korrekt):
theta1 = 3 * m.sin(3.141*time + 1.1)
theta2 = 2 * m.sin(3.141*time + 0.4)

om1 = 9.42*m.cos(3.14*time + 1.1)
om2 = 6.28*m.cos(3.14*time + 0.4) + 9.42*m.cos(3.14*time + 1.1)

acc1 = -29.5788*m.sin(3.14*time + 1.1)
acc2 = -19.7192*m.sin(3.14*time + 0.4) - 29.5788*m.sin(3.14*time + 1.1)

print("θ1 og θ2 ved t=0.1s:")
print(f"θ1 = {theta1:.6f} rad")
print(f"θ2 = {theta2:.6f} rad")

# Trigonometri
c1 = m.cos(theta1)
s1 = m.sin(theta1)
phi = theta1 + theta2
c12 = m.cos(phi)
s12 = m.sin(phi)

# ======================
# 2. ACCELERATIONER MED KORREKTE AKSER
# ======================
# C1 acceleration - RIGTIG FORMEL:
# Position: r_C1 = (l1/2 cosθ1, l1/2 sinθ1)
# Hastighed: v_C1 = (-l1/2 sinθ1 * ω1, l1/2 cosθ1 * ω1)
# Acceleration: a_C1 = (-l1/2 sinθ1 * α1 - l1/2 cosθ1 * ω1², 
#                       l1/2 cosθ1 * α1 - l1/2 sinθ1 * ω1²)

a_C1x = -0.5*l1*s1*acc1 - 0.5*l1*c1*om1**2
a_C1y = 0.5*l1*c1*acc1 - 0.5*l1*s1*om1**2

print(f"\na_C1 = [{a_C1x:.6f}, {a_C1y:.6f}] m/s²")

# A acceleration
a_Ax = -l1*s1*acc1 - l1*c1*om1**2
a_Ay = l1*c1*acc1 - l1*s1*om1**2

# C2 acceleration
omega_AB = om2
alpha_AB = acc2

# Position: r_C2 = (l1 cosθ1 + l2/2 cos(θ1+θ2), l1 sinθ1 + l2/2 sin(θ1+θ2))
# Acceleration fra kinematik:
a_C2x = a_Ax - 0.5*l2*s12*alpha_AB - 0.5*l2*c12*omega_AB**2
a_C2y = a_Ay + 0.5*l2*c12*alpha_AB - 0.5*l2*s12*omega_AB**2

print(f"a_C2 = [{a_C2x:.6f}, {a_C2y:.6f}] m/s²")

# ======================
# 3. NEWTON'S 2. LOV MED KORREKTE AKSER
# ======================
# Lad os definere akserne:
# x: horisontal, positiv til højre
# y: vertikal, positiv opad
# g = -9.8ĵ (tyngdekraft nedad)

# 6 ligninger, 6 ubekendte: Ox, Oy, Ax, Ay, τ1, τ2
A = np.zeros((6, 6))
b = np.zeros(6)

# Ligning 1: OA translation x: m1*a_C1x = Ox - Ax
A[0, 0] = 1   # Ox
A[0, 2] = -1  # -Ax
b[0] = m1 * a_C1x

# Ligning 2: OA translation y: m1*a_C1y = Oy - Ay - m1*(-g) 
# Tyngdekraft er m1*g NEDAD, dvs. -m1*g i y-retning
# Så: m1*a_C1y = Oy - Ay - (-m1*g) = Oy - Ay + m1*g
# Alternativt: m1*a_C1y = Oy - Ay + m1*g
# Eller: m1*a_C1y - m1*g = Oy - Ay
# Fasit har: Oy - Ay - 29.4 = -5.496, så 29.4 = m1*g = 3*9.8
# Derfor: Oy - Ay = -5.496 + 29.4 = 23.904

# Så den rigtige er: m1*a_C1y = Oy - Ay + m1*g
# Eller: m1*a_C1y - m1*g = Oy - Ay
A[1, 1] = 1   # Oy
A[1, 3] = -1  # -Ay
b[1] = m1 * a_C1y - m1 * g

# Ligning 3: AB translation x: m2*a_C2x = Ax
A[2, 2] = 1   # Ax
b[2] = m2 * a_C2x

# Ligning 4: AB translation y: m2*a_C2y = Ay - m2*(-g) = Ay + m2*g
# Eller: m2*a_C2y - m2*g = Ay
A[3, 3] = 1   # Ay
b[3] = m2 * a_C2y - m2 * g

# Ligning 5: OA rotation om C1 
# Moment af kraft (Fx, Fy) på punkt (x,y): M = x*Fy - y*Fx
# For kraft i O på punkt O relativt til C1:
# r_O/C1 = (-l1/2 cosθ1, -l1/2 sinθ1)
# M_O = (-l1/2 cosθ1)*Oy - (-l1/2 sinθ1)*Ox = -l1/2 cosθ1*Oy + l1/2 sinθ1*Ox

# For kraft -A på punkt A relativt til C1:
# r_A/C1 = (l1/2 cosθ1, l1/2 sinθ1)
# M_A = (l1/2 cosθ1)*(-Ay) - (l1/2 sinθ1)*(-Ax) = -l1/2 cosθ1*Ay + l1/2 sinθ1*Ax

# Så: I1*α1 = τ1 - τ2 + M_O + M_A
#     = τ1 - τ2 + l1/2 sinθ1*Ox - l1/2 cosθ1*Oy + l1/2 sinθ1*Ax - l1/2 cosθ1*Ay

A[4, 0] = 0.5*l1 * s1    # Ox: (l1/2)sinθ1
A[4, 1] = -0.5*l1 * c1   # Oy: -(l1/2)cosθ1  
A[4, 2] = 0.5*l1 * s1    # Ax: (l1/2)sinθ1
A[4, 3] = -0.5*l1 * c1   # Ay: -(l1/2)cosθ1
A[4, 4] = 1              # τ1
A[4, 5] = -1             # -τ2
b[4] = I1 * acc1

# Ligning 6: AB rotation om C2
# r_A/C2 = (-l2/2 cos(θ1+θ2), -l2/2 sin(θ1+θ2))
# Moment af A om C2: M = r_x*Fy - r_y*Fx
# = (-l2/2 c12)*Ay - (-l2/2 s12)*Ax = -l2/2 c12*Ay + l2/2 s12*Ax

A[5, 2] = 0.5*l2 * s12   # Ax: (l2/2)sin(θ1+θ2)
A[5, 3] = -0.5*l2 * c12  # Ay: -(l2/2)cos(θ1+θ2)
A[5, 5] = 1              # τ2
b[5] = I2 * acc2

# ======================
# 4. LØS SYSTEMET
# ======================
x = np.linalg.solve(A, b)
Ox, Oy, Ax, Ay, tau1, tau2 = x

print(f"\n{'='*60}")
print("OPGAVE 3 & 4 - KORREKTE RESULTATER:")
print(f"{'='*60}")
print(f"Ox = {Ox:.4f} N")
print(f"Oy = {Oy:.4f} N")
print(f"Ax = {Ax:.4f} N")
print(f"Ay = {Ay:.4f} N")
print(f"τ1 = {tau1:.4f} N·m")
print(f"τ2 = {tau2:.4f} N·m")

# ======================
# 5. SAMMENLIGN MED FASIT
# ======================
print(f"\n{'='*60}")
print("SAMMENLIGN MED FASIT:")
print(f"{'='*60}")
print(f"1. Ox - Ax = {Ox - Ax:.4f} (fasit: 21.30)")
print(f"2. Oy - Ay - {m1*g:.1f} = {Oy - Ay - m1*g:.4f} (fasit: -5.496)")

# Beregn fasit-ligningerne
fasit_eq3 = tau1 - tau2 + 0.4920*Ax - 0.08880*Ay - 1.3055
fasit_eq4 = tau2 + 0.06380*Ax + 0.1358*Ay

print(f"3. τ1 - τ2 + 0.4920*Ax - 0.08880*Ay - 1.3055 = {fasit_eq3:.4f} (fasit: -20.11)")
print(f"4. τ2 + 0.06380*Ax + 0.1358*Ay = {fasit_eq4:.4f} (fasit: -12.65)")

# ======================
# 6. OPSKRIV LIGNINGER SOM I FASIT
# ======================
print(f"\n{'='*60}")
print("OPGAVE 3 - LIGNINGER AF BEVÆGELSE (Newton's 2. lov):")
print(f"{'='*60}")

print(f"Link OA:")
print(f"  m1*a_C1x = Ox - Ax")
print(f"  {m1} * {a_C1x:.4f} = Ox - Ax")
print(f"  => Ox - Ax = {m1*a_C1x:.4f}")

print(f"\n  m1*a_C1y = Oy - Ay + m1*g")
print(f"  {m1} * {a_C1y:.4f} = Oy - Ay + {m1*g:.1f}")
print(f"  => Oy - Ay = {m1*a_C1y:.4f} - {m1*g:.1f} = {m1*a_C1y - m1*g:.4f}")

print(f"\n  I1*α1 = τ1 - τ2 + (l1/2)[sinθ1*(Ox+Ax) - cosθ1*(Oy+Ay)]")
print(f"  {I1} * {acc1:.4f} = τ1 - τ2 + {0.5*l1*s1:.4f}*(Ox+Ax) - {0.5*l1*c1:.4f}*(Oy+Ay)")

print(f"\nLink AB:")
print(f"  m2*a_C2x = Ax")
print(f"  {m2} * {a_C2x:.4f} = Ax")
print(f"  => Ax = {m2*a_C2x:.4f}")

print(f"\n  m2*a_C2y = Ay + m2*g")
print(f"  {m2} * {a_C2y:.4f} = Ay + {m2*g:.1f}")
print(f"  => Ay = {m2*a_C2y:.4f} - {m2*g:.1f} = {m2*a_C2y - m2*g:.4f}")

print(f"\n  I2*α_AB = τ2 + (l2/2)[sin(θ1+θ2)*Ax - cos(θ1+θ2)*Ay]")
print(f"  {I2} * {acc2:.4f} = τ2 + {0.5*l2*s12:.4f}*Ax - {0.5*l2*c12:.4f}*Ay")

# ======================
# 7. KONTROLLER TRIG-VÆRDIER
# ======================
print(f"\n{'='*60}")
print("TRIG-VÆRDIER TIL FASIT:")
print(f"{'='*60}")
print(f"(l1/2)sinθ1 = {0.5*l1*s1:.6f} ≈ 0.4920?")
print(f"-(l1/2)cosθ1 = {-0.5*l1*c1:.6f} ≈ -0.08880?")
print(f"(l2/2)sin(θ1+θ2) = {0.5*l2*s12:.6f} ≈ 0.06380?")
print(f"-(l2/2)cos(θ1+θ2) = {-0.5*l2*c12:.6f} ≈ 0.1358?")
print("\nNB: Fasit har fortegn omvendt i nogle tilfælde - måske pga deres moment-definition.")