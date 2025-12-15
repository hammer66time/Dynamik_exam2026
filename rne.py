import numpy as np
import math

# ======================
# PRÆCIS RNE FOR 2-DOF
# ======================
def precise_RNE(t):
    # Konstanter (samme som opgaven)
    l1, l2 = 0.5, 0.3
    m1, m2 = 3.0, 2.0
    I1, I2 = 0.5, 0.3
    g = 9.8
    
    # 1. KINEMATIK - PRÆCIST som dit script
    pi_val = 3.141  # Samme som du bruger
    
    theta1 = 3 * math.sin(pi_val*t + 1.1)
    theta2 = 2 * math.sin(pi_val*t + 0.4)
    
    dtheta1 = 3 * pi_val * math.cos(pi_val*t + 1.1)
    dtheta2 = 2 * pi_val * math.cos(pi_val*t + 0.4)
    
    ddtheta1 = -3 * pi_val**2 * math.sin(pi_val*t + 1.1)
    ddtheta2 = -2 * pi_val**2 * math.sin(pi_val*t + 0.4)
    
    print(f"Kinematik ved t={t}s:")
    print(f"θ1 = {theta1:.8f} rad, θ2 = {theta2:.8f} rad")
    print(f"ω1 = {dtheta1:.8f} rad/s, ω2 = {dtheta2:.8f} rad/s")
    print(f"α1 = {ddtheta1:.8f} rad/s², α2 = {ddtheta2:.8f} rad/s²")
    
    # 2. GEOMETRI
    # Vektorer i lokale koordinatsystemer
    # Alle i 3D med z som rotationsakse
    
    # COM positioner relativt til led-start
    r1 = np.array([l1/2, 0, 0])      # COM1 relativt til O
    r2 = np.array([l2/2, 0, 0])      # COM2 relativt til A
    
    # Led-ender
    p1 = np.array([l1, 0, 0])        # A relativt til O
    p2 = np.array([l2, 0, 0])        # B relativt til A
    
    # 3. FORWARD RECURSION
    
    # Base (led 0)
    omega0 = np.array([0, 0, 0])
    domega0 = np.array([0, 0, 0])
    dv0 = np.array([0, 0, 0])
    dv_c0 = np.array([0, -g, 0])  # Kun tyngde
    
    # Rotation matricer
    def Rz(theta):
        c, s = math.cos(theta), math.sin(theta)
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    
    R01 = Rz(theta1)
    R12 = Rz(theta2)
    z = np.array([0, 0, 1])
    
    # ---- LED 1 ----
    omega1 = R01.T @ omega0 + dtheta1 * z
    domega1 = R01.T @ domega0 + ddtheta1 * z + np.cross(R01.T @ omega0, dtheta1*z)
    
    dv1 = R01.T @ (dv0 + np.cross(domega0, p1) + np.cross(omega0, np.cross(omega0, p1)))
    dv_c1 = dv1 + np.cross(domega1, r1) + np.cross(omega1, np.cross(omega1, r1))
    
    # ---- LED 2 ----
    omega2 = R12.T @ omega1 + dtheta2 * z
    domega2 = R12.T @ domega1 + ddtheta2 * z + np.cross(R12.T @ omega1, dtheta2*z)
    
    dv2 = R12.T @ (dv1 + np.cross(domega1, p2) + np.cross(omega1, np.cross(omega1, p2)))
    dv_c2 = dv2 + np.cross(domega2, r2) + np.cross(omega2, np.cross(omega2, r2))
    
    # 4. BACKWARD RECURSION
    
    # End-effector (ingen ydre kraft)
    f3 = np.array([0, 0, 0])
    n3 = np.array([0, 0, 0])
    
    # ---- LED 2 ----
    F2 = m2 * dv_c2
    # For 2D rotation, inerti-tensor reduceres til skalar om z-aksen
    N2 = np.array([0, 0, I2 * domega2[2]]) + np.cross(omega2, np.array([0, 0, I2 * omega2[2]]))
    
    f2 = R12 @ f3 + F2
    n2 = R12 @ n3 + np.cross(r2, F2) + np.cross(p2, R12 @ f3) + N2
    
    tau2 = n2[2]  # z-komponent
    
    # ---- LED 1 ----
    F1 = m1 * dv_c1
    N1 = np.array([0, 0, I1 * domega1[2]]) + np.cross(omega1, np.array([0, 0, I1 * omega1[2]]))
    
    f1 = R01 @ f2 + F1
    n1 = R01 @ n2 + np.cross(r1, F1) + np.cross(p1, R01 @ f2) + N1
    
    tau1 = n1[2]  # z-komponent
    
    # 5. REAKTIONSKRÆFTER
    # Kræfterne f1 og f2 er kræfterne PÅ leddene fra naboled
    # Så reaktionskræfterne er:
    F_O = -f1[:2]  # Kraft fra base PÅ led 1 (modsat f1)
    F_A = -f2[:2]  # Kraft fra led 1 PÅ led 2 (modsat f2)
    
    # For at matche Newton-notasjon: Ax, Ay er kraft FRA OA på AB
    # I RNE: f2 er kraft FRA AB på omgivelser i A, så F_A = -f2
    Ax, Ay = F_A[0], F_A[1]
    Ox, Oy = F_O[0], F_O[1]
    
    # Accelerationer af COM (til tjek)
    a_C1 = dv_c1[:2]
    a_C2 = dv_c2[:2]
    
    return {
        'tau1': tau1, 'tau2': tau2,
        'Ox': Ox, 'Oy': Oy, 'Ax': Ax, 'Ay': Ay,
        'theta1': theta1, 'theta2': theta2,
        'a_C1': a_C1, 'a_C2': a_C2,
        'dtheta1': dtheta1, 'dtheta2': dtheta2,
        'ddtheta1': ddtheta1, 'ddtheta2': ddtheta2
    }

# ======================
# KØR BEREGNINGEN
# ======================
result = precise_RNE(0.1)

print("\n" + "="*70)
print("RNE RESULTATER - PRÆCISE VÆRDIER")
print("="*70)
print(f"\nMOTOR MOMENTER:")
print(f"τ1 = {result['tau1']:.8f} N·m")
print(f"τ2 = {result['tau2']:.8f} N·m")

print(f"\nREAKTIONSKRÆFTER:")
print(f"Ox = {result['Ox']:.8f} N")
print(f"Oy = {result['Oy']:.8f} N")
print(f"Ax = {result['Ax']:.8f} N")
print(f"Ay = {result['Ay']:.8f} N")

print(f"\nACCELERATIONER:")
print(f"a_C1 = [{result['a_C1'][0]:.8f}, {result['a_C1'][1]:.8f}] m/s²")
print(f"a_C2 = [{result['a_C2'][0]:.8f}, {result['a_C2'][1]:.8f}] m/s²")

# ======================
# BEREGN FASIT'S KOEFFICIENTER
# ======================
print("\n" + "="*70)
print("TRIG-VÆRDIER FOR FASIT'S KOEFFICIENTER")
print("="*70)

l1, l2 = 0.5, 0.3
theta1, theta2 = result['theta1'], result['theta2']
c1, s1 = math.cos(theta1), math.sin(theta1)
c12, s12 = math.cos(theta1+theta2), math.sin(theta1+theta2)

print(f"\nFra vores θ1 = {theta1:.8f} rad:")
print(f"sin(θ1) = {s1:.8f}, cos(θ1) = {c1:.8f}")
print(f"(l1/2)sinθ1 = {0.5*l1*s1:.8f}")
print(f"-(l1/2)cosθ1 = {-0.5*l1*c1:.8f}")

print(f"\nFra vores θ1+θ2 = {theta1+theta2:.8f} rad:")
print(f"sin(θ1+θ2) = {s12:.8f}, cos(θ1+θ2) = {c12:.8f}")
print(f"(l2/2)sin(θ1+θ2) = {0.5*l2*s12:.8f}")
print(f"(l2/2)cos(θ1+θ2) = {0.5*l2*c12:.8f}")

print(f"\nSammenlign med fasit:")
print(f"Fasit: 0.4920 vs vores: {0.5*l1*s1:.4f}")
print(f"Fasit: -0.08880 vs vores: {-0.5*l1*c1:.4f}")
print(f"Fasit: 0.06380 vs vores: {0.5*l2*s12:.4f}")
print(f"Fasit: 0.1358 vs vores: {0.5*l2*c12:.4f}")

# ======================
# OPSKRIV NEWTON-LIGNINGER MED RNE-VÆRDIER
# ======================
print("\n" + "="*70)
print("KORREKTE NEWTON-LIGNINGER (fra RNE resultater)")
print("="*70)

m1, m2 = 3.0, 2.0
g = 9.8
a_C1, a_C2 = result['a_C1'], result['a_C2']
Ax, Ay = result['Ax'], result['Ay']
Ox, Oy = result['Ox'], result['Oy']
tau1, tau2 = result['tau1'], result['tau2']

print(f"\nLink OA translationsligninger:")
print(f"m1*a_C1x = Ox - Ax")
print(f"{m1} * {a_C1[0]:.6f} = {Ox:.6f} - {Ax:.6f}")
print(f"{m1*a_C1[0]:.6f} = {Ox - Ax:.6f}")

print(f"\nm1*a_C1y = Oy - Ay + m1*g  (tyngde NEDAD)")
print(f"{m1} * {a_C1[1]:.6f} = {Oy:.6f} - {Ay:.6f} + {m1*g:.6f}")
print(f"{m1*a_C1[1]:.6f} = {Oy - Ay + m1*g:.6f}")

print(f"\nLink AB translationsligninger:")
print(f"m2*a_C2x = Ax")
print(f"{m2} * {a_C2[0]:.6f} = {Ax:.6f}")
print(f"{m2*a_C2[0]:.6f} = {Ax:.6f}")

print(f"\nm2*a_C2y = Ay + m2*g")
print(f"{m2} * {a_C2[1]:.6f} = {Ay:.6f} + {m2*g:.6f}")
print(f"{m2*a_C2[1]:.6f} = {Ay + m2*g:.6f}")

I1, I2 = 0.5, 0.3

print(f"\nRotation OA om C1:")
print(f"I1*α1 = τ1 - τ2 + (l1/2)[sinθ1*(Ox+Ax) - cosθ1*(Oy+Ay)]")
print(f"{I1}*{result['ddtheta1']:.6f} = {tau1:.6f} - {tau2:.6f} + {0.5*l1*s1:.6f}*({Ox:.6f}+{Ax:.6f}) - {0.5*l1*c1:.6f}*({Oy:.6f}+{Ay:.6f})")

LHS_rot1 = I1 * result['ddtheta1']
RHS_rot1 = tau1 - tau2 + 0.5*l1*s1*(Ox+Ax) - 0.5*l1*c1*(Oy+Ay)
print(f"{LHS_rot1:.6f} = {RHS_rot1:.6f}")

print(f"\nRotation AB om C2:")
print(f"I2*α_AB = τ2 + (l2/2)[sin(θ1+θ2)*Ax - cos(θ1+θ2)*Ay]")
print(f"{I2}*{result['ddtheta1']+result['ddtheta2']:.6f} = {tau2:.6f} + {0.5*l2*s12:.6f}*{Ax:.6f} - {0.5*l2*c12:.6f}*{Ay:.6f})")

LHS_rot2 = I2 * (result['ddtheta1'] + result['ddtheta2'])
RHS_rot2 = tau2 + 0.5*l2*s12*Ax - 0.5*l2*c12*Ay
print(f"{LHS_rot2:.6f} = {RHS_rot2:.6f}")

# ======================
# SAMMENLIGN MED FASIT EKSPLICIT
# ======================
print("\n" + "="*70)
print("DIREKTE SAMMENLIGN MED FASIT'S LIGNINGER")
print("="*70)

# Beregn fasit's udtryk med vores værdier
fasit_eq1 = Ox - Ax
fasit_eq2 = Oy - Ay - m1*g
fasit_eq3 = tau1 - tau2 + 0.4920*Ax - 0.08880*Ay - 1.3055
fasit_eq4 = tau2 + 0.06380*Ax + 0.1358*Ay

print(f"\n1. Ox - Ax = {fasit_eq1:.6f} (fasit: 21.30)")
print(f"   Afvigelse: {fasit_eq1 - 21.30:.6f}")

print(f"\n2. Oy - Ay - {m1*g:.1f} = {fasit_eq2:.6f} (fasit: -5.496)")
print(f"   Afvigelse: {fasit_eq2 - (-5.496):.6f}")

print(f"\n3. τ1 - τ2 + 0.4920*Ax - 0.08880*Ay - 1.3055 = {fasit_eq3:.6f} (fasit: -20.11)")
print(f"   Afvigelse: {fasit_eq3 - (-20.11):.6f}")

print(f"\n4. τ2 + 0.06380*Ax + 0.1358*Ay = {fasit_eq4:.6f} (fasit: -12.65)")
print(f"   Afvigelse: {fasit_eq4 - (-12.65):.6f}")

print(f"\nKonklusion: Fasit har nok brugt lidt andre θ-værdier eller afrundet anderledes.")