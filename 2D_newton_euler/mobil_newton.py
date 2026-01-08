# ============================================================================
# MOBILE ROBOT – NEWTON-EULER LØSNING (SYMPY)
# ============================================================================
# Opgave:
# Beregn friktionskraften ved forhjul B samt normal-kræfterne
# under hjul A og B for en mobil robot.
#
# Metode:
# - Newtons 2. lov (sum af kræfter)
# - Momentligning omkring tyngdepunktet G
# ============================================================================

import sympy as sp

# ----------------------------------------------------------------------------
# SYMBOLER (ukendte kræfter)
# ----------------------------------------------------------------------------
F_friction, F_A, F_B = sp.symbols("F_friction F_A F_B")

# ----------------------------------------------------------------------------
# KONSTANTER (givne værdier fra opgaven)
# ----------------------------------------------------------------------------
m = 30          # masse [kg]
a = -1.5         # acceleration [m/s^2]
g = 9.81        # tyngdeacceleration [m/s^2]
F_ext = 5       # ekstern kraft [N]
theta = sp.rad(30)  # kraftens vinkel [rad]

# Geometri (afstande fra tyngdepunkt G)
x_A = 0.4       # afstand til hjul A [m]
x_B = 0.3       # afstand til hjul B [m]
h = 0.3         # højde til kraftens angrebspunkt [m]

# ----------------------------------------------------------------------------
# LIGNING 1: SUM AF KRÆFTER I x-RETNING
# ΣFx = m·ax
#
# Retning:
# - Robotten accelererer mod venstre ⇒ ax = -a
# ----------------------------------------------------------------------------
eq_x = (
    -F_friction
    + F_ext * sp.cos(theta)
    + m * a
)

# ----------------------------------------------------------------------------
# LIGNING 2: SUM AF KRÆFTER I y-RETNING
# ΣFy = 0 (ingen vertikal acceleration)
# ----------------------------------------------------------------------------
eq_y = (
    F_A
    + F_B
    + F_ext * sp.sin(theta)
    - m * g
)

# ----------------------------------------------------------------------------
# LIGNING 3: MOMENT OM TYNGDEPUNKT G
# ΣMG = 0
#
# Positiv retning: mod uret
# ----------------------------------------------------------------------------
eq_M = (
    x_A * F_A
    - x_B * F_B
    - h * F_ext * sp.cos(theta)
)

# ----------------------------------------------------------------------------
# LØS LIGNINGSSYSTEMET
# ----------------------------------------------------------------------------
solution = sp.solve(
    [eq_x, eq_y, eq_M],
    [F_friction, F_A, F_B],
    dict=True
)[0]

# ----------------------------------------------------------------------------
# RESULTATER
# ----------------------------------------------------------------------------
print("\nRESULTATER (numerisk):\n")

for key, value in solution.items():
    print(f"{key} = {float(value):.2f} N")

# ----------------------------------------------------------------------------
# PÆN SYMBOLSK UDSKRIFT (til rapport/eksamen)
# ----------------------------------------------------------------------------
print("\nSYMBOLSK FORM:\n")
sp.pprint(solution)
