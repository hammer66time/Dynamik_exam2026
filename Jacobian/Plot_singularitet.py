# plot_singularities.py
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

# --------------------------
# Symbolic setup (optional, for clarity / exact det(J))
# --------------------------
theta1, theta2, s = sp.symbols("theta1 theta2 s")
y1 = 5 * sp.cos(theta1) + 3 * sp.cos(theta1 + theta2)
y2 = 5 * sp.sin(theta1) + 3 * sp.sin(theta1 + theta2)
y3 = 0.5 * s + 3

y = sp.Matrix([y1, y2, y3])
q = sp.Matrix([theta1, theta2, s])
J_sym = y.jacobian(q)
detJ_sym = sp.simplify(J_sym.det())  # should simplify to 7.5*sin(theta2)

# If you want the LaTeX string for reports:
# detJ_latex = sp.latex(detJ_sym)

# --------------------------
# Numeric / plotting utilities
# --------------------------
L1 = 5.0
L2 = 3.0

def end_effector(theta1_val, theta2_val):
    """Return (x,y) end-effector coordinates for given angles"""
    x = L1 * np.cos(theta1_val) + L2 * np.cos(theta1_val + theta2_val)
    y = L1 * np.sin(theta1_val) + L2 * np.sin(theta1_val + theta2_val)
    return x, y

def plot_singular_loci(save_path=None, show=True):
    """Plot loci for theta2=0 and theta2=pi and sample robot poses."""
    theta1_vals = np.linspace(-np.pi, np.pi, 400)
    theta2_cases = [0.0, np.pi]

    fig, ax = plt.subplots(figsize=(6,6))

    # loci curves
    for t2 in theta2_cases:
        xs = L1 * np.cos(theta1_vals) + L2 * np.cos(theta1_vals + t2)
        ys = L1 * np.sin(theta1_vals) + L2 * np.sin(theta1_vals + t2)
        label = f"θ₂ = {t2:.2f} rad (radius ~ {abs(L1 + (L2 if t2==0 else -L2)):.1f})"
        ax.plot(xs, ys, label=label, linewidth=2)

    # base marker
    ax.scatter([0], [0], zorder=5)
    ax.text(0, 0, " base (0,0)", va='bottom', ha='left')

    # sample poses to show collinearity
    sample_theta1 = [-1.0, -0.2, 0.5, 1.2]
    for t2 in theta2_cases:
        for th1 in sample_theta1:
            x0, y0 = 0.0, 0.0
            x1 = L1 * np.cos(th1)
            y1 = L1 * np.sin(th1)
            x2 = x1 + L2 * np.cos(th1 + t2)
            y2 = y1 + L2 * np.sin(th1 + t2)
            ax.plot([x0, x1, x2], [y0, y1, y2], marker='o', linewidth=1, markersize=4)

    ax.set_aspect('equal', 'box')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('End-effector loci & sample poses for singular configurations\n(θ₂ = 0 -> stretched, θ₂ = π -> folded)')
    ax.legend()
    ax.grid(True)

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

def plot_determinant(save_path=None, show=True):
    """Plot det(J) = 7.5*sin(theta2) across a range and mark zeros."""
    # Use the simplified symbolic determinant if desired
    # detJ_sym is a sympy expression; convert to numeric function:
    detJ_func = sp.lambdify(theta2, detJ_sym, "numpy")

    theta2_plot = np.linspace(-2*np.pi, 2*np.pi, 800)
    det_vals = detJ_func(theta2_plot)

    fig, ax = plt.subplots(figsize=(8,4))
    ax.plot(theta2_plot, det_vals, linewidth=2)
    ax.axhline(0, linestyle='--')
    zeros = np.arange(-2*np.pi, 2.001*np.pi, np.pi)
    ax.scatter(zeros, np.zeros_like(zeros), marker='x', zorder=6)

    ax.set_xlabel('θ₂ (rad)')
    ax.set_ylabel('det(J)')
    ax.set_title(f"det(J) = {sp.srepr(detJ_sym)}  -> zeros at θ₂ = k·π")
    ax.grid(True)

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    plt.close(fig)

# --------------------------
# Example main - call from your project
# --------------------------
if __name__ == "__main__":
    # Show plots interactively:
    plot_singular_loci(save_path="singular_loci.png", show=True)
    plot_determinant(save_path="determinant_plot.png", show=True)

    # Or, if you want to only save and not show, call:
    # plot_singular_loci(save_path="singular_loci.png", show=False)
    # plot_determinant(save_path="determinant_plot.png", show=False)
