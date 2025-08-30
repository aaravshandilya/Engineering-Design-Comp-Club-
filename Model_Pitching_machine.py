import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import scrolledtext
import threading

# ---------------- Ball parameters ----------------
ball_diameter = 0.0508      # meters
ball_radius = ball_diameter / 2
ball_mass = 0.00454         # kg
I_ball = (2/5) * ball_mass * ball_radius**2   # Moment of inertia for solid sphere

# ---------------- Wheel parameters ----------------
wheel_diameter = 0.0635     # meters
wheel_radius = wheel_diameter / 2
rpm_wheel = 200
omega_wheel = (rpm_wheel * 2 * np.pi) / 60   # rad/s
v_surface = omega_wheel * wheel_radius       # surface speed of wheel

# ---------------- Contact mechanics ----------------
mu = 0.7                 # friction coefficient
contact_time = 0.02      # seconds
k_contact = 1e4          # N/m spring constant
c_contact = 0.01         # damping coefficient
n_hertz = 1.5            # Hertzian contact exponent (nonlinear foam contact)

# ---------------- Air properties ----------------
rho_air = 1.225          # kg/m^3
Cd = 0.47                # drag coefficient for sphere
A_cross = np.pi * ball_radius**2
Cl_spin = 0.2            # lift coefficient due to spin (Magnus effect)
spin_decay = 0.1         # approximate decay of spin per second

# ---------------- Simulation ----------------
dt = 1e-4
t_max = 5.0

# ---------------- Contact Dynamics ----------------
def contact_dynamics():
    """
    Simulate the short-time contact of the ball with the wheel.
    Uses a nonlinear spring-damper model to compute the exit velocity and spin.
    """
    x = 0.0
    v = 0.0
    omega = 0.0
    dt_local = 1e-5
    t_local = 0.0
    while t_local < contact_time:
        # Nonlinear spring-damper force (Hertzian contact for foam)
        F_spring = k_contact * x**n_hertz
        F_damp = c_contact * v
        F_normal = F_spring + F_damp

        # Coulomb friction to build spin
        F_friction = np.clip(mu * F_normal, -ball_mass*100, ball_mass*100)

        # Ball acceleration
        a = F_friction / ball_mass
        v += a * dt_local
        x += v * dt_local

        # Torque and angular acceleration
        torque = F_friction * ball_radius
        omega += torque / I_ball * dt_local

        t_local += dt_local
    return v, omega

v_exit, spin_exit = contact_dynamics()

# ---------------- Forces on ball ----------------
def forces(vx, vy, omega):
    """
    Computes net forces on the ball including:
    - Air drag
    - Magnus lift (spin)
    - Gravity
    """
    v = np.sqrt(vx**2 + vy**2)
    if v == 0:
        return 0, -ball_mass * 9.81

    # Drag
    Fd = 0.5 * rho_air * Cd * A_cross * v**2
    Fx_drag = -Fd * (vx / v)
    Fy_drag = -Fd * (vy / v)

    # Magnus effect (spin lift)
    Fm = 0.5 * rho_air * A_cross * Cl_spin * v**2
    Fx_magnus = -Fm * (vy / v)
    Fy_magnus = Fm * (vx / v)

    # Gravity
    Fy_gravity = -ball_mass * 9.81

    Fx_total = Fx_drag + Fx_magnus
    Fy_total = Fy_drag + Fy_magnus + Fy_gravity

    return Fx_total, Fy_total

# ---------------- RK4 Integration ----------------
def rk4_step(x, y, vx, vy, omega, dt):
    """
    RK4 numerical integration for projectile motion with drag and spin lift
    """
    def dv(vx, vy):
        Fx, Fy = forces(vx, vy, omega)
        return Fx/ball_mass, Fy/ball_mass

    k1_vx, k1_vy = dv(vx, vy)
    k1_x, k1_y = vx, vy

    k2_vx, k2_vy = dv(vx + 0.5*k1_vx*dt, vy + 0.5*k1_vy*dt)
    k2_x, k2_y = vx + 0.5*k1_vx*dt, vy + 0.5*k1_vy*dt

    k3_vx, k3_vy = dv(vx + 0.5*k2_vx*dt, vy + 0.5*k2_vy*dt)
    k3_x, k3_y = vx + 0.5*k2_vx*dt, vy + 0.5*k2_vy*dt

    k4_vx, k4_vy = dv(vx + k3_vx*dt, vy + k3_vy*dt)
    k4_x, k4_y = vx + k3_vx*dt, vy + k3_vy*dt

    vx_new = vx + (dt/6)*(k1_vx + 2*k2_vx + 2*k3_vx + k4_vx)
    vy_new = vy + (dt/6)*(k1_vy + 2*k2_vy + 2*k3_vy + k4_vy)
    x_new = x + (dt/6)*(k1_x + 2*k2_x + 2*k3_x + k4_x)
    y_new = y + (dt/6)*(k1_y + 2*k2_y + 2*k3_y + k4_y)
    return x_new, y_new, vx_new, vy_new

# ---------------- Launch Settings ----------------
angle_deg = 20
angle = np.radians(angle_deg)
vx0 = v_exit * np.cos(angle)
vy0 = v_exit * np.sin(angle)
x, y = 0, 1.0
vx, vy = vx0, vy0
omega = spin_exit

# ---------------- Tkinter Window Function ----------------
def show_formulas():
    """
    Shows a detailed Tkinter window with all physics concepts and formulas
    """
    window = tk.Tk()
    window.title("Physics Formulas and Concepts Used")
    window.geometry("800x600")
    
    text = scrolledtext.ScrolledText(window, wrap=tk.WORD, font=("Courier", 11))
    text.pack(expand=True, fill='both')

    content = """
DETAILED PHYSICS CONCEPTS AND FORMULAS:

1. Ball-Wheel Contact Mechanics
   - Nonlinear spring-damper (Hertzian contact for foam):
        F_normal = k * x^n + c * v
     where x = compression, v = compression velocity, n ~ 1.5 for foam
   - Friction builds spin:
        F_friction = mu * F_normal
   - Torque on ball:
        tau = F_friction * radius
        angular acceleration: alpha = tau / I
   - Exit linear velocity and spin are computed over contact_time (~0.02 s)

2. Drag Force
   - Opposes motion through air:
        F_drag = 0.5 * rho_air * Cd * A * v^2
   - Components:
        Fx_drag = -F_drag * vx/v
        Fy_drag = -F_drag * vy/v

3. Magnus Effect (Spin Lift)
   - Lift force perpendicular to velocity due to spin:
        F_magnus = 0.5 * rho_air * A * Cl_spin * v^2
   - Direction depends on rotation:
        Fx_magnus = -F_magnus * vy/v
        Fy_magnus = F_magnus * vx/v
   - Spin decay applied to simulate energy loss

4. Gravity
   - Constant downward force:
        F_gravity = m * g, g = 9.81 m/s^2

5. Motion Integration
   - Runge-Kutta 4th order (RK4) for accuracy in non-linear forces
   - Updates: x, y, vx, vy each time step dt
   - Stable for small dt (~1e-4 s)

6. Simulation Notes
   - Initial launch angle: angle_deg
   - Contact exit velocity: v_exit
   - Magnus spin: omega_exit
   - All forces combined to get realistic projectile trajectory
   - No live plotting: final trajectory plotted after simulation

7. Other Concepts
   - Hertzian contact exponent n ~ 1.5 models foam deformation
   - Moment of inertia for solid sphere: I = 2/5 m r^2
   - Rolling vs sliding friction approximated via Coulomb limit
   - Air density, drag, and lift coefficients calibrated for small foam balls
"""
    text.insert(tk.END, content)
    text.configure(state='disabled')
    window.mainloop()

# ---------------- Run Tkinter in a Separate Thread ----------------
tk_thread = threading.Thread(target=show_formulas)
tk_thread.daemon = True
tk_thread.start()

# ---------------- SIMULATE TRAJECTORY ----------------
xs, ys = [x], [y]
t = 0.0

while y >= 0 and t < t_max:
    x, y, vx, vy = rk4_step(x, y, vx, vy, omega, dt)
    omega = max(0, omega - spin_decay*dt)
    xs.append(x)
    ys.append(y)
    t += dt

# ---------------- PLOT STATIC FINAL TRAJECTORY ----------------
plt.figure(figsize=(8,5))
plt.plot(xs, ys, 'b-', label=f"{v_exit:.2f} m/s, {angle_deg}° launch")
plt.scatter([xs[0]], [ys[0]], color='red', label='Launch point', zorder=5)
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Foam Ball Trajectory with Drag + Magnus Effect")
plt.legend()
plt.grid(True)
plt.show()

print(f"Flight time: {t:.3f} s, Range: {xs[-1]:.3f} m")
