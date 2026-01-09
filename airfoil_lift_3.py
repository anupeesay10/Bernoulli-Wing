import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# Physical constants
# -------------------------------
rho = 1.225  # air density (kg/m^3)
p0 = 101325  # atmospheric pressure (Pa)

# -------------------------------
# Airfoil geometry
# -------------------------------
def airfoil_camber(x, camber=0.05):
    """Simple cambered airfoil shape"""
    return camber * np.sin(np.pi * x)

# -------------------------------
# Velocity field (Bernoulli-based)
# -------------------------------
def velocity_field(x, y, V=30.0, camber=0.05, alpha_deg=5):
    """
    Returns local horizontal velocity (u_top, u_bottom) over airfoil.
    """
    alpha = np.deg2rad(alpha_deg)
    U_inf = V * np.cos(alpha)
    influence = camber * np.exp(-((y - airfoil_camber(x)) ** 2) / 0.02)
    u_top = U_inf * (1 + influence)
    u_bottom = U_inf * (1 - influence)
    return u_top, u_bottom

# -------------------------------
# Bernoulli pressure
# -------------------------------
def pressure_from_velocity(v):
    """Compute pressure from velocity using Bernoulli"""
    return p0 - 0.5 * rho * v ** 2

# -------------------------------
# Lift calculation (manual trapezoid integration)
# -------------------------------
def compute_lift(x, p_top, p_bottom):
    y = p_bottom - p_top
    dx = np.diff(x)
    return np.sum((y[:-1] + y[1:]) / 2 * dx)

# -------------------------------
# Streamlit App
# -------------------------------
st.title("✈️ 2D Airfoil Lift and Takeoff Visualization")

# -------------------------------
# Plane selection
# -------------------------------
plane_options = {
    "Small Plane": 1500,   # Required lift in N/m
    "Medium Plane": 3000,
    "Large Plane": 5000
}
plane_choice = st.selectbox("Select a plane", list(plane_options.keys()))
L_required = plane_options[plane_choice]

# -------------------------------
# User sliders
# -------------------------------
V = st.slider("Freestream velocity (m/s)", 10, 100, 30)
camber = st.slider("Airfoil camber", 0.0, 0.2, 0.05)
alpha = st.slider("Angle of attack (deg)", -10, 20, 5)

# -------------------------------
# Discretize chord
# -------------------------------
x = np.linspace(0, 1, 400)
y_top = airfoil_camber(x)
y_bottom = -airfoil_camber(x)

# -------------------------------
# Velocity fields
# -------------------------------
u_top, u_bottom = velocity_field(x, y_top, V, camber, alpha)

# -------------------------------
# Pressure fields
# -------------------------------
p_top = pressure_from_velocity(u_top)
p_bottom = pressure_from_velocity(u_bottom)

# -------------------------------
# Lift calculation
# -------------------------------
lift = compute_lift(x, p_top, p_bottom)
st.write(f"Estimated Lift per unit span: {lift:.2f} N/m")

# -------------------------------
# Pressure distribution plot
# -------------------------------
fig1, ax1 = plt.subplots(figsize=(10, 4))
ax1.plot(x, p_top, label="Top Surface Pressure")
ax1.plot(x, p_bottom, label="Bottom Surface Pressure")
ax1.invert_yaxis()
ax1.set_xlabel("Chord Position")
ax1.set_ylabel("Pressure (Pa)")
ax1.set_title("Pressure Distribution Over Airfoil")
ax1.legend()
ax1.grid(True)
st.pyplot(fig1)

# -------------------------------
# Lift vs Airspeed plot
# -------------------------------
V_range = np.linspace(1, 100, 200)
lift_vals = []

for V_test in V_range:
    u_top, u_bottom = velocity_field(x, y_top, V_test, camber, alpha)
    p_top = pressure_from_velocity(u_top)
    p_bottom = pressure_from_velocity(u_bottom)
    lift_vals.append(compute_lift(x, p_top, p_bottom))

# Compute takeoff speed for selected plane
lift_vals = np.array(lift_vals)
V_takeoff_idx = np.where(lift_vals >= L_required)[0]
if len(V_takeoff_idx) > 0:
    V_takeoff = V_range[V_takeoff_idx[0]]
else:
    V_takeoff = np.nan

fig2, ax2 = plt.subplots(figsize=(10, 4))
ax2.plot(V_range, lift_vals, color="blue", label="Lift Curve")
if not np.isnan(V_takeoff):
    ax2.scatter(V_takeoff, L_required, color="red", s=100, zorder=5, label="Takeoff Point")
    ax2.axhline(L_required, color="red", linestyle="--", alpha=0.5)
    ax2.axvline(V_takeoff, color="red", linestyle="--", alpha=0.5)
    ax2.text(V_takeoff+1, L_required+50, f"V_takeoff={V_takeoff:.1f} m/s", color="red")
ax2.set_xlabel("Airflow Speed (m/s)")
ax2.set_ylabel("Lift per unit span (N/m)")
ax2.set_title(f"Lift vs Airspeed ({plane_choice})")
ax2.grid(True)
ax2.legend()
st.pyplot(fig2)

