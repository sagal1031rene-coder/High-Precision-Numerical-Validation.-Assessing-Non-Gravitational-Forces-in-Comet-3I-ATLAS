# ==============================================================================
# 3I/ATLAS HIGH-PRECISION ORBITAL SIMULATION CODE (RK45)
# ==============================================================================
# Built for Google Colab (2026) to validate the 3iatlas orbital platform.
# This code integrates the perturbed equations of motion for the hyperbolic 
# comet 3I/ATLAS, incorporating solar gravitation, solar radiation pressure,
# and asymmetric outgassing with an 8-degree jet misalignment.
# ==============================================================================

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# 1. PHYSICAL AND ORBITAL CONSTANTS
# ------------------------------------------------------------------------------
G_SUN = 1.32712440018e11      # Parameter G*M for the Sun (km^3/s^2)
AU_TO_KM = 149597870.7        # Astronomical Unit in km
C_LIGHT = 299792.458          # Speed of light in km/s

# ------------------------------------------------------------------------------
# 2. 3I/ATLAS PHYSICAL PARAMETERS (CONSENSUAL VALUES)
# ------------------------------------------------------------------------------
RADIUS = 1300.0               # Radius of nucleus (meters, 1.3 km)
DENSITY = 550.0               # Bulk density (kg/m^3)
VOLUME = (4.0 / 3.0) * np.pi * (RADIUS ** 3) # Volume in m^3 (~9.20e9 m^3)
MASS_INITIAL = DENSITY * VOLUME  # Initial mass in kg (~5.06e12 kg)
AREA = np.pi * (RADIUS ** 2)  # Cross-sectional area in m^2 (~5.31e6 m^2)
SOLAR_CONSTANT = 1361.0       # Solar constant at 1 AU (W/m^2)

# ------------------------------------------------------------------------------
# 3. EQUATIONS OF MOTION WITH PERTURBATIONS
# ------------------------------------------------------------------------------
def equations_of_motion(t, y, outgassing_active=True):
    """
    Computes the derivatives of the state vector y = [x, y, z, vx, vy, vz].
    Coordinates are in km and km/s, centered on the Sun.
    """
    r_vec = y[:3]
    v_vec = y[3:]
    r_mag = np.linalg.norm(r_vec)
    r_unit = r_vec / r_mag
    
    # 1. Gravitational Acceleration (Monolithic Sun)
    a_grav = -G_SUN * r_vec / (r_mag ** 3)
    
    # 2. Solar Radiation Pressure Acceleration
    r_au = r_mag / AU_TO_KM
    # Force in Newtons: F = (S_0 / c) * A / d_au^2
    p_rad_force = (SOLAR_CONSTANT / (C_LIGHT * 1000.0)) * (AREA / (r_au ** 2))
    # Convert acceleration to km/s^2 (F_N / m_kg * 1e-3)
    a_rad_mag = (p_rad_force / MASS_INITIAL) * 1e-3
    a_rad = a_rad_mag * r_unit
    
    # 3. Asymmetric Outgassing Thrust Acceleration
    # Activated when inside 2.5 AU, with an 8-degree jet misalignment
    if outgassing_active and r_au < 2.5:
        # Thrust magnitude scales directly with solar flux intensity (20% of solar radiation force)
        thrust_magnitude = a_rad_mag * 0.2
        
        # Rotation angle for the jet misalignment (8 degrees)
        theta = np.radians(8)
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        
        # Outgassing jet operates in the orbital plane (xy rotation of r_unit vector)
        a_thrust_x = thrust_magnitude * (cos_t * r_unit[0] - sin_t * r_unit[1])
        a_thrust_y = thrust_magnitude * (sin_t * r_unit[0] + cos_t * r_unit[1])
        a_thrust = np.array([a_thrust_x, a_thrust_y, 0.0])
    else:
        a_thrust = np.zeros(3)
        
    # Total acceleration
    dvdt = a_grav + a_rad + a_thrust
    return np.concatenate([v_vec, dvdt])

# ------------------------------------------------------------------------------
# 4. SIMULATION EXECUTION AND PLOTTING
# ------------------------------------------------------------------------------
def run_simulation():
    # Initial conditions based on NASA hyperbolic excess speed of 64.04 km/s
    # State vector: [x, y, z, vx, vy, vz] in km and km/s
    y0 = [1.5e8, 0.0, 0.0, -20.0, 60.8, 0.0]
    t_span = (0, 1.2e7)  # ~138 days simulation span
    
    print("Integrating equations of motion using adaptive RK45 solver (rtol=1e-9)...")
    sol = solve_ivp(equations_of_motion, t_span, y0, method='RK45', rtol=1e-9, atol=1e-9)
    print("Simulation completed successfully.")
    
    # Extract data
    times = sol.t
    positions = sol.y[:3, :].T
    velocities = sol.y[3:, :].T
    distances = np.linalg.norm(positions, axis=1)
    
    # Minimum heliocentric distance (Perihelion)
    min_dist_idx = np.argmin(distances)
    min_dist = distances[min_dist_idx]
    
    # Maximum Non-Gravitational Acceleration calculation
    a_non_grav_max = 0.0
    for r in distances:
        r_au = r / AU_TO_KM
        if r_au < 2.5:
            # Solar radiation pressure acceleration at distance r
            p_rad_force = (SOLAR_CONSTANT / (C_LIGHT * 1000.0)) * (AREA / (r_au ** 2))
            a_rad_mag = (p_rad_force / MASS_INITIAL) * 1e-3
            # Combined acceleration of radiation and 20% outgassing thrust
            a_total = a_rad_mag * 1.2
            if a_total > a_non_grav_max:
                a_non_grav_max = a_total
                
    print(f"Perihelion (Minimum heliocentric distance): {min_dist:,.2f} km ({min_dist/AU_TO_KM:.4f} AU)")
    print(f"Maximum Non-Gravitational Acceleration: {a_non_grav_max:.2e} km/s^2")
    
    # Generate Plots
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
    
    # Plot 1: Heliocentric Trajectory
    plt.figure(figsize=(7, 7))
    plt.plot(positions[:, 0] / AU_TO_KM, positions[:, 1] / AU_TO_KM, 'b-', label='3I/ATLAS Trajectory', linewidth=2)
    plt.plot(0, 0, 'yo', markersize=12, label='Sun')
    # Draw reference planetary orbits
    theta_circle = np.linspace(0, 2*np.pi, 200)
    plt.plot(np.cos(theta_circle), np.sin(theta_circle), 'g--', label='Earth Orbit (1.0 AU)', alpha=0.6)
    plt.plot(1.524*np.cos(theta_circle), 1.524*np.sin(theta_circle), 'r--', label='Mars Orbit (1.52 AU)', alpha=0.6)
    plt.xlim(-1.8, 1.8)
    plt.ylim(-1.8, 1.8)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('X (AU)', fontsize=12, fontweight='bold')
    plt.ylabel('Y (AU)', fontsize=12, fontweight='bold')
    plt.title('Heliocentric Orbit of 3I/ATLAS', fontsize=14, fontweight='bold', pad=15)
    plt.legend(frameon=True, facecolor='white', framealpha=0.9, edgecolor='gray', loc='upper right')
    plt.tight_layout()
    plt.savefig('trajectory_3iatlas.png', dpi=300)
    plt.close()
    
    # Plot 2: Velocity Magnitude vs Time
    v_mags = np.linalg.norm(velocities, axis=1)
    plt.figure(figsize=(8, 4))
    plt.plot(times / 86400, v_mags, 'darkorange', linewidth=2, label='Orbital Speed')
    plt.axhline(64.04, color='red', linestyle=':', label='Hyperbolic Excess Speed ($v_\infty = 64.04$ km/s)', linewidth=1.5)
    plt.xlabel('Time (Days)', fontsize=11, fontweight='bold')
    plt.ylabel('Velocity (km/s)', fontsize=11, fontweight='bold')
    plt.title('Magnitude of Velocity vs. Time', fontsize=13, fontweight='bold')
    plt.legend(frameon=True, facecolor='white', framealpha=0.9)
    plt.tight_layout()
    plt.savefig('velocity_3iatlas.png', dpi=300)
    plt.close()
    
    # Plot 3: Specific Orbital Energy vs Time
    energies = 0.5 * (v_mags**2) - (G_SUN / distances)
    plt.figure(figsize=(8, 4))
    plt.plot(times / 86400, energies, 'purple', linewidth=2, label='Specific Energy ($\epsilon$)')
    plt.xlabel('Time (Days)', fontsize=11, fontweight='bold')
    plt.ylabel('Specific Energy (km^2/s^2)', fontsize=11, fontweight='bold')
    plt.title('Total Specific Orbital Energy vs. Time', fontsize=13, fontweight='bold')
    plt.legend(frameon=True, facecolor='white', framealpha=0.9)
    plt.tight_layout()
    plt.savefig('energy_3iatlas.png', dpi=300)
    plt.close()
    
    print("All three figures saved successfully.")

if __name__ == "__main__":
    run_simulation()
