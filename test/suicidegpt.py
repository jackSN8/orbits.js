import numpy as np
from math import sqrt, log

# =====================================================
# USER INPUT SECTION
# =====================================================

# Spacecraft parameters
m0     = 20000      # wet mass (kg)
mdry   = 12000      # dry mass (kg)
Isp    = 350        # seconds
thrust = 200000     # Newtons

# Target body parameters
mu      = 4.9048695e12   # GM of the Moon (m^3/s^2) example
R_body  = 1737.4e3        # radius of the Moon (m)

# Orbit parameters
h_initial = 100e3   # 100 km
h_final   = -10e3   # -10 km (i.e., periapsis 10 km below surface)
# =====================================================


# -----------------------------------------------------
# Helper functions
# -----------------------------------------------------

def vis_viva(mu, r, a):
    """Orbital velocity using the vis-viva equation."""
    return sqrt(mu * (2/r - 1/a))

def delta_v_from_mass(m0, mf, Isp):
    """Δv from Tsiolkovsky."""
    g0 = 9.80665
    return Isp * g0 * log(m0 / mf)

# -----------------------------------------------------
# Suicide burn solver
# -----------------------------------------------------

def suicide_burn_dv(mu, R_body, h_initial, h_final,
                    m0, mdry, Isp, thrust):

    r0 = R_body + h_initial
    rp = R_body + h_final
    ra = r0
    a  = 0.5*(rp + ra)

    # Velocity at pericenter
    v_p = vis_viva(mu, rp, a)

    # Local gravity at pericenter
    g_p = mu / (rp**2)

    # Effective deceleration
    def a_eff(m):
        return thrust/m - g_p

    # Find altitude where burn must begin
    altitudes = np.linspace(h_initial, h_final, 2000)
    m = m0

    for h in altitudes[::-1]:
        distance_left = h - h_final
        if a_eff(m) <= 0:
            continue
        if v_p**2 / (2 * a_eff(m)) < distance_left:
            h_start = h
            break

    # Δv ≈ pericenter velocity + gravity losses
    dv_total = v_p + 5.0  # crude 5 m/s gravity-loss correction

    return {
        "pericenter_velocity": v_p,
        "burn_start_altitude_m": h_start,
        "dv_required_m_s": dv_total
    }


# =====================================================
# RUN COMPUTATION
# =====================================================

results = suicide_burn_dv(
    mu, R_body,
    h_initial, h_final,
    m0, mdry, Isp, thrust
)

print("\n--- Suicide Burn Calculation ---")
print(f"Velocity at pericenter:        {results['pericenter_velocity']:.2f} m/s")
print(f"Start suicide burn at:         {results['burn_start_altitude_m']:.1f} m")
print(f"Required Δv:                    {results['dv_required_m_s']:.2f} m/s")
