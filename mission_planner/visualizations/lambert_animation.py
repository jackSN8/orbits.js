"""
Lambert's Problem Visualization - Animated Transfer Orbits

Demonstrates how different transfer times produce different orbital arcs
connecting two points in space (Lambert's problem).

Based on the geometric construction from Wikipedia's Lambert's problem article.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle, FancyArrowPatch, Wedge
import sys
import os


class LambertVisualizer:
    """Visualize Lambert arc solutions with proper geometric construction."""

    def __init__(self, r1_km, r2_km, mu=398600.4418):
        """
        Initialize Lambert problem visualization.

        Args:
            r1_km: Departure position vector [x, y, z] (km)
            r2_km: Arrival position vector [x, y, z] (km)
            mu: Gravitational parameter (km^3/s^2), default is Earth
        """
        self.r1 = np.array(r1_km[:2])  # Use only x,y for 2D
        self.r2 = np.array(r2_km[:2])
        self.mu = mu

        # Calculate geometric parameters
        self.r1_mag = np.linalg.norm(self.r1)
        self.r2_mag = np.linalg.norm(self.r2)
        self.c = np.linalg.norm(self.r2 - self.r1)  # Chord length
        self.s = (self.r1_mag + self.r2_mag + self.c) / 2  # Semi-perimeter

        # Calculate transfer angle
        cos_dtheta = np.dot(self.r1, self.r2) / (self.r1_mag * self.r2_mag)
        self.dtheta = np.arccos(np.clip(cos_dtheta, -1, 1))

        # Determine direction (short way vs long way)
        cross_z = self.r1[0] * self.r2[1] - self.r1[1] * self.r2[0]
        self.short_way = True  # Use short way for visualization

        if cross_z < 0:
            self.dtheta = 2 * np.pi - self.dtheta

    def parabolic_tof(self):
        """Calculate time of flight for parabolic (minimum energy) transfer."""
        # For parabolic orbit, a → infinity
        # This gives the minimum time for the given geometry
        a_min = self.s / 2
        tof = np.sqrt(2 * a_min**3 / self.mu) * (1 - (self.s - self.c) / (self.s))
        return tof

    def calculate_lambert_arc(self, tof_seconds):
        """
        Calculate transfer orbit using Universal Variable formulation.

        This follows the geometric construction more closely.
        """
        # Use Battin's method with universal variables
        # Simplified for visualization

        # Initial guess for semi-major axis
        a = self.s / 2  # Minimum energy orbit

        # Iterate to find correct a for given tof
        for iteration in range(100):
            # Calculate alpha and beta parameters
            sin_half_sum = np.sqrt(self.s / (2 * a))
            sin_half_diff = np.sqrt((self.s - self.c) / (2 * a))

            # Clamp to valid range
            sin_half_sum = np.clip(sin_half_sum, 0, 1)
            sin_half_diff = np.clip(sin_half_diff, 0, 1)

            alpha = 2 * np.arcsin(sin_half_sum)
            beta = 2 * np.arcsin(sin_half_diff)

            # Time of flight for this semi-major axis
            tof_calc = np.sqrt(a**3 / self.mu) * (
                (alpha - np.sin(alpha)) - (beta - np.sin(beta))
            )

            error = tof_calc - tof_seconds

            if abs(error) < 0.1:  # Converged
                break

            # Newton-Raphson update
            dtof_da = (3/2) * np.sqrt(a / self.mu) * (
                (alpha - np.sin(alpha)) - (beta - np.sin(beta))
            )

            if abs(dtof_da) > 1e-10:
                a = a - error / dtof_da

            # Keep a positive
            if a <= 0:
                a = self.s / 4

        # Calculate eccentricity
        p = (4 * a * (self.s - self.r1_mag) * (self.s - self.r2_mag) / self.c**2) * np.sin((alpha + beta)/2)**2
        e = np.sqrt(max(0, 1 - p / a))

        # Generate the actual orbit points that pass through r1 and r2
        orbit_points = self.generate_orbit_through_points(a, e, alpha, beta)

        # Calculate velocities using vis-viva
        v1_mag = np.sqrt(self.mu * (2/self.r1_mag - 1/a))
        v2_mag = np.sqrt(self.mu * (2/self.r2_mag - 1/a))

        return {
            'a': a,
            'e': e,
            'orbit_points': orbit_points,
            'v1': v1_mag,
            'v2': v2_mag,
            'alpha': alpha,
            'beta': beta,
            'tof': tof_calc
        }

    def generate_orbit_through_points(self, a, e, alpha, beta):
        """
        Generate orbit points that actually pass through r1 and r2.

        Uses proper orbital geometry with focus at origin.
        """
        # The orbit has one focus at origin (central body)
        # We need to orient it so it passes through r1 and r2

        # Calculate true anomalies at r1 and r2
        # Using the relationship: r = a(1-e²)/(1+e*cos(θ))

        # True anomaly at r1
        cos_nu1 = (a * (1 - e**2) / self.r1_mag - 1) / e if e > 1e-6 else 0
        cos_nu1 = np.clip(cos_nu1, -1, 1)
        nu1 = np.arccos(cos_nu1)

        # True anomaly at r2
        cos_nu2 = (a * (1 - e**2) / self.r2_mag - 1) / e if e > 1e-6 else 0
        cos_nu2 = np.clip(cos_nu2, -1, 1)
        nu2 = np.arccos(cos_nu2)

        # Adjust nu2 to be ahead of nu1
        if self.short_way:
            if nu2 < nu1:
                nu2 = 2 * np.pi - nu2

        # Generate true anomaly values from nu1 to nu2
        n_points = 100
        if self.short_way:
            nu_vals = np.linspace(nu1, nu2, n_points)
        else:
            nu_vals = np.linspace(nu1, nu1 + 2*np.pi - (nu2 - nu1), n_points)

        # Calculate r(θ) for each true anomaly
        r_vals = a * (1 - e**2) / (1 + e * np.cos(nu_vals))

        # Convert to Cartesian coordinates in orbital plane
        x_orbit = r_vals * np.cos(nu_vals)
        y_orbit = r_vals * np.sin(nu_vals)

        # Rotate to match actual r1 orientation
        theta_r1 = np.arctan2(self.r1[1], self.r1[0])
        rotation_angle = theta_r1 - nu1

        cos_rot = np.cos(rotation_angle)
        sin_rot = np.sin(rotation_angle)

        x_final = x_orbit * cos_rot - y_orbit * sin_rot
        y_final = x_orbit * sin_rot + y_orbit * cos_rot

        return np.column_stack([x_final, y_final])


def create_lambert_animation(output_file='lambert_arc_animation.gif', fps=20, duration=10):
    """Create animated visualization of Lambert's problem."""

    # Setup: Two points on circular orbits around Earth
    # Use different altitudes and a non-degenerate transfer angle
    r_earth = 6378  # km
    h1 = 300  # Lower orbit altitude (km)
    h2 = 800  # Higher orbit altitude (km)

    # Departure at 0 degrees, arrival at 120 degrees (avoids degenerate 180° case)
    angle1 = 0
    angle2 = 120 * np.pi / 180

    r1 = (r_earth + h1) * np.array([np.cos(angle1), np.sin(angle1)])
    r2 = (r_earth + h2) * np.array([np.cos(angle2), np.sin(angle2)])

    visualizer = LambertVisualizer(r1, r2)

    # Time range: from 0.3x to 3x the parabolic time
    tof_parabolic = visualizer.parabolic_tof()
    tof_range = np.linspace(tof_parabolic * 0.3, tof_parabolic * 2.5, fps * duration)

    # Setup figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle('Lambert\'s Problem: Transfer Orbit Solutions', fontsize=16, fontweight='bold')

    # Left plot: Orbital geometry
    ax1.set_aspect('equal')
    ax1.grid(True, alpha=0.3, linestyle='--')
    ax1.set_xlabel('X (km)', fontsize=12)
    ax1.set_ylabel('Y (km)', fontsize=12)
    ax1.set_title('Transfer Orbits Connecting Two Points')

    # Draw Earth
    earth = Circle((0, 0), r_earth, color='dodgerblue', alpha=0.4, label='Earth', zorder=1)
    ax1.add_patch(earth)

    # Draw initial circular orbits
    theta_circ = np.linspace(0, 2*np.pi, 100)
    orbit1_x = (r_earth + h1) * np.cos(theta_circ)
    orbit1_y = (r_earth + h1) * np.sin(theta_circ)
    orbit2_x = (r_earth + h2) * np.cos(theta_circ)
    orbit2_y = (r_earth + h2) * np.sin(theta_circ)
    ax1.plot(orbit1_x, orbit1_y, 'g--', alpha=0.3, linewidth=1, label='Departure orbit')
    ax1.plot(orbit2_x, orbit2_y, 'r--', alpha=0.3, linewidth=1, label='Arrival orbit')

    # Mark departure and arrival
    ax1.plot(*r1, 'go', markersize=14, label='Departure', zorder=5, markeredgecolor='darkgreen', markeredgewidth=2)
    ax1.plot(*r2, 'ro', markersize=14, label='Arrival', zorder=5, markeredgecolor='darkred', markeredgewidth=2)

    # Draw chord
    ax1.plot([r1[0], r2[0]], [r1[1], r2[1]], 'gray', linestyle=':',
             linewidth=2, alpha=0.5, label='Chord')

    # Transfer orbit line
    orbit_line, = ax1.plot([], [], 'r-', linewidth=3, label='Transfer Orbit', zorder=3)

    # Velocity arrows (will be updated)
    v1_arrow = None
    v2_arrow = None

    # Right plot: Parameters
    ax2.set_xlabel('Time of Flight (hours)', fontsize=12)
    ax2.set_ylabel('Semi-major Axis (km)', fontsize=12)
    ax2.set_title('Orbit Energy vs Transfer Time')
    ax2.grid(True, alpha=0.3)

    tof_history = []
    a_history = []

    param_line, = ax2.plot([], [], 'b-', linewidth=2.5, label='Semi-major axis')
    current_point, = ax2.plot([], [], 'ro', markersize=12, markeredgecolor='darkred',
                              markeredgewidth=2, label='Current', zorder=5)

    # Add reference line for parabolic orbit
    ax2.axhline(y=visualizer.s/2, color='green', linestyle='--',
                alpha=0.5, label='Minimum energy')

    # Info box
    info_text = ax1.text(0.02, 0.98, '', transform=ax1.transAxes,
                         verticalalignment='top', fontsize=10,
                         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9),
                         family='monospace')

    ax1.legend(loc='upper right', fontsize=9, framealpha=0.9)
    ax2.legend(loc='upper left', fontsize=10, framealpha=0.9)

    # Set limits
    max_coord = max(np.linalg.norm(r1), np.linalg.norm(r2)) * 1.6
    ax1.set_xlim(-max_coord, max_coord)
    ax1.set_ylim(-max_coord, max_coord)

    def init():
        orbit_line.set_data([], [])
        param_line.set_data([], [])
        current_point.set_data([], [])
        info_text.set_text('')
        return orbit_line, param_line, current_point, info_text

    def animate(frame):
        nonlocal v1_arrow, v2_arrow

        tof = tof_range[frame]

        try:
            result = visualizer.calculate_lambert_arc(tof)

            # Update orbit
            orbit_points = result['orbit_points']
            orbit_line.set_data(orbit_points[:, 0], orbit_points[:, 1])

            # Remove old arrows
            if v1_arrow is not None:
                v1_arrow.remove()
                v2_arrow.remove()

            # Calculate velocity directions
            v_scale = 800
            v1_dir = orbit_points[1] - orbit_points[0]
            v1_dir = v1_dir / np.linalg.norm(v1_dir) * v_scale

            v2_dir = orbit_points[-1] - orbit_points[-2]
            v2_dir = v2_dir / np.linalg.norm(v2_dir) * v_scale

            # Draw velocity vectors
            v1_arrow = FancyArrowPatch(r1, r1 + v1_dir,
                                      arrowstyle='->', mutation_scale=25,
                                      linewidth=2.5, color='green', zorder=4)
            v2_arrow = FancyArrowPatch(r2, r2 + v2_dir,
                                      arrowstyle='->', mutation_scale=25,
                                      linewidth=2.5, color='red', zorder=4)
            ax1.add_patch(v1_arrow)
            ax1.add_patch(v2_arrow)

            # Update info
            orbit_type = 'Ellipse' if result['e'] < 1 else 'Hyperbola'
            energy = -visualizer.mu / (2 * result['a']) if result['a'] > 0 else float('inf')

            info_text.set_text(
                f"Transfer Time: {tof/3600:.2f} hours\n"
                f"Semi-major axis: {result['a']:.0f} km\n"
                f"Eccentricity: {result['e']:.4f}\n"
                f"Orbit: {orbit_type}\n"
                f"v₁: {result['v1']:.1f} m/s\n"
                f"v₂: {result['v2']:.1f} m/s"
            )

            # Update parameter history
            tof_history.append(tof / 3600)
            a_history.append(result['a'])

            param_line.set_data(tof_history, a_history)
            current_point.set_data([tof / 3600], [result['a']])

            # Auto-scale
            if len(tof_history) > 1:
                ax2.set_xlim(min(tof_history) * 0.95, max(tof_history) * 1.05)
                ax2.set_ylim(min(a_history) * 0.8, max(a_history) * 1.2)

        except Exception as e:
            print(f"Frame {frame} error: {e}")

        return orbit_line, param_line, current_point, info_text

    # Create animation
    print(f"Creating {len(tof_range)} frame animation...")
    anim = FuncAnimation(fig, animate, init_func=init,
                        frames=len(tof_range), interval=1000/fps,
                        blit=False, repeat=True)

    # Save
    writer = PillowWriter(fps=fps)
    anim.save(output_file, writer=writer)
    print(f"✓ Animation saved: {output_file}")

    plt.close()


if __name__ == "__main__":
    print("=" * 70)
    print("Lambert's Problem: Geometric Transfer Orbit Visualization")
    print("=" * 70)
    print("\nShows how transfer orbits connect two fixed points in space")
    print("as a function of transfer time.\n")

    create_lambert_animation('lambert_arc_animation.gif', fps=20, duration=10)

    print("\n✓ Complete! The animation demonstrates:")
    print("  - How orbit shape changes with transfer time")
    print("  - Departure and arrival velocity vectors")
    print("  - Semi-major axis evolution")
    print("\nThis is the foundation of porkchop plot calculations!")
