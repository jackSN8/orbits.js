"""
Payload calculation engine for interplanetary missions.

This module implements the CRITICAL MISSING PIECE: reverse calculation from
required delta-V to maximum deliverable payload using binary search on the
Tsiolkovsky rocket equation.
"""
import numpy as np
from typing import Optional
import sys
import os

# Add parent directory to path for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data.rockets import Rocket, G0, EFFICIENCY_FACTOR, HEO_BONUS_MPS


def calculate_delta_v(rocket: Rocket, payload_kg: float) -> float:
    """
    Calculate achievable delta-V for a given payload mass.

    This is the FORWARD calculation from starship_science_graphs.ipynb:649ca5cb

    Args:
        rocket: Rocket specification
        payload_kg: Payload mass (kg)

    Returns:
        Achievable delta-V in m/s (0 if infeasible)
    """
    # Available propellant (can't exceed tank capacity or LEO payload budget)
    m_prop_available = min(
        rocket.leo_payload_kg - payload_kg - rocket.dry_mass_kg,
        rocket.wet_mass_kg - rocket.dry_mass_kg
    )

    if m_prop_available <= 0:
        return 0.0

    # Tsiolkovsky equation with efficiency factor
    m_wet = rocket.dry_mass_kg + payload_kg + m_prop_available
    m_dry = rocket.dry_mass_kg + payload_kg

    delta_v = EFFICIENCY_FACTOR * rocket.isp_s * G0 * np.log(m_wet / m_dry)

    # Add HEO bonus if applicable (from starship_science_graphs.ipynb)
    if rocket.starting_orbit == 'HEO':
        delta_v += HEO_BONUS_MPS

    return delta_v


def calculate_max_payload(rocket: Rocket, delta_v_required_mps: float,
                         tolerance_kg: float = 1.0) -> float:
    """
    Binary search to find maximum payload deliverable for required delta-V.

    This is the INVERSE calculation - the critical missing piece needed for
    the web application.

    Args:
        rocket: Rocket with wet_mass_kg, dry_mass_kg, isp_s, leo_payload_kg
        delta_v_required_mps: Total delta-V needed (m/s)
        tolerance_kg: Binary search tolerance (default 1 kg)

    Returns:
        Maximum payload mass in kg (0 if trajectory is infeasible)

    Algorithm:
        1. Binary search between 0 kg and max LEO payload
        2. For each payload guess, calculate achievable delta-V
        3. If achievable >= required, try higher payload (can carry more)
        4. If achievable < required, try lower payload (need to reduce)
        5. Converge within tolerance
    """
    # Handle edge cases
    if delta_v_required_mps <= 0:
        return rocket.leo_payload_kg  # No delta-V needed, full payload

    # Binary search bounds
    m_payload_min, m_payload_max = 0.0, rocket.leo_payload_kg

    # Check if trajectory is even feasible with zero payload
    max_dv = calculate_delta_v(rocket, 0.0)
    if max_dv < delta_v_required_mps:
        return 0.0  # Trajectory is infeasible even with no payload

    # Binary search
    while m_payload_max - m_payload_min > tolerance_kg:
        m_payload = (m_payload_min + m_payload_max) / 2

        # Calculate achievable delta-V with this payload
        delta_v_achievable = calculate_delta_v(rocket, m_payload)

        # Binary search update
        if delta_v_achievable >= delta_v_required_mps:
            m_payload_min = m_payload  # Can carry more
        else:
            m_payload_max = m_payload  # Need to reduce payload

    return (m_payload_min + m_payload_max) / 2


def calculate_payload_matrix(rocket: Rocket, dv_matrix: np.ndarray) -> np.ndarray:
    """
    Vectorized payload calculation for entire porkchop plot matrix.

    Args:
        rocket: Rocket specification
        dv_matrix: 2D array of delta-V values (m/s)

    Returns:
        2D array of payload masses (kg), same shape as dv_matrix
    """
    payload_matrix = np.zeros_like(dv_matrix, dtype=float)

    for i in range(dv_matrix.shape[0]):
        for j in range(dv_matrix.shape[1]):
            dv = dv_matrix[i, j]
            if dv > 0:  # Skip invalid entries
                payload_matrix[i, j] = calculate_max_payload(rocket, dv)

    return payload_matrix


def heliocentric_dv_to_departure_dv(v_infinity_mps: float) -> float:
    """
    Convert heliocentric delta-V (v-infinity) to departure delta-V from LEO.

    The porkchop data shows heliocentric delta-V (Sun reference frame), but
    we need the actual delta-V burn from LEO to achieve that trajectory.

    Args:
        v_infinity_mps: Hyperbolic excess velocity / heliocentric delta-V (m/s)

    Returns:
        Delta-V required from LEO (m/s)
    """
    v_leo = 7780  # m/s (circular velocity at ~200 km altitude)
    v_escape = 10930  # m/s (escape velocity at same altitude)

    # v_infinity = hyperbolic excess velocity
    # To achieve this, departure velocity must be:
    # v_departure^2 = v_infinity^2 + v_escape^2
    v_departure = np.sqrt(v_infinity_mps**2 + v_escape**2)

    # Delta-V from circular LEO
    delta_v = v_departure - v_leo

    return delta_v


def departure_dv_to_heliocentric_dv(delta_v_from_leo: float) -> float:
    """
    Convert departure delta-V from LEO to heliocentric delta-V (v-infinity).

    Inverse of heliocentric_dv_to_departure_dv.

    Args:
        delta_v_from_leo: Delta-V from LEO (m/s)

    Returns:
        Hyperbolic excess velocity (m/s)
    """
    v_leo = 7780  # m/s
    v_escape = 10930  # m/s

    v_departure = v_leo + delta_v_from_leo
    v_infinity = np.sqrt(v_departure**2 - v_escape**2)

    return v_infinity


def calculate_c3_from_delta_v(delta_v_mps: float) -> float:
    """
    Calculate characteristic energy C3 from delta-V from LEO.

    C3 is used to characterize interplanetary trajectories.
    From starship_science_graphs.ipynb:649ca5cb

    Args:
        delta_v_mps: Delta-V from LEO in m/s

    Returns:
        C3 in km^2/s^2
    """
    v_leo = 7780  # m/s (circular velocity at ~200 km)
    v_escape = 10930  # m/s (escape velocity at same altitude)

    c3_m2s2 = (v_leo + delta_v_mps)**2 - v_escape**2
    c3_km2s2 = c3_m2s2 / 1e6  # Convert to km^2/s^2

    return c3_km2s2


def calculate_c3_from_heliocentric_dv(v_infinity_mps: float) -> float:
    """
    Calculate C3 directly from heliocentric delta-V (v-infinity).

    C3 = v_infinity^2

    Args:
        v_infinity_mps: Hyperbolic excess velocity (m/s)

    Returns:
        C3 in km^2/s^2
    """
    c3_m2s2 = v_infinity_mps**2
    c3_km2s2 = c3_m2s2 / 1e6

    return c3_km2s2


def validate_round_trip(rocket: Rocket, delta_v_test: float,
                       tolerance_percent: float = 1.0) -> bool:
    """
    Validate payload calculation by round-trip test.

    Test: payload → delta-V → payload should recover within tolerance

    Args:
        rocket: Rocket specification
        delta_v_test: Test delta-V value (m/s)
        tolerance_percent: Acceptable error percentage (default 1%)

    Returns:
        True if round-trip is within tolerance
    """
    # Forward: delta-V → payload
    payload1 = calculate_max_payload(rocket, delta_v_test)

    # Backward: payload → delta-V
    dv_calc = calculate_delta_v(rocket, payload1)

    # Forward again: delta-V → payload
    payload2 = calculate_max_payload(rocket, dv_calc)

    # Check if we recovered the payload within tolerance
    error_percent = abs(payload1 - payload2) / payload1 * 100 if payload1 > 0 else 0

    return error_percent <= tolerance_percent


if __name__ == "__main__":
    # Quick validation tests
    from data.rockets import ROCKETS

    print("=== Payload Calculation Engine Validation ===\n")

    # Test 1: Starship refueled in LEO
    starship = ROCKETS['Starship refuelled in LEO']
    print(f"Rocket: {starship.name}")
    print(f"LEO Payload Capacity: {starship.leo_payload_kg/1000:.1f} tons\n")

    # Test different delta-V requirements
    test_dvs = [3000, 5000, 7000, 9000]  # m/s
    print("Delta-V Required (m/s) | Max Payload (tons) | C3 (km²/s²)")
    print("-" * 60)
    for dv in test_dvs:
        payload_kg = calculate_max_payload(starship, dv)
        c3 = calculate_c3_from_delta_v(dv)
        print(f"{dv:>22} | {payload_kg/1000:>18.1f} | {c3:>12.1f}")

    # Test 2: Round-trip validation
    print("\n=== Round-Trip Validation ===")
    for rocket_name in ['Atlas V 551', 'Falcon 9 expended', 'Starship refuelled in LEO']:
        rocket = ROCKETS[rocket_name]
        dv_test = 5000  # m/s
        is_valid = validate_round_trip(rocket, dv_test)
        status = "✓ PASS" if is_valid else "✗ FAIL"
        print(f"{rocket_name:30} {status}")

    # Test 3: Edge cases
    print("\n=== Edge Case Tests ===")

    # Zero delta-V (should return full payload)
    payload_zero_dv = calculate_max_payload(starship, 0.0)
    print(f"Zero delta-V payload: {payload_zero_dv/1000:.1f} tons (expected {starship.leo_payload_kg/1000:.1f})")

    # Infeasible trajectory (very high delta-V)
    payload_infeasible = calculate_max_payload(starship, 20000)
    print(f"Infeasible trajectory payload: {payload_infeasible:.1f} kg (expected 0)")

    # HEO vs LEO Starship
    starship_heo = ROCKETS['Starship refuelled in HEO']
    dv_compare = 7000
    payload_leo = calculate_max_payload(starship, dv_compare)
    payload_heo = calculate_max_payload(starship_heo, dv_compare)
    print(f"\nFor {dv_compare} m/s delta-V:")
    print(f"  LEO Starship: {payload_leo/1000:.1f} tons")
    print(f"  HEO Starship: {payload_heo/1000:.1f} tons")
    print(f"  HEO advantage: {(payload_heo - payload_leo)/1000:.1f} tons")
