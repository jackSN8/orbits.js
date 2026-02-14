"""
Porkchop plot data loader for pre-computed trajectory data.

Loads CSV files from porkchop-data-generator.ipynb outputs.
"""
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from pathlib import Path
from typing import Tuple, Optional
import os


# Base data directory (relative to this file)
BASE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = BASE_DIR / "data"


# Pre-computed porkchop data file locations
PORKCHOP_FILES = {
    'Mars': DATA_DIR / 'mars_dvs.csv',
    'Jupiter': DATA_DIR / 'jupiter_dvs_paper.csv',
    'Saturn': DATA_DIR / 'saturn_dvs.csv',
    'Venus': DATA_DIR / 'venus_dvs.csv',
}


# Metadata for pre-computed datasets (from porkchop-data-generator.ipynb)
PORKCHOP_METADATA = {
    'Mars': {
        'start_date': datetime(2025, 2, 17),  # Approximate from notebooks
        'tej_res_days': 5,  # Launch date resolution
        'tof_res_days': 5,  # Time-of-flight resolution
        'max_tof_years': 5,
    },
    'Jupiter': {
        'start_date': datetime(2030, 1, 1),
        'tej_res_days': 5,
        'tof_res_days': 5,
        'max_tof_years': 3,
    },
    'Saturn': {
        'start_date': datetime(2030, 1, 1),
        'tej_res_days': 5,
        'tof_res_days': 5,
        'max_tof_years': 3,
    },
    'Venus': {
        'start_date': datetime(2025, 1, 1),
        'tej_res_days': 5,
        'tof_res_days': 5,
        'max_tof_years': 3,
    },
}


def load_porkchop_csv(destination: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load pre-computed porkchop data from CSV file.

    Args:
        destination: Planet name ('Mars', 'Jupiter', 'Saturn', 'Venus')

    Returns:
        Tuple of:
        - dv_matrix: 2D array of delta-V values (m/s), shape (num_launch_dates, num_tofs)
        - launch_dates: 1D array of launch dates (datetime objects)
        - tof_days: 1D array of time-of-flight values (days)

    Raises:
        FileNotFoundError: If porkchop file doesn't exist
        KeyError: If destination is not recognized
    """
    if destination not in PORKCHOP_FILES:
        available = ', '.join(PORKCHOP_FILES.keys())
        raise KeyError(f"Destination '{destination}' not recognized. Available: {available}")

    filepath = PORKCHOP_FILES[destination]

    if not filepath.exists():
        raise FileNotFoundError(
            f"Porkchop data file not found: {filepath}\n"
            f"Run porkchop-data-generator.ipynb to generate data."
        )

    # Load CSV (no headers, just numeric matrix)
    dv_matrix = np.loadtxt(filepath, delimiter=',')

    # Generate launch date array from metadata
    metadata = PORKCHOP_METADATA[destination]
    num_rows = dv_matrix.shape[0]
    launch_dates = np.array([
        metadata['start_date'] + timedelta(days=i * metadata['tej_res_days'])
        for i in range(num_rows)
    ])

    # Generate TOF array
    num_cols = dv_matrix.shape[1]
    tof_days = np.array([
        j * metadata['tof_res_days']
        for j in range(num_cols)
    ])

    return dv_matrix, launch_dates, tof_days


def filter_porkchop_data(
    dv_matrix: np.ndarray,
    launch_dates: np.ndarray,
    tof_days: np.ndarray,
    date_range: Optional[Tuple[datetime, datetime]] = None,
    tof_range_days: Optional[Tuple[int, int]] = None
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Filter porkchop data by date and TOF ranges.

    Args:
        dv_matrix: Full delta-V matrix
        launch_dates: Full launch date array
        tof_days: Full TOF array
        date_range: (start_date, end_date) tuple to filter launch dates
        tof_range_days: (min_tof, max_tof) tuple to filter TOF

    Returns:
        Filtered (dv_matrix, launch_dates, tof_days)
    """
    # Filter by launch date
    if date_range is not None:
        start_date, end_date = date_range
        date_mask = (launch_dates >= start_date) & (launch_dates <= end_date)
        dv_matrix = dv_matrix[date_mask, :]
        launch_dates = launch_dates[date_mask]

    # Filter by TOF
    if tof_range_days is not None:
        min_tof, max_tof = tof_range_days
        tof_mask = (tof_days >= min_tof) & (tof_days <= max_tof)
        dv_matrix = dv_matrix[:, tof_mask]
        tof_days = tof_days[tof_mask]

    return dv_matrix, launch_dates, tof_days


def get_delta_v_at_point(
    dv_matrix: np.ndarray,
    launch_dates: np.ndarray,
    tof_days: np.ndarray,
    launch_date: datetime,
    tof: int
) -> float:
    """
    Get delta-V at a specific trajectory point.

    Args:
        dv_matrix: Delta-V matrix
        launch_dates: Launch date array
        tof_days: TOF array
        launch_date: Desired launch date
        tof: Desired time-of-flight (days)

    Returns:
        Delta-V in m/s (interpolated if necessary)
    """
    # Find nearest launch date index
    date_idx = np.argmin(np.abs(launch_dates - launch_date))

    # Find nearest TOF index
    tof_idx = np.argmin(np.abs(tof_days - tof))

    return dv_matrix[date_idx, tof_idx]


def get_available_destinations() -> list[str]:
    """Return list of destinations with pre-computed data."""
    return list(PORKCHOP_FILES.keys())


def check_data_availability(destination: str) -> bool:
    """Check if porkchop data exists for destination."""
    if destination not in PORKCHOP_FILES:
        return False
    return PORKCHOP_FILES[destination].exists()


if __name__ == "__main__":
    # Quick test of data loading
    print("=== Porkchop Data Loader Test ===\n")

    for dest in get_available_destinations():
        print(f"Loading {dest} data...")

        if not check_data_availability(dest):
            print(f"  ✗ Data file not found: {PORKCHOP_FILES[dest]}")
            continue

        try:
            dv_matrix, launch_dates, tof_days = load_porkchop_csv(dest)
            print(f"  ✓ Loaded successfully")
            print(f"    Shape: {dv_matrix.shape} (launch dates × TOF)")
            print(f"    Launch date range: {launch_dates[0].date()} to {launch_dates[-1].date()}")
            print(f"    TOF range: {tof_days[0]:.0f} to {tof_days[-1]:.0f} days")
            print(f"    Min delta-V: {np.min(dv_matrix[dv_matrix > 0]):.0f} m/s")
            print(f"    Max delta-V: {np.max(dv_matrix):.0f} m/s")

            # Test filtering
            date_range = (datetime(2030, 1, 1), datetime(2035, 1, 1))
            tof_range = (100, 300)
            filtered_dv, filtered_dates, filtered_tof = filter_porkchop_data(
                dv_matrix, launch_dates, tof_days,
                date_range=date_range,
                tof_range_days=tof_range
            )
            print(f"    Filtered shape (2030-2035, 100-300 days): {filtered_dv.shape}")

        except Exception as e:
            print(f"  ✗ Error loading: {e}")

        print()
