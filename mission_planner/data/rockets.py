"""
Rocket database for interplanetary mission planning.

Extracted from starship_science_graphs.ipynb cell 'da7de3c0'.
Contains specifications for various launch vehicles.
"""
from dataclasses import dataclass
from typing import Optional


@dataclass
class Rocket:
    """
    Rocket vehicle specification.

    Attributes:
        name: Vehicle name
        wet_mass_kg: Fully fueled mass including propellant (kg)
        dry_mass_kg: Empty mass of upper stage without propellant or payload (kg)
        isp_s: Specific impulse in vacuum (seconds)
        leo_payload_kg: Maximum payload capacity to LEO (kg)
        starting_orbit: 'LEO' or 'HEO' - affects delta-V calculations
    """
    name: str
    wet_mass_kg: float
    dry_mass_kg: float
    isp_s: float
    leo_payload_kg: float
    starting_orbit: str = 'LEO'  # 'LEO' or 'HEO'

    def __str__(self) -> str:
        return f"{self.name} (LEO payload: {self.leo_payload_kg/1000:.1f} tons, Isp: {self.isp_s:.1f}s)"


# Rocket database extracted from starship_science_graphs.ipynb
# performance dictionary cell 'da7de3c0'
# [wet_mass_kg, dry_mass_kg, Isp_s, LEO_payload_kg]

ROCKETS = {
    'Atlas V 551': Rocket(
        name='Atlas V 551',
        wet_mass_kg=23077,
        dry_mass_kg=2247,
        isp_s=449.7,
        leo_payload_kg=18850,
        starting_orbit='LEO'
    ),
    'Falcon 9 RTLS': Rocket(
        name='Falcon 9 RTLS',
        wet_mass_kg=111500,
        dry_mass_kg=4200,
        isp_s=348,
        leo_payload_kg=15000,
        starting_orbit='LEO'
    ),
    'Falcon 9 expended': Rocket(
        name='Falcon 9 expended',
        wet_mass_kg=111500,
        dry_mass_kg=4200,
        isp_s=348,
        leo_payload_kg=22800,
        starting_orbit='LEO'
    ),
    'Falcon Heavy expended': Rocket(
        name='Falcon Heavy expended',
        wet_mass_kg=111500,
        dry_mass_kg=4200,
        isp_s=348,
        leo_payload_kg=63800,
        starting_orbit='LEO'
    ),
    'Delta IV Heavy': Rocket(
        name='Delta IV Heavy',
        wet_mass_kg=30700,
        dry_mass_kg=3480,
        isp_s=465.5,
        leo_payload_kg=28790,
        starting_orbit='LEO'
    ),
    'SLS Block 1B': Rocket(
        name='SLS Block 1B',
        wet_mass_kg=143000,
        dry_mass_kg=14000,
        isp_s=460,
        leo_payload_kg=105000,
        starting_orbit='LEO'
    ),
    'Saturn V': Rocket(
        name='Saturn V',
        wet_mass_kg=123000,
        dry_mass_kg=13500,
        isp_s=420,
        leo_payload_kg=145000,
        starting_orbit='LEO'
    ),
    'Unrefuelled Starship': Rocket(
        name='Unrefuelled Starship',
        wet_mass_kg=1300000,
        dry_mass_kg=120000,
        isp_s=380,
        leo_payload_kg=250000,
        starting_orbit='LEO'
    ),
    'Starship refuelled in LEO': Rocket(
        name='Starship refuelled in LEO',
        wet_mass_kg=1300000,
        dry_mass_kg=120000,
        isp_s=380,
        leo_payload_kg=1450000,
        starting_orbit='LEO'
    ),
    'Starship refuelled in HEO': Rocket(
        name='Starship refuelled in HEO',
        wet_mass_kg=1300000,
        dry_mass_kg=120000,
        isp_s=380,
        leo_payload_kg=1450000,
        starting_orbit='HEO'
    ),
    'Light Starship refuelled in HEO': Rocket(
        name='Light Starship refuelled in HEO',
        wet_mass_kg=1300000,
        dry_mass_kg=80000,
        isp_s=380,
        leo_payload_kg=1450000,
        starting_orbit='HEO'
    ),
}


def get_all_rockets() -> list[Rocket]:
    """Return list of all available rockets."""
    return list(ROCKETS.values())


def get_rocket(name: str) -> Optional[Rocket]:
    """Get rocket by name."""
    return ROCKETS.get(name)


def get_rocket_names() -> list[str]:
    """Return list of all rocket names."""
    return list(ROCKETS.keys())


# Physical constants
G0 = 9.80665  # Standard gravity (m/s^2)
EFFICIENCY_FACTOR = 0.87  # Delta-V efficiency factor from starship_science_graphs.ipynb
HEO_BONUS_MPS = 2000  # Extra delta-V for HEO-starting Starship (m/s)
