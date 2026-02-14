# Interplanetary Mission Planner

Interactive web application for planning interplanetary missions using Lambert arc solutions.

## Features

- 🚀 **Rocket Database**: 10+ launch vehicles including Starship, Falcon Heavy, SLS, Atlas V
- 📊 **Interactive Porkchop Plots**: Visualize delta-V requirements across launch windows
- 📦 **Payload Calculations**: Automatically calculate deliverable payload for any trajectory
- 🔄 **Rocket Comparison**: Compare multiple vehicles for the same mission
- 🌍 **Multiple Destinations**: Pre-computed data for Mars, Jupiter, Saturn, Venus

## Installation

```bash
# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`.

## How It Works

### Payload Calculation Algorithm

The critical feature is calculating **maximum payload** given a required delta-V:

1. **Input**: Rocket specs (wet/dry mass, Isp) + Required delta-V
2. **Binary Search**: Find payload mass that achieves exactly the required delta-V
3. **Output**: Maximum deliverable payload in kg/tons

This is the reverse of the Tsiolkovsky equation, solved numerically.

### Data Sources

- **Rocket specs**: Extracted from [starship_science_graphs.ipynb](../starship_science_graphs.ipynb)
- **Porkchop data**: Pre-computed using [porkchop-data-generator.ipynb](../porkchop-data-generator.ipynb)
- **Orbital mechanics**: Poliastro library (Lambert solver, ephemeris)

## Project Structure

```
mission_planner/
├── app.py                      # Main Streamlit application
├── data/
│   ├── rockets.py              # Rocket database (10+ vehicles)
│   ├── porkchop_loader.py      # CSV data loading
│   └── precomputed/            # Pre-computed porkchop CSVs
├── calculations/
│   ├── payload_calc.py         # Payload calculation engine
│   ├── lambert.py              # Lambert solver wrapper (TODO)
│   └── trajectory.py           # Trajectory utilities (TODO)
├── ui/
│   ├── porkchop_plot.py        # Interactive plot component (TODO)
│   └── rocket_selector.py     # Rocket selection widget (TODO)
└── tests/                      # Unit tests (TODO)
```

## Usage Examples

### Select Mission Parameters

1. Choose rocket (e.g., "Starship refuelled in LEO")
2. Select destination (Mars, Jupiter, Saturn, Venus)
3. Pick launch window (e.g., 2030-2040)
4. Set time-of-flight range (e.g., 100-300 days)

### View Porkchop Plot

- **Color**: Delta-V requirement (red = low energy, purple/black = high energy)
- **Hover**: Shows launch date, TOF, delta-V, and **calculated payload**
- **Click**: View detailed trajectory information

### Compare Rockets

Select multiple rockets to see side-by-side payload capacity for the same trajectory.

## Technical Details

### Rocket Database

From [starship_science_graphs.ipynb](../starship_science_graphs.ipynb):

- **Efficiency Factor**: 0.87 (applied to delta-V calculations)
- **HEO Bonus**: +2000 m/s for Starship starting from high Earth orbit
- **LEO Payload**: Maximum mass deliverable to low Earth orbit

### Porkchop Plot Format

CSV files with:
- **Rows**: Launch dates (5-day spacing)
- **Columns**: Time-of-flight (5-day increments)
- **Values**: Delta-V in m/s (Sun reference frame)

## Development

### Run Tests

```bash
# Test payload calculation
python calculations/payload_calc.py

# Test data loader
python data/porkchop_loader.py
```

### Add New Rocket

Edit `data/rockets.py` and add to the `ROCKETS` dictionary:

```python
'My Rocket': Rocket(
    name='My Rocket',
    wet_mass_kg=100000,
    dry_mass_kg=10000,
    isp_s=350,
    leo_payload_kg=5000,
    starting_orbit='LEO'
)
```

## Future Enhancements

- [ ] Launch window optimizer (find best opportunities across decades)
- [ ] Gravity assist trajectory planning
- [ ] 3D trajectory visualization
- [ ] Cost optimization ($/kg)
- [ ] User accounts and mission saving
- [ ] PDF mission reports

## References

- **Poliastro**: https://docs.poliastro.space/
- **Lambert's Problem**: https://en.wikipedia.org/wiki/Lambert%27s_problem
- **Porkchop Plot**: https://en.wikipedia.org/wiki/Porkchop_plot

## License

MIT License - see parent project
