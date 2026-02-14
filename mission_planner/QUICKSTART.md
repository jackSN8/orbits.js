# Quick Start Guide

Get your Interplanetary Mission Planner running in 5 minutes!

## Prerequisites

- Python 3.9 or higher
- pip package manager

## Installation

### 1. Install Dependencies

```bash
cd mission_planner
pip install -r requirements.txt
```

This will install:
- Streamlit (web framework)
- Poliastro (orbital mechanics)
- Plotly (interactive plots)
- NumPy, Pandas, and other scientific computing libraries

### 2. Verify Installation

Test the payload calculation engine:

```bash
python calculations/payload_calc.py
```

You should see validation output showing payload calculations for different delta-V values.

### 3. Check Data Availability

Test the data loader:

```bash
python data/porkchop_loader.py
```

This will show which pre-computed porkchop datasets are available. Expected files:
- `../data/mars_dvs.csv` (if exists)
- `../jupiter_dvs_paper.csv` (if exists)
- `../saturn_dvs.csv` (if exists)
- `../venus_dvs.csv` (if exists)

**Note**: If these files don't exist, you'll need to generate them by running the parent directory's `porkchop-data-generator.ipynb` notebook first.

## Running the App

### Local Development

```bash
streamlit run app.py
```

The app will automatically open in your browser at `http://localhost:8501`.

### Docker Deployment

```bash
# Build the Docker image
docker build -t mission-planner .

# Run the container
docker run -p 8501:8501 mission-planner
```

Access at `http://localhost:8501`.

## Using the App

### Mission Planner Page

1. **Select Rocket**: Choose from 10+ launch vehicles (sidebar)
2. **Select Destination**: Mars, Jupiter, Saturn, or Venus
3. **Set Launch Window**: Pick date range (e.g., 2030-2040)
4. **Set Time of Flight**: Adjust TOF slider (e.g., 100-300 days)
5. **View Porkchop Plot**: Interactive heatmap shows delta-V requirements
6. **Hover for Details**: See launch date, TOF, delta-V, and **calculated payload**
7. **Find Optimal**: Check the Mission Statistics panel for the best trajectory

### Rocket Comparison Page

1. Navigate to "Rocket Comparison" in the sidebar
2. Select 2-6 rockets to compare
3. Choose comparison mode:
   - **Fixed ΔV**: See max payload each rocket can deliver
   - **Fixed Payload**: See max ΔV each rocket can achieve
4. View bar chart and detailed comparison table

## Example Missions

### Mars 2033 Opposition

- **Rocket**: Starship refuelled in LEO
- **Destination**: Mars
- **Launch Window**: 2032-2034
- **TOF**: 150-250 days
- **Expected**: ~200-300 tons payload with ~5,000 m/s delta-V

### Jupiter Hohmann Transfer

- **Rocket**: Starship refuelled in HEO
- **Destination**: Jupiter
- **Launch Window**: 2030-2040
- **TOF**: 600-1,000 days
- **Expected**: ~30-100 tons payload with ~8,000-9,000 m/s delta-V

## Troubleshooting

### "Pre-computed data not found"

Run the porkchop data generator from the parent directory:

```bash
cd ..
jupyter notebook porkchop-data-generator.ipynb
```

Generate data for your desired destinations (Mars, Jupiter, etc.) and save as CSV files.

### "ModuleNotFoundError"

Make sure you've installed all dependencies:

```bash
pip install -r requirements.txt
```

### App won't start

Check that port 8501 is not already in use:

```bash
lsof -i :8501
```

If occupied, you can specify a different port:

```bash
streamlit run app.py --server.port=8502
```

### Porkchop plot is empty

- Check that your date range overlaps with the pre-computed data
- Expand the TOF range slider
- Verify CSV files exist and contain data

## Next Steps

- **Customize Rockets**: Edit `data/rockets.py` to add your own launch vehicles
- **Generate More Data**: Run porkchop generator for more destinations or date ranges
- **Deploy to VPS**: See [README.md](README.md) for production deployment instructions
- **Explore Code**: Check out the payload calculation algorithm in `calculations/payload_calc.py`

## Support

For issues or questions:
1. Check the main [README.md](README.md)
2. Review the implementation plan at `~/.claude/plans/harmonic-forging-kite.md`
3. Open an issue in the parent GitHub repository

## Performance Tips

- Use date ranges of 1-5 years for best plot performance
- TOF ranges of 200-500 days load fastest
- First load may be slow as Streamlit builds cache
- Subsequent loads are much faster thanks to `@st.cache_data`

Enjoy planning your interplanetary missions! 🚀
