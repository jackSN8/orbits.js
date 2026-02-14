# Lambert's Problem Visualization

An animated demonstration of how different transfer orbits connect two points in space.

## What is Lambert's Problem?

Lambert's problem is fundamental to orbital mechanics: **Given two positions and a transfer time, find the orbit that connects them.**

This is used for:
- **Interplanetary mission planning** (your porkchop plots!)
- **Rendezvous calculations** (docking spacecraft)
- **Orbital intercepts** (satellite servicing)

## The Mathematics

For two position vectors **r₁** and **r₂** separated by time **t**, Lambert's problem finds the unique conic section (ellipse, parabola, or hyperbola) connecting them.

### Key Parameters:
- **c** = |r₂ - r₁| (chord length)
- **s** = (r₁ + r₂ + c) / 2 (semi-perimeter)
- **a** = semi-major axis (determined by transfer time)

### The Solution:
The transfer time **t** is related to the geometry by:
```
t = √(a³/μ) × [α - β - (sin α - sin β)]
```

Where:
- **α** and **β** are geometric angles derived from the semi-perimeter
- **μ** is the gravitational parameter

## Running the Animation

```bash
cd lambert_animation
python3 lambert_animation.py
```

This creates `lambert_arc_animation.gif` showing:

### Left Panel: Orbital Geometry
- **Blue circle**: Earth (or central body)
- **Green dot**: Departure point
- **Red dot**: Arrival point
- **Red curve**: Transfer orbit
- **Green arrow**: Departure velocity vector
- **Red arrow**: Arrival velocity vector

### Right Panel: Parameter Evolution
- Shows how semi-major axis changes with transfer time
- **Red dot**: Current transfer time

## What the Animation Shows

As the animation progresses through different transfer times:

1. **Fast transfers** (short time):
   - High energy orbits
   - Large semi-major axis
   - May be hyperbolic (e > 1)
   - High ΔV required

2. **Minimum energy transfer** (Hohmann-like):
   - Optimal balance of time and energy
   - Smallest semi-major axis
   - Elliptical orbit
   - Minimum ΔV

3. **Slow transfers** (long time):
   - Lower energy but longer duration
   - Larger semi-major axis again
   - Very elongated ellipse
   - Moderate ΔV

## Connection to Your Porkchop Plots

Each point on your porkchop plot represents a different Lambert arc solution:
- **X-axis (launch date)**: Determines r₁ (Earth position) and r₂ (Mars position)
- **Y-axis (TOF)**: The transfer time **t**
- **Color (payload)**: Result of the ΔV required by the Lambert solution

The **red regions** on porkchop plots correspond to low-energy Lambert arcs (near minimum energy transfer), while **black regions** are high-energy arcs requiring more ΔV than the rocket can provide.

## Customization

Edit the script to explore different scenarios:

```python
# Change departure and arrival positions
r1 = np.array([7000, 0])      # Departure (km)
r2 = np.array([-8000, 0])     # Arrival (km)

# Adjust animation parameters
create_lambert_animation(
    output_file='my_animation.gif',
    fps=30,        # Frames per second
    duration=15    # Total duration (seconds)
)
```

## Technical Details

### Iteration Method
The script uses Newton-Raphson iteration to solve for the semi-major axis **a** that produces the desired transfer time **t**.

### Geometric Construction
1. Calculate chord **c** and semi-perimeter **s**
2. Iterate to find **a** matching the transfer time
3. Calculate eccentricity **e** from geometry
4. Generate orbit points using polar equation: r = a(1-e²)/(1+e cos θ)
5. Rotate to match the actual r₁, r₂ orientation

### Limitations
- Simplified 2D coplanar orbits (not 3D)
- Assumes point-mass gravity (no perturbations)
- Short-way transfers only (could be extended for long-way)
- May have numerical issues for very short or very long transfers

## Further Reading

- **Wikipedia**: [Lambert's problem](https://en.wikipedia.org/wiki/Lambert%27s_problem)
- **Poliastro docs**: [Lambert solver](https://docs.poliastro.space/)
- **Battin**: *An Introduction to the Mathematics and Methods of Astrodynamics*
- **Vallado**: *Fundamentals of Astrodynamics and Applications*

## Generated Output

The script creates `lambert_arc_animation.gif` - a 12-second looping animation at 15 fps showing 180 different transfer orbits.

File size: ~2-5 MB depending on complexity

## Requirements

```bash
pip install numpy matplotlib
# Optional: pip install poliastro  # For enhanced Lambert solver
```

---

**Tip**: The animation makes it visually clear why mission planners care so much about launch windows - only certain geometric configurations allow efficient (low ΔV) transfers!
