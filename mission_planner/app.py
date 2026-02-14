"""
Interplanetary Mission Planner - Main Streamlit Application

Web interface for planning interplanetary missions with interactive
porkchop plots and payload calculations.
"""
import streamlit as st
import plotly.graph_objects as go
import numpy as np
from datetime import datetime, date
import sys
import os

# Add current directory to path for imports
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from data.rockets import get_all_rockets, get_rocket_names, get_rocket
from data.porkchop_loader import (
    load_porkchop_csv,
    filter_porkchop_data,
    get_available_destinations,
    check_data_availability
)
from calculations.payload_calc import (
    calculate_payload_matrix,
    calculate_c3_from_delta_v,
    heliocentric_dv_to_departure_dv
)
from ui.rocket_comparison import render_comparison_page


# Page configuration
st.set_page_config(
    page_title="Interplanetary Mission Planner",
    page_icon="🚀",
    layout="wide",
    initial_sidebar_state="expanded"
)


# Custom CSS for better styling
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        margin-bottom: 2rem;
    }
    .metric-container {
        background-color: #f0f2f6;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
    }
</style>
""", unsafe_allow_html=True)


def create_porkchop_plot(rocket, destination, dv_matrix, launch_dates, tof_days):
    """
    Create interactive Plotly heatmap for porkchop plot with payload as color.

    Args:
        rocket: Rocket object
        destination: Destination planet name
        dv_matrix: 2D array of delta-V values (m/s)
        launch_dates: Array of launch dates
        tof_days: Array of time-of-flight values (days)

    Returns:
        Plotly Figure object
    """
    # Convert heliocentric delta-V to departure delta-V from LEO
    # Porkchop data shows Sun-frame delta-V, but payload calc needs LEO departure delta-V
    departure_dv_matrix = np.vectorize(heliocentric_dv_to_departure_dv)(dv_matrix)

    # Calculate payload for each trajectory point
    payload_matrix = calculate_payload_matrix(rocket, departure_dv_matrix)

    # Convert payload to tons for display
    payload_tons = payload_matrix / 1000

    # Convert launch dates to strings for display
    launch_date_strings = [d.strftime('%Y-%m-%d') for d in launch_dates]

    # Create custom hover template
    # Note: x=launch dates, y=TOF after transposing
    hover_template = (
        "<b>Launch:</b> %{x}<br>"
        "<b>TOF:</b> %{y} days<br>"
        "<b>Payload:</b> %{z:.1f} tons<br>"
        "<b>ΔV:</b> %{customdata:,.0f} m/s<br>"
        "<extra></extra>"
    )

    # Transpose the data (swap axes)
    payload_tons_T = payload_tons.T
    dv_matrix_T = dv_matrix.T

    # Create year markers for x-axis
    years_set = sorted(set(d.year for d in launch_dates))
    year_positions = []
    year_labels = []
    for year in years_set:
        # Find first occurrence of this year
        for i, d in enumerate(launch_dates):
            if d.year == year:
                year_positions.append(i)
                year_labels.append(str(year))
                break

    # Create heatmap with payload as color (black = 0, bright colors = high payload)
    fig = go.Figure(data=go.Heatmap(
        z=payload_tons_T,  # Transposed: rows=TOF, cols=launch dates
        x=launch_date_strings,  # Launch dates on X-axis
        y=tof_days,  # TOF on Y-axis
        customdata=dv_matrix_T,  # Delta-V shown in hover
        hovertemplate=hover_template,
        colorscale=[
            [0.0, 'black'],      # Zero payload = black
            [0.1, 'darkblue'],
            [0.3, 'blue'],
            [0.5, 'cyan'],
            [0.7, 'yellow'],
            [0.9, 'orange'],
            [1.0, 'red']         # Max payload = red
        ],
        colorbar=dict(
            title="Payload<br>(tons)",
            tickformat=",.0f",
        ),
        zmin=0,  # Start at zero
        zmax=np.percentile(payload_tons[payload_tons > 0], 98) if np.any(payload_tons > 0) else 1
    ))

    # Update layout with swapped axes
    fig.update_layout(
        title=f"Porkchop Plot: Earth to {destination} ({rocket.name})",
        xaxis_title="Launch Date",
        yaxis_title="Time of Flight (days)",
        height=600,
        hovermode='closest',
        font=dict(size=12),
        xaxis=dict(
            tickmode='array',
            tickvals=[launch_date_strings[i] for i in year_positions],
            ticktext=year_labels,
            tickangle=-45
        ),
        yaxis=dict(
            tickformat='d'
        )
    )

    return fig


def main():
    """Main application."""

    # Header
    st.markdown('<p class="main-header">🚀 Interplanetary Mission Planner</p>', unsafe_allow_html=True)
    st.markdown(
        '<p class="sub-header">Plan interplanetary missions with interactive porkchop plots and payload calculations</p>',
        unsafe_allow_html=True
    )

    # Page navigation
    page = st.sidebar.radio(
        "Navigation",
        ["Mission Planner", "Rocket Comparison"],
        index=0
    )

    st.sidebar.markdown("---")

    # Route to appropriate page
    if page == "Rocket Comparison":
        render_comparison_page()
        return

    # Sidebar: Mission Parameters
    st.sidebar.title("Mission Parameters")

    # Rocket selection
    rocket_names = get_rocket_names()
    selected_rocket_name = st.sidebar.selectbox(
        "🚀 Rocket Type",
        rocket_names,
        index=rocket_names.index('Starship refuelled in LEO') if 'Starship refuelled in LEO' in rocket_names else 0
    )
    rocket = get_rocket(selected_rocket_name)

    # Show rocket specs
    with st.sidebar.expander("📊 Rocket Specifications"):
        st.write(f"**Wet Mass:** {rocket.wet_mass_kg/1000:.1f} tons")
        st.write(f"**Dry Mass:** {rocket.dry_mass_kg/1000:.1f} tons")
        st.write(f"**Specific Impulse:** {rocket.isp_s:.1f} s")
        st.write(f"**LEO Payload:** {rocket.leo_payload_kg/1000:.1f} tons")
        st.write(f"**Starting Orbit:** {rocket.starting_orbit}")

    # Destination selection
    available_destinations = get_available_destinations()
    destination = st.sidebar.selectbox(
        "🌍 Destination",
        available_destinations,
        index=0 if 'Mars' not in available_destinations else available_destinations.index('Mars')
    )

    # Check data availability
    if not check_data_availability(destination):
        st.error(
            f"⚠️ Pre-computed data for {destination} not found!\n\n"
            f"Run `porkchop-data-generator.ipynb` to generate data."
        )
        st.stop()

    # Load porkchop data
    with st.spinner(f"Loading {destination} trajectory data..."):
        try:
            dv_matrix, launch_dates, tof_days = load_porkchop_csv(destination)
        except Exception as e:
            st.error(f"Error loading data: {e}")
            st.stop()

    # Date range filter
    st.sidebar.subheader("🗓️ Launch Window")
    min_date = launch_dates[0].date()
    max_date = launch_dates[-1].date()

    date_range = st.sidebar.date_input(
        "Date Range",
        value=(min_date, min_date + (max_date - min_date) // 10),
        min_value=min_date,
        max_value=max_date
    )

    # TOF range filter
    st.sidebar.subheader("⏱️ Time of Flight")
    min_tof = int(tof_days[0])
    max_tof = int(tof_days[-1])

    tof_range = st.sidebar.slider(
        "TOF Range (days)",
        min_value=min_tof,
        max_value=max_tof,
        value=(100, 300)
    )

    # Filter data
    if len(date_range) == 2:
        start_date = datetime.combine(date_range[0], datetime.min.time())
        end_date = datetime.combine(date_range[1], datetime.min.time())
        dv_filtered, dates_filtered, tof_filtered = filter_porkchop_data(
            dv_matrix, launch_dates, tof_days,
            date_range=(start_date, end_date),
            tof_range_days=tof_range
        )
    else:
        dv_filtered, dates_filtered, tof_filtered = dv_matrix, launch_dates, tof_days

    # Main content area
    col1, col2 = st.columns([3, 1])

    with col1:
        st.subheader(f"📊 Porkchop Plot: Earth → {destination}")

        if dv_filtered.size == 0:
            st.warning("No data available for selected date/TOF range. Try expanding the filters.")
        else:
            # Create and display plot
            fig = create_porkchop_plot(rocket, destination, dv_filtered, dates_filtered, tof_filtered)
            st.plotly_chart(fig, use_container_width=True)

            # Info about the plot
            with st.expander("ℹ️ How to Read This Plot"):
                st.markdown("""
                - **Color Scale**: Black (zero payload) → Blue → Cyan → Yellow → Orange → Red (maximum payload)
                - **Hover**: Shows launch date, time of flight, deliverable payload, and required ΔV
                - **Bright Red/Orange Regions**: Optimal launch windows with highest payload capacity
                - **Black Regions**: Infeasible trajectories (trajectory requires more ΔV than rocket can provide)
                - **Blue/Cyan Regions**: Marginal trajectories (low payload capacity)

                **Tip**: Look for the brightest red/orange areas to find launch windows that maximize payload delivery!
                """)

    with col2:
        st.subheader("📈 Mission Statistics")

        # Find optimal trajectory
        valid_dv = dv_filtered[dv_filtered > 0]
        if valid_dv.size > 0:
            min_dv = np.min(valid_dv)
            min_idx = np.unravel_index(np.argmin(dv_filtered + (dv_filtered == 0) * 1e10), dv_filtered.shape)
            optimal_launch = dates_filtered[min_idx[0]]
            optimal_tof = tof_filtered[min_idx[1]]

            # Calculate payload for optimal trajectory
            from calculations.payload_calc import calculate_max_payload
            optimal_payload = calculate_max_payload(rocket, min_dv)
            optimal_c3 = calculate_c3_from_delta_v(min_dv)

            # Display metrics
            st.metric("🎯 Minimum ΔV", f"{min_dv:,.0f} m/s")
            st.metric("📦 Max Payload (@ min ΔV)", f"{optimal_payload/1000:.1f} tons")
            st.metric("🌟 Optimal Launch", optimal_launch.strftime('%Y-%m-%d'))
            st.metric("⏱️ Flight Time", f"{optimal_tof:.0f} days")
            st.metric("🔋 C3 Energy", f"{optimal_c3:.1f} km²/s²")

            # Trajectory details
            st.markdown("---")
            st.markdown("**Trajectory Details:**")
            st.markdown(f"- Arrival: {(optimal_launch + np.timedelta64(int(optimal_tof), 'D')).strftime('%Y-%m-%d')}")
            st.markdown(f"- Total ΔV: {min_dv:,.0f} m/s")
            st.markdown(f"- Payload: {optimal_payload:,.0f} kg")

        else:
            st.info("Adjust filters to see trajectory statistics")

        # Data summary
        st.markdown("---")
        st.markdown("**Dataset Info:**")
        st.markdown(f"- Launch dates: {len(dates_filtered)}")
        st.markdown(f"- TOF points: {len(tof_filtered)}")
        st.markdown(f"- Total trajectories: {dv_filtered.size:,}")


if __name__ == "__main__":
    main()
