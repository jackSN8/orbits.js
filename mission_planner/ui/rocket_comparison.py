"""
Rocket comparison tool for side-by-side payload analysis.
"""
import streamlit as st
import plotly.graph_objects as go
import numpy as np
from typing import List

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from data.rockets import Rocket
from calculations.payload_calc import calculate_max_payload, calculate_delta_v


def create_comparison_chart(rockets: List[Rocket], delta_v_required: float) -> go.Figure:
    """
    Create bar chart comparing payload capacity across multiple rockets.

    Args:
        rockets: List of Rocket objects to compare
        delta_v_required: Required delta-V in m/s

    Returns:
        Plotly Figure object
    """
    rocket_names = []
    payloads_tons = []
    colors = []

    for rocket in rockets:
        payload_kg = calculate_max_payload(rocket, delta_v_required)
        payload_tons = payload_kg / 1000

        rocket_names.append(rocket.name)
        payloads_tons.append(payload_tons)

        # Color code based on payload
        if payload_tons > 100:
            colors.append('green')
        elif payload_tons > 10:
            colors.append('orange')
        elif payload_tons > 0:
            colors.append('red')
        else:
            colors.append('gray')

    fig = go.Figure(data=[
        go.Bar(
            x=rocket_names,
            y=payloads_tons,
            marker_color=colors,
            text=[f"{p:.1f}t" for p in payloads_tons],
            textposition='outside',
        )
    ])

    fig.update_layout(
        title=f"Payload Comparison for ΔV = {delta_v_required:,.0f} m/s",
        xaxis_title="Rocket",
        yaxis_title="Maximum Payload (tons)",
        height=500,
        showlegend=False,
        xaxis=dict(tickangle=-45)
    )

    return fig


def create_delta_v_chart(rockets: List[Rocket], payload_kg: float) -> go.Figure:
    """
    Create bar chart showing achievable delta-V for fixed payload.

    Args:
        rockets: List of Rocket objects to compare
        payload_kg: Fixed payload mass in kg

    Returns:
        Plotly Figure object
    """
    rocket_names = []
    delta_vs = []

    for rocket in rockets:
        # Check if this payload is even possible
        if payload_kg > rocket.leo_payload_kg:
            delta_v = 0
        else:
            delta_v = calculate_delta_v(rocket, payload_kg)

        rocket_names.append(rocket.name)
        delta_vs.append(delta_v)

    fig = go.Figure(data=[
        go.Bar(
            x=rocket_names,
            y=delta_vs,
            marker_color='steelblue',
            text=[f"{dv:,.0f}" for dv in delta_vs],
            textposition='outside',
        )
    ])

    fig.update_layout(
        title=f"ΔV Comparison for {payload_kg/1000:.1f} ton Payload",
        xaxis_title="Rocket",
        yaxis_title="Achievable ΔV (m/s)",
        height=500,
        showlegend=False,
        xaxis=dict(tickangle=-45)
    )

    return fig


def render_comparison_page():
    """Render the rocket comparison page in Streamlit."""
    from data.rockets import get_all_rockets

    st.title("🔄 Rocket Comparison Tool")
    st.markdown("Compare payload capacity across different launch vehicles for the same mission.")

    # Rocket selection
    st.subheader("Select Rockets to Compare")
    all_rockets = get_all_rockets()
    rocket_names = [r.name for r in all_rockets]

    # Default selection: Starship variants + common rockets
    default_selections = [
        'Starship refuelled in LEO',
        'Starship refuelled in HEO',
        'Falcon 9 expended',
        'Atlas V 551'
    ]
    default_indices = [i for i, name in enumerate(rocket_names) if name in default_selections]

    selected_indices = st.multiselect(
        "Choose rockets (select 2-6 for best comparison)",
        range(len(all_rockets)),
        default=default_indices[:4],
        format_func=lambda i: all_rockets[i].name
    )

    if len(selected_indices) < 2:
        st.warning("Please select at least 2 rockets to compare.")
        st.stop()

    selected_rockets = [all_rockets[i] for i in selected_indices]

    # Comparison mode
    st.subheader("Comparison Mode")
    mode = st.radio(
        "What would you like to compare?",
        ["Fixed ΔV (find max payload)", "Fixed Payload (find max ΔV)"]
    )

    col1, col2 = st.columns(2)

    if mode == "Fixed ΔV (find max payload)":
        # User specifies delta-V, we calculate payload
        with col1:
            delta_v = st.slider(
                "Required ΔV (m/s)",
                min_value=1000,
                max_value=15000,
                value=5000,
                step=100
            )

            # Reference points
            st.markdown("**Reference missions:**")
            st.markdown("- Mars: ~3,500 - 5,000 m/s")
            st.markdown("- Jupiter: ~6,000 - 9,000 m/s")
            st.markdown("- Saturn: ~8,000 - 11,000 m/s")

        with col2:
            fig = create_comparison_chart(selected_rockets, delta_v)
            st.plotly_chart(fig, use_container_width=True)

        # Detailed comparison table
        st.subheader("Detailed Comparison")
        comparison_data = []
        for rocket in selected_rockets:
            payload_kg = calculate_max_payload(rocket, delta_v)
            comparison_data.append({
                'Rocket': rocket.name,
                'Max Payload (tons)': f"{payload_kg/1000:.2f}",
                'LEO Capacity (tons)': f"{rocket.leo_payload_kg/1000:.1f}",
                'Isp (s)': f"{rocket.isp_s:.1f}",
                'Starting Orbit': rocket.starting_orbit
            })

        import pandas as pd
        st.dataframe(pd.DataFrame(comparison_data), use_container_width=True)

    else:
        # User specifies payload, we calculate delta-V
        with col1:
            payload_tons = st.slider(
                "Payload Mass (tons)",
                min_value=1.0,
                max_value=200.0,
                value=50.0,
                step=1.0
            )
            payload_kg = payload_tons * 1000

        with col2:
            fig = create_delta_v_chart(selected_rockets, payload_kg)
            st.plotly_chart(fig, use_container_width=True)

        # Detailed comparison table
        st.subheader("Detailed Comparison")
        comparison_data = []
        for rocket in selected_rockets:
            if payload_kg > rocket.leo_payload_kg:
                delta_v = 0
                status = "❌ Exceeds LEO capacity"
            else:
                delta_v = calculate_delta_v(rocket, payload_kg)
                status = "✓"

            comparison_data.append({
                'Rocket': rocket.name,
                'ΔV (m/s)': f"{delta_v:,.0f}",
                'Status': status,
                'LEO Capacity (tons)': f"{rocket.leo_payload_kg/1000:.1f}",
                'Payload Fraction': f"{(payload_kg/rocket.leo_payload_kg*100):.1f}%"
            })

        import pandas as pd
        st.dataframe(pd.DataFrame(comparison_data), use_container_width=True)


if __name__ == "__main__":
    # Quick test
    from data.rockets import ROCKETS

    rockets = [
        ROCKETS['Starship refuelled in LEO'],
        ROCKETS['Falcon Heavy expended'],
        ROCKETS['Atlas V 551']
    ]

    fig = create_comparison_chart(rockets, 5000)
    print("Comparison chart created successfully!")
