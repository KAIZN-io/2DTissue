import pandas as pd
import numpy as np
from numpy.polynomial import Polynomial
from pathlib import Path
import os
import plotly.graph_objs as go
from plotly.subplots import make_subplots

current_path = Path(os.getcwd())
data = pd.read_csv(current_path.parents[1] / "assets" / 'v_order_analysis_low_noise.csv', sep=';')
data = data[data['particles'] % 100 == 0]  # drop all rows where particles number is not Mod 100
# data = pd.read_csv(current_path.parents[1] / "assets" / 'v_order_data_0_to_80_noise.csv', sep=',')

# Extract x and y values
x = data['particles']
y = data['order_parameter']

# Fit the data with a 3rd-degree polynomial
poly_fit = Polynomial.fit(x, y, 3)

# Create a more detailed x range for a smoother polynomial fit
x_fit = np.arange(100, 1001, 1)

# Create a subplot with Plotly
fig = make_subplots(specs=[[{"secondary_y": True}]])

# Add the box plot to the subplot
fig.add_trace(
    go.Box(
        x=data['particles'],
        y=data['order_parameter'],
        name='Order Parameter',
        fillcolor='rgba(230, 230, 250, 1)',
        boxmean=True,  # Show the mean as a dashed line
        marker_color='rgba(0, 76, 153, 0.5)',  # Change the box plot marker color
        line=dict(color='rgba(0, 76, 153, 1)')
    ),
    secondary_y=False,
)

# Add the polynomial fit to the subplot
fig.add_trace(
    go.Scatter(
        x=x_fit,
        y=poly_fit(x_fit),
        mode='lines',
        line=dict(color='red', width=2),
        name='Polynomial Fit (3rd degree)',
    ),
    secondary_y=False,
)

# Our design for the x and y axes
axis_design = dict(
    title_standoff=25,
    title_font=dict(size=32, color='black'),
    tickfont=dict(size=26),
    ticks='outside',  # Draw ticks outside the plot area
    ticklen=5,  # Length of the ticks
    tickwidth=2,  # Width of the ticks
    tickcolor='black',  # Color of the ticks
    gridcolor='rgba(211, 211, 211, 0.5)',  # light gray with alpha 0.5
    zerolinecolor='lightgray',
    linecolor='black',  # Set the line color for the x-axis
    linewidth=2
)

fig.update_xaxes(title_text="Particles Number", **axis_design)
fig.update_yaxes(title_text="Order Parameter", **axis_design)

# Hide the secondary y-axis
fig.update_yaxes(showticklabels=False, secondary_y=True)

# Customize the plot appearance
fig.update_layout(
    # title="Order Parameter vs Particles Number",
    # title_font=dict(size=36),
    hovermode="x unified",
    plot_bgcolor="white",
    legend=dict(
        font=dict(size=32),
        yanchor="top",
        y=0.99,
        xanchor="right",
        x=1.40
    ),
    margin=dict(l=80, r=50, t=80, b=80)
)

# Display the plot
fig.show()
