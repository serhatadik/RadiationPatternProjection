import numpy as np
import plotly.express as px

def plotter(value_map, title, UTM_long, UTM_lat, colorbar_title):
    value_map = value_map[::3, ::3]
    fig = px.imshow(value_map[1:value_map.shape[0] - 1, 1:value_map.shape[1] - 1], origin='lower',title=title)
    fig.update_layout(coloraxis_colorbar_x=0.61)
    fig.update_layout(coloraxis_colorbar_y=0.5)
    fig.update_layout(height=500)
    UTM_long = np.squeeze(UTM_long)
    UTM_lat = np.squeeze(UTM_lat)

    interval = max(1, len(UTM_lat) // 10)  # interval for ticks to avoid overlap, you can adjust if needed
    fig.update_xaxes(
        tickvals=list(range(0, len(UTM_lat), interval)),
        ticktext=UTM_lat[::interval],
        title_text="UTM<sub>E</sub> [m]",
        title_standoff=10,
        tickfont=dict(size=14),
        side='bottom'  # to place the xlabel at the bottom of the plot
    )
    fig.update_layout(
        title={
            'y': 1.0,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top',
            'font': dict(
                family="Times New Roman",  # Set font family to Times New Roman
                size=25
            )
        }
    )

    interval = max(1, len(UTM_long) // 10)  # interval for ticks to avoid overlap, you can adjust if needed
    fig.update_yaxes(
        tickvals=list(range(0, len(UTM_long), interval)),
        ticktext=UTM_long[::interval],
        title_text="UTM<sub>N</sub> [m]",
        title_standoff=10,
        tickfont=dict(size=14)
    )
    fig.update_layout(
        coloraxis_colorbar_x=0.64,
        coloraxis_colorbar_y=0.5,
        coloraxis_colorbar=dict(
            title=colorbar_title, # Add label to the colorbar
            title_side = "right"
        ),

        height=500
    )

    fig.show()
    return None
