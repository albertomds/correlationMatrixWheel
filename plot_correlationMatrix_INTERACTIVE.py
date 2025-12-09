import numpy as np
import plotly.graph_objects as go
import math
import dash
from dash import dcc, html
from dash.dependencies import Input, Output

#lsof -i:8050 & kill -9 to kill previous jobs

filenames = ["correl_01_9mku.dat", "correl_01_9mku.dat", "correl_01_9mku.dat", "correl_01_9mku.dat" ]  # Max 4
labels = ["FILE1", "FILE2", "FILE3", "FILE4"]

max_arcs = 2500           # Max arcs drawn per region pair
min_corr_threshold = 0.1  # Only draw arcs above/below this |value|
index_correction = +0     # Shift residue numbers shown in labels (This just affects the final visualization)
label_radius = 1.02       # Radius for residue-number labels
radius = 1.0              # Base circle radius

showIndices = False       # Show residue index numbers on outer ring
showColorBar = False      # Show colorbar 

#The Regions' number is as appear in your input, not affected by index_correction
regions = [
    (953, 1000),   # Region 1: BH-H1 
    (1009, 1021),  # Region 2: Lid
    (1296, 1303),  # Region 3: DNA2
    (1304, 1311)   # Region 4: DNA
]
region_names  = ['BH-H1', 'Lid', 'DNA2', 'DNA']
region_colors = ["#e181b0", "#00ff22", "#f01be2", "#f0db1b"]


#######################################################################

n_regions = len(regions)
angle_per_region = 360 / n_regions


def residue_angles(region_idx):
    start, end = regions[region_idx]
    n_res = end - start + 1
    rel_pos = np.linspace(0, 1, n_res, endpoint=False) + 0.5 / n_res
    base_angle = angle_per_region * region_idx
    angles_deg = base_angle + rel_pos * angle_per_region
    return np.deg2rad(angles_deg)


def interpolate_color(corr):
    corr = max(-1, min(1, corr))
    
    if corr < 0:
        t = (corr + 1)
        r = int(t * 255)
        g = int(t * 255)
        b = 255
    else:
        t = corr
        r = 255
        g = int((1 - t) * 255)
        b = int((1 - t) * 255)

    return f'rgba({r},{g},{b},0.5)'

def generate_figure(data, title, show_colorbar=False):
    fig = go.Figure()
    region_angles_arr = [residue_angles(i) for i in range(n_regions)]
    drawn_pairs = set()
    shapes = []

    for i in range(n_regions):
        start_i, end_i = regions[i]
        angles_i = region_angles_arr[i]
        indices = np.arange(start_i, end_i + 1) + index_correction

        x_lbl = label_radius * radius * np.cos(angles_i)
        y_lbl = label_radius * radius * np.sin(angles_i)

        fig.add_trace(go.Scatter(
            x=x_lbl,
            y=y_lbl,
            mode='text',
            text=[str(idx) for idx in indices] if showIndices else None,
            textposition='top center',
            textfont=dict(size=8, color='black'),
            hoverinfo='none',
            showlegend=False
        ))

        start_angle = math.radians(angle_per_region * i)
        end_angle = math.radians(angle_per_region * (i + 1))
        arc_x = 1.15 * radius * np.cos(np.linspace(start_angle, end_angle, 30))
        arc_y = 1.15 * radius * np.sin(np.linspace(start_angle, end_angle, 30))

        fig.add_trace(go.Scatter(
            x=arc_x,
            y=arc_y,
            mode='lines',
            line=dict(width=6, color=region_colors[i]),
            showlegend=False,
            hoverinfo='skip'
        ))

        mid_angle = (start_angle + end_angle) / 2
        x_txt = 1.5 * radius * math.cos(mid_angle)
        y_txt = 1.5 * radius * math.sin(mid_angle)

        fig.add_trace(go.Scatter(
            x=[x_txt],
            y=[y_txt],
            mode='text',
            text=[region_names[i]],
            textfont=dict(size=18, color=region_colors[i]),
            showlegend=False
        ))

    for i in range(n_regions):
        start_i, end_i = regions[i]
        angles_i = region_angles_arr[i]

        for j in range(i + 1, n_regions):
            start_j, end_j = regions[j]
            angles_j = region_angles_arr[j]

            submatrix = data[start_i:end_i + 1, start_j:end_j + 1]
            mask = np.abs(submatrix) >= min_corr_threshold
            if not np.any(mask):
                continue

            r_idx, c_idx = np.where(mask)
            corrs = submatrix[r_idx, c_idx]

            if len(corrs) > max_arcs:
                sorted_idx = np.argsort(np.abs(corrs))[::-1][:max_arcs]
                r_idx = r_idx[sorted_idx]
                c_idx = c_idx[sorted_idx]
                corrs = corrs[sorted_idx]

            x1 = radius * np.cos(angles_i[r_idx])
            y1 = radius * np.sin(angles_i[r_idx])
            x2 = radius * np.cos(angles_j[c_idx])
            y2 = radius * np.sin(angles_j[c_idx])

            for xi1, yi1, xi2, yi2, corr, ri, ci in zip(x1, y1, x2, y2, corrs, r_idx, c_idx):
                res_i = start_i + ri + index_correction
                res_j = start_j + ci + index_correction

                pair_key = tuple(sorted((res_i, res_j)))
                if pair_key in drawn_pairs:
                    continue

                drawn_pairs.add(pair_key)

                color = interpolate_color(corr)
                width = 0.5 + abs(corr) * 3.5
                path = f'M {xi1},{yi1} Q 0,0 {xi2},{yi2}'

                shapes.append(dict(
                    type='path',
                    path=path,
                    line=dict(color=color, width=width),
                    layer='below',
                    name=f'{res_i},{res_j}'
                ))

                all_arcs[filenames.index(title)].append((i, j, res_i, res_j, corr))

    fig.update_layout(
        width=450,
        height=400,
        plot_bgcolor='white',
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        margin=dict(l=10, r=70, t=60, b=20),
        shapes=shapes
    )

    if show_colorbar:
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode='markers',
            marker=dict(
                colorscale=[[0, 'blue'], [0.5, 'white'], [1, 'red']],
                cmin=-1,
                cmax=1,
                colorbar=dict(
                    title="Correlation",
                    tickvals=[-1, -0.5, 0, 0.5, 1],
                    ticktext=["-1", "-0.5", "0", "0.5", "1"],
                    len=0.6,
                    x=1.2,
                ),
                showscale=True,
                size=0
            ),
            hoverinfo='none',
            showlegend=False
        ))

    fig.update_yaxes(scaleanchor='x', scaleratio=1)
    return fig



all_arcs = [[] for _ in range(len(filenames))]

figures = [
    generate_figure(
        np.loadtxt(fname, comments='#'),
        fname,
        show_colorbar=(showColorBar and i == len(filenames) - 1)
    )
    for i, fname in enumerate(filenames)
]

# --- Dash ---
app = dash.Dash(__name__)

app.layout = html.Div([
    html.H2("Filter Residues"),
    dcc.Input(
        id='residue-input',
        type='text',
        placeholder='e.g., 45,102,200',
        style={'width': '600px'}
    ),
    html.Div([
        html.Button('Highlight', id='submit-button', n_clicks=0),
        html.Button('Table', id='table-button', n_clicks=0, style={'marginLeft': '10px'})
    ]),
    html.Div("Leave empty and press highlight to reset.", style={'marginBottom': '10px'}),

    html.Div(id='graph-container', children=[
        html.Div([
            html.H4(labels[i], style={'textAlign': 'center'}),
            dcc.Graph(id=f'graph-{i}', figure=figures[i])
        ], style={'width': '23%', 'minWidth': '300px', 'margin': '10px'})
        for i in range(len(filenames))
    ], style={'display': 'flex', 'flexWrap': 'wrap', 'justifyContent': 'center'}),

    html.Div([
        html.H4("Highlighted Correlations Table"),
        html.Pre(id='table-output', style={
            'whiteSpace': 'pre-wrap',
            'fontFamily': 'monospace',
            'border': '1px solid #ccc',
            'padding': '10px',
            'marginTop': '20px',
            'maxHeight': '300px',
            'overflowY': 'auto',
            'width': '50%'
        })
    ])
])


@app.callback(
    [Output(f'graph-{i}', 'figure') for i in range(len(filenames))],
    Input('submit-button', 'n_clicks'),
    Input('residue-input', 'value')
)
def update_all_graphs(n_clicks, input_value):
    selected = set()
    if input_value:
        parts = input_value.split(',')
        for part in parts:
            part = part.strip()
            if '-' in part:
                try:
                    start, end = map(int, part.split('-'))
                    selected.update(str(i) for i in range(start, end + 1))
                except ValueError:
                    continue
            else:
                if part.isdigit():
                    selected.add(part)

    output_figs = []
    for fig in figures:
        new_fig = go.Figure(fig)
        for shape in new_fig.layout.shapes:
            if 'name' not in shape:
                continue
            res_pair = shape['name'].split(',')
            shape.visible = (
                True if not selected or any(r in selected for r in res_pair) else False
            )
        output_figs.append(new_fig)
    return output_figs


@app.callback(
    Output('table-output', 'children'),
    Input('table-button', 'n_clicks'),
    Input('residue-input', 'value'),
    prevent_initial_call=True
)
def print_table(n_clicks_table, input_value):
    selected = set()
    if input_value:
        parts = input_value.split(',')
        for part in parts:
            part = part.strip()
            if '-' in part:
                try:
                    start, end = map(int, part.split('-'))
                    selected.update(range(start, end + 1))
                except ValueError:
                    continue
            elif part.isdigit():
                selected.add(int(part))

    output = ""
    for idx, arcs in enumerate(all_arcs):
        filtered = [
            (region1, region2, res1, res2, corr)
            for region1, region2, res1, res2, corr in arcs
            if not selected or res1 in selected or res2 in selected
        ]

        if not filtered:
            continue

        filtered.sort(key=lambda x: x[4], reverse=True)

        output += f"\n=== {labels[idx]} ===\n"
        output += f"{'Region1':<10} {'Region2':<10} {'Residue1':<10} {'Residue2':<10} {'Correlation':>10}\n"

        for region1, region2, res1, res2, corr in filtered:
            output += f"{region_names[region1]:<10} {region_names[region2]:<10} {res1:<10} {res2:<10} {corr:10.3f}\n"

    return output or "No correlations matched the filter."


if __name__ == '__main__':
    app.run(debug=True)
