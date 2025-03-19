import plotly.graph_objs as go
import json
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Generate CNV plots from JSON data.')
parser.add_argument('json_file', type=str, help='Path to the JSON file containing SV data')
args = parser.parse_args()

# Load your JSON data
with open(args.json_file) as f:
    sv_data = json.load(f)

# Loop through each SV (assuming your JSON contains multiple SVs)
for sv in sv_data:
    print(type(sv))

    # Extract data for plotting
    positions_before = sv['before_sv']['positions']
    b_allele_freq_before = sv['before_sv']['b_allele_freq']
    positions_after = sv['after_sv']['positions']
    b_allele_freq_after = sv['after_sv']['b_allele_freq']
    
    # Generate hover text (optional, can be customized)
    hover_text_before = [f"Position: {pos}, BAF: {baf}" for pos, baf in zip(positions_before, b_allele_freq_before)]
    hover_text_after = [f"Position: {pos}, BAF: {baf}" for pos, baf in zip(positions_after, b_allele_freq_after)]
    
    # Plotting data for 'before_sv' and 'after_sv'
    baf_trace_before = go.Scatter(
        x=positions_before,
        y=b_allele_freq_before,
        mode="markers+lines",
        name="B-Allele Frequency (Before SV)",
        text=hover_text_before,
        marker=dict(
            color='blue',
            size=10
        ),
        line=dict(
            color="black",
            width=0
        ),
        showlegend=False
    )

    baf_trace_after = go.Scatter(
        x=positions_after,
        y=b_allele_freq_after,
        mode="markers+lines",
        name="B-Allele Frequency (After SV)",
        text=hover_text_after,
        marker=dict(
            color='red',
            size=10
        ),
        line=dict(
            color="black",
            width=0
        ),
        showlegend=False
    )
    
    # Create layout for the plot
    layout = go.Layout(
        title=f"SV Plot: {sv['chromosome']} {sv['start']}-{sv['end']} ({sv['sv_type']})",
        xaxis=dict(title="Position"),
        yaxis=dict(title="B-Allele Frequency"),
        hovermode='closest'
    )
    
    # Create figure with data and layout
    fig = go.Figure(data=[baf_trace_before, baf_trace_after], layout=layout)
    
    # Save the plot to an HTML file (use a unique filename per SV)
    file_name = f"output/SV_{sv['chromosome']}_{sv['start']}_{sv['end']}.html"
    fig.write_html(file_name)

    print(f"Plot saved as {file_name}")
