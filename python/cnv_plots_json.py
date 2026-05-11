#!/usr/bin/env python3
import argparse
import json
import os

DEFAULT_MIN_SV_LENGTH = 50000
MARKER_SIZE = 8
STATIC_FORMATS = {"svg", "pdf", "png", "jpg", "jpeg", "webp", "eps"}
ALLOWED_FORMATS = STATIC_FORMATS.union({"html"})

STATE_COLORS = {
    "1": "darkred",
    "2": "red",
    "3": "gray",
    "4": "green",
    "5": "blue",
    "6": "darkblue",
}

SV_TYPE_LABELS = {
    "DEL": "Deletion",
    "DUP": "Duplication",
    "INV": "Inversion",
}

REQUIRED_SECTION_KEYS = {
    "positions",
    "b_allele_freq",
    "population_freq",
    "log2_ratio",
    "is_snp",
}

REQUIRED_SV_KEYS = {
    "chromosome",
    "start",
    "end",
    "sv_type",
    "size",
    "before_sv",
    "sv",
    "after_sv",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Generate CNV plots from JSON data.")
    parser.add_argument("json_file", type=str, help="Path to the JSON file containing SV data")
    parser.add_argument(
        "chromosome",
        type=str,
        nargs="?",
        default=None,
        help="Chromosome to filter SVs by (e.g., chr3)",
    )
    parser.add_argument(
        "--formats",
        type=str,
        default="html",
        help="Comma-separated output formats (e.g., html,svg,pdf,png)",
    )
    parser.add_argument("--width", type=int, default=1200, help="Figure width in pixels for static exports")
    parser.add_argument("--height", type=int, default=800, help="Figure height in pixels for static exports")
    parser.add_argument("--scale", type=float, default=2.0, help="Scale factor for raster exports")
    parser.add_argument(
        "--min-sv-length",
        type=int,
        default=DEFAULT_MIN_SV_LENGTH,
        help="Minimum SV length in base pairs to plot",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Directory to save output plots (default: ./CNV_Plots)",
    )
    return parser.parse_args()


def parse_formats(formats_arg):
    formats = [fmt.strip().lower() for fmt in formats_arg.split(",") if fmt.strip()]
    invalid = [fmt for fmt in formats if fmt not in ALLOWED_FORMATS]
    if invalid:
        allowed = ", ".join(sorted(ALLOWED_FORMATS))
        bad = ", ".join(invalid)
        raise ValueError(f"Unsupported format(s): {bad}. Allowed formats are: {allowed}")
    return formats


def load_json_records(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Input JSON file not found: {path}")

    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)

    if not isinstance(data, list):
        raise ValueError("Input JSON must contain a top-level list of SV records.")

    return data


def validate_record(record, idx):
    missing_sv_keys = REQUIRED_SV_KEYS - set(record.keys())
    if missing_sv_keys:
        missing = ", ".join(sorted(missing_sv_keys))
        raise ValueError(f"Record {idx} missing required top-level key(s): {missing}")

    for section in ["before_sv", "sv", "after_sv"]:
        missing_section_keys = REQUIRED_SECTION_KEYS - set(record[section].keys())
        if missing_section_keys:
            missing = ", ".join(sorted(missing_section_keys))
            raise ValueError(f"Record {idx}, section {section} missing key(s): {missing}")


def build_hover_text(section, positions, states, log2_ratio, is_snp, b_allele_freq, population_freq):
    hover_text = []
    for i, position in enumerate(positions):
        if section == "sv":
            hover_text.append(
                f"Position: {position}<br>"
                f"State: {states[i]}<br>"
                f"Log2 Ratio: {log2_ratio[i]}<br>"
                f"SNP: {is_snp[i]}<br>"
                f"BAF: {b_allele_freq[i]}<br>"
                f"Population Frequency: {population_freq[i]}<br>"
            )
        else:
            hover_text.append(
                f"Position: {position}<br>"
                f"Log2 Ratio: {log2_ratio[i]}<br>"
                f"BAF: {b_allele_freq[i]}<br>"
                f"Population Frequency: {population_freq[i]}<br>"
            )
    return hover_text


def add_section_traces(fig, record, section, start, end):
    positions = record[section]["positions"]
    b_allele_freq = record[section]["b_allele_freq"]
    population_freq = record[section]["population_freq"]
    log2_ratio = record[section]["log2_ratio"]
    is_snp = record[section]["is_snp"]

    b_allele_freq = [freq if snp_flag else float("nan") for freq, snp_flag in zip(b_allele_freq, is_snp)]

    if section == "sv":
        states = record[section].get("states", ["NA"] * len(positions))
        state_colors = [STATE_COLORS.get(str(state), "black") for state in states]
        marker_symbols = ["circle" if snp_flag else "circle-open" for snp_flag in is_snp]
    else:
        states = ["NA"] * len(positions)
        state_colors = ["black"] * len(positions)
        marker_symbols = ["circle" if snp_flag else "circle-open" for snp_flag in is_snp]

    hover_text = build_hover_text(section, positions, states, log2_ratio, is_snp, b_allele_freq, population_freq)

    import plotly

    log2_trace = plotly.graph_objs.Scatter(
        x=positions,
        y=log2_ratio,
        mode="markers+lines",
        name="Log2 Ratio",
        text=hover_text,
        hoverinfo="text",
        marker={"color": state_colors, "size": MARKER_SIZE, "symbol": marker_symbols},
        line={"color": "black", "width": 0},
        showlegend=False,
    )

    baf_trace = plotly.graph_objs.Scatter(
        x=positions,
        y=b_allele_freq,
        mode="markers+lines",
        name="B-Allele Frequency",
        text=hover_text,
        hoverinfo="text",
        marker={"color": state_colors, "size": MARKER_SIZE, "symbol": marker_symbols},
        line={"color": "black", "width": 0},
        showlegend=False,
    )

    if section == "sv":
        fig.add_vrect(x0=start, x1=end, fillcolor="black", layer="below", line_width=0, opacity=0.1)
        fig.add_vline(x=start, line_width=2, line_color="black", layer="below")
        fig.add_vline(x=end, line_width=2, line_color="black", layer="below")

    fig.append_trace(log2_trace, row=1, col=1)
    fig.append_trace(baf_trace, row=2, col=1)


def build_figure(record, width, height):
    from plotly.subplots import make_subplots

    chromosome = record["chromosome"]
    start = record["start"]
    end = record["end"]
    sv_type = record["sv_type"]
    sv_length = record["size"]

    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.05,
        subplot_titles=("SNP Log2 Ratio", "SNP B-Allele Frequency"),
    )

    for section in ["before_sv", "sv", "after_sv"]:
        add_section_traces(fig, record, section, start, end)

    fig.update_xaxes(title_text="Chromosome Position", row=2, col=1)
    fig.update_yaxes(title_text="Log2 Ratio", range=[-2.0, 2.0], row=1, col=1)
    fig.update_yaxes(title_text="B-Allele Frequency", range=[-0.2, 1.2], row=2, col=1)

    title_label = SV_TYPE_LABELS.get(sv_type, sv_type)
    fig.update_layout(
        title_text=f"{title_label} at {chromosome}:{start}-{end} ({sv_length} bp)",
        title_x=0.5,
        showlegend=False,
        template="simple_white",
        font={"family": "Arial", "size": 20, "color": "black"},
        width=width,
        height=height,
        margin={"l": 100, "r": 30, "t": 120, "b": 90},
    )
    fig.update_xaxes(showline=True, linewidth=2, linecolor="black", mirror=True, ticks="outside")
    fig.update_yaxes(showline=True, linewidth=2, linecolor="black", mirror=True, ticks="outside")
    return fig


def write_outputs(fig, base_name, output_dir, formats, width, height, scale):
    if "html" in formats:
        html_path = os.path.join(output_dir, f"{base_name}.html")
        fig.write_html(html_path)
        print(f"Plot saved as {html_path}")

    requested_static_formats = [fmt for fmt in formats if fmt in STATIC_FORMATS]
    if requested_static_formats:
        try:
            for fmt in requested_static_formats:
                out_path = os.path.join(output_dir, f"{base_name}.{fmt}")
                fig.write_image(out_path, format=fmt, width=width, height=height, scale=scale)
                print(f"Plot saved as {out_path}")
        except ValueError as err:
            print("Static image export requires Kaleido. Install with: pip install -U kaleido")
            raise err


def main():
    args = parse_args()
    try:
        import plotly  # noqa: F401
    except ImportError as err:
        raise ImportError(
            "Missing required dependency 'plotly'. Install with: conda install -c conda-forge plotly"
        ) from err

    formats = parse_formats(args.formats)
    records = load_json_records(args.json_file)

    output_dir = args.output_dir if args.output_dir else os.path.join(os.getcwd(), "CNV_Plots")
    os.makedirs(output_dir, exist_ok=True)

    skip_count_chrom = 0
    skip_count_length = 0
    save_count = 0

    for idx, record in enumerate(records, start=1):
        validate_record(record, idx)

        if args.chromosome and record["chromosome"] != args.chromosome:
            skip_count_chrom += 1
            continue

        if abs(record["size"]) < args.min_sv_length:
            skip_count_length += 1
            continue

        fig = build_figure(record, args.width, args.height)

        sv_length = record["size"]
        svlen_kb = sv_length // 1000
        base_name = (
            f"SV_{record['chromosome']}_{record['start']}_{record['end']}_"
            f"{record['sv_type']}_{svlen_kb}kb"
        )

        write_outputs(fig, base_name, output_dir, formats, args.width, args.height, args.scale)
        save_count += 1

    print(
        f"Finished processing {save_count} SVs. "
        f"Skipped {skip_count_chrom} SVs due to chromosome filter and "
        f"{skip_count_length} SVs due to length filter."
    )


if __name__ == "__main__":
    main()
