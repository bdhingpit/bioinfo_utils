#!/home/bdhingpit/miniconda3/envs/my_bioinfo_env/bin/python

"""
For practice only.

CLI-based utilty to show visualize quality profile of a FASTQ file (i.e. similar to FastQC).
"""

import plotly.graph_objects as go
import colorlover as cl
import skbio
from skbio.sequence import DNA
import pandas as pd
import numpy as np
import argparse
import time
from warnings import simplefilter


def define_layout(plot_title, x_title, y_title, y_range, show_legend):
    """
    Define layout for plot

    Module: per_base_qc_plot
    Module: fastq_stat
    """
    layout = go.Layout(
        title=dict(text=plot_title, x=0.5),
        template="plotly_white",
        yaxis=dict(title=y_title, showgrid=True, gridcolor="#d9d4d3", zerolinecolor="#d9d4d3", range=y_range),
        xaxis=dict(
            title=x_title,
        ),
        font=dict(family="Open Sans", size=16, color="#2e1c18"),
        showlegend=show_legend,
    )

    return layout


def module_per_base_qc_plot(args):
    """
    Display per-base quality plot

    Module: per_base_qc_plot
    """
    fastq_df = tabulate_fastq_info(args.fastq, args.encoding, args.num_seqs)
    show_per_base_qc_plot(fastq_df, None)


def tabulate_fastq_info(fastq, encoding, num_seqs):
    """
    Store in a dataframe the per-sequence quality information of the FASTQ library

    Module: per_base_qc_plot
    """
    imported_fastq = list(skbio.io.read(fastq, format='fastq', verify=False, variant=encoding))

    fastq_df = pd.DataFrame()

    for count, seq in enumerate(imported_fastq):
        fastq_df[count] = seq.positional_metadata.quality

        if count == num_seqs - 1:
            break

    fastq_df.index += 1

    return fastq_df


def define_color_scale():
    """
    Define color scale for quality scores

    Module: per_base_qc_plot
    """
    col_scales = cl.scales['4']['div']['RdYlGn']
    col_scales_40 = cl.interp(col_scales, 40)

    return col_scales_40


def create_boxes(fastq_df, col_scales):
    """
    Create boxplots per base position

    Module: per_base_qc_plot
    """
    traces = []

    for base in range(len(fastq_df)):
        traces.append(
            go.Box(
                # y=fastq_df.iloc[base].values,
                name="Base Position Quality",
                x=[base + 1],
                boxpoints=False,
                whiskerwidth=0.5,
                marker=dict(size=0.1, color=col_scales[int(round(fastq_df.iloc[base].mean(), 0))]),
                line=dict(width=1),
                q1=[np.nanpercentile(fastq_df.iloc[base].values, 25)],
                q3=[np.nanpercentile(fastq_df.iloc[base].values, 75)],
                median=[np.nanpercentile(fastq_df.iloc[base].values, 50)],
                lowerfence=[np.nanpercentile(fastq_df.iloc[base].values, 10)],
                upperfence=[np.nanpercentile(fastq_df.iloc[base].values, 90)],
                hoverlabel=dict(namelength=-1, align="left"),
            )
        )

    return traces


def show_per_base_qc_plot(fastq_df, output):
    """
    Display plot

    Module: per_base_qc_plot
    """
    col_scales_40 = define_color_scale()
    traces = create_boxes(fastq_df, col_scales_40)
    layout = define_layout(
        "Per-Base Quality Score", "Base Position", "Quality Score", [0, 41], False  # Max possible Q-score = 41
    )

    fig = go.Figure(data=traces, layout=layout)
    # fig.update_yaxes(ticksuffix = "    ")

    fig.show()


def module_fastq_stat(args):
    """
    Show basic FASTQ library stats

    Module: fastq_stat
    """
    # TODO: `fastq_stat` module is unfinished
    imported_fastq = list(skbio.io.read(args.fastq, format='fastq', verify=False, variant=args.encoding))

    show_fastq_stat_plots(imported_fastq)
    # get_gc_content(imported_fastq)


def show_fastq_stat_plots(fastq):
    """
    Display distribution of sequence lengths

    Module: fastq_stat
    """
    hist_trace = create_seq_length_hist(fastq)

    fig = go.Figure(data=hist_trace)
    fig.show()


def create_seq_length_hist(fastq):
    """
    Create histogram for sequence length

    Module: fastq_stat
    """
    # TODO: Probably just merge this with `show_fastq_stat_plots`
    fastq_seq_lengths = [len(seq) for seq in fastq]

    hist_trace = go.Histogram(x=fastq_seq_lengths)

    return hist_trace


def create_gc_content(fastq):
    """
    Estimate GC content of FASTQ library

    Module: fastq_stat
    """
    # TODO: Unfinished
    concat_seq = ""

    for seq in fastq:
        concat_seq += str(seq)

    gc_content = DNA(concat_seq).gc_content()

    return gc_content


def get_num_seqs(fastq):
    """
    Get number of sequences in FASTQ library

    Module: fastq_stat
    """
    # TODO: Unfinished
    len(fastq)


def main(args):
    # Ignore pandas warnings
    simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

    args.func(args)


if __name__ == "__main__":
    parent_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter, description="Inspect FASTQ files interactively", add_help=False
    )

    # Define arguments always required
    required_args = parent_parser.add_argument_group('Required arguments')

    required_args.add_argument("--fastq", dest="fastq", type=str, required=True, metavar="STRING", help="Path to FASTQ file")

    required_args.add_argument(
        "--encoding",
        dest="encoding",
        type=str,
        required=False,
        choices=["sanger", "illumina1.3", "illumina1.8", "solexa"],
        default="illumina1.8",
        metavar="STRING",
        help="Encoding system used in FASTQ file. Default: \"illumina1.8\" || Choices: [\"sanger\", \"illumina1.3\", \"illumina1.8\", \"solexa\"]",
    )

    # Define main parser
    main_parser = argparse.ArgumentParser()

    subparser = main_parser.add_subparsers(title="Analysis type", dest="analysis_type", description="Available subcommands")

    # Parser for per_base_qc_plot module
    per_base_qc_plot_subparser = subparser.add_parser("per_base_qc_plot", parents=[required_args])

    per_base_qc_plot_subparser.add_argument(
        "--num_seqs",
        dest="num_seqs",
        type=int,
        required=False,
        default=10000,
        metavar="INTEGER",
        help="Number of sequences to subsample for per-base QC plot. Default: 10000",
    )

    per_base_qc_plot_subparser.set_defaults(func=module_per_base_qc_plot)

    # Parser for fastq_stat module
    fastq_stat_subparser = subparser.add_parser("fastq_stat", parents=[required_args])

    fastq_stat_subparser.set_defaults(func=module_fastq_stat)

    args = main_parser.parse_args()

    main(args)
