#!/home/bdhingpit/miniconda3/envs/my_bioinfo_env/bin/python

"""
For practice only. Visualize as a stacked bar plot a profile (e.g. taxonomic profile)

CLI-based utility to
"""


import plotly.express as px
import pandas as pd
import argparse
import colorlover as cl
import random
from warnings import simplefilter


# def check_metadata_var_in_profile_table(profile_df, metadata_df):
#     # Check if the metadata variables are present in the profile table
#     # Metadata variables are used to group together the samples accordingly
#     pass


def normalize_table(profile_df, metadata_df):
    """Normalize profile counts using sampling depth"""
    num_vars = len(metadata_df.columns[1:])

    # Convert variable type
    feature_cols = profile_df.iloc[:, 1:-num_vars].columns
    convert_dict = dict.fromkeys(feature_cols, "float64")
    profile_df = profile_df.astype(convert_dict)

    # Get depth per sample and normalize
    per_sample_depth = profile_df.iloc[:, 1:-num_vars].sum(axis=1)
    profile_df.iloc[:, 1:-num_vars] = profile_df.iloc[:, 1:-num_vars].div(per_sample_depth, axis=0) * 100

    return profile_df


def melt_feature_table(profile_df, metadata_df):
    """Melt profile using sample names/columns as IDs"""
    num_vars = len(metadata_df.columns[1:])

    sample_col = profile_df.columns[0]
    # metadata_col = list(profile_df.columns[-num_vars:-1])
    features_col = list(profile_df.columns[1:-num_vars])

    melt_profile_df = pd.melt(profile_df, id_vars=[sample_col], value_vars=features_col, var_name="Features")

    return melt_profile_df


def define_layout(fig):
    """Define plot layout"""
    # Define font formats
    title_font_format = dict(size=20, family="Open Sans", color="#2e1c18")
    axis_font_format = dict(size=12, family="Open Sans", color="#2e1c18")
    legend_font_format = dict(size=8, family="Open Sans", color="#2e1c18")

    # Define layout
    fig.update_layout(
        template="plotly_white",
        title=dict(text="Taxonomy Relative Abundance Barplot", x=0.5, font=title_font_format),
        yaxis_title="% Abundance",
        yaxis_title_font=axis_font_format,
        xaxis_title="Index",
        xaxis_title_font=axis_font_format,
        legend=dict(font=legend_font_format),
    )

    return fig


def create_stacked_barplot(profile_df):
    """Create the stacked barplot"""
    # Define color scale for bars
    col_scale = cl.scales["10"]["div"]["Spectral"]
    col_scale_features = cl.interp(col_scale, len(profile_df["Features"].unique()))
    random.Random(4).shuffle(col_scale_features)

    fig = px.bar(profile_df, x="index", y="value", color="Features", color_discrete_sequence=col_scale_features)
    fig = define_layout(fig)
    fig.show()

    return fig


def main(args):
    simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

    profile_df = pd.read_csv(args.profile, sep="\t")
    metadata_df = pd.read_csv(args.metadata, sep="\t", comment="#")

    norm_profile_df = normalize_table(profile_df, metadata_df)
    melt_profile_df = melt_feature_table(norm_profile_df, metadata_df)
    fig = create_stacked_barplot(melt_profile_df)

    fig.write_html("./stacked_barplot.html")


if __name__ == "__main__":
    parent_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Produce interactive stacked barplot from a feature table profile",
        add_help=True,
    )

    # Define arguments always required
    required_args = parent_parser.add_argument_group('Required arguments')

    required_args.add_argument(
        "--profile", dest="profile", type=str, required=True, metavar="STRING", help="Path to profile table in TSV format"
    )

    required_args.add_argument(
        "--format",
        dest="format",
        type=str,
        required=True,
        choices=["QIIME2"],
        metavar="STRING",
        help="Format of the profile table. Default: \"QIIME2\" || Choices: [\"QIIME2\"]",
    )

    required_args.add_argument(
        "--metadata", dest="metadata", type=str, required=True, metavar="STRING", help="Path to metadata file in TSV format"
    )

    args = parent_parser.parse_args()

    main(args)
