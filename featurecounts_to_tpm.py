#!/home/bdhingpit/miniconda3/bin/python3

"""
CLI-based utility to normalize gene counts inferred by featureCounts to transcripts-per-million (TPM)
"""

import argparse
import pandas as pd


def convert_featurecounts_to_tpm(featurecounts_table):
    """Normalize gene counts to TPM"""
    featurecounts_df = pd.read_csv(featurecounts_table, sep="\t", header=1)

    # Divide by 1000 so gene length is expressed as per kbp
    featurecounts_df.Length = featurecounts_df.Length / 1000

    # Divide read counts by gene length = reads per kilobase (RPK)
    featurecounts_df.iloc[:, 6:] = featurecounts_df.iloc[:, 6:].div(featurecounts_df.Length, axis=0)

    # Divide RPK by sample/column sum of RPK
    featurecounts_df.iloc[:, 6:] = featurecounts_df.iloc[:, 6:].div(featurecounts_df.iloc[:, 6:].sum(axis=0), axis=1)

    # Multiply by 1e6 to get TPM
    featurecounts_df.iloc[:, 6:] = featurecounts_df.iloc[:, 6:] * 1000000

    # Revert length to original value
    featurecounts_df.Length = featurecounts_df.Length * 1000

    return featurecounts_df


def main(featurecounts_table, tpm_table):
    converted_fc_table = convert_featurecounts_to_tpm(featurecounts_table)
    converted_fc_table.to_csv(tpm_table, sep="\t", index=False)


if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        description='Convert featureCounts abundance table to TPM',
    )

    parser.add_argument(
        '--fc_table_in',
        dest='featurecounts_table',
        type=str,
        required=True,
        metavar='PATH',
        help='Path to the featureCounts input table',
    )

    parser.add_argument(
        '--tpm_table_out', dest='tpm_table', type=str, required=True, metavar='PATH', help='Path to the TPM output table'
    )

    args = parser.parse_args()

    featurecounts_table = args.featurecounts_table
    tpm_table = args.tpm_table

    main(featurecounts_table, tpm_table)
