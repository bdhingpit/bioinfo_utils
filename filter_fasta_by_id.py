#!/home/bdhingpit/miniconda3/bin/python3

"""
CLI-based utility to filter a FASTA file based on a list of FASTA headers
"""
# TODO: Instead of indicating filt_file_type, just use a parameter that indicates the column number

import argparse
import pandas as pd


def get_contig_ids_from_dvf(filt_file, conf_thresh=0.9, pval_thresh=0.01):
    """Filter DeepVirFinder output based on confidence and p-values"""
    dvf_df = pd.read_csv(filt_file, sep='\t')

    filt_dvf_df = dvf_df.loc[(dvf_df['score'] >= conf_thresh) & (dvf_df['pvalue'] <= pval_thresh)].reset_index()

    return list(filt_dvf_df['name'])


def get_contig_ids(filt_file):
    """Get the IDs listed in the MMSeqs2 output tsv file"""
    filt_file_df = pd.read_csv(filt_file, sep='\t', header=None)

    return list(filt_file_df[0])


def fasta_to_dict(fasta_file):
    """Represent FASTA file as a dictionary"""
    # TODO: Perhaps using Biopython would be more efficient here

    # Initiate variables
    fasta_dict = dict()
    fasta_header_old = ''  # previous fasta header
    fasta_header_cur = ''  # current fasta header
    sequence = ''  # fasta sequence

    with open(fasta_file, 'r') as fasta:
        # lines = fasta.readlines()
        # last = lines[-1]

        for line in fasta:
            if line.startswith('>'):
                # For the very first line
                if fasta_header_old == '':
                    fasta_header_old = line.strip()
                    fasta_header_cur = line.strip()

                # If a new fasta header is encountered
                else:
                    fasta_dict['{}'.format(fasta_header_cur)] = sequence
                    fasta_header_cur = line.strip()

            else:
                # If still in the same fasta entry
                if fasta_header_cur == fasta_header_old:
                    sequence += line

                # If in new fasta entry
                else:
                    sequence = line
                    fasta_header_old = fasta_header_cur

        # For last fasta entry; placed outside of loop which means loop has ended already
        fasta_dict['{}'.format(fasta_header_cur)] = sequence

    return fasta_dict


def dict_to_fasta(fasta_dict):
    """Convert id-sequence dictionary to FASTA"""
    fasta_string = ''

    for fasta_header in fasta_dict.keys():
        fasta_string += fasta_header + '\n'
        fasta_string += fasta_dict['{}'.format(fasta_header)]

    return fasta_string


def filter_fasta(fasta_dict, ids_to_inc):
    """Filter FASTA file"""
    filt_fasta_dict = dict()

    for header in fasta_dict.keys():
        if header.strip('>') in ids_to_inc:
            filt_fasta_dict['{}'.format(header)] = fasta_dict['{}'.format(header)]

    return filt_fasta_dict


def main(filt_file, fasta_file, filt_type, output_file):
    if filt_type == "DVF":
        filt_contig_ids = get_contig_ids_from_dvf(filt_file)
    elif filt_type == "GENERAL":
        filt_contig_ids = get_contig_ids(filt_file)

    fasta_dict = fasta_to_dict(fasta_file)
    filt_fasta_dict = filter_fasta(fasta_dict, filt_contig_ids)
    fasta_string = dict_to_fasta(filt_fasta_dict)

    output_fasta_file = open('{}'.format(output_file), 'w')
    output_fasta_file.write(fasta_string)
    output_fasta_file.close()


if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        usage=argparse.SUPPRESS,
        description='Filter FASTA file based on DeepVirFinder-predicted viral contigs',
    )

    parser.add_argument(
        '--filt_file_in',
        dest='filt_file',
        type=str,
        required=True,
        metavar='PATH',
        help='Path to the file to be used for filtering',
    )
    parser.add_argument(
        '--fasta_in', dest='fasta_file', type=str, required=True, metavar='PATH', help='Path to FASTA file to filter'
    )
    parser.add_argument(
        '--filt_file_type',
        dest='filt_type',
        type=str,
        required=True,
        choices=['DVF', 'GENERAL'],
        metavar='TEXT',
        help='Type of filter file. For GENERAL, contig IDs should be listed in the first column of the tsv file. Choices: DVF|GENERAL',
    )
    parser.add_argument(
        '--fasta_out', dest='output_file', required=True, metavar='PATH', help='Path to output filtered FASTA file'
    )

    args = parser.parse_args()

    filt_file = args.filt_file
    fasta_file = args.fasta_file
    filt_type = args.filt_type
    output_file = args.output_file

    main(filt_file, fasta_file, filt_type, output_file)
