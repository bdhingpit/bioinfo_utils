#!/home/bdhingpit/miniconda3/bin/python3

"""
CLI-based utility to create to add feature abundance information based on gene or contig abundance
"""

import argparse
import pandas as pd


def merge_contigs_abund_and_annot_tables(args):
    """Function connected with argparse `contigs` subcommand"""
    merged_df = merge_contig_tables(args.abund_table, args.annot_table, args.split_lineage)

    return merged_df


def merge_genes_abund_and_contig_tables(args):
    """Function connected with argparse `genes` subcommand"""
    merged_df = merge_gene_tables(args.abund_table, args.annot_table, args.annot_mode)

    return merged_df


def parse_featureCounts_prodigal_table(annot_df):
    """Create new columns for gene-related properties"""
    # Split Geneid column and append to Chr (contig name) the gene number
    annot_df[["Geneid", "Gene_num"]] = annot_df["Geneid"].str.split("_", expand=True)
    annot_df["Gene_name"] = annot_df["Chr"] + "_" + annot_df["Gene_num"].astype(str)

    # Select relevant columns
    annot_df = annot_df.iloc[:, [-1] + list(range(6, annot_df.shape[1] - 2))]

    return annot_df


def merge_contig_tables(abund_table, annot_table, will_split_lineage):
    """Merge contig abundance and gene annotation tables"""
    # Load
    abund_df = pd.read_csv(abund_table, sep="\t")
    annot_df = pd.read_csv(annot_table, sep="\t", header=None)

    # Merge
    merged_df = abund_df.merge(annot_df.iloc[:, list(range(0, 4)) + [-1]], left_on='Contig', right_on=0)

    # Rearrange the columns
    merged_df = merged_df.iloc[
        :, [0] + list(range(abund_df.shape[1] + 1, merged_df.shape[1])) + list(range(1, abund_df.shape[1]))
    ]

    # Rename headers of columns 2-5
    merged_df = merged_df.rename(columns={1: 'TaxonID', 2: 'TaxonLevel', 3: 'TaxonName', 8: 'TaxonLineage'})

    if will_split_lineage:
        merged_df = split_mmseqs_lineage(merged_df)

    return merged_df


def split_mmseqs_lineage(merged_df):
    """Create new columns with separated taxonomic lineage"""
    # Define the regex expressions
    superkingdom = r'd_(?P<superkingdom>[\w\s]+)'
    kingdom = r'\bk_(?P<kingdom>[\w\s]+)'
    phylum = r'\bp_(?P<phylum>[\w\s]+)'
    class_ = r'\bc_(?P<class>[\w\s]+)'
    order = r'\bo_(?P<order>[\w\s]+)'
    family = r'\bf_(?P<family>[\w\s]+)'
    genus = r'\bg_(?P<genus>[\w\s]+)'
    species = r'\bs_(?P<species>[\w\s]+)'

    ranks_regex_dict = dict(
        superkingdom=superkingdom,
        kingdom=kingdom,
        phylum=phylum,
        class_=class_,
        order=order,
        family=family,
        genus=genus,
        species=species,
    )

    split_lineage_df = pd.DataFrame()
    level = 0

    for rank in ranks_regex_dict.keys():
        rank_names = merged_df['TaxonLineage'].str.extract(ranks_regex_dict[rank])
        split_lineage_df[rank] = rank_names

        level += 1

    split_lineage_df = split_lineage_df.T.ffill().T

    merged_split_df = pd.concat([split_lineage_df, merged_df], axis=1)
    merged_split_df = merged_split_df.drop(['TaxonID', 'TaxonLevel', 'TaxonName', 'TaxonLineage'], axis=1)

    return merged_split_df


def merge_gene_tables(abund_table, annot_table, annot_mode):
    """Merge gene abundance table and contigs"""
    # Load
    abund_df = pd.read_csv(abund_table, sep="\t")

    if annot_mode == "eggnog":
        annot_df = parse_eggnog_table(annot_table)
    elif annot_mode == "kofamkoala":
        annot_df = parse_kofamkoala_table(annot_table)

    # Clean abund_df from featureCounts
    abund_df[["Geneid", "Gene_num"]] = abund_df["Geneid"].str.split("_", expand=True)
    abund_df["Gene_name"] = abund_df["Chr"] + "_" + abund_df["Gene_num"].astype(str)
    abund_df = abund_df.iloc[:, [-1] + list(range(6, abund_df.shape[1] - 2))]

    # Merge then rearrange
    merged_df = abund_df.merge(annot_df, on="Gene_name")
    merged_df = merged_df.iloc[:, [-1, 0] + list(range(1, merged_df.shape[1] - 1))]

    return merged_df


def parse_eggnog_table(annot_table):
    """Transform annotation table result from eggNOG"""
    # TODO: Add other options which gene groups (KO, COG, etc) to use

    annot_df = pd.read_csv(annot_table, sep="\t", header=4)

    annot_df = annot_df[["#query", "KEGG_ko"]]

    # If multiple KOs are present, take the first one only
    annot_df["KEGG_ko"] = annot_df["KEGG_ko"].str.split(",", expand=True)[0]

    # Remove rows with missing KO annotations
    annot_df = annot_df[annot_df["KEGG_ko"].str.contains("-") == False].reset_index(drop=True)

    annot_df["KEGG_ko"] = annot_df["KEGG_ko"].str[3:]
    annot_df.columns = ["Gene_name", "Annotation"]

    return annot_df


def parse_kofamkoala_table(annot_df):
    # TODO: Finish function
    pass


def main(args):
    merged_df = args.func(args)
    merged_df.to_csv(args.out, sep='\t', index=False)

    return None


if __name__ == '__main__':
    # Define parent parser
    parent_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description="Merge together an abundance table (of contigs or genes) and the corresponding feature's annotation",
        add_help=False,
    )

    # Define arguments always required
    required_args = parent_parser.add_argument_group('Required arguments')

    required_args.add_argument(
        '--abund_table', dest='abund_table', type=str, required=True, metavar='PATH', help='Path to the abundance table input'
    )

    required_args.add_argument(
        '--annot_table', dest='annot_table', type=str, required=True, metavar='PATH', help='Path to the annotation table input'
    )

    required_args.add_argument(
        '--out', dest='out', type=str, required=True, metavar='PATH', help='Path to output merged table'
    )

    # Define main parser
    main_parser = argparse.ArgumentParser()

    # Define subcommands
    subparsers = main_parser.add_subparsers(title='Feature type', dest='feature_type', description='Available subcommands')

    # Define contigs_subparser arguments
    contigs_subparser = subparsers.add_parser("contigs", parents=[required_args])
    contigs_subparser.add_argument(
        '--split_lineage',
        dest='split_lineage',
        type=str,
        required=False,
        metavar='BOOL',
        default=True,
        help='Whether or not to split the taxonomic lineage generated by MMSeqs2',
    )
    contigs_subparser.set_defaults(func=merge_contigs_abund_and_annot_tables)

    # Define genes_subparser arguments
    genes_subparser = subparsers.add_parser("genes", parents=[required_args])
    genes_subparser.add_argument(
        '--annot_mode',
        dest='annot_mode',
        type=str,
        required=True,
        choices=['eggnog', 'kofamkoala'],
        metavar='BOOL',
        default=True,
        help='Gene annotation method used. Choices ["eggnog"]',
    )
    genes_subparser.set_defaults(func=merge_genes_abund_and_contig_tables)

    args = main_parser.parse_args()

    main(args)
