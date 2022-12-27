#!/usr/bin/env python3

'''
Definition of terms:
> ATTRIBUTE or ATTR - refer to any fields (mostly the LOCUS_TAG and ID) in the 9th column of the GFF file
> LOCUS_TAG - the locus tag field in the 9th column of the GFF file
> ID - the id field in the 9th column of the GFF file
> PRODUCT_NAME - the assigned product name as listed in the *.product_name file
'''

import sys
import argparse
import pandas as pd

pd.set_option('display.max_colwidth', None)
pd.set_option('display.max_columns', None)

###################################
#
# GENERAL FILE PROCESSING FUNCTIONS
#
###################################

#Load GFF file
def load_gff(gff_file):
	return pd.read_csv(gff_file, sep='\t', header=None)

#Load *.product_name file
def load_prod_name(prod_names_file):
	return pd.read_csv(prod_names_file, sep='\t', header=None)

#Check *.product_namees and *.gff file if there are duplicate LOCUS TAGs
def check_dup_locus_tags(prod_names_loc_tag, gff_loc_tag):
	prod_names_loc_tag_dup = prod_names_loc_tag[prod_names_loc_tag.duplicated(cols=0)]

	pass

#Gets the attribute "locus_tag=" in col9 of GFF file
def get_gff_loc_tags(gff_df):
	gff_col9_df = gff_df.iloc[:, 8].str.split(';', expand=True) #split using ; as separator
	loc_tag_bool_df = gff_col9_df.apply(lambda row: row.astype(str).str.contains('locus_tag='), axis=1) #make truth table of columns that contain "locus_tags"
	loc_tag_series = gff_col9_df[loc_tag_bool_df].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1) #get columns that are non-NaN based on truth table
	loc_tag_only = loc_tag_series.str.replace('locus_tag=', '', regex=True) #removed "locus_tag="

	return loc_tag_only

#Gets the attribute "ID=" in col9 of GFF file
def get_gff_ids(gff_df):
	gff_col9_df = gff_df.iloc[:, 8].str.split(';', expand=True) #split using ; as separator
	id_bool_df = gff_col9_df.apply(lambda row: row.astype(str).str.contains('ID='), axis=1) #make truth table of columns that contain "locus_tags"
	id_series = gff_col9_df[id_bool_df].apply(lambda x: ','.join(x.dropna().astype(str)), axis=1) #get columns that are non-NaN based on truth table
	id_only = id_series.str.replace('[ID=|.]', '', regex=True)

	return id_only

#Check if all ID field in the GFF file is the first ATTRIBUTE in the 9th column
def is_id_field_set_first(gff_df):
	col9_df = gff_df[8].str.split(';', expand=True)

	id_field_bool_df = col9_df[0].str.contains('ID=')

	idx_not_id_first_list = id_field_bool_df[~id_field_bool_df].index.values

	if len(idx_not_id_first_list) == 0:
		return True

	else:
		return False


##################
#
# SUBCOMMAND MODES
#
##################

#Converts the -1, 1 in col7 to - and +, respectively
def parse_col7(args):
	#load gff file
	gff_df = load_gff(args.gff_file)

	#in col7, if -1 replace by -, if +1, replace by +, if neither, replace by .
	gff_df.loc[gff_df[6] == -1, 6] = '-'
	gff_df.loc[gff_df[6] == 1, 6] = '+'
	gff_df.loc[(gff_df[6] != '+') & (gff_df[6] != '-'), 6] = '.'

	print(gff_df)

	#Save file
	gff_df.to_csv(args.output_prefix+'_PARSE_COL7.gff', sep='\t', header=False, index=False)
	return gff_df


#Return the GFF file with the PRODUCT_NAME added to the ID field in GFF's 9th col (for now, this fxn uses the LOCUS_TAG to do the matching to PRODUCT_NAME)
def add_prod_name_to_id(args):
	#load needed files
	gff_df = load_gff(args.gff_file)
	prod_names_df = load_prod_name(args.product_name_file)

	#Check first if ID fields are first ATTRIBUTES
	if is_id_field_set_first(gff_df) == False:
		print('Exiting - Not all ID fields are set as first ATTRIBUTE in column 9. Row number(s):')
		print(idx_not_id_first_list+1)

		return

	loc_tags = get_gff_loc_tags(gff_df) #get the LOCUS_TAG ATTRIBUTE
	gff_w_loc_tags = pd.concat([gff_df, loc_tags.to_frame(name='loc_tags')], axis=1) #concat to the orig GFF the LOCUS_TAGs

	#merge GFF and product names file based on LOCUS_TAGs
	merged_df = gff_w_loc_tags.merge(prod_names_df, left_on='loc_tags', right_on=0).iloc[:, 0:12]
	merged_df.columns = list(range(len(merged_df.columns)))

	#Split col9 using ; as separator
	splt_col9_df = merged_df[8].str.split(';', expand=True)

	#Create a new ID string (merged_df[10] is the LOCUS_TAG and merged_df[11] is the PRODUCT_NAME)
	new_id_field_ser = 'Name=' + merged_df[10] + '_' + merged_df[11]

	#Replace the orig ID string by the new ID string
	splt_col9_df[0] = new_id_field_ser.to_frame(name=0)

	#Rejoin together col9 
	new_col9_ser = splt_col9_df.apply(lambda x: ';'.join(x.values.astype(str)), axis=1)

	#Remove quotation marks in some of the col9 rows
	new_col9_ser = new_col9_ser.map(lambda x: x.lstrip('"').rstrip('"'))

	#Cleanup; remove ";None*" at the end part of the entire col9 string
	new_col9_ser = new_col9_ser.str.replace(';None*', '', regex=True)

	#Apply new col9 to GFF file
	merged_df[8] = new_col9_ser.to_frame(name=8)

	#Drop the cols past the ATTRIBUTE cols (col9)
	merged_df.drop(columns=[9,10,11], inplace=True)

	print(merged_df)

	#Save file
	merged_df.to_csv(args.output_prefix+'_ID_w_PROD_NAME.gff', sep='\t', header=False, index=False)
	return gff_df


#Adds an attribute to col9 of GFF file
def add_attribute(args):
	gff_df = load_gff(args.gff_file)
	name_attr_map_df = pd.read_csv(args.locus_attr_map_file, sep='\t', header=None)

	#Aggregate ATTR_VALUE based on same GENE_IDS
	new_name_attr_map_df = name_attr_map_df[[0, 1]].drop_duplicates().reset_index(drop=True)
	new_name_attr_map_df[1] = new_name_attr_map_df.groupby([0])[1].transform(lambda x: ','.join(x))
	new_name_attr_map_df.drop_duplicates(inplace=True)

	#Get LOCUS_TAGS
	loc_tags = get_gff_loc_tags(gff_df) #get the LOCUS_TAG ATTRIBUTE
	gff_w_loc_tags = pd.concat([gff_df, loc_tags.to_frame(name='loc_tags')], axis=1) #concat to the orig GFF the LOCUS_TAGs

	#merge GFF and ATTR map file file based on LOCUS_TAGs
	merged_df = gff_w_loc_tags.merge(new_name_attr_map_df, left_on='loc_tags', right_on=0, how='left').iloc[:, 0:12]
	merged_df.columns = list(range(len(merged_df.columns)))

	#Add ATTR_CLASS to col9 where applicable
	idx_with_match = merged_df[10][~merged_df[10].isnull()].index.tolist() #Get all index in GFF file with a match in the map file
	idx_wo_attr_class = merged_df[8][~merged_df[8].str.contains(args.attr_class)].index.tolist() #Get all index in GFF col9 without the ATTR_CLASS yet
	
	intersect = list(set(idx_wo_attr_class) & set(idx_with_match)) #Get the intersection of the 2 lists above
	intersect.sort() #Sort the list

	merged_df.iloc[intersect, 8] = merged_df.iloc[intersect, 8] + args.attr_class + '=' #Add the ATTR_CLASS to the rows without the ATTR_CLASS yet

	#Split the rows of interest in col9 into different ATTR
	splt_df = merged_df.iloc[idx_with_match, 8].str.split(';', expand=True)

	#concat the ATTR_VAL and splt_df
	splt_df = pd.concat([splt_df, merged_df.iloc[idx_with_match, -1]], axis=1)
	splt_df.columns = list(range(len(splt_df.columns)))

	#Merge column containing ATTR_CLASS and ATTR_VAL
	splt_df_bool = splt_df.apply(lambda x: x.astype(str).str.contains(args.attr_class)) #Find where the ATTR_CLASS are

	attr_class_val_ser = pd.Series(0, index=idx_with_match) #Instantiate a series
	exstng_attr_class = splt_df[splt_df_bool].apply(lambda x: ''.join(x.dropna().astype(str)), axis=1) #Create a Series for the existing and newly-created attr_class

	for i in idx_with_match:
		if exstng_attr_class[i].endswith('{}='.format(args.attr_class)):
			attr_class_val_ser[i] = exstng_attr_class[i] + splt_df.loc[i, splt_df.columns[-1]] + ';'

		else:
			attr_class_val_ser[i] = exstng_attr_class[i] + ',' + splt_df.loc[i, splt_df.columns[-1]]

	for i in splt_df.index.tolist():
		splt_df.loc[i, :].replace(to_replace=r'^{}.*'.format(args.attr_class), regex=True, value=attr_class_val_ser[i], inplace=True) #***DOUBLE CHECK THIS ONE, MIGHT HAVE A BUG WHEN ADDING MORE ATTR_VALS TO THE SAME ATTR_CLASS

	#Rejoin splt_df
	col9_mod = splt_df.iloc[:, 0:-1].apply(lambda x: ';'.join(x.dropna().astype(str)), axis=1)

	#Replace orig col9 in gff_df by modified col9
	gff_df.loc[col9_mod.index, 8] = col9_mod

	print(gff_df)

	#Export final gff file
	gff_df.to_csv(args.output_prefix+'_w_ADDED_ATTR.gff', sep='\t', header=False, index=False)
	

######
#
# MAIN
#
######

def main():
	#Argument parser
	parser = argparse.ArgumentParser(
		prog='gff_parser.py',
		description='Perform different processes to GFF files'
		)

	subparsers = parser.add_subparsers()
	subparsers.metavar = 'Sub-commands:'

	#1st subcommand
	parser_fxn1 = subparsers.add_parser('parse_col7', help='Changes -1 and +1 in 7th column to - and +, respectively')
	parser_fxn1.add_argument('gff_file', help='Path to GFF file')
	parser_fxn1.add_argument('output_prefix', help='Prefix of the output reformatted GFF file')
	parser_fxn1.set_defaults(func=parse_col7)

	#2nd subcommand
	parser_fxn2 = subparsers.add_parser('add_prod_name_to_id', help='Adds the PRODUCT_NAME to the ID field in the 9th column of the GFF file')
	parser_fxn2.add_argument('gff_file', help='Path to GFF file')
	parser_fxn2.add_argument('product_name_file', help='Path to *.product_name file')
	parser_fxn2.add_argument('output_prefix', help='Prefix of the output reformatted GFF file')
	parser_fxn2.set_defaults(func=add_prod_name_to_id)

	#3rd subcommand
	parser_fxn3 = subparsers.add_parser('add_attribute', help='Add ATTRIBUTES to col9 given a map of LOCUS_TAG to ATTRIBUTE_VALUE')
	parser_fxn3.add_argument('gff_file', help='Path to GFF file')
	parser_fxn3.add_argument('locus_attr_map_file', help='A TSV file containing the LOCUS_TAG in the first column and the ATTRIBUTE_VALUE in the second column')
	parser_fxn3.add_argument('attr_class', help='Name of attribute class (e.g. ')
	parser_fxn3.add_argument('output_prefix', help='Prefix of the output reformatted GFF file')
	parser_fxn3.set_defaults(func=add_attribute)


	args = parser.parse_args()
	args.func(args)


if __name__ == '__main__':
	main()