"""
- grabs the smallest p-value and the fifth-percentile p-value (only 5% of values are smaller than it) for each model type
- outputs a single file for the smallest and the 5th-ptile values for each model
"""

### imports ==================================================================
from pprint import pprint as pprint 	# QOL
import os
import csv
import sys
import math


### variables ==================================================================
pval_dir = "/home/ualcpr/QTL/DLinkMaP/parallelizer/out_data/p-vals"
pval_ptile_dir = "/home/ualcpr/QTL/DLinkMaP/parallelizer/out_data/percentile_data"
pval_ptile_0 = "{0}/{1}".format( pval_ptile_dir, "pval_ptile_0.csv" )
pval_ptile_5 = "{0}/{1}".format( pval_ptile_dir, "pval_ptile_5.csv" )


### functions ==================================================================

def dir_getter( path, out_type="full" ):
	"""
	returns the full path of directories, or only the names of the directories
	"""
	import os
	directory_contents = os.listdir(path)

	full_paths = [ "{0}/{1}".format( path, subdir ) for subdir in directory_contents ]

	if out_type != "full":
		return( directory_contents )

	elif out_type == "full":
		return( full_paths )

### code ==================================================================


pval_files = dir_getter( pval_dir )
pval_files.sort()

header_row = [
	"run_number",
	"P-Val - Additive",
	"P-Val - Dominant",
	"P-Val - Full",
	"P-Val - Main",
	"P-Val - Add-D",
	"P-Val - Dom-D",
	"P-Val - Full-D",
	"P-Val - Diet"
]

with open( pval_ptile_0, "w" ) as file0, open( pval_ptile_5, "w" ) as file5:
	print( ",".join( header_row ), file=file0 )
	print( ",".join( header_row ), file=file5 )

	for pval_file in pval_files:

		to_print_file0 = []
		to_print_file5 = []

		permutation = pval_file[-8:-4]
		to_print_file0.append( permutation )
		to_print_file5.append( permutation )

		pval_add = []
		pval_dom = []
		pval_full = []
		pval_main = []
		pval_addD = []
		pval_domD = []
		pval_fullD = []
		pval_diet = []

		open_file = open( pval_file, "r" )
		reader = csv.reader( open_file )
		for row in reader:
			if "LR - Additive" not in row:
				pval_add.append( row[-8] )
				pval_dom.append( row[-7] )
				pval_full.append( row[-6] )
				pval_main.append( row[-5] )
				pval_addD.append( row[-4] )
				pval_domD.append( row[-3] )
				pval_fullD.append( row[-2] )
				pval_diet.append( row[-1] )
		open_file.close()

		pval_add = [ float(item) for item in pval_add if item != "NA" ]
		pval_dom = [ float(item) for item in pval_dom if item != "NA" ]
		pval_full = [ float(item) for item in pval_full if item != "NA" ]
		pval_main = [ float(item) for item in pval_main if item != "NA" ]
		pval_addD = [ float(item) for item in pval_addD if item != "NA" ]
		pval_domD = [ float(item) for item in pval_domD if item != "NA" ]
		pval_fullD = [ float(item) for item in pval_fullD if item != "NA" ]
		pval_diet = [ float(item) for item in pval_diet if item != "NA" ]

		pval_add.sort()
		pval_dom.sort()
		pval_full.sort()
		pval_main.sort()
		pval_addD.sort()
		pval_domD.sort()
		pval_fullD.sort()
		pval_diet.sort()

		ptile_5_index = math.floor( len( pval_add )*95 / 100 )

		to_print_file0.append( str(pval_add[-1]) )
		to_print_file0.append( str(pval_dom[-1]) )
		to_print_file0.append( str(pval_full[-1]) )
		to_print_file0.append( str(pval_main[-1]) )
		to_print_file0.append( str(pval_addD[-1]) )
		to_print_file0.append( str(pval_domD[-1]) )
		to_print_file0.append( str(pval_fullD[-1]) )
		to_print_file0.append( str(pval_diet[-1]) )

		to_print_file5.append( str(pval_add[ptile_5_index]) )
		to_print_file5.append( str(pval_dom[ptile_5_index]) )
		to_print_file5.append( str(pval_full[ptile_5_index]) )
		to_print_file5.append( str(pval_main[ptile_5_index]) )
		to_print_file5.append( str(pval_addD[ptile_5_index]) )
		to_print_file5.append( str(pval_domD[ptile_5_index]) )
		to_print_file5.append( str(pval_fullD[ptile_5_index]) )
		to_print_file5.append( str(pval_diet[ptile_5_index]) )

		print( "{0:>04} / {1:>04}".format( pval_files.index( pval_file )+1, len( pval_files ) ).center( 100, " " ), end="\r" ) # remove if enabling tqdm

		print( ",".join( to_print_file0 ), file=file0 )
		print( ",".join( to_print_file5 ), file=file5 )



# file = "/home/ualcpr/QTL/DLinkMaP/parallelizer/out_data/p-vals/male_pval_0001.csv"
# opened_csv_file = open(file, 'r')
# reader = csv.reader(opened_csv_file)
#
# outer_fruits_list = []
# for row in reader:
#         outer_fruits_list.append(row)
#
# opened_csv_file.close()
#
# pprint(outer_fruits_list[:3])
