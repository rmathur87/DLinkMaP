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
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


### variables ==================================================================
pval_dir = "/Users/rele.c/Downloads/DLinkMaP/parallelizer/out_data/p-vals"
pval_ptile_dir = "/Users/rele.c/Downloads/DLinkMaP/parallelizer/out_data/percentile_data"
pval_ptile_0 = "{0}/{1}".format( pval_ptile_dir, "pval_ptile_0.csv" )
pval_ptile_5 = "{0}/{1}".format( pval_ptile_dir, "pval_ptile_5.csv" )
histogram_dir = "/Users/rele.c/Downloads/DLinkMaP/parallelizer/out_data/percentile_data/histograms"
hist_0_dir = "{0}/{1}".format( histogram_dir, "hist_0" )
hist_5_dir = "{0}/{1}".format( histogram_dir, "hist_5" )

os.system( "mkdir {0}".format( pval_ptile_dir ) )


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

def make_hist( data_frame, header_name, image_save_dir ):
	"""
	Creates a histogram and saves it to the particular location along with proper names
	"""
	import seaborn as sns
	import matplotlib.pyplot as plt
	import pandas as pd
	import numpy as np

	data_frame[ header_name ] = pd.to_numeric( data_frame[ header_name ],errors='coerce')

	f, ( ax_box, ax_hist ) = plt.subplots( 2, sharex=True, gridspec_kw={"height_ratios": (.1, .9)} )

	# adding n_bins using the Freedman-Diaconis rule.
	Q1 = np.percentile(data_frame[ header_name ], 25, interpolation = 'midpoint')
	Q3 = np.percentile(data_frame[ header_name ], 75, interpolation = 'midpoint')
	IQR = Q3-Q1

	n_bins = int(round(( max(data_frame[ header_name ]) - min(data_frame[ header_name ]) )/( (2*IQR)/len(data_frame[ header_name ])**(1./3.) ), 0))

	print((n_bins))
	print(type(n_bins))

	sns.boxplot( data_frame[header_name], ax=ax_box )
	sns.histplot( data=data_frame, x=header_name, ax=ax_hist, bins=n_bins )

	ax_box.set(xlabel='')
	# plt.show()
	plt.savefig( "{0}/{1}.png".format( image_save_dir, header_name.lower().replace(" ", "_") ) )



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

histogram_data_0 = []
histogram_data_5 = []

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

		histogram_data_0.append( to_print_file0 )
		histogram_data_5.append( to_print_file5 )

		print( "run_{2} = {0:>04} / {1:>04}".format( pval_files.index( pval_file )+1, len( pval_files ), permutation ).center( 100, " " ), end="\r" ) # remove if enabling tqdm

		print( ",".join( to_print_file0 ), file=file0 )
		print( ",".join( to_print_file5 ), file=file5 )

# make histograms
hist_0_df = pd.DataFrame( histogram_data_0 , columns = header_row )
hist_5_df = pd.DataFrame( histogram_data_5 , columns = header_row )

for feature in header_row:
	make_hist( hist_0_df, feature, hist_0_dir )
	make_hist( hist_5_df, feature, hist_5_dir )
