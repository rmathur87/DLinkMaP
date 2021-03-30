"""
copies all the data for each permutation into the out_data directory

Method:
- read in all data
- find which ones have been run (non-zero size)
- for those, remane the files, and copy the files into the appropriate out_data subdirectories
"""
### imports ==================================================================
from pprint import pprint as pprint
import os
import glob
import shutil
import time

### variables ==================================================================

run_directory = "/scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts"
loglike_dir = "/home/ualcpr/QTL/DLinkMaP/parallelizer/out_data/logLike"
pval_dir = "/home/ualcpr/QTL/DLinkMaP/parallelizer/out_data/p-vals"

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

### code ==================================================================""


run_dirs = dir_getter( run_directory )
n_runs = len( run_dirs )
run_dirs.sort()

finished_runs = [] # permutations that have not been run

for dir in run_dirs:
	maleWt_csv = [ "{0}/maleWt/{1}".format( dir, f ) for f in os.listdir("{0}/maleWt".format( dir )) if f.endswith('.csv') ]

	total_file_size = 0

	for file in maleWt_csv:
		total_file_size += os.path.getsize( file )

	if total_file_size != 0:
		finished_runs.append( dir[-4:] )

	print( "        Checking run_{0}: {1:04}/{2:04}".format( dir[-4:], run_dirs.index( dir ), len( run_dirs ) ), end="\r" )

print()

# copy files

for run in finished_runs:

	# scratch_dir = "/scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0}".format( run )
	# home_dir = "/home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0}".format( run )
	#
	# if os.path.exists( home_dir ):
	# 	shutil.rmtree( home_dir )
	#
	# shutil.copytree( scratch_dir, home_dir )
	# time.sleep( 1 )

	# copy p-val file
	pval_src = "/scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0}/maleWt/p-value_Male_avgbyvial.csv".format( run )
	pval_dst = "/home/ualcpr/QTL/DLinkMaP/parallelizer/out_data/p-vals/male_pval_{0}.csv".format( run )
	shutil.copyfile( pval_src, pval_dst )

	# copy logLike file
	logLike_src = "/scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0}/maleWt/logLike_Male_avgbyvial.csv".format( run )
	logLike_dst = "/home/ualcpr/QTL/DLinkMaP/parallelizer/out_data/logLike/male_logLike_{0}.csv".format( run )
	shutil.copyfile( logLike_src, logLike_dst )

	print( "Finished copying run_{0}: {1:04}/{2:04}".format( run, finished_runs.index( run ), len( finished_runs ) ), end="\r" )
