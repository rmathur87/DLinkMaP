"""
Identifies which runs "failed", or otherwise did not run.
Returns a list of numbers, so that those specific runs can be run again manually.

If the number of runs is too large, we may have to make it create a bash file to run (likely course of action).

- If any of the file sizes are 0, it is likely that it did not run
"""

### imports ==================================================================
from pprint import pprint as pprint
import os
import glob

### variables ==================================================================

run_directory = "/home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts"
rerun_script = "/scratch/ualcpr/QTL/DLinkMaP/parallelizer/code/rerun.sh"
scratch_copy = "/home/ualcpr/QTL/DLinkMaP/parallelizer/code/copy_to_scratch.sh"

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

not_run = [] # permutations that have not been run

for dir in run_dirs:
	maleWt_csv = [ "{0}/maleWt/{1}".format( dir, f ) for f in os.listdir("{0}/maleWt".format( dir )) if f.endswith('.csv') ]

	total_file_size = 0

	for file in maleWt_csv:
		total_file_size += os.path.getsize( file )

	if total_file_size == 0:
		not_run.append( dir[-4:] )

with open( rerun_script, "w" ) as scfile:
	print( "# to rerun only models that have not been previously run/failed to run before", file=scfile)
	print( "# {1} models not/failed to run: {0}".format( "\t".join(not_run), len(not_run) ), file=scfile )
	print( file=scfile )
	for run in not_run:

		print( "cd /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0}".format( run ), file=scfile )
		print( "rm perm{0}shSCRIPT.*".format( run ), file=scfile )
		print( "rm maleWt/*.csv".format( run ), file=scfile )
		# print( "mv perm{0}.sh temp.sh".format( run ), file=scfile )
		# print( "sed 's/home/scratch/g' temp.sh > perm{0}.sh".format( run ), file=scfile )
		print( "chmod +x perm{0}.sh".format( run ), file=scfile )
		print( "run_script perm{0}.sh".format( run ), file=scfile )
		print( "sleep 10".format( run ), file=scfile )
		print( file=scfile )

		print( "Made script for run_{0}: {1:04}/{2:04} = {3:.3f}%".format( run, not_run.index( run ), len( not_run ), (not_run.index(run)/len(not_run))*100 ), end="\r" )


# with open( scratch_copy, "w" ) as scratch_file:
# 	print( "# to copy uncompleted models from home to scratch", file=scratch_file )
# 	print( "# copying models: {0}".format( "\t".join(not_run) ), file=scratch_file )
# 	print( file=scratch_file )
# 	for run in not_run:
# 		home_dir = "/home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0}".format( run )
# 		scratch_dir = "/scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0}".format( run )
#
# 		print( "cp -r {0} {1}".format( home_dir, scratch_dir ), file=scratch_file )
# 		print( "sleep 1", file=scratch_file )
# 		print( file=scratch_file )
