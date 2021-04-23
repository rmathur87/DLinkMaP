"""
- creates multiple scripts for QTL analysis to run on ASC
"""

### FUNCTIONS ==================================================================

def split_pi( pi_file, n_runs ):
	"""
	- splits pi into n sizes to get "random" seeds for the QTL runs
	"""

	ret_lst_full = []

	with open( pi_file, "r" ) as infile:
		for line in infile:
			temp = line.strip()
			ret_lst_full = [(temp[i:i+5]) for i in range(0, len(temp), 5)]

	return( ret_lst_full[:n_runs] )


#
#    ,,
#    db                                                    mm
#                                                          MM
#  `7MM  `7MMpMMMb.pMMMb.  `7MMpdMAo.  ,pW"Wq.  `7Mb,od8 mmMMmm  ,pP"Ybd
#    MM    MM    MM    MM    MM   `Wb 6W'   `Wb   MM' "'   MM    8I   `"
#    MM    MM    MM    MM    MM    M8 8M     M8   MM       MM    `YMMMa.
#    MM    MM    MM    MM    MM   ,AP YA.   ,A9   MM       MM    L.   I8
#  .JMML..JMML  JMML  JMML.  MMbmmd'   `Ybmd9'  .JMML.     `Mbmo M9mmmP'
#                            MM
#                          .JMML.

import sys
from pprint import pprint as pprint
import os
import shutil
import time
from pathlib import Path 				# to import files
# from tqdm import tqdm as tqdm

#
#                                 ,,            ,,          ,,
#                                 db           *MM        `7MM
#                                               MM          MM
#  `7M'   `MF' ,6"Yb.  `7Mb,od8 `7MM   ,6"Yb.   MM,dMMb.    MM   .gP"Ya  ,pP"Ybd
#    VA   ,V  8)   MM    MM' "'   MM  8)   MM   MM    `Mb   MM  ,M'   Yb 8I   `"
#     VA ,V    ,pm9MM    MM       MM   ,pm9MM   MM     M8   MM  8M////// `YMMMa.
#      VVV    8M   MM    MM       MM  8M   MM   MM.   ,M9   MM  YM.    , L.   I8
#       W     `Moo9^Yo..JMML.   .JMML.`Moo9^Yo. P^YbmdP'  .JMML. `Mbmmd' M9mmmP'
#

# n_runs = 10	00
n_runs = 10
run_scripts_path = "/scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts"
mini_script_path = "/scratch/ualcpr/QTL/DLinkMaP/parallelizer/code/script.sh"

#         ,,
#       `7MM             mm
#         MM             MM
#    ,M""bMM   ,6"Yb.  mmMMmm   ,6"Yb.
#  ,AP    MM  8)   MM    MM    8)   MM
#  8MI    MM   ,pm9MM    MM     ,pm9MM
#  `Mb    MM  8M   MM    MM    8M   MM
#   `Wbmd"MML.`Moo9^Yo.  `Mbmo `Moo9^Yo.
#

pi_file = "/scratch/ualcpr/QTL/DLinkMaP/parallelizer/stock_data/pi.digits"
# pi_file = "/Users/rele.c/Downloads/DLinkMaP/parallelizer/stock_data/pi.digits"
params_stock = "/scratch/ualcpr/QTL/DLinkMaP/parallelizer/stock_data/parameters.csv"
script_stock = "/scratch/ualcpr/QTL/DLinkMaP/parallelizer/stock_data/script_example.sh"

#
#                            ,,
#                          `7MM
#                            MM
#   ,p6"bo   ,pW"Wq.    ,M""bMM   .gP"Ya
#  6M'  OO  6W'   `Wb ,AP    MM  ,M'   Yb
#  8M       8M     M8 8MI    MM  8M//////
#  YM.    , YA.   ,A9 `Mb    MM  YM.    ,
#   YMbmd'   `Ybmd9'   `Wbmd"MML. `Mbmmd'
#
#

pi_seeds = split_pi( pi_file, n_runs )

# print( "".center( 100, "=" ) )
# print( " BEGINNING RUN FOR {:04} RUNS ".format( n_runs ).center( 100, "=" ) )
# print( "".center( 100, "=" ) )
# print()


try:
	shutil.rmtree( run_scripts_path )
	time.sleep( 1 )
except OSError as e:
	print("Error: %s : %s" % (run_scripts_path, e.strerror))

try:
	os.mkdir( run_scripts_path )
except Exception as e:
	print( e )
	sys.exit( "Could not create {} directory.".format( run_scripts_path ) )

os.system( "rm {0}".format( mini_script_path ) )
os.system( "touch {0}".format( mini_script_path ) )

for i in range( n_runs ):
# for i in tqdm(range( n_runs ), ascii=True, desc="Creating Run Files"): # comment this out if it does not work

	# print( "SEED = {0} ".format( pi_seeds[i] ).ljust( 50, "=" ), end="" )
	# print( " RUN {0:04}".format( i+1 ).rjust( 50, "=" ) )

	# print( "Run {0:04}\tSeed = {1}".format( i+1, pi_seeds[i] ) )

	# make directory
	dir_name = run_scripts_path + "/run_{:04}".format(i+1)
	try:
		os.mkdir( dir_name )
	except Exception as e:
		print( "Could not create directory for run {:04}.".format( i+1 ) )
		sys.exit()

	# print( "{:04} - Copying femaleWt ".format( i+1 ).ljust( 30, "." ), end = "" )
	# shutil.copytree(
	# 	"/scratch/ualcpr/QTL/DLinkMaP/results/femaleWt",
	# 	dir_name + "/femaleWt"
	# )
	# print( " DONE" )
	# time.sleep( 0.1 )

	# print( "{:04} - Copying maleWt ".format( i+1 ).ljust( 30, "." ), end = "" )
	shutil.copytree(
		"/scratch/ualcpr/QTL/DLinkMaP/results/TG",
		dir_name + "/TG"
	)
	# print( " DONE" )
	time.sleep( 0.1 )

	# print( "{:04} - Copying Trehalose ".format( i+1 ).ljust( 30, "." ), end = "" )
	# shutil.copytree(
	# 	"/scratch/ualcpr/QTL/DLinkMaP/results/Trehalose",
	# 	dir_name + "/Trehalose"
	# )
	# print( " DONE" )
	# time.sleep( 0.1 )

	# copy parameters
	params_file = Path( params_stock ).read_text().split("\n")
	param_file_for_run = dir_name + "/params_{:04}.csv".format(i+1)

	with open( param_file_for_run, "w" ) as pfile:
		for line in params_file:
			add_line = line
			if line.split(",")[0] == "numPermutations":
				temp_line = [
					line.split(",")[0],
					line.split(",")[1],
					"1"
				]
				add_line = ",".join(temp_line)
			if line.split(",")[0] == "seed":
				temp_line = [
					line.split(",")[0],
					line.split(",")[1],
					str( pi_seeds[i] )
				]
				add_line = ",".join(temp_line)
			if line.split(",")[0] == "commDir":
				temp_line = [
					line.split(",")[0],
					line.split(",")[1],
					"/scratch/ualcpr/QTL/DLinkMaP/mapping/"
				]
				add_line = ",".join(temp_line)
			if line.split(",")[0] == "nullDir":
				temp_line = [
					line.split(",")[0],
					line.split(",")[1],
					"/scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{:04}/TG".format( i+1 )
				]
				add_line = ",".join(temp_line)
			if line.split(",")[0] == "outDir":
				temp_line = [
					line.split(",")[0],
					line.split(",")[1],
					"/scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{:04}/TG".format( i+1 )
				]
				add_line = ",".join(temp_line)
			if line.split(",")[0] == "outDir":
				temp_line = [
					line.split(",")[0],
					line.split(",")[1],
					"/scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{:04}/TG".format( i+1 )
				]
				add_line = ",".join(temp_line)
			if line.split(",")[0] == "phenotype":
				temp_line = [
					line.split(",")[0],
					line.split(",")[1],
					"TG"
				]
				add_line = ",".join(temp_line)
			if line.split(",")[0] == "fileName":
				temp_line = [
					line.split(",")[0],
					line.split(",")[1],
					"/scratch/ualcpr/QTL/DLinkMaP/phenotypes/TG/TG_data_updated_5_8_13_formatted.csv"
				]
				add_line = ",".join(temp_line)
			print( add_line, file=pfile )
	# print( "params_{:04}.csv created".format(i+1).rjust( 60, "." ) )
	time.sleep( 0.1 )

	# copy_shell
	script_file = Path( script_stock ).read_text().split("\n")
	script_file_for_run = dir_name + "/perm{:04}.sh".format(i+1)
	copy_home_script = dir_name + "/copy{:04}.sh".format(i+1)

	with open( script_file_for_run, "w" ) as sfile:

		print( "#!/bin/bash", file=sfile )
		print( file=sfile )
		print( "# load R", file=sfile )
		print( "module purge", file=sfile )
		print( "sleep 5", file=sfile )
		print( "source /opt/asn/etc/asn-bash-profiles-special/modules.sh", file=sfile )
		print( "sleep 5", file=sfile )
		print( "module load R/4.0.5_scratch", file=sfile )
		print( "sleep 5", file=sfile )
		print( file=sfile )
		print( "# R CMD INSTALL ~/QTL/DLinkMap/DSPRqtl_2.0-5.tar.gz", file=sfile )
		print( "rm /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/*RData".format( i+1), file=sfile )
		print( "sleep 5", file=sfile )
		print( "Rscript /scratch/ualcpr/QTL/DLinkMaP/mapping/NullSetUpBla_TG_Trehalose_Weights_permutation.R /scratch/ualcpr/QTL/DLinkMaP/mapping/ {1} /scratch/ualcpr/QTL/DLinkMaP/phenotypes/TG/TG_data_updated_5_8_13_formatted.csv tg /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG 1".format( i+1, pi_seeds[i] ), file=sfile ) # I might have to change the `tg` in this line to `TG` if it fails to run
		print( "mv /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/Data_tg_permutation1.csv /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/readable_data_{0:04}.txt".format( i+1), file=sfile ) # may have to change `tg` here to `TG` depending on previous line
		print( "mv /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/*RData /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/Null_Male_tg_avgbyvial.RData".format( i+1), file=sfile ) # may have to change `tg` here to `TG` depending on previous lines
		print( file=sfile )
		print( "# Run R Script", file=sfile )
		print( "Rscript /scratch/ualcpr/QTL/DLinkMaP/mapping/MAP_general.R /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/params_{0:04}.csv".format( i+1 ), file=sfile )
		print( file=sfile )
		print( "# run copy script", file=sfile )
		print( "cd {0}".format( dir_name ), file=sfile )
		print( "chmod +x {0}".format( copy_home_script ), file=sfile )
		print( "bash copy{:04}.sh".format(i+1), file=sfile )
		print( file=sfile )

	with open( copy_home_script, "w" ) as copyfile:
		# need to run uniq fromwithin the data to get uniqed lines
		print( "# uniqing {0}".format( i+1 ), file=copyfile )
		print( "cd {0}".format( dir_name ), file=copyfile )
		print( "cp /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/p-value_Male_avgbyvial.csv /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/temp-pval.csv".format( i+1 ), file=copyfile )
		print( "cat /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/temp-pval.csv | uniq > /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/p-value_Male_avgbyvial.csv".format( i+1 ), file=copyfile )
		print( "cp /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/p-value_Male_avgbyvial.csv /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/p-value_Male_avgbyvial.csv".format( i+1 ), file=copyfile )
		print( "cp /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/logLike_Male_avgbyvial.csv /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/temp-log_like.csv".format( i+1 ), file=copyfile )
		print( "cat /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/temp-log_like.csv | uniq > /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/logLike_Male_avgbyvial.csv".format( i+1 ), file=copyfile )
		print( "cp /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/logLike_Male_avgbyvial.csv /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/TG/logLike_Male_avgbyvial.csv".format( i+1 ), file=copyfile )

	# print( "perm{:04}.sh created".format(i+1).rjust( 60, "." ) )
	time.sleep( 0.1 ) # just in case we need it

	# os.system( "pwd" )
	# os.system( "ls -ltrh" )
	#
	# print( "Ready to run. Sleep 1".rjust( 90, "." ) )
	# time.sleep( 1 )
	#
	# print( "Changing Directory. Sleep 1" )
	# os.system( "cd ./../run_scripts/run_{:04}".format( i+1 ) )
	# time.sleep(1)
	#
	# os.system( "pwd" )
	# os.system( "ls -ltrh" )
	#
	# print( "Assigning Run Permission. Sleep 1".rjust( 90, "." ) )
	# os.system( "chmod +x perm{:04}.sh".format( i+1 ) )
	# time.sleep( 1 )
	#
	# print( "Attempting Run".rjust( 90, "." ) )
	# run_command = "run_script perm{:04}.sh".format( i+1 )
	# print( "Running: ", run_command )
	# os.system( run_command )
	# print( "Finished Running. Sleep 5".rjust( 90, "." ) )
	# time.sleep( 5 )
	#
	# print( "Changing Directory back. Sleep 5" )
	# os.system( "cd /scratch/ualcpr/QTL/Downloads/DLinkMaP/parallelizer/code" )
	# time.sleep( 5 )
	#
	# os.system( "pwd" )
	# os.system( "ls -ltrh" )

	with open( mini_script_path, "a" ) as outfile:
		print( "cd /scratch/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{:04}".format( i+1 ), file=outfile )
		print( "chmod +x perm{:04}.sh".format( i+1 ), file=outfile )
		print( "run_script perm{:04}.sh".format( i+1 ), file=outfile )
		print( "sleep 10", file=outfile )
		print( file=outfile )

	print( "{0:>04} / {1}".format( i+1, n_runs ).center( 100, " " ), end="\r" ) # remove if enabling tqdm

	# print()
	# print()
	# print()

os.system( "chmod +x {0}".format( mini_script_path ) )

# to split bash file
# lines_per_file = 500
# smallfile = None
# with open( mini_script_path, "r" ) as bigfile:
# 	for lineno, line in enumerate(bigfile):
# 		if lineno % lines_per_file == 0:
# 			if smallfile:
# 				smallfile.close()
# 			small_filename = 'mini_script-{}.sh'.format(lineno + lines_per_file)
# 			smallfile = open(small_filename, "w")
# 		smallfile.write(line)
# 	if smallfile:
# 		smallfile.close()

# print( "".center( 100, "=" ) )
# print( " FINISHED RUNNING ALL ".center( 100, "=" ) )
# print( "".center( 100, "=" ) )

# os.system( "tree ../" )
# sys.exit( "Creation" )


# Rscript /scratch/ualcpr/QTL/DLinkMaP/mapping/NullSetUpBla_TG_Trehalose_Weights_permutation.R /scratch/ualcpr/QTL/DLinkMaP/mapping/ 12345 /scratch/ualcpr/QTL/DLinkMaP/phenotypes/MaleWt/male_Weight_Formatted_AvgByVial.txt Weight /scratch/ualcpr/QTL/DLinkMaP/mapping 1
