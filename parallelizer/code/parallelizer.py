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

# n_runs = 1000
n_runs = 12
run_scripts_path = "/home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts"
mini_script_path = "/home/ualcpr/QTL/DLinkMaP/parallelizer/code/script.sh"

#         ,,
#       `7MM             mm
#         MM             MM
#    ,M""bMM   ,6"Yb.  mmMMmm   ,6"Yb.
#  ,AP    MM  8)   MM    MM    8)   MM
#  8MI    MM   ,pm9MM    MM     ,pm9MM
#  `Mb    MM  8M   MM    MM    8M   MM
#   `Wbmd"MML.`Moo9^Yo.  `Mbmo `Moo9^Yo.
#

pi_file = "/home/ualcpr/QTL/DLinkMaP/parallelizer/stock_data/pi.digits"
# pi_file = "/Users/rele.c/Downloads/DLinkMaP/parallelizer/stock_data/pi.digits"
params_stock = "/home/ualcpr/QTL/DLinkMaP/parallelizer/stock_data/parameters.csv"
script_stock = "/home/ualcpr/QTL/DLinkMaP/parallelizer/stock_data/script_example.sh"

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

print( "".center( 100, "=" ) )
print( " BEGINNING RUN FOR {:04} RUNS ".format( n_runs ).center( 100, "=" ) )
print( "".center( 100, "=" ) )
print()


try:
    shutil.rmtree( run_scripts_path )
except OSError as e:
    print("Error: %s : %s" % (run_scripts_path, e.strerror))

try:
	os.mkdir( run_scripts_path )
except Exception as e:
	print( "Could not create {} directory.".format( run_scripts_path ) )

os.system( "rm {0}".format( mini_script_path ) )
os.system( "touch {0}".format( mini_script_path ) )

for i in range( n_runs ):

	print( "SEED = {0} ".format( pi_seeds[i] ).ljust( 50, "=" ), end="" )
	print( " RUN {0:04}".format( i+1 ).rjust( 50, "=" ) )

	# print( "Run {0:04}\tSeed = {1}".format( i+1, pi_seeds[i] ) )

	# make directory
	dir_name = run_scripts_path + "/run_{:04}".format(i+1)
	try:
		os.mkdir( dir_name )
	except Exception as e:
		print( "Could not create directory for run {:04}.".format( i+1 ) )
		continue

	print( "{:04} - Copying femaleWt ".format( i+1 ).ljust( 30, "." ), end = "" )
	shutil.copytree(
		"/home/ualcpr/QTL/DLinkMaP/results/femaleWt",
		dir_name + "/femaleWt"
	)
	print( " DONE" )
	time.sleep( 0.1 )

	print( "{:04} - Copying maleWt ".format( i+1 ).ljust( 30, "." ), end = "" )
	shutil.copytree(
		"/home/ualcpr/QTL/DLinkMaP/results/maleWt",
		dir_name + "/maleWt"
	)
	print( " DONE" )
	time.sleep( 0.1 )

	print( "{:04} - Copying Trehalose ".format( i+1 ).ljust( 30, "." ), end = "" )
	shutil.copytree(
		"/home/ualcpr/QTL/DLinkMaP/results/Trehalose",
		dir_name + "/Trehalose"
	)
	print( " DONE" )
	time.sleep( 0.1 )

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
					"/home/ualcpr/QTL/DLinkMaP/mapping/"
				]
				add_line = ",".join(temp_line)
			if line.split(",")[0] == "nullDir":
				temp_line = [
					line.split(",")[0],
					line.split(",")[1],
					"/home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{:04}/maleWt".format( i+1 )
				]
				add_line = ",".join(temp_line)
			if line.split(",")[0] == "outDir":
				temp_line = [
					line.split(",")[0],
					line.split(",")[1],
					"/home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{:04}/maleWt".format( i+1 )
				]
				add_line = ",".join(temp_line)
			if line.split(",")[0] == "fileName":
				temp_line = [
					line.split(",")[0],
					line.split(",")[1],
					"/home/ualcpr/QTL/DLinkMaP/phenotypes/MaleWt/male_Weight_Formatted_AvgByVial.txt"
				]
				add_line = ",".join(temp_line)
			print( add_line, file=pfile )
	print( "params_{:04}.csv created".format(i+1).rjust( 60, "." ) )
	time.sleep( 0.5 )

	# copy_shell
	script_file = Path( script_stock ).read_text().split("\n")
	script_file_for_run = dir_name + "/perm{:04}_script.sh".format(i+1)

	with open( script_file_for_run, "w" ) as sfile:
		for line in script_file:
			add_line = line.strip()
			if "Rscript" in line:
				temp = [
					line.split()[0],
					line.split()[1],
					"/home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{0:04}/params_{0:04}.csv".format( i+1 )
				]
				add_line = " ".join( temp )
			print( add_line, file=sfile )
	print( "perm{:04}_script.sh created".format(i+1).rjust( 60, "." ) )

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
	# os.system( "chmod +x perm{:04}_script.sh".format( i+1 ) )
	# time.sleep( 1 )
	#
	# print( "Attempting Run".rjust( 90, "." ) )
	# run_command = "run_script perm{:04}_script.sh".format( i+1 )
	# print( "Running: ", run_command )
	# os.system( run_command )
	# print( "Finished Running. Sleep 5".rjust( 90, "." ) )
	# time.sleep( 5 )
	#
	# print( "Changing Directory back. Sleep 5" )
	# os.system( "cd /home/ualcpr/QTL/Downloads/DLinkMaP/parallelizer/code" )
	# time.sleep( 5 )
	#
	# os.system( "pwd" )
	# os.system( "ls -ltrh" )

	with open( mini_script_path, "a" ) as outfile:
		print( "cd /home/ualcpr/QTL/DLinkMaP/parallelizer/run_scripts/run_{:04}".format( i+1 ), file=outfile )
		print( "chmod +x perm{:04}_script.sh".format( i+1 ), file=outfile )
		print( "run_script perm{:04}_script.sh".format( i+1 ), file=outfile )
		print( "sleep 10", file=outfile )
		print()

	print()
	print()
	print()


print( "".center( 100, "=" ) )
print( " FINISHED RUNNING ALL ".center( 100, "=" ) )
print( "".center( 100, "=" ) )

# os.system( "tree ../" )
# sys.exit( "Creation" )
