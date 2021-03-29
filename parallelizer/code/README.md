# Paralellizer

This folder parallelizes the script.

### `paralellizer.py`

- This script paralellizes the code.
- It creates multiple mini-scripts that subset the run into 10 runs of 100 permutations each.
	- Between each permutation, remember to delete the `lmod` database from ASC using: `rm /apps/db/lmod-ualcpr.db`. This resets the database so all imports will work properly (most of the time).
- It is only intended to be used once.

### `fail_run_finder.py`

- This script finds runs that failed due to failing to load the module.
- It does this by looking at the sizes of the p-value and logLike files.
	- If the size is 0 bytes, then it failed to run.
	- If the size is not 0 bytes, it ran.
- It then finds the run number, and creates a shell file with only the runs that failed called `rerun.sh`.

### `data_copy.py`

- Copies the p-value and the logLike data into the `out_data` directory for storage.
- It renames the files as well to reflect the permutation number.
- These can then be parsed by other programs since they are in a single place.

### `percentiler.py`

- Finds the 100th percentile and 95th percentile p-values are creates the data in the `percentile_data` directory.


## Run in following order:
1. HOME: fail_run_finder
2. HOME: copy_to_scratch.sh
3. HOME: rerun.sh
4. HOME: data_copy
