### How to RUN

```
python3 paralellizer.py
bash script.sh
```

### How it works

- Creates files and directory structure required for run on ASC
- Creates a script file that can then be used to start jobs
- Uses digits of pi as a random seed
- Intended to run 1000 jobs of 1 permutation each

### Some Models not running due to inability to load R

There are some models that do not run due to inability to load R.
- [ ] I will have to create a script to check which models were not run (by searching for the temp.csv files, and then rerun only those instead of rerunning everything).
	- [ ] I may have to also use the size of th two csv files (p-val and log_like) to check if the jobs have been run. If `size == 0`, means not run, and needs to be run again
- [ ] This script would also have to copy the files for each permutation elsewhere to save (for those that running was possible for).
- [ ] When making this script, use it to alter the shell file that is finally run, do not delete any data, just replace the `script.sh` shell with a smaller version for only those jobs that have not been run.

### Additional Scripts

There are additional scripts that perform cursory tasks, and they are explained in the README.md within the `code` directory
