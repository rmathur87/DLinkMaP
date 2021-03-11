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

#### Changing number of permutations

- change `n_runs` in `parallelizer.py` to be the number of permutations required.
