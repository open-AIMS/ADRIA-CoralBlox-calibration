# AIMS HPC job submission guide

See setup instructions in the README.md.

Run `sbatch` from the **repository root** (`~/ADRIA-CoralBlox-calibration/`), not from
`scripts/` — all paths below are relative to that directory:

```bash
cd ~/ADRIA-CoralBlox-calibration/

sbatch --job-name=cb-calib --chdir ~/ADRIA-CoralBlox-calibration/ scripts/submit.slurm scripts/run_loc_calib.jl --output scripts/calib_output.log
```

A log of outputs from the optimizer will be saved in `scripts/calib_output.log`.
