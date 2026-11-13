# SiPM test analysis

Compiled ROOT analysis of Janus spectroscopy CSV files.

Requirements: CMake 3.16+, a C++17 compiler, and ROOT with Spectrum and Minuit.
The ROOT environment must provide `root-config` on `PATH`.

## Run

```bash
./run_analysis.sh
# Or specify an input CSV and an optional output directory:
./run_analysis.sh /path/to/Run94_list.csv /path/to/plots
```

The run script calls `build.sh` automatically. Its default input is
`data/Run94_list.csv`; output goes to `plots/`. Build and analysis messages
overwrite `log/<input-stem>.log` on each run.

The program reads the fixed Janus CSV columns and plots TStamp, PHA_LG,
and PHA_HG. ChannelMask and DataType are ignored. Peak candidates must have
more than 100 counts in their raw histogram bin. Gaussian fits use TMinuit,
with half-height pre-fits followed by final fits within the pre-fit mean
plus or minus two pre-fit standard deviations. Only converged fits with
reliable covariance matrices are returned.

Mean and sigma graphs use the parameter uncertainties as vertical error
bars. The HG mean graph shows the first six retained peaks and fits points
1 through 5 with a weighted straight line. Peak indices are sequential;
missing or duplicate candidates require care when interpreting them as p.e.
counts. Numerical convergence alone does not imply good fit quality.

CSV data, build artifacts, plots, and logs are excluded from version control.
