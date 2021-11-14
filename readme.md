<img src="alpha-art.png" width="600" />

# `plot`

[plot](https://github.com/michael-isichenko/plot) is a python program for command-line plotting numeric data from files or stdin.  A simple command-line plotting tool usable with unix pipes, rather than interactive 'ecosystems' like R, jupyter or mathematica notebook, or gnuplot, is something which is unfortunately missing from standard shell-based tools.  `plot` is an attempt to fill this void.

The program name is admittedly too generic but is in line with the convention of most unix shell commands being (very) short.

The source also includes functions callable from user applications after `import plot`.

## Supported inputs

CSV or whitespace-separated utf-8 text.  Compressed files are handled transparently.  The input is normally expected in a dataframe-like format starting with a header line with one or more field names followed by data rows with as many numeric fields as there are field names.  `plot` will plot any number of 'y' columns vs one 'x' column in the dataframe, and more (see below).  Column names (or zero-based column numbers) to plot are supplied by the user on comannd line. Datetime 'x' data is supported.  Plotting from multiple files is also supported.  Most useful for a unix environment is plotting from standard input, so one can plot data passed through `grep` or any other program.  `plot -t` prints columns of random data usable for self-testing and demo plots.

## Supported graphics backends

* `pyplot`.  Uses `matplotlib.pyplot` supplied with python installation.  Requires PLOT_DRIVER=pyplot environment variable or `--driver pyplot` option.
* `gnuplot`. Requires a [gnuplot5](http://www.gnuplot.info/docs_5.0/gnuplot.pdf) installation.  This is plot's default driver.

`plot` can be used on unix-based systems with a graphics layer such as linux (X11) or OSX (Quartz, Cocoa).  Either graphics driver generates popup windows supported by forked background processes (`python` for `pyplot` and `gnuplot_qt` gnuplot). When no loner needed, the graphics windows need to be closed one by one by hand or all by the companion `unplot` script.  Both gnuplot and pyplot drivers support graphics-free plotting to a pdf file (`--output` option).

## Supported plot types

* Scatter plots
* 2D graphs for one or more y(x)
* 1D histograms of column(s) in a dataframe
* regressograms for one or more y(x)

2D graphs, histograms, and regressograms can be smoothed by local linear regression using user-supplied bandwidth (`--llr` option).  Additional data manipulation include splines supported by gnuplot (`--smooth` option), cumulation (`--cumsum`), differencing (`--diff`), subsetting (`--start`, `--end`), and printing basic statistics (`--stats`).

## Regressogram

Regressogram (`--rgram` option) is a useful exploratory data analysis (EDA) tool for visual detection of relationships, including nonlinear, among data fields. Regressogram was the primary reason for writing the `plot` tool.

[Regressogram](https://bookdown.org/egarpor/NP-UC3M/kre-i-reg.html) is a refinement of scatter plot and offers better visualization of high-noise data.  Given data arrays `(x[],y[])`, `x`-bins `B[]` are chosen either uniform in `x` or to contain samples of equal count (weight). An regressogram displays mean and optionally standard
deviation of the samples `y[i]` where `x[i]` fall in bin `B` vs the bin position.  This representation of `y(x)` is similar to a [KNN](https://en.wikipedia.org/wiki/K-nearest_neighbors_algorithm) ML data learner whose complexity increases with finer bins.  `plot` supports weighted regressograms by using a weight data column when computing the mean and the standard deviation.

Command lines and resulting plots below demonstrate different views of the same noisy data (generated in CSV format by the same tool):

### Scatter plot: `plot -t | plot 1-4 -p2 plt`
<img src="scatter.png" width="600" />

This view is not particularly telling.


### Regressogram with errorbars: `plot -t | plot 1-4 -rERW 5 -B 60`
<img src="rgram-errorbars.png" width="600" />

When error bars are large, the dependence of Mean(y) of x is not very visible.  Larger bins will generate smaller error bars at the price of fewer details of the y(x) dependency (bias-variance tradeoff).


### Smoothed regressogram: `plot -t | plot 1-4 -zERW 5 -B 60 -L 0.4`
<img src="rgram-smooth.png" width="600" />

Smoothing removes some noise and gives a better idea whether and how y depends on x.


### Smoothed histogram: `plot -t | plot 2-4 -HzL 0.5`
<img src="hgram-smooth.png" width="600" />

Multiple histograms can be displayed simultaneously.


## Usage

`plot` and the companion script `unplot` have help/examples options:

```
$ plot -h
usage: plot [-h] [-X] [-Y] [-T TITLE] [-q] [-z] [-w] [-Q] [-P] [-n] [-c]
            [-O SMOOTH] [-S STATS] [-s START] [-e END] [-x] [-L bandwidth]
            [-d] [-f FILTERY] [-p] [-D] [-l LOGSCALE] [-F] [-H] [-R] [-M] [-a]
            [-y YSCALE] [-m XTRIM] [-B NBINS] [-W WTS_COL] [-r] [-E]
            [-o OUTPUT] [-v] [-K] [-t] [-3] [-C] [-1 SEED] [-2 DRIVER]
            [-g GNUPLOT] [-N NTESTS]
            [files_and_columns ...]

Plot columnar data in files or stdin, normally separated by comma or whitespace and with a header line
The columns to plot can be specified by name, int index, inclusive range of ints, or a mix thereof
Column examples: foo | 1,0,5 | foo,bar,4-8,2 | 1,3- (trailing dash means through the end of available columns)
Default action: Print columns in the input available for plotting

positional arguments:
  files_and_columns

optional arguments:
  -h, --help            show this help message and exit
  -X, --examples        Print use examples
  -Y, --show_linetypes  Print help on supported line types
  -T TITLE, --title TITLE
                        Use this window title [filename]
  -q, --nolegend        Do not mark plotted line in legend
  -z, --xzeroaxis       Add horizontal zero line
  -w, --whitespace      Split input by whitespace [comma]
  -Q, --seqnum          Use sequence number for x data [first column]
  -P, --separate        Make plot for each data column separate [together]
  -n, --noheader        Indicate that data has no header.  Generate header names F0, F1, etc
  -c, --cumsum          Plot cumulative sums
  -O SMOOTH, --smooth SMOOTH
                        Apply gnuplot smoothing, supported types: s?bezier|(c|mc|ac)splines|(f|c)normal|kdensity <bwidth>
  -S STATS, --stats STATS
                        Add stats (before any cumsum or diff) to plots: c=count, m=mean, s=sdev, a=sharpe, d=Drawdown, D=DrawdownToSdev, A=all, P=statPrintStdout
  -s START, --start START
                        Skip data for x < start
  -e END, --end END     Skip data for x >= end
  -x, --noaxes          Generate plot without axes or border
  -L bandwidth, --llr bandwidth
                        Smooth data using local linear regression (LLR) over this bandwidth. Weird (buggy?) with --equal_wt option
  -d, --diff            Plot differences
  -f FILTERY, --filtery FILTERY
                        Filter (transform) y data through an expression to be applied via eval().  Examples: "y-y**3", "np.sin(y)"
  -p, --points          Plot with points (scatter plot) [with lines]
  -D, --dots            Plot with small dots (scatter plot) [with lines]
  -l LOGSCALE, --logscale LOGSCALE
                        Use log scale for x|y|xy.  Uppercase means select automatically based on data.  Useful with --separate option
  -F, --ols_fit         Add an OLS linear fit to each plotted y(x) line
  -H, --hgram           Plot histogram of indicated data column(s)
  -R, --rgram           Plot regressogram of data columns: first treated as x, rest as y(s)
  -M, --esm             For regressogram with yerr: plot error bars of sample mean (x 1/sqrt(Neff)) rather than yerr
  -a, --swapxy          For rgram: plot first column (x turned y) vs other indicated columns (yy turned xx) separately
  -y YSCALE, --yscale YSCALE
                        Scale y values for rgram and dataplot by this factor
  -m XTRIM, --xtrim XTRIM
                        For hgrm or rgram: trim this percentage of the x distribution tails on both sides [0.001]
  -B NBINS, --nbins NBINS
                        For histogram: use this many bins: sqrt: sqrt(size)/2, qbrt: size**(1/3), or <int> [sqrt]
  -W WTS_COL, --wts_col WTS_COL
                        For hgram or rgram: use this column for weights
  -r, --yerr            For regressogram: plot with yerrorbars
  -E, --equal_wt        Use equal-weight (histo|regresso)gram bins. Implies a density plot for histogram [equal-size].  Can look weird for highly peaked histogram
  -o OUTPUT, --output OUTPUT
                        Plot to this pdf file [qt window].  Seems to work better with -D plt option
  -v, --verbose         Print gnuplot command line
  -K, --keeptmp         Keep tmp file for gnuplot (helpful with -v option) [delete]
  -t, --test_csv        Print a csv stream for testing
  -3, --test_parametric
                        For test_csv: Print data for parametric plots
  -C, --test_clean      For test_csv: Print simple clean data for plots
  -1 SEED, --seed SEED  Use this seed for random test_csv data
  -2 DRIVER, --driver DRIVER
                        Use this backend graphics driver: gnuplot (gnu) or pyplot (plt) [$PLOT_DRIVER (if set) or gnuplot]
  -g GNUPLOT, --gnuplot GNUPLOT
                        Use this gnuplot binary [gnuplot]
  -N NTESTS, --ntests NTESTS
                        test_csv: print this many columns [4]


$ plot -X
# optionally: export PLOT_DRIVER=plt
  ls -l /usr/bin|grep -v total|plot -cnQw 4 # cumsum of file sizes (input without header)
  plot -t | plot 0-4                        # plot columns 1-4 vs column 0
  plot -t | plot 0-4 -cS msaD -s 500 -e 900 # cumulative plots for x in [500,900) with statistics in legend
  for i in 3 4 5; do plot -t1 $i > tmp$i.csv; done; plot 0,2-3 tmp{3,4,5}.csv -zc # cumulative plots from multiple files
  plot -t | plot 0-4 -zO bezier             # bezier-smoothed plots (gnuplot driver only)
  plot -t | plot 0,1,bar,baz,fee -zL 10     # same plots, mixed column spec, LLR smoothed
  plot -t | plot 1-4 -p                     # scatter plots
  plot -t | plot 1   -H                     # histograms of column 1 data
  plot -t | plot 2-4 -HL 0.5                # histograms of multiple columns, LLR-smoothed
  plot -t | plot 1-4 -Rz                    # regressograms of columns 2-4 (y) vs column 1 (x)
  plot -t | plot 1-4 -RzL 0.4               # regressograms smoothed with bandwidth 0.4
  plot -t | plot 1-4 -PRzW 5                # weighted regressograms in separate windows
  plot -t | plot 1-4 -rERW 5 -B 60 -L 0.4   # weighted regressograms with 60 equal-weight bins and errorbars
  plot -tN 500 -B 100|plot 1-500 -Rq        # spaghetti art
  plot -t3N 200|plot 1-199 -qxO bezier      # alpha art
  plot -X | head -16 | bash                 # run all of the above
  unplot                                    # kill all active gnuplot windows
```

## Implementation notes

`plot` loads input data into a pandas dataframe using `pandas.read_csv`.  The reason for using `pandas` is that parsing and tokenizing by `read_csv` is faster than in pure python.  Support for alignment/merging of data from multiple files is also a plus.  When the gnuplot backend is used, the data (after any smoothing or other processing) is saved to a `tempfile` passed to gnuplot executed via `subprocess`.  The gnuplot script can be viewed using `--verbose` option, and the temporary file can be kept (`--keeptmp`) for debugging.


## Installation

The script plot.py starts with the hashbang
`#!/usr/local/miniconda3/bin/python`.  If this is not where your python is, a shell script ~/bin/plot can be created to contain something like
```
#!/bin/bash
exec /your/python /path/to/plot.py $@
```
The script is best symlinked to ~/bin/plot or /usr/local/bin/plot.  All examples assume you have `plot` in your $PATH.


## Version info

* Initial release by Michael Isichenko
* Tested with python 3.7.6 and gnuplot 5.4 on mac OSX and linux


## TODO

* More robust Datetime data support
* Input filtering (a.k.a. 'where' clause) based on data values
* `eval`-based support for arithmetic and boolean expressions with data columns
* Maybe: support surface/contour/heatmap plots while keeping a clean CLI
* Maybe: add [GP](https://en.wikipedia.org/wiki/Gaussian_process) in addition to LLR and splines
