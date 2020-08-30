# `plot`

`plot` is python program for command-line plotting numeric data in
files or stdin and is convenient to use with unix pipes.  A
command-line plotting tool, rather than 'ecosystems' like R, jupyter,
or gnuplot, is something which is missing in standard shell-based
tools.  `plot` is an attempt to fill this void.

The source also includes functions callable from user applications
after `import plot`.

## Supported inputs

CSV or whitespace-separated utf-8 text.  The input is expected in the
form a dataframe with a header line with one or more field names
followed by data rows with as many numeric fields as there are field
names.  `plot` will plot any number of 'y' columns vs one 'x' column
in the dataframe, and more (see below).  Zero-based column numbers to
plot are supplied by user.

## Supported graphic backends

* `gnuplot`. Requires a [gnuplot5](http://www.gnuplot.info/docs_5.0/gnuplot.pdf) installation.  Generates cleaner graphics in a popup window without blocking command line or calling application.
* `pyplot`.  Uses matplotlib.pyplot supplied with python istallation.

## Supported plot types

* 2D graphs for one or more y(x)
* Scatter plots
* 1D histogram of a column in a dataframe
* xy-histograms of one or more y(x)

# xy-histogram

is a useful exploratory data analysis (EDA) tool for visual detection
of relationships, including nonlinear, among data fields.
xy-histogram is the primary reason for writing the `plot` tool.

An xy-histogram is a refinement of scatter plot and offers better
visualization of high-noise data.  Given data arrays $x_i$, $y_i$,
$i<N$, the bins $B_j$, $j<K\ll N$, are chosen either uniformly in $x$
or to contain equal count (weight) samples, xy-histogram displays
mean, and optionally standfard deviation of $\{y_i: x_i\in B_j\}$ vs
the $x$-bin $B_j$.  `plot` supports weighted xy-histograms by using a
weight data column for computing the mean and standard deviiation.

Command lines creating plots below demonstrate different view of the
same self-generated synthetic data:

### data:

`plot -t | head`
idx,foo,bar,boo,baz,wt
0,-2.7252320310513154,0.7515473227059941,0.22724419051121147,0.5009288010301421,716.1529207087109
1,-1.4503227442225184,0.9606778632248073,1.344122923318801,-1.333776384728359,605.2745872681331
2,-0.7418797345577915,-0.3479186065443103,-0.25589206096413675,0.9414321759657861,861.1656707423165
3,-0.815781602897208,0.5577812012125812,0.15225066771275708,-0.059347183822418326,624.2262162421391
4,-0.7435488913386459,1.2189909563056986,0.3467639364580495,-1.2059644351468646,620.2993697018401
5,2.522194888685291,-0.45901722833801906,-0.04625315311327499,1.9144753824473426,377.1847966081002
6,-1.7608271675310445,0.5679662325404847,1.2224698144847723,-2.2232181858855893,938.35563950581
7,-1.475958595171189,1.656643800125099,-1.1206972897423884,0.1774228658598056,565.4256143984067
8,0.6414610595883813,-1.871405226037516,1.431733155166075,-0.7128272011603145,452.54319181621423

### scatter plot:

`plot -t | plot 1-4 -zp`
<img src="xy-scatter.png" width="800" />

### xy-histogram with errorbars:

`plot -t | plot 1-4 -zeEHW 5`
@import "xy-hist-errorbars.png"

### xy-histogram bezier-smoothed:

`plot -t | plot 1-4 -zsEHW 5`
@import "xy-hist-smooth.png"


# Usage

`plot` and the companion script `unplot` have help/examples options:

```
$ plot -h
usage: plot [-h] [-T TITLE] [-z] [-w] [-l] [-n] [-c] [-s] [-d] [-p] [-H]
            [-B NBINS] [-W WTS_COL] [-e] [-E] [-g GREP_KEY] [-D DRIVER] [-v]
            [-t]
            [files_and_columns [files_and_columns ...]]

Plot columnar data in files or stdin. Default action: Print columns available for plotting

positional arguments:
  files_and_columns

optional arguments:
  -h, --help            show this help message and exit
  -T TITLE, --title TITLE
                        Use this window title [filename]
  -z, --xzeroaxis       Add horizontal zero line
  -w, --whitespace      Split input by whitespace [comma]
  -l, --list            Use sequence number for x data [first column]
  -n, --noheader        Indicate that data has no header.  Generate header names F0, F1, etc
  -c, --cumsum          Plot cumulative sums
  -s, --smooth          Plot bezier-smooth data
  -d, --diff            Plot differences
  -p, --points          Plot with points (e.g. for a scatter plot) [with lines]
  -H, --hist            Plot histogram: regular for single data column or xy histogram(s) for multiple columns (first treated as x)
  -B NBINS, --nbins NBINS
                        For histogram: use this many bins: sqrt: size**(1/2), qbrt: size**(1/3), or <int> [sqrt]
  -W WTS_COL, --wts_col WTS_COL
                        For histogram: use this column for weights
  -e, --yerr            For xy-historgram: plot with yerrorbars
  -E, --equal_wt        Use equal-weight histogram bins. Implies a density plot for x-histogram [equal-size]
  -g GREP_KEY, --grep_key GREP_KEY
                        Skip input lines without this word
  -D DRIVER, --driver DRIVER
                        Use this backend graphics driver: gnuplot or pyplot [gnuplot]
  -v, --verbose         Print gnuplot command line
  -t, --test_csv        Print a csv stream for testing

Examples:
  ls -l /usr/bin|grep -v total|plot -clw 4 # cumsum of file sizes
  plot -t | plot 0-4               # plot columns 1-4 vs column 0
  plot -t | plot 0-4 -zc           # cumulative plots
  plot -t | plot 0-4 -zs           # smoothed plots
  plot -t | plot 1-4 -p            # scatter plots
  plot -t | plot 1   -Hz           # histograms of column 1 data
  plot -t | plot 1-4 -Hz           # xy-histograms of columns 2-4 (y) vs column 1 (x)
  plot -t | plot 1-4 -HzW 5        # weighted xy-histograms
  plot -t | plot 1-4 -Hzs          # smoothed xy-histograms
  plot -t | plot 1-3 -eEHW 5 -B 10 # xy-histograms with equal-weight bins and errorbars
  unplot                           # kill all active gnuplot windows
```

# Version info

* Initial release V. 0.5 by Michael Isichenko
* Tested with python 3.7.6 and gnuplot 5.4

# TODO

* Date/time data support
* Handling of missing values
* Maybe: support surface/contour/heatmap plots while keeping a clean CLI
