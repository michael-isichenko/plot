#!/usr/local/miniconda3/bin/python

# Copyright (c) 2020 Michael Isichenko. All rights reserved
# This code is licensed under the terms of the MIT license

import sys, os, subprocess, argparse, tempfile
import numpy as np

#
# Constants
#
VERSION            = '0.5 as of 2020-08-30'
DEFAULT_DRIVER     = 'gnuplot'
CHUNK              = 1024
HIST_DEFAULT_NBINS = 'sqrt'
HIST_EQUAL_WT      = 1 << 0
HIST_DENSITY       = 1 << 1
HIST_YERR          = 1 << 2

#
# Generic utilities
#
def MeanStd(x, wts=None):
    """Compute weighted mean and std"""
    if wts is None:
        return np.mean(x), np.std(x)
    assert len(x) == len(wts)
    if len(x) == 0:
        return np.nan, np.nan
    mean = np.average(x, weights=wts)
    var  = np.average((x - mean)**2, weights=wts)
    return mean, np.sqrt(max(var, 0.0))

def ParseColumns(columns, seq_list, hist=False):
    cols = []
    groups = columns.split(',')
    for group in groups:
        if '-' in group:
            beg, end = group.split('-')
            for c in range(int(beg), int(end) + 1):
                cols.append(c)
        else:
            cols.append(int(group))
    if hist:                         return cols
    elif seq_list or len(cols) == 1: return -1, cols
    else:                            return cols[0], cols[1:]

def PrintColumns(fname, opt):
    if os.path.isfile(fname):
        with open(fname) as input:
            line = input.readline()
    else:
        line = sys.stdin.readline()
    cols = line.rstrip().split(None if opt.whitespace else ',')
    for c, col in enumerate(cols):
        print('{:3d} {}'.format(c, col))

def PrintTestDataframe():
    np.random.seed(0) # for reproducible output
    N = 1000
    K = 4
    mean = np.zeros(K)
    cov = np.ones((K,K))
    for i in range(K):
        for j in range(i+1, K):
            cov[i,j] = 1.0*(j - i)/K - 0.5
            cov[j,i] = cov[i,j]
    print('idx,foo,bar,boo,baz,wt')
    X = np.random.multivariate_normal(mean, cov, size=N)
    for idx in range(N):
        print(f'{idx},{",".join([str(x) for x in X[idx]])},{np.random.uniform(100,1000)}')

def ReadDataframe(fnames, cols, opt):
    data = np.zeros((CHUNK, len(cols)))
    data_names = []
    wts = None if opt.wts_col is None else np.zeros(CHUNK)
    wts_name = None
    size = 0
    with (open(fnames[0], encoding='utf-8') if fnames[0] != '-' else sys.stdin) as input:
        for line in input:
            if line.startswith('#'):
                continue
            if opt.grep_key and len(header) > 0 and opt.grep_key not in line:
                continue
            tokens = line.rstrip().split() if opt.whitespace else line.rstrip().split(',')
            if len(data_names) == 0:
                if opt.noheader:
                    data_names = ['seqno' if c == -1 else f'F{c}' for c in cols]
                else:
                    data_names = [('seqno' if c == -1 else tokens[c]) for c in cols]
                if wts is not None:
                    wts_name = tokens[opt.wts_col]
                if not opt.noheader:
                    continue
            if size >= data.shape[0]:
                data.resize(size + CHUNK, len(cols))
                if wts is not None:
                    wts.resize(size + CHUNK)
            data[size] = [float(size if c == -1 else tokens[c]) for c in cols]
            if wts is not None:
                wts[size] = float(tokens[opt.wts_col])
            size += 1
    data.resize(size, data.shape[1])
    if wts is not None:
        wts.resize(size)
    return data, data_names, wts, wts_name

#
# Histograms
#
def _Parse_nbins(nbins, nx):
    if   isinstance(nbins, int):  return nbins
    elif isinstance(nbins, str):
        if   nbins == 'sqrt':     return int(np.ceil(np.sqrt(nx)))
        elif nbins == 'qbrt':     return int(np.ceil(nx**(1.0/3)))
        else:                     return int(nbins) # will raise ValueError for bad string
    raise ValueError(f'unsupported nbins spec "{nbins}"')

def ComputeBinEdges(x, wts, nbins, flags):
    if (flags & HIST_EQUAL_WT): # O(size*log(size))
        perm = np.argsort(x)
        sorted_x = x[perm]
        sorted_wts = wts[perm] if wts is not None else np.ones(len(x))
        cum_wts = np.cumsum(sorted_wts)
        total_wt = cum_wts[-1]
        wt_levels = np.arange(0, total_wt*0.999999, total_wt/nbins)
        assert len(wt_levels) == nbins
        indices = np.searchsorted(cum_wts, wt_levels, 'left')
        assert indices[0] == 0
        bin_edges = np.append(sorted_x[indices], sorted_x[-1])
    else: # O(size)
        xmin, xmax = np.min(x), np.max(x)
        step = 1.0*(xmax - xmin)/nbins
        bin_edges = np.arange(xmin, xmax + 1e-6*step, step)
        assert len(bin_edges) == nbins + 1
    return bin_edges

def ComputeHistogram(X, wts, nbins_spec, flags):
    """
    Bin by values of the first row of X.
    If X has more rows, they are interpreted as Ys for computing XY histgograms.
    If X.ndim == 1, a row of 1's is implicitly added for regular (X only) histogram.
    Supported nbins_specs:
      'sqrt':  sqrt(len(x))
      'qbrt':  len(x)**(1./3)
      25:      25
    flags is optional bitmask which can include these bits:
       HIST_EQUAL_WT - split data into equal-weigh bins (or compute densities for unweighted x-histogram)
       HIST_YERR     - return both y.mean and y.std for each bin
    Return: bin_edges, (ymeans, yerr), kind('x' or 'xy')
    """
    assert X.ndim <= 2, f'unsupported X shape {X.shape}'
    xyhist = X.ndim == 2 and X.shape[0] > 1
    x = X[0] if xyhist else X.flatten()
    nx = len(x)
    #print(f'XXX X: {X.shape} nx={nx}')
    nbins = _Parse_nbins(nbins_spec, nx)
    bin_edges = ComputeBinEdges(x, wts, nbins, flags)
    if xyhist: return bin_edges, ComputeXYHist(x, X[1:], wts, bin_edges, flags), 'xy'
    else:      return bin_edges, ComputeXHist(x,         wts, bin_edges, flags), 'x'

def ComputeXYHist(x, yy, wts, bin_edges, flags):
    nx = x.shape[0]
    ny = yy.shape[0]
    assert nx > 0 and ny > 0, f'invalid shapes: x: {x.shape}, yy: {yy.shape}'
    nbins = len(bin_edges) - 1
    ymeans = np.zeros((ny, nbins))
    yerr = np.zeros((ny, nbins)) if (flags & HIST_YERR) else None
    for bin in range(nbins):
        condition = (bin_edges[bin] <= x) & (x < bin_edges[bin + 1])
        for iy in range(ny):
            y = yy[iy]
            bin_y = y[condition]
            bin_wts = wts[condition] if wts is not None else None
            if (flags & HIST_YERR):
                ymeans[iy, bin], yerr[iy, bin] = MeanStd(bin_y, bin_wts)
            else:
                ymeans[iy, bin] = np.average(bin_y, weights=bin_wts) if len(bin_y) else np.nan
    return ymeans, yerr

def ComputeXHist(x, wts, bin_edges, flags): # use stdlib
    hist, _ = np.histogram(x, bins=bin_edges, weights=wts, density=bool(flags & HIST_EQUAL_WT))
    return hist.reshape(1, len(bin_edges) - 1), None

def RunHistogram(fnames, columns, opt):
    """Plot x-histogram for a single column of input or
       xy-histogram(s) for of y column(s) vs x column
    """
    assert columns, 'culumns required for a histogram'
    cols = ParseColumns(columns, opt.list, hist=True)
    assert cols[0] >= 0, f'invalid columns {cols}'
    data, names, wts, wts_name = ReadDataframe(fnames, cols, opt)
    flags = 0
    if opt.equal_wt:
        flags |= HIST_EQUAL_WT
    if opt.yerr:
        flags |= HIST_YERR
    bin_edges, (ymeans, yerr), kind = ComputeHistogram(data.transpose(), wts, opt.nbins, flags)
    if kind == 'x':
        if wts is None:    names.append('count')
        elif opt.equal_wt: names.append('density')
        else:              names.append('weight')
    if not opt.title:
        opt.title = 'xyhist' if len(cols) > 1 else 'hist'
    if wts is not None:
        opt.title += f' wts={wts_name}'
    ny, nbins = ymeans.shape
    if yerr is None:
        X = np.zeros((2*len(bin_edges), 1+ny)) # dataframe for a ladder function
        for bin, bin_edge in enumerate(bin_edges):
            X[2*bin,   0] = bin_edge
            X[2*bin+1, 0] = bin_edge
            for iy in range(ny):
                X[2*bin,   1+iy] = ymeans[iy, bin-1] if bin > 0     else np.nan
                X[2*bin+1, 1+iy] = ymeans[iy, bin  ] if bin < nbins else np.nan
        errorbars = None
    else:
        assert kind == 'xy'
        assert yerr.shape == (ny, nbins)
        X = np.zeros((nbins, 1 + ny)) # dataframe for centers of the bins
        for bin in range(nbins):
            X[bin, 0] = 0.5*(bin_edges[bin] + bin_edges[bin + 1])
        X[:,1:] = ymeans.transpose()
        errorbars = yerr.transpose()
    PlotDataframe(X, names, errorbars, **vars(opt))

#
# Generic plotting frontend
#
def PlotDataframe(X, names, errorbars, **kwarg):
    """
    Plot columns of dataframe X.  Default kwarg['cols']=range(X.shape[1]).
    Use x=X[:,cols[0]] and one or more Y=X[:,cols[1:]]
    The function writes a tmp datafile and runs a gnuplot script on it
    Supported kwargs:
    how:     ['with lines']
    cumsum:  [False] -- X is modified in-place
    diff:    [False] -- X is modified in-place
    smooth:  [False] -- use gnuplot 'smooth bezier' option
    title:   ['']
    driver:  'gnuplot' or 'pyplot' [plot.DEFAULT_DRIVER]
    gnuplot options:
      gnuplot: ['gnuplot'] -- gnuplot executable
      term:    ['qt']      -- gnuplot terminal
      verbose: [False]     -- print gnuplot command line to stdout
    """
    assert X.ndim == 2
    nrows = X.shape[0]
    ncols = X.shape[1]
    assert len(names) == ncols, f'mismatched names={names} vs ncols={ncols}'
    cols = kwarg.get('cols', range(ncols))
    assert max(cols) < ncols, f'cols={cols} out of range {ncols}'
    xcol = cols[0]
    ycols = cols[1:]
    if kwarg.get('cumsum', False):
        for ycol in ycols:
            X[:,ycol] = np.cumsum(X[:,ycol])
    elif kwarg.get('diff', False):
        for ycol in ycols:
            X[1:,ycol] = np.diff(X[:,ycol])
            X[0,ycol] = 0.0
    title   = kwarg.get('title', '')
    if title is None:
        title = ''
    if errorbars is not None:
        assert errorbars.shape == (nrows, len(ycols))
    cumsum  = kwarg.get('cumsum', False)
    diff    = kwarg.get('diff', False)
    smooth  = kwarg.get('smooth', False)
    verbose = kwarg.get('verbose', False)
    if cumsum: title += ' cumsum'
    if diff:   title += ' diff'
    if smooth: title += ' smooth'
    if   kwarg.get('driver', DEFAULT_DRIVER) == 'gnuplot':  ExecuteGnuplot(X, xcol, ycols, names, title, errorbars, kwarg)
    elif kwarg.get('driver', DEFAULT_DRIVER) == 'pyplot':   ExecutePyplot(X, xcol, ycols, names, title, errorbars, kwarg)
    else: assert False, 'no plotting driver set'

#
# Gnuplot driver
#
def ExecuteGnuplot(X, xcol, ycols, names, title, errorbars, kwarg):
    nrows = X.shape[0]
    ncols = X.shape[1]
    ny = ncols - 1
    verbose = kwarg.get('verbose', False)
    cols = [xcol] + list(ycols)
    script = []
    #script.append(f"set datafile separator ','")
    script.append(f"set xlabel '{names[xcol]}'")
    script.append(f"set key autotitle columnhead")
    script.append(f"set term {kwarg.get('term', 'qt')} title '{title}'")
    if kwarg.get('xzeroaxis', False):
        script.append('set xzeroaxis')
    how   = kwarg.get('how', 'with lines')
    extra = 'smooth bezier' if kwarg.get('smooth', False) else ''
    with tempfile.NamedTemporaryFile(delete=True) as tmp:
        line = ' '.join((names[c] for c in cols)) + '\n'
        tmp.write(line.encode())
        for row in range(nrows):
            tokens = [str(v) for v in X[row,:]]
            if errorbars is not None:
                tokens += (str(v) for v in errorbars[row,:])
            line = ' '.join(tokens) + '\n'
            tmp.write(line.encode())
        tmp.flush()
        items= []
        for yidx, ycol in enumerate(ycols):
            f = tmp.name if 0 == yidx else ''
            items.append(f"'{f}' using {xcol + 1}:{ycol + 1} {extra} {how}")
            if errorbars is not None:
                items.append(f"'{f}' using {xcol + 1}:{ycol + 1}:{ycol + 1 + ny} with errorbars notitle")
        all_items = ', '.join(items)
        script.append(f"plot {all_items}")
        cmd_line = [kwarg.get('gnuplot', 'gnuplot'), '-persist', '-e', '; '.join(script)]
        if verbose:
            print(f'# cmd_line: {cmd_line}')
        subprocess.run(cmd_line, shell=False, stderr=None if verbose else subprocess.DEVNULL)

def ExecutePyplot(X, xcol, ycols, names, title, errorbars, kwarg):
    import matplotlib.pyplot as plt
    nrows = X.shape[0]
    ncols = X.shape[1]
    ny = ncols - 1
    verbose = kwarg.get('verbose', False)
    cols = [xcol] + list(ycols)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    lw = 1 # default (2?) seems thick
    ls = 'solid'
    for iy in range(ny):
        x = X[:, 0]
        y = X[:, 1 + iy]
        label = names[1 + iy]
        if errorbars is not None:
            plt.errorbar(x, y, errorbars[:, iy], label=label, linewidth=lw, linestyle=ls)
        elif kwarg.get('points', False):
            plt.scatter(x, y, label=label)
        else:
            plt.plot(x, y, label=label, linewidth=lw, linestyle=ls)
    ax.set_xlabel(names[0])
    plt.legend()
    plt.show(block=True)

#
# 2D data plotter
#   
def RunDataplot(fnames, columns, opt):
    if columns is None:             xcol, ycols = 0,[1]
    elif isinstance(columns, list): xcol, *ycols = columns
    else:                           xcol, ycols = ParseColumns(columns, opt.list)
    cols = [xcol] + ycols
    data, names, wts, wts_name = ReadDataframe(fnames, cols, opt)
    PlotDataframe(data, names, None, **vars(opt))
    
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='Plot columnar data in files or stdin. Default action: Print columns available for plotting',
                                     epilog='Examples:\n'
                                     + '  ls -l /usr/bin|grep -v total|plot -clw 4 # cumsum of file sizes\n'
                                     + '  plot -t | plot 0-4               # plot columns 1-4 vs column 0\n'
                                     + '  plot -t | plot 0-4 -zc           # cumulative plots\n'
                                     + '  plot -t | plot 0-4 -zs           # smoothed plots\n'
                                     + '  plot -t | plot 1-4 -p            # scatter plots\n'
                                     + '  plot -t | plot 1   -Hz           # histograms of column 1 data\n'
                                     + '  plot -t | plot 1-4 -Hz           # xy-histograms of columns 2-4 (y) vs column 1 (x)\n'
                                     + '  plot -t | plot 1-4 -HzW 5        # weighted xy-histograms\n'
                                     + '  plot -t | plot 1-4 -Hzs          # smoothed xy-histograms\n'
                                     + '  plot -t | plot 1-3 -eEHW 5 -B 10 # xy-histograms with equal-weight bins and errorbars\n'
                                     + '  unplot                           # kill all active gnuplot windows'
    )
    parser.add_argument('-T', '--title',                           help='Use this window title [filename]')
    parser.add_argument('-z', '--xzeroaxis',  action='store_true', help='Add horizontal zero line')
    parser.add_argument('-w', '--whitespace', action='store_true', help='Split input by whitespace [comma]')
    parser.add_argument('-l', '--list',       action='store_true', help='Use sequence number for x data [first column]')
    parser.add_argument('-n', '--noheader',   action='store_true', help='Indicate that data has no header.  Generate header names F0, F1, etc')
    parser.add_argument('-c', '--cumsum',     action='store_true', help='Plot cumulative sums')
    parser.add_argument('-s', '--smooth',     action='store_true', help='Plot bezier-smooth data')
    parser.add_argument('-d', '--diff',       action='store_true', help='Plot differences')
    parser.add_argument('-p', '--points',     action='store_true', help='Plot with points (e.g. for a scatter plot) [with lines]')
    parser.add_argument('-H', '--hist',       action='store_true', help='Plot histogram: regular for single data column or xy histogram(s) for multiple columns (first treated as x)')
    parser.add_argument('-B', '--nbins',                           help=f'For histogram: use this many bins: sqrt: size**(1/2), qbrt: size**(1/3), or <int> [{HIST_DEFAULT_NBINS}]', default=HIST_DEFAULT_NBINS)
    parser.add_argument('-W', '--wts_col',    type=int,            help='For histogram: use this column for weights')
    parser.add_argument('-e', '--yerr',       action='store_true', help='For xy-historgram: plot with yerrorbars')
    parser.add_argument('-E', '--equal_wt',   action='store_true', help='Use equal-weight histogram bins. Implies a density plot for x-histogram [equal-size]')
    parser.add_argument('-g', '--grep_key',                        help='Skip input lines without this word')
   #parser.add_argument('-M', '--multikey_col', nargs=1,           help='Plot separate lines for for each value in this column (TODO)')
    parser.add_argument('-D', '--driver',                          help=f'Use this backend graphics driver: gnuplot or pyplot [{DEFAULT_DRIVER}]', default=DEFAULT_DRIVER)
    parser.add_argument('-v', '--verbose',    action='store_true', help='Print gnuplot command line')
    parser.add_argument('-t', '--test_csv',   action='store_true', help='Print a csv stream for testing')
    parser.add_argument('files_and_columns',  nargs='*')
    args = parser.parse_args()

    #assert not args.multikey_col, 'multikey_col not supported (todo)'
    if args.test_csv:
        PrintTestDataframe()
    else:
        args.how = 'with points' if args.points else 'with lines'
        fnames = []
        columns = None
        for arg in args.files_and_columns:
            if os.path.isfile(arg) or arg == '-':
                fnames.append(arg)
            else:
                assert columns is None, "duplicare columns specs not supported"
                columns = arg
        if not fnames:
            fnames.append('-')
        if columns is None:
            PrintColumns(fnames[0], args)
        elif args.hist:
            RunHistogram(fnames, columns, args)
        else:
            RunDataplot(fnames, columns, args)

if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        pass
