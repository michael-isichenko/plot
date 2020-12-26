#!/usr/local/miniconda3/bin/python

# Copyright (c) 2020 Michael Isichenko. All rights reserved
# This code is licensed under the terms of the MIT license

import sys, os, subprocess, argparse, tempfile, re
from itertools import permutations 
import numpy as np
import pandas as pd

#
# Constants
#
VERSION            = '0.5 as of 2020-08-30'
DEFAULT_DRIVER     = 'gnuplot'
CHUNK              = 1024
DEFAULT_NTESTS     = 4
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

def ParseColumns(cols_spec):
    cols = []
    groups = cols_spec.split(',')
    for group in groups:
        if '-' in group:
            beg, end = group.split('-')
            for c in range(int(beg), int(end) + 1):
                cols.append(c)
        else:
            cols.append(int(group))
    return cols

def PrintColumns(fname, opt):
    if os.path.isfile(fname):
        with open(fname) as input:
            line = input.readline()
    else:
        line = sys.stdin.readline()
    cols = line.rstrip().split(None if opt.whitespace else ',')
    for c, col in enumerate(cols):
        print('{:3d} {}'.format(c, col))

def ToeplitzCov(K, mincor):
    def Elem(i, j): return np.exp(np.log(mincor)*np.abs(i - j)/K)
    return np.fromfunction(Elem, (K,K))

def NaiveCov(K):
    cov  = np.ones((K,K))
    for i in range(K):
        for j in range(i+1, K):
            cov[i,j] = 1.0*(j - i)/K - 0.5
            cov[j,i] = cov[i,j]
    return cov

def Profiles(N, K):
    def Profile(n, k): return 0.1*((0.7 - n/N)**4 + 5*(0.5 - n/N)**5)*np.sin(1.0*k/K)
    return np.fromfunction(Profile, (N, K))

def PrintTestDataframe(opt):
    np.random.seed(opt.seed) # for reproducible output
    N = 50 if opt.test_parametric else 1000
    K = max(opt.ntests, 4)
    mean = np.zeros(K)
    #cov = NaiveCov(K)
    cov = ToeplitzCov(K, 0.01)
    header = 'idx'
    if opt.ntests < 9:
        short_names = ('foo', 'bar', 'baz', 'fee', 'dee', 'doo', 'zee', 'zoo')
        for i in range(opt.ntests):
            header += f',{short_names[i]}'
    else:
        for i, p in enumerate(permutations('abcdefghijklmnopqrstuvwxyz', 3)):
            header += f',{"".join(p)}'
            if i == opt.ntests - 1:
                break
    header += ',wt'
    print(header)
    X  = np.random.multivariate_normal(mean, cov, size=N) # N K-dimensional points drawn from N(0, cov)
    if opt.test_parametric:
        tt = np.linspace(0, 2*np.pi, N, endpoint=False)
        X[:, 0] = np.cos(tt) - 0.1*tt
        for k in range(1, K):
            #X[:, k] = 0.01*X[:, k] + np.cos(np.sqrt(4.5*np.pi*tt))
            X[:, k] = 0.08*X[:, k] + (1 + np.exp(-tt))*np.cos(1.55*(1 + 0.05*(1.0*k/K - 0.5))*tt)
        for idx, t in enumerate(tt):
            wt = np.random.uniform(10,100)
            print(f'{t},{",".join([str(x) for x in X[idx]])},{wt}')
    else:
        X += Profiles(N, K)
        for idx in range(N):
            wt = np.random.uniform(10,100)
            print(f'{idx},{",".join([str(x) for x in X[idx]])},{wt}')

def ReadDataframe(fnames, cols, merge_by_x, opt):
    """
    Read data from one or more compatible files into pandas dataframe.
    The reasons for using pandas are: (1) fast read_csv() and (2) merge() functionality for multiple files/
    For multiple files, generate a wide dataframe with columns like <fname>:<column_name>
    Any weight columns must be included in the dataframe
    """
    sep = r'\s+' if opt.whitespace else ','
    if opt.noheader:
        header = None                 
        names  = [f'C{c}' for c in cols]
    else:
         header =  0
         names = None
    if   opt.seqnum: index_col = False # generate index 0,1,...
    elif merge_by_x: index_col = 0     # index contains x data index_col is per usecols, not the whole input
    else:            index_col = None  # no index
    for i, fname in enumerate(fnames):
        #print(f'XXX1 fname={fname} sep="{sep}" names={names} usecols={cols} index_col={index_col}')
        file_df = pd.read_csv(sys.stdin if fname == '-' else fname, sep=sep, header=header, names=names, usecols=cols, index_col=index_col)
        if len(fnames) > 1:
            file_df.columns = [fname + ':' + c for c in file_df.columns]
        if i == 0:
            df = file_df
        else:
            df = pd.merge(df, file_df, how='outer', left_index=True, right_index=True)
    #print(f'XXX2 df columns={df.columns.tolist()}:\ndf.index:\n{df.index}\ndf:\n{df}')
    return df

def ReadDataframeNoPandas(fname, cols, opt):
    data = np.zeros((CHUNK, len(cols)))
    data_names = []
    wts = None if opt.wts_col is None else np.zeros(CHUNK)
    wts_name = None
    size = 0
    with (open(fname, encoding='utf-8') if fname != '-' else sys.stdin) as input:
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
    #print(f'XXX3 x: {type(x)} x={x} wts={wts}')
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
        #print(f'XXX4 x={x} perm={perm} sorted_x={sorted_x}')
        bin_edges = np.append(sorted_x[indices], sorted_x[-1])
    else: # O(size)
        xmin, xmax = np.min(x), np.max(x)
        step = 1.0*(xmax - xmin)/nbins
        bin_edges = np.arange(xmin, xmax + 1e-6*step, step)
        assert len(bin_edges) == nbins + 1
    return bin_edges

def ComputeHistogramOld(X, wts, nbins_spec, rgram, flags):
    """
    Generate a histgram for each row of X, or regressograms (if rgram is True) for x,y[:] = X[0],X[1:]
    Supported nbins_specs:
      'sqrt':  sqrt(len(x))
      'qbrt':  len(x)**(1./3)
      25:      25
    flags is optional bitmask which can include these bits:
       HIST_EQUAL_WT - split data into equal-weigh bins (or compute densities for unweighted x-histogram)
       HIST_YERR     - return both y.mean and y.std for each bin
    Return: bin_edges, (ymeans, yerr), kind('x' or 'xy')
    """
    assert X.ndim == 2 and X.shape[0] > 1, f'unsupported X shape {X.shape}'
    x = X[0] if rgram else X.flatten()
    nbins = _Parse_nbins(nbins_spec, len(x))
    bin_edges = ComputeBinEdges(x, wts, nbins, flags)
    if rgram: return bin_edges, ComputeXYHist(x, X[1:], wts, bin_edges, flags)
    else:     return bin_edges, ComputeXHist(X,         wts, bin_edges, flags)

def ComputeHistogramFromDataframe(df, nbins_spec, opt, flags):
    """
    Generate a histgram for each column of df, or regressograms (if rgram is True) for x,y[:] = df[:,0],df[:,1:]
    Supported nbins_specs:
      'sqrt':  sqrt(len(x))
      'qbrt':  len(x)**(1./3)
      25:      25
    flags is optional bitmask which can include these bits:
       HIST_EQUAL_WT - split data into equal-weigh bins (or compute densities for unweighted x-histogram)
       HIST_YERR     - return both y.mean and y.std for each bin
    Return: bin_edges, (ymeans, yerr), kind('x' or 'xy')
    """
    if opt.wts_col is not None:
        end = len(df.columns) -1
        wts = df.iloc[:, -1].to_numpy()
        x = df.iloc[:,0].to_numpy() if opt.rgram else df.iloc[:, 0:end].to_numpy().flatten()
    else:
        end = len(df.columns)
        wts = None
        x = df.iloc[:,0].to_numpy() if opt.rgram else df.to_numpy().flatten()
    nbins = _Parse_nbins(nbins_spec, len(x))
    bin_edges = ComputeBinEdges(x, wts, nbins, flags)
    #print(f'XXX5 bin_edges={bin_edges}')
    if opt.rgram: return bin_edges, ComputeXYHist(x, df.iloc[:, 1:end].to_numpy().transpose(), wts, bin_edges, flags)
    else:         return bin_edges, ComputeXHist(    df.iloc[:, 0:end].to_numpy().transpose(), wts, bin_edges, flags)

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

def ComputeXHist(X, wts, bin_edges, flags): # use stdlib
    nx = X.shape[0]
    nbins = len(bin_edges) - 1
    hists = np.zeros((nx, nbins))
    for i in range(nx):
        hists[i, :], _ = np.histogram(X[i, :], bins=bin_edges, weights=wts, density=bool(flags & HIST_EQUAL_WT))
    #return hist.reshape(1, len(bin_edges) - 1), None
    return hists, None

def MakeLadderDataframe(bin_edges, ymeans):
    ny, nbins = ymeans.shape
    X = np.empty((2*len(bin_edges), 1+ny))
    for bin, bin_edge in enumerate(bin_edges):
        X[2*bin,   0] = bin_edge
        X[2*bin+1, 0] = bin_edge
        for iy in range(ny):
            X[2*bin,   1+iy] = ymeans[iy, bin-1] if bin > 0     else np.nan
            X[2*bin+1, 1+iy] = ymeans[iy, bin  ] if bin < nbins else np.nan
    return X

def MakeCenterBinDataframe(bin_edges, ymeans):
    ny, nbins = ymeans.shape
    X = np.zeros((nbins, 1 + ny))
    for bin in range(nbins):
        X[bin, 0] = 0.5*(bin_edges[bin] + bin_edges[bin + 1])
    X[:,1:] = ymeans.transpose()
    return X

def RunHistogram(fnames, cols_spec, opt):
    """Plot histogram or regressogram for one or more data columns
    """
    cols = ParseColumns(cols_spec)
    if False:
        data, names, wts, wts_name = ReadDataframeNoPandas(fnames[0], cols, opt)
    if True:
        if opt.wts_col is not None:
            cols.append(opt.wts_col) # weights are always in the last column
        df = ReadDataframe(fnames, cols, False, opt)
        data, names = DfAdaptor(df, False)
        #print(f'XXX6 names={names}')
    flags = 0
    if opt.equal_wt:
        flags |= HIST_EQUAL_WT
    if opt.yerr:
        flags |= HIST_YERR
    #bin_edges, (ymeans, yerr) = ComputeHistogramOld(data.transpose(), wts, opt.nbins, opt.rgram, flags) # ymeans can contain nans
    bin_edges, (ymeans, yerr) = ComputeHistogramFromDataframe(df, opt.nbins, opt, flags) # ymeans can contain nans
    ny, nbins = ymeans.shape
    if not opt.title:
        opt.title = f"{'rgram' if opt.rgram else 'hgram'} ({nbins} {'EW ' if opt.equal_wt else ''}bins)"
    if opt.hgram:
        names.insert(0, 'range')
    if opt.wts_col is not None:
        opt.title += f' wts={names[-1]}'
        del names[-1]
        del cols[-1]
    #print(f'XXX7 ymeans: {ymeans.shape}:\n{ymeans}')
    X = MakeCenterBinDataframe(bin_edges, ymeans) # bezier-smoothing a ladder is no good
    if yerr is None:
        #X = MakeLadderDataframe(bin_edges, ymeans)
        errorbars = None
    else:
        assert opt.rgram
        assert yerr.shape == (ny, nbins)
        errorbars = yerr.transpose()
    #print(f'XXX8 cols={cols} names={names} X: {X.shape}:\n{X}')
    PlotDataframe(X, names, errorbars, **vars(opt))

#
# Generic plotting frontend
#
def PlotDataframe(X, names, errorbars, **kwargs):
    """
    Plot columns of dataframe X.  Default kwargs['cols']=range(X.shape[1]).
    Use x=X[:,cols[0]] and one or more Y=X[:,cols[1:]]
    The function writes a tmp datafile and runs a gnuplot script on it
    Supported kwargs:
    how:      ['with lines']
    cumsum:   [False] -- X is modified in-place
    diff:     [False] -- X is modified in-place
    bezier:   [False] -- use gnuplot 'smooth bezier' option
    llr:      [None]  -- use LLR kernel smoothing over this bandwidth
    title:    ['']
    logscale: [''] -- x, y, or xy
    filtery:  [''] -- to be applied via eval to y data, X modified in-place
    driver:  'gnuplot' or 'pyplot' (or 'plt') [plot.DEFAULT_DRIVER]
    output:   [''] -- save graphics to this pdf file
    gnuplot-specific options:
      gnuplot: ['gnuplot'] -- gnuplot executable
      term:    ['qt']      -- gnuplot terminal
      verbose: [False]     -- print gnuplot command line to stdout
    """
    assert X.ndim == 2
    nrows = X.shape[0]
    ncols = X.shape[1]
    assert len(names) == ncols, f'mismatched names={names} vs ncols={ncols}'
    cols = range(ncols)
    xcol = cols[0]
    ycols = cols[1:]
    filtery = kwargs.get('filtery', '') # like y-y**3
    if kwargs.get('nozeros', False):
        for ycol in ycols:
            y = X[:, ycol]
            X[:, ycol] = np.where(y == 0, np.nan, y)
    if filtery:
        for ycol in ycols:
            X[:, ycol] = eval(re.sub(r'\by\b', 'X[:, ycol]', filtery))
    if kwargs.get('cumsum', False):
        for ycol in ycols:
            X[:, ycol] = np.nancumsum(X[:, ycol])
    elif kwargs.get('diff', False):
        for ycol in ycols:
            X[1:, ycol] = np.diff(X[:, ycol])
            X[0, ycol] = 0.0
    title = kwargs.get('title', '')
    if errorbars is not None:
        assert errorbars.shape == (nrows, len(ycols))
    cumsum  = kwargs.get('cumsum', False)
    diff    = kwargs.get('diff', False)
    bezier  = kwargs.get('bezier', False)
    verbose = kwargs.get('verbose', False)
    if filtery: title += f' ({filtery})'
    if cumsum:  title += ' cumsum'
    if diff:    title += ' diff'
    if bezier:  title += ' bezier'
    llr_bandwidth = kwargs.get('llr', 0)
    if llr_bandwidth:
        SmoothY(X, xcol, ycols, llr_bandwidth)
        title += f' LLR({llr_bandwidth})'
    #print(f'XXX9 cols={cols} xcol={xcol} ycols={ycols} names={names}') # names[0] == None ??
    driver = kwargs.get('driver', DEFAULT_DRIVER)
    if driver in ('gnuplot', 'gnu'):  ExecuteGnuplot(X, xcol, ycols, names, title, errorbars, kwargs)
    elif driver in ('pyplot', 'plt'): ExecutePyplot(X, xcol, ycols, names, title, errorbars, kwargs)
    else:                             assert False, f'unsupported plotting driver "{driver}"'

def DfAdaptor(df, x_in_index):
    if x_in_index:
        x = df.index.to_numpy()
        Y = df.to_numpy()
        xname = df.index.name
        if xname is None:
            xname = 'seqnum'
        ynames = list(df.columns)
        return np.hstack([x.reshape(len(df), 1), Y]), [xname] + ynames
    else:
        return df.to_numpy(), list(df.columns)

#
# Gnuplot driver
#
def ExecuteGnuplot(X, xcol, ycols, names, title, errorbars, kwargs):
    nrows = X.shape[0]
    ncols = X.shape[1]
    ny = ncols - 1
    verbose = kwargs.get('verbose', False)
    logscale = kwargs.get('logscale', '')
    cols = [xcol] + list(ycols)
    script = []
    if False: # maybe play with palette
        script.append('set style function pm3d')
        script.append('set palette color')
        script.append('f(x)=(x+10)/20')
        script.append('set cbrange [f(-10):f(10)]') # [0:1]
    #script.append(f"set datafile separator ','")
    script.append(f"set datafile missing 'nan'")
    script.append(f"set xlabel '{names[xcol]}'")
    if kwargs.get('nolegend', False):
        script.append('unset key')
    else:
        script.append('set key autotitle columnhead')
    if logscale:
        script.append(f'set logscale {logscale}')
    if kwargs.get('xzeroaxis', False):
        script.append('set xzeroaxis')
    how   = kwargs.get('how', 'with lines')
    extra = 'smooth bezier' if kwargs.get('bezier', False) else ''
    output = kwargs.get('output', None)
    if output is not None:
        script.append('set terminal pdfcairo')
        #script.append('set terminal pdfcairo enhanced color notransparent font "Arial,10"')
        script.append(f"set output '{output}'")
        #script.append(f"set title '{title}'")
    else:
        script.append(f"set term qt title '{title}'")
    if kwargs.get('noaxes', False):
        script.append('unset xtics')
        script.append('unset ytics')
        script.append('unset xlabel')
        script.append('unset border')
        script.append('unset title')
    with tempfile.NamedTemporaryFile(delete=not kwargs.get('keeptmp', False)) as tmp:
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
            items.append(f"'{f}' using {xcol + 1}:(${ycol + 1}) {extra} {how}") # 1:($2) will make line breaks on missing data
            if errorbars is not None:
                items.append(f"'{f}' using {xcol + 1}:{ycol + 1}:{ycol + 1 + ny} with errorbars notitle")
        all_items = ', '.join(items)
        script.append(f"plot {all_items}")
        cmd_line = [kwargs.get('gnuplot', 'gnuplot'), '-persist', '-e', '; '.join(script)]
        if verbose:
            print(f'# cmd_line: {cmd_line}')
        subprocess.run(cmd_line, shell=False, stderr=None if verbose else subprocess.DEVNULL)

def EpanechnikovKernel(x): # a reasonable default kernel, numpy-vectorized
    return 0.75*np.where(np.abs(x) < 1, 1 - x**2, 0)

def LLR_values(xx, yy, h, at_x_points=None, kernel=EpanechnikovKernel, weights=None):
    """
    Compute Local Linear Regression values at specified x points [xx by default].
    """
    nx = len(at_x_points) if at_x_points else len(xx)
    yy_smooth = np.empty(nx)
    yy_smooth[:] = np.nan
    for i in range(nx):
        x = at_x_points[i] if at_x_points else xx[i]
        kk = kernel((xx - x)/h) if weights is None else weights*kernel((xx - x)/h)
        M1 = np.dot(kk, xx - x)
        M2 = np.dot(kk, (xx - x)**2)
        s  = kk*M2 - kk*(xx - x)*M1
        denom = np.sum(s)
        if denom:
            s /= denom
        #print(f'XXX11 s={s}\nyy={yy}\nnansum={np.nansum(s*yy)}\ndot={np.dot(s,yy)}')
        if np.any(s):
            yy_smooth[i] = np.nansum(s*yy) # np.dot(s, yy)
    return yy_smooth

def SmoothY(X, xcol, ycols, llr_bandwidth): # works with nans
    for ycol in ycols:
        #print(f'XXX10 y[{ycol}] before:\n{X[:, ycol]}')
        X[:, ycol] = LLR_values(X[:, xcol], X[:, ycol], llr_bandwidth)
        #print(f'XXX10 y[{ycol}] after:\n{X[:, ycol]}')

def ExecutePyplot(X, xcol, ycols, names, title, errorbars, kwargs):
    import matplotlib.pyplot as plt
    nrows = X.shape[0]
    ncols = X.shape[1]
    ny = ncols - 1
    verbose = kwargs.get('verbose', False)
    logscale = kwargs.get('logscale', '')
    cols = [xcol] + list(ycols)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    lw = 1 # default (2?) seems thick
    ls = 'solid'
    if logscale and 'x' in logscale: plt.xscale('log')
    if logscale and 'y' in logscale: plt.yscale('log')
    for iy in range(ny):
        x = X[:, 0]
        y = X[:, 1 + iy]
        if kwargs.get('bezier', False):
            assert False, 'bezier not supported by pyplot driver'
            if 0 == iy:
                from scipy.interpolate import interp1d
            y = interp1d(x, y, kind='cubic')(x)
        if kwargs.get('nolegend', False):
            label = None
        else:
            label = names[1 + iy]
        if errorbars is not None:
            plt.errorbar(x, y, errorbars[:, iy], label=label, linewidth=lw, linestyle=ls)
        elif kwargs.get('points', False):
            plt.scatter(x, y, label=label)
        else:
            plt.plot(x, y, label=label, linewidth=lw, linestyle=ls)
    if kwargs.get('noaxes', False):
        plt.axis('off')
    else:
        ax.set_xlabel(names[0])
        if not kwargs.get('nolegend', False):
            plt.legend()
    output = kwargs.get('output', None)
    if output is not None:
        plt.show(block=False)
        fig.savefig(output, bbox_inches='tight')
    else:
        plt.show(block=True)

#
# 2D data plotter
#   
def RunDataplot(fnames, cols_spec, opt):
    cols = ParseColumns(cols_spec)
    df = ReadDataframe(fnames, cols, True, opt)
    opt_dict = vars(opt)
    if 'title' not in opt_dict or opt.title is None:
        opt.title = ' '.join(fnames)
    #print(f'XXX12 df:\n{df}')
    X, names = DfAdaptor(df, True)
    PlotDataframe(X, names, None, **vars(opt))

def PrintExamples():
    print('  ls -l /usr/bin|grep -v total|plot -cnQw 4 # cumsum of file sizes (input without header)\n'
          + '  plot -t | plot 0-4                        # plot columns 1-4 vs column 0\n'
          + '  for i in 3 4 5; do plot -ts $i > tmp$i.csv; done; plot 0,2-3 tmp{3,4,5}.csv -zc # cumulative plots from multiple files\n'
          + '  plot -t | plot 0-4 -zb                    # bezier-smoothed plots (gnuplot driver only)\n'
          + '  plot -t | plot 0-4 -zL 20                 # same plots, LLR smoothed\n'
          + '  plot -t | plot 1-4 -p                     # scatter plots\n'
          + '  plot -t | plot 1   -Hz                    # histograms of column 1 data\n'
          + '  plot -t | plot 2-4 -HzL 0.5               # histograms of multiple columns, LLR-smoothed\n'
          + '  plot -t | plot 1-4 -Rz                    # regressograms of columns 2-4 (y) vs column 1 (x)\n'
          + '  plot -t | plot 1-4 -RzL 0.4               # regressograms smoothed with bandwidth 0.4\n'
          + '  plot -t | plot 1-4 -RzW 5                 # weighted regressograms\n'
          + '  plot -t | plot 1-4 -eERW 5 -B 60 -L 0.4   # weighted regressograms with 60 equal-weight bins and errorbars\n'
          + '  plot -tN 500 -B 100|plot 1-500 -RqL 0.3   # spaghetti art\n'
          + '  plot -tPN 200|plot 1-199 -bqx             # alpha art\n'
          + '  unplot                                    # kill all active gnuplot windows\n'
          + '  plot -X | head -14 | bash                 # run all of the above'
)
    
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='Plot columnar data in files or stdin. Default action: Print columns available for plotting')
    parser.add_argument('-X', '--examples',   action='store_true', help='Print use examples')
    parser.add_argument('-T', '--title',                           help='Use this window title [filename]')
    parser.add_argument('-q', '--nolegend',   action='store_true', help='Do not mark plotted line in legend')
    parser.add_argument('-z', '--xzeroaxis',  action='store_true', help='Add horizontal zero line')
    parser.add_argument('-w', '--whitespace', action='store_true', help='Split input by whitespace [comma]')
    parser.add_argument('-Q', '--seqnum',     action='store_true', help='Use sequence number for x data [first column]')
    parser.add_argument('-n', '--noheader',   action='store_true', help='Indicate that data has no header.  Generate header names F0, F1, etc')
    parser.add_argument('-c', '--cumsum',     action='store_true', help='Plot cumulative sums')
    parser.add_argument('-b', '--bezier',     action='store_true', help='Plot bezier-smooth data (gnuplot driver only)')
    parser.add_argument('-x', '--noaxes',     action='store_true', help='Generate plot without axes or border')
    parser.add_argument('-L', '--llr', metavar='bandwidth', type=float, help='Smooth data using local linear regression (LLR) over this bandwidth')
    parser.add_argument('-d', '--diff',       action='store_true', help='Plot differences')
    parser.add_argument('-f', '--filtery',    type=str,            help='Filter (transform) y data through an expression to be applied via eval().  Examples: "y-y**3", "np.sin(y)"')
    parser.add_argument('-Z', '--nozeros',    action='store_true', help='Filter out zero y data by replacing 0 -> nan')
    parser.add_argument('-p', '--points',     action='store_true', help='Plot with points (e.g. for a scatter plot) [with lines]')
    parser.add_argument('-l', '--logscale',   type=str,            help='Use log scale for x|y|xy')
    parser.add_argument('-H', '--hgram',      action='store_true', help='Plot histogram of indicated data column(s)')
    parser.add_argument('-R', '--rgram',      action='store_true', help='Plot regressogram of data columns: first treated as x, rest as y(s)')
    parser.add_argument('-B', '--nbins',                           help=f'For histogram: use this many bins: sqrt: size**(1/2), qbrt: size**(1/3), or <int> [{HIST_DEFAULT_NBINS}]', default=HIST_DEFAULT_NBINS)
    parser.add_argument('-W', '--wts_col',    type=int,            help='For histogram: use this column for weights')
    parser.add_argument('-e', '--yerr',       action='store_true', help='For regressogram: plot with yerrorbars')
    parser.add_argument('-E', '--equal_wt',   action='store_true', help='Use equal-weight (histo|regresso)gram bins. Implies a density plot for histogram [equal-size]')
   #parser.add_argument('-g', '--grep_key',                        help='Skip input lines without this word')
   #parser.add_argument('-M', '--multikey_col', nargs=1,           help='Plot separate lines for for each value in this column (TODO)')
    parser.add_argument('-D', '--driver',                          help=f'Use this backend graphics driver: gnuplot or pyplot [{DEFAULT_DRIVER}]', default=DEFAULT_DRIVER)
    parser.add_argument('-o', '--output',                          help='Plot to this pdf file [qt window].  Seems to work better with -D plt option')
    parser.add_argument('-v', '--verbose',    action='store_true', help='Print gnuplot command line')
    parser.add_argument('-K', '--keeptmp',    action='store_true', help='Keep tmp file for gnuplot (helpful with -v option) [delete]')
    parser.add_argument('-t', '--test_csv',   action='store_true', help='Print a csv stream for testing')
    parser.add_argument('-P', '--test_parametric', action='store_true', help='For test_csv: Print data for parametric plots')
    parser.add_argument('-s', '--seed',       type=int,            help='Use this seed for random test_csv data', default=0)
    parser.add_argument('-N', '--ntests',     type=int,            help=f'test_csv: print this many columns [{DEFAULT_NTESTS}]', default=DEFAULT_NTESTS)
    parser.add_argument('files_and_columns',  nargs='*')
    args = parser.parse_args()

    if args.examples:
        PrintExamples()
    elif args.test_csv:
        PrintTestDataframe(args)
    else:
        args.how = 'with points' if args.points else 'with lines'
        fnames = []
        columns = None
        for arg in args.files_and_columns:
            if os.path.isfile(arg) or arg == '-':
                fnames.append(arg)
            elif not re.search(r'[_a-zA-Z]', arg): # assume this is a columns spec (numbers, commas, dashes)
                assert columns is None, f'multiple columns specs ({arg} after {columns}) not supported'
                columns = arg
            else:
                assert False, f'file {arg} not found'
        if not fnames:
            fnames.append('-')
        if columns is None:            PrintColumns(fnames[0], args)
        elif args.hgram or args.rgram: RunHistogram(fnames, columns, args)
        else:                          RunDataplot(fnames, columns, args)

if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        pass
