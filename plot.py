#!/usr/local/miniconda3/bin/python

# Copyright (c) 2020 Michael Isichenko. All rights reserved
# This code is licensed under the terms of the MIT license

import sys, os, io, subprocess, argparse, tempfile, re
from itertools import permutations 
import numpy as np
import pandas as pd
import lzma
from scipy import stats

#
# Constants
#
VERSION            = '0.5 as of 2020-08-30'
_DEFAULT_DRIVER    = 'gnuplot'
DEFAULT_DRIVER     = os.environ['PLOT_DRIVER'] if 'PLOT_DRIVER' in os.environ else _DEFAULT_DRIVER
CHUNK              = 1024
DEFAULT_NTESTS     = 4
HIST_DEFAULT_NBINS = 'sqrt'

_linetype_table = (
    # code gnuplot             pyplot
    ('l', 'with lines',       'solid'),
    ('p', 'with points',      'scatter'),
    ('d', 'with linespoints', 'dashed'),
    ('D', 'with dots',        'dotted'),
    ('i', 'with impulses',     None),
    )

def SupportedLinetype(how):
    for code, _, _ in _linetype_table:
        if code == how:
            return True
    return False

#
# Generic utilities
#
def ComputeMI(xx, yy, wts):
    assert len(xx) == len(yy)
    nbins = int(np.sqrt(len(xx)))
    nbins = min(500, max(10, nbins))
    if np.isnan(xx).any() or np.isnan(yy).any(): # scipy does not support skipping nans
        good = ~np.isnan(xx) & ~np.isnan(yy)
        xx = xx[good]
        yy = yy[good]
        wts = wts[good]
    if np.sum(wts) == 0:
        return np.nan
    wxy, ingnore_x_edge, ignore_y_edge, ignore_bin_number = stats.binned_statistic_2d(xx, yy, wts, 'sum', nbins)
    wx, wy, w = np.sum(wxy, axis=1), np.sum(wxy, axis=0), np.sum(wxy)
    good = np.where(wxy > 0)
    wxy = wxy[good]
    wxwy = np.outer(wx, wy)[good]
    mi = np.sum(wxy/w*np.log(wxy*w/wxwy))
    return mi

def Neff(wts): # weighted laws of large numbers (for esm - error of sampled mean)
    return np.sum(wts)**2/np.sum(wts**2)

def WeightedMeanStd(x, wts):
    """Compute weighted mean and std"""
    if len(x) == 0:
        return np.nan, np.nan
    assert not np.isnan(wts).any(), f'detected nans in wts={wts}'
    if np.isnan(x).any():
        good = ~np.isnan(x)
        x = x[good]
        wts = wts[good]
    mean = np.average(x, weights=wts)
    var  = max(np.average((x - mean)**2, weights=wts), 0.0)
    return mean, np.sqrt(max(var, 0.0))

def ParseColumns(columns, opt):
    cols = []
    hows = []
    groups = columns.split(',')
    for group in groups:
        if '-' in group: # 3-10 or 3-10:p, only int ranges are supported
            beg, end = group.split('-')
            if len(end) == 0:
                end = None
                how = opt.how
            elif ':' in end:
                end, how = end.split(':')
                assert SupportedLinetype(how), f'unsupported linetype "{how}" in {columns}'
            else:
                how = opt.how
            if end is None:
                cols.append(int(beg))
                cols.append(-1) # open range
                hows.append(how)
                hows.append(how)
            else:
                for col in range(int(beg), int(end) + 1):
                    cols.append(col)
                    hows.append(how)
        else: # 5 or 5:p
            assert len(cols) == 0 or cols[-1] != -1, f'invalid use of open range for cols={cols}'
            if ':' in group:
                group, how = group.split(':')
                assert SupportedLinetype(how), f'unsupported linetype "{how}" in {columns}'
            else:
                how = opt.how
            try: # support for mixed int or str columns:
                col = int(group)
            except ValueError:
                col = group
            cols.append(col)
            hows.append(how)
    assert len(cols) == len(hows) # even for x column
    return cols, hows

def PrintColumns(fname, opt):
    if os.path.isfile(fname):
        if(fname.endswith(".xz")):
            with lzma.open(fname) as input:
                while True:
                    line = input.readline().decode()
                    if len(line) and line[0] != '#':
                        break
        else:
            with open(fname) as input:
                while True:
                    line = input.readline()
                    if len(line) and line[0] != '#':
                        break
    else:
        while True:
            line = sys.stdin.readline()
            if len(line) and line[0] != '#':
                break
    cols = line.rstrip().split(None if opt.whitespace else ',')
    for c, col in enumerate(cols):
        print('{:3d} {}'.format(c, col))

def IsMixed(cols):
    have_int = False
    have_str = False
    for col in cols:
        if   isinstance(col, int): have_int = True
        elif isinstance(col, str): have_str = True
        if have_int and have_str:
            return True
    return False

def NamesToNumbersFromLine(cols, line, sep):
    #print(f'XXX line: {line}')
    new_cols = []
    tokens = line.rstrip().split(sep)
    for col in cols:
        if isinstance(col, int):
            if col >= 0:
                new_cols.append(col)
            else: # open range
                assert len(cols) > 1
                assert col == cols[-1]
                for c in range(cols[-2] + 1, len(tokens)):
                    new_cols.append(c)
        else: # not isinstance(col, str)
            #print(f'XXX type({col})={type(col)} -- not an int')
            colon = col.find(':')
            new_cols.append(tokens.index(col if colon < 0 else col[0:colon]))
    new_names = [tokens[c] for c in new_cols]
    return new_cols, new_names
    
def NamesToNumbersFromBuffer(cols, header_buffer, sep):
    for line in header_buffer.decode().split('\n'):
        if line.startswith('#'):
            continue
        return NamesToNumbersFromLine(cols, line, sep)
    return None

def NamesToNumbersFromFile(cols, fname, sep):
    if False and fname == '-': # nice try
        # io.DEFAULT_BUFFER_SIZE = 8192 is not enough for input with a very wide header
        with os.fdopen(sys.stdin.fileno(), 'r', 2*io.DEFAULT_BUFFER_SIZE) as newin:
            for line in newin:
                if line.startswith('#'):
                    continue
                return NamesToNumbersFromLine(cols, line, sep)
    else:
        with (lzma.open(fname) if fname.endswith('.xz') else open(fname)) as file:
            #print(f'XXX {file}')
            for line in file:
                if line.startswith('#'):
                    continue
                #print(f'XXX "{line}"')
                return NamesToNumbersFromLine(cols, line, sep)
    return None

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
    #N = 40 if opt.test_parametric else 1000 # 60 ok
    N = 200 if opt.test_parametric else 1000
    if opt.test_clean:
        print('x,true_y,noisy_y')
        for idx in range(N):
            x = 2*np.pi*idx/N
            true_y = np.sin(x)
            noisy_y = true_y + np.random.normal()
            print(f'{x},{true_y},{noisy_y}')
        return
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

def Within(x, beg, end): return (beg <= x) & (x < end)
def AreSorted(xx):
    for i in range(len(xx) - 1):
        if xx[i] >= xx[i + 1]:
            return False
    return True

def AllInts(xx):
    for x in xx:
        if not isinstance(x, int):
            return False
    return True

def GrepInput(input, grep_key):
    if type(input) is str:
        pass

def ReadDataframe(fnames, cols, wts_col, merge_by_x, opt):
    """
    read_csv wrapper supporting mixed int/str columns, one or multiple files including stdin:
    Read data from one or more compatible files into pandas dataframe.
    The reasons for using pandas are: (1) fast pd.read_csv() and (2) pd.merge() functionality for multiple files.
    For multiple files, generate a wide dataframe with columns like <fname>:<column_name>
    Any weight columns must be included in the dataframe
    """
    sep = r'\s+' if opt.whitespace else ','
    separate_wts_col = wts_col is not None and wts_col not in cols
    all_cols = cols + wts_col if separate_wts_col else cols
    #print(f'XXX cols={cols} all_cols={all_cols}')
    if opt.noheader:
        assert all_cols[-1] != -1, 'open range requires a header (TODO: support this)'
        header = None                 
        names  = [f'C{c}' for c in all_cols]
    else:
         header = 0 # same as 'infer'
         names = None
    if   opt.seqnum: index_col = False # generate index 0,1,...
    elif merge_by_x: index_col = 0     # index contains x data.  NOTE: index_col is per usecols, not the whole input
    else:            index_col = None  # no index
    skiprows = None
    if False:
        # XXX pd.read_csv skiprows feature is buggy: header is incorrectly taken from last skipped row
        if   opt.start is not None and opt.end is not None: skiprows = lambda x: x < float(opt.start) or x >= float(opt.end)
        elif opt.start is not None:                         skiprows = lambda x: x < float(opt.start)
        elif opt.end   is not None:                         skiprows = lambda x: x >= float(opt.end)
    keep_rows = None # less efficient workaround: subsetting dataframe with keep_rows
    #if opt.start is not None: print(f'XXX opt.start={opt.start}')
    #if opt.end is not None: print(f'XXX opt.end={opt.end}')
    if   opt.start is not None and opt.end is not None: keep_rows = lambda df: Within(df.index,   float(opt.start), float(opt.end))
    elif opt.start is not None:                         keep_rows = lambda df:        df.index >= float(opt.start) # float(opt.start)
    elif opt.end   is not None:                         keep_rows = lambda df:        df.index <  float(opt.end) # float(opt.end)
    for i, fname in enumerate(fnames):
        needs_reorder = opt.swapxy or AllInts(all_cols) and not AreSorted(all_cols)
        if (IsMixed(all_cols) or # read_csv does not support mixed int/str usecols, a workaround uses peaking into input
            needs_reorder):  # read_csv generates a col-sorted columns regardless of usecols order.  df can be reordered only using string names
            if fname == '-': # io.DEFAULT_BUFFER_SIZE is 8192.  if header of stdin is bigger, you are out of luck: read file or use numeric columns then
                all_cols, all_col_names = NamesToNumbersFromBuffer(all_cols, sys.stdin.buffer.peek(io.DEFAULT_BUFFER_SIZE), sep)
            else:
                all_cols, all_col_names = NamesToNumbersFromFile(all_cols, fname, sep)
            #print(f'XXX all_cols={all_cols}')
        file_df = pd.read_csv(sys.stdin if fname == '-' else fname,
                              sep=sep, header=header, parse_dates=not opt.seqnum, skiprows=skiprows,
                              names=names, usecols=all_cols, index_col=index_col, dtype=float, comment='#')
        #print(f'XXX df.index.name={file_df.index.name} df.columns={file_df.columns} all_cols={all_cols}')
        if needs_reorder:
            file_df = file_df[all_col_names[1:] if merge_by_x else all_col_names] # abscissa x is in df.index
        if len(fnames) > 1:
            file_df.columns = [fname + ':' + c for c in file_df.columns]
        if keep_rows is not None:
            file_df = file_df.loc[keep_rows, :]
        if i == 0:
            df = file_df
        else:
            df = pd.merge(df, file_df, how='outer', left_index=True, right_index=True)
    #print(f'XXX2 df columns={df.columns.tolist()}:\ndf.index:\n{df.index}\ndf:\n{df}')
    if separate_wts_col:
        cols = all_cols[0:-1]
        opt.wts_name = df.columns[-1]
    else:
        cols = all_cols
        opt.wts_name = None
    return df, cols # df contains any wts data, cols does not reference it

def ReadDataframeNoPandas(fname, cols, merge_by_x, opt): # XXX unused
    data = np.zeros((CHUNK, len(cols)))
    names = []
    size = 0
    with (open(fname, encoding='utf-8') if fname != '-' else sys.stdin) as input:
        for line in input:
            if line.startswith('#'):
                continue
            if opt.grep_key and not line.startswith(opt.grep_key):
                continue
            tokens = line.rstrip().split() if opt.whitespace else line.rstrip().split(',')
            if len(names) == 0:
                if opt.noheader:
                    names = ['seqno' if c == -1 else f'F{c}' for c in cols]
                else:
                    names = [tokens[c] for c in cols]
                continue
            if size >= data.shape[0]:
                data.resize(size + CHUNK, len(cols))
            data[size] = [float(size if c == -1 else tokens[c]) for c in cols]
            size += 1
    data.resize(size, data.shape[1])
    return data, names

#
# Histograms
#
def _Parse_nbins(nbins, nx):
    if   isinstance(nbins, int):  return nbins
    elif isinstance(nbins, str):
        if   nbins == 'sqrt':     return max(int(np.ceil(np.sqrt(nx)/2)), 20)
        elif nbins == 'cbrt':     return max(int(np.ceil(nx**(1.0/3))), 20)
        else:                     return int(nbins) # will raise ValueError for bad string
    raise ValueError(f'unsupported nbins spec "{nbins}"')

def ComputeBinEdges(x, wts, nbins, opt):
    if opt.xtrim: assert 0 < opt.xtrim and opt.xtrim < 0.5
    good_data = ~np.isnan(x) & (wts > 0)
    x = x[good_data]
    wts = wts[good_data]
    if np.sum(wts) == 0:
        return None
    if opt.equal_wt:
        perm = np.argsort(x)
        xmin, xmax = x[perm[0]], x[perm[-1]]
        cum_wts = np.cumsum(wts[perm])
        total = cum_wts[-1]
        wt_levels = np.linspace(total*opt.xtrim, total*(1 - opt.xtrim), num=nbins + 1)
        indices = np.searchsorted(cum_wts, wt_levels, 'left')
        bin_edges = np.array([x[perm[idx]] for idx in indices])
        bin_edges[-1] += (bin_edges[-1] - bin_edges[0])*1e-6 # extend right edge a little
    else:
        xmin, xmax = np.min(x), np.max(x)
        bin_edges = np.linspace(xmin + opt.xtrim*(xmax - xmin), xmax - opt.xtrim*(xmax - xmin), num=nbins + 1) # max included
    return bin_edges

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

def MakeCenterBinDataframe(bin_edges, ymeans, ynames, yerrs=None):
    #print(f'XXX bin_edges={bin_edges} ymeans={ymeans}')
    ny, nbins = ymeans.shape
    x = np.zeros(nbins)
    for bin in range(nbins):
        x[bin] = 0.5*(bin_edges[bin] + bin_edges[bin + 1])
    if yerrs is None:
        return pd.DataFrame(data=ymeans.transpose(), index=x, columns=ynames)
    else:
        return pd.DataFrame(data=np.hstack((ymeans.transpose(), yerrs.transpose())), index=x, columns=ynames + [yname + '_err' for yname in ynames])

def SetCustomWts(df, wts_col):
    wts_field, wts_suffix = wts_col.split(':') # ex: lmdv:5, lmdv:all
    wts = df[wts_field].to_numpy()
    if wts_suffix != 'all':
        level = int(wts_suffix)
        wts = np.where(np.floor(wts) == level, 1, 0)
    return wts
    
def RunHRgram(fnames, columns, opt): # actually hgram or rgram
    """Plot histogram or regressogram for one or more data columns
    """
    cols, hows = ParseColumns(columns, opt)
    if opt.wts_col is not None:
        wts_col = opt.wts_col.split(':')[0] if ':' in opt.wts_col else opt.wts_col
        wts_col, _ = ParseColumns(wts_col, opt)
    else:
        wts_col = None
    dont_merge_by_x = False
    df, cols = ReadDataframe(fnames, cols, wts_col, dont_merge_by_x, opt)
    #print(f'XXX cols={cols} df.columns={df.columns} wts_col={wts_col} opt.wts_name={opt.wts_name}')
    #print(f'XXX df.head:\n{df.head()}')
    if opt.wts_col is None:
        wts = np.ones(len(df))
    else:
        if ':' in opt.wts_col:
            wts = SetCustomWts(df, opt.wts_col)
        else:
            wts = df.iloc[:, len(df.columns) - 1].to_numpy()
        del df[opt.wts_name]
    while len(hows) < len(cols):
        hows.append(hows[-1])
    if len(fnames) > 1:
        orig_hows = hows
        for _ in range(1, len(fnames)):
            hows.extend(orig_hows)
    data_width = len(df.columns)
    if opt.rgram:
        xcol = 0
        ycols = list(range(1, data_width))
        if opt.yscale != 1:
            #print(f'XXX before applying yscale={opt.yscale}:\n{df.head()}')
            if opt.swapxy:
                df.iloc[:, xcol] *= opt.yscale
            else:
                for ycol in ycols:
                    df.iloc[:, ycol] *= opt.yscale
            #print(f'XXX after:\n{df.head()}')
    elif opt.hgram:
        xcol = None
        ycols = list(range(0, data_width))
    #print(f'XXX RunHRgram: xcol={xcol} ({df.columns[xcol]}) ycols={list(ycols)} ({[df.columns[ycol] for ycol in ycols]})')
    if opt.rgram and opt.swapxy:
        opt.nplots = len(ycols)
        for ycol in ycols:
            RunHRgramOnDataframe(df, ['l'], ycol, [xcol], wts, opt)
    elif opt.separate:
        opt.nplots = len(ycols)
        for ycol in ycols:
            RunHRgramOnDataframe(df, hows, xcol, [ycol], wts, opt)
    else:
        RunHRgramOnDataframe(df, hows, xcol, ycols, wts, opt)

def RunHRgramOnDataframe(df, hows, xcol, ycols, wts, opt): # actually hgram or rgram
    ny = len(ycols) # for the output plot
    # assemble out_df as follows:
    # hgram: x: merged bins of each df[ycol], yy: obs counts in the bins
    # rgram: x: bins of df[xcol], yy: the means of each df[ycol] in the bins
    out_df = None # modestly sized (nbins x ny) and cheap to copy
    cms = np.zeros((ny, 3)) # count, mean, sdev
    if opt.hgram: # generally multiple yy to histogram and to merge
        for iy, ycol in enumerate(ycols):
            y = df.iloc[:, ycol].to_numpy().flatten()
            good = ~np.isnan(y) & (wts > 0)
            good_y = y[good]
            if len(good_y) == 0:
                print(f"# no good {df.columns[ycol]} data: all {len(y)} are missing")
                continue
            good_wts = wts[good]
            cms[iy, 0] = len(good_y)
            cms[iy, 1] = np.average(good_y, weights=good_wts)
            cms[iy, 2] = np.sqrt(np.average((good_y - cms[iy, 1])**2, weights=good_wts))
            nbins = _Parse_nbins(opt.nbins, len(good_y))
            if 0 == np.sum(good_wts):
                print(f'# no good data for y {df.columns[ycol]}')
                continue
            bin_edges = ComputeBinEdges(good_y, good_wts, nbins, opt)
            assert len(bin_edges) == nbins + 1
            ymeans, ignore_bin_edges = np.histogram(good_y, bins=bin_edges, weights=good_wts, density=opt.equal_wt)
            hgram_df = MakeCenterBinDataframe(bin_edges, ymeans.reshape(1, nbins), [df.columns[ycol]]) # bezier-smoothing a ladder is no good
            hgram_df.index.name = ''
            if 0 == iy:
                out_df = hgram_df
            else:
                out_df = pd.merge(out_df, hgram_df, how='outer', left_index=True, right_index=True)
        yerr_cols = None
    elif opt.rgram: # single x, generally multiple yy
        x  = df.iloc[:, xcol ].to_numpy()
        yy = df.iloc[:, ycols].to_numpy()
        if np.isnan(x).all():
            print(f"# no good {df.columns[xcol]} data: all {len(x)} are missing")
            return
        nbins = _Parse_nbins(opt.nbins, len(x))
        bin_edges = ComputeBinEdges(x, wts, nbins, opt)
        if bin_edges is None:
            print(f"# no good weights for {df.columns[xcol]} data")
            return            
        ymeans = np.zeros((ny, nbins))
        ynames = [df.columns[ycol] for ycol in ycols]
        yerrs = np.zeros((ny, nbins)) if opt.yerr else None
        mis = np.zeros(ny)
        neff, fmean, fstd = np.zeros(ny), np.zeros(ny), np.zeros(ny)
        time_rgram = False
        if time_rgram:
            import time
            tbeg = time.time()
        good_x_wts = ~np.isnan(x) & (wts > 0)
        if True: # direct rgram
            in_bins = np.zeros((nbins, len(x)), dtype=bool)
            for bin in range(nbins):
                in_bins[bin, :] = good_x_wts & (bin_edges[bin] <= x) & (x < bin_edges[bin + 1])
            for iy in range(ny):
                good_y = ~np.isnan(yy[:, iy])
                all_good = good_x_wts & good_y & (x >= bin_edges[0]) & (x < bin_edges[-1])
                mis[iy] = ComputeMI(x[all_good], yy[all_good, iy], wts[all_good])
                total_wts_in_bin = np.zeros(nbins)
                for bin in range(nbins):
                    in_bin = in_bins[bin, :] & good_y
                    wts_in_bin = wts[in_bin]
                    total_wts_in_bin[bin] = np.sum(wts_in_bin)
                    if 0 == total_wts_in_bin[bin]:
                        ymeans[iy, bin] = np.nan
                        if opt.yerr:
                            yerrs[iy, bin] = np.nan
                    else:
                        if opt.yerr:
                            ymeans[iy, bin], yerrs[iy, bin] = WeightedMeanStd(yy[in_bin, iy], wts_in_bin)
                            if opt.esm and not np.isnan(yerrs[iy, bin]):
                                yerrs[iy, bin] /= np.sqrt(Neff(wts_in_bin))
                        else:
                            ymeans[iy, bin] = np.average(yy[in_bin, iy], weights=wts_in_bin)
                neff[iy] = Neff(wts[all_good])
                fmean[iy], fstd[iy] = WeightedMeanStd(ymeans[iy, :], total_wts_in_bin)
        else: # rgram using binned_statistic, seems slightly faster
            with np.errstate(invalid='ignore', divide='ignore'):
                for iy in range(ny):
                    all_good = good_x_wts & ~np.isnan(yy[:, iy]) & (x >= bin_edges[0]) & (x < bin_edges[-1])
                    x_good = x[all_good]
                    y_good = yy[all_good, iy]
                    wts_good = wts[all_good]
                    mis[iy] = ComputeMI(x_good, y_good, wts_good)
                    wy, _, _ =  stats.binned_statistic(x_good, wts_good*y_good, statistic='sum', bins=bin_edges)
                    w,  _, _ =  stats.binned_statistic(x_good, wts_good,        statistic='sum', bins=bin_edges)
                    ymeans[iy, :] = wy/w # contains nans for bins with no good observations
                    neff[iy] = Neff(wts_good)
                    fmean[iy], fstd[iy] = WeightedMeanStd(ymeans[iy, :], w)
                    if opt.yerr:
                        wyy, _, _ =  stats.binned_statistic(x_good, wts_good*y_good**2, statistic='sum', bins=bin_edges)
                        #yerrs[iy, :] = np.where(w > 0, wyy/w - (wy/w)**2, 0)
                        yerrs[iy, :] = np.sqrt(wyy/w - (wy/w)**2) # contains nans
                        #print(f'XXX w={w}\nwy={wy}\nwyy={wyy}\nymeans={ymeans[iy, :]}\nyerrs={yerrs[iy, :]}')
                        if opt.esm:
                            ww,  _, _ = stats.binned_statistic(x_good, wts_good**2, statistic='sum', bins=bin_edges)
                            neff_by_bin = np.where(np.isnan(ww), 1, w**2/ww) # by bin
                            yerrs[iy, :] /= np.sqrt(neff_by_bin)
        if time_rgram: print(f'# rgram time: {time.time() - tbeg}')
        out_df = MakeCenterBinDataframe(bin_edges, ymeans, ynames, yerrs=yerrs) # bezier-smoothing a ladder is no good
        #print(f'XXX xcol={xcol} ycols={ycols} df.index.name={df.index.name} df.columns={df.columns}')
        out_df.index.name = df.columns[xcol] # and not df.index.name
        yerr_cols = list(range(ny, 2*ny)) if opt.yerr else None
    PlotDataframe(out_df, list(range(ny)), yerr_cols, cms if opt.hgram else (mis, neff, fmean, fstd), hows, opt)

def Drawdown(yy):
    lastPeak = 0.0
    drawdown = 0.0
    for cumPnl in np.nancumsum(yy):
        fromLastPeak = cumPnl - lastPeak;
        if fromLastPeak > 0:
            lastPeak = cumPnl
        elif fromLastPeak < drawdown:
            drawdown = fromLastPeak;
    assert drawdown <= 0
    return -drawdown

def ComputeStats(stats_spec, yy):
    if not stats_spec:
        return None
    ny = yy.shape[1]
    ystats = [{} for iy in range(ny)]
    all = 'A' in stats_spec
    for iy in range(ny):
        y = yy[:, iy]
        mean = np.nanmean(y)
        sdev = np.nanstd(y)
        if all or 'c' in stats_spec: ystats[iy]['count' ] = np.count_nonzero(~np.isnan(y)) # len(y)
        if all or 'm' in stats_spec: ystats[iy]['mean'  ] = mean
        if all or 's' in stats_spec: ystats[iy]['sdev'  ] = sdev
        if all or 'a' in stats_spec: ystats[iy]['sharpe'] = mean/sdev*np.sqrt(252.0)
        if all or 'd' in stats_spec: ystats[iy]['dd'    ] = Drawdown(y)
        if all or 'D' in stats_spec: ystats[iy]['dd2s'  ] = Drawdown(y)/sdev if sdev > 0 else np.nan
    return ystats

def Digits(x):
    #print(f'XXX Digits({x})')
    rc = 0 if x == 0 else np.max(int(3 - np.log10(np.abs(x))), 0)
    #print(f'XXX Digits({x})={rc}')
    return rc

def FormatNum(x):
    if np.isnan(x):
        return 'nan'
    if np.abs(x) >= 1e6:
        return f'{x:.2e}'
    else:
        return f'{round(x, Digits(x))}'

def FormatInt(x): return f'{round(x):,}'
        
_stats_fields = 'count mean sdev sharpe dd dd2s'.split()
def PrintStats(ystats, names):
    fields = [] # same for all ycols
    for c, ystat in enumerate(ystats):
        if not ystat:
            continue
        if not fields:
            fields = [f for f in _stats_fields if f in ystat]
            print(','.join(['name'] + fields))
        print(','.join([names[c]] + [f'{FormatNum(ystat[f])}' for f in fields]))

def FormatYStats(ystats, iy):
    ystat = ystats[iy]
    return '(' + ' '.join([f'{f}={FormatNum(ystat[f])}' for f in _stats_fields if f in ystat]) + ')' if ystat else ''

def FormatHgramStats(cms, iy):
    count, mean, sdev = cms[iy]
    return f'(nob={FormatInt(count)}: {FormatNum(mean)} +/- {FormatNum(sdev)})'

def FormatRgramStats(cms, iy):
    mi, neff, fmean, fstd = cms[0][iy], cms[1][iy], cms[2][iy], cms[3][iy]
    return f'(MI={FormatNum(mi)} neff={FormatNum(neff)} pred={FormatNum(fmean)} +/- {FormatNum(fstd)})'
    
def IsDate(name):
    return name.lower() in ('date', 'time')

def ComputeOlsFit(x, yy):
    nob, ny = yy.shape
    icepts = np.zeros(ny)
    betas = np.zeros(ny)
    xmean = np.mean(x)
    denom = np.dot(x - xmean, x - xmean)
    assert denom != np.nan and denom > 0
    if denom != 0:
        for iy in range(ny):
            y = yy[:, iy]
            ymean = np.nanmean(y)
            #num = np.dot(x - xmean, y - ymean) # yy can have nans
            num = np.nansum((x - xmean)*(y - ymean))
            betas[iy] = num/denom
            icepts[iy] = ymean - betas[iy]*xmean
    return icepts, betas

#
# Generic plotting frontend
#
global_plot_idx = 0
def PlotDataframe(df, ycols, yerr_cols, cms, hows, opt):
    """
    Plot ycols of dataframe vs x contained in its index, with optional y data manipulation
    """
    # exporting df to np arrays is convenient for pre-plot manipulations
    #print(f'XXX df.index.name={df.index.name} df.columns={df.columns} ycols={ycols}')
    x         = df.index.to_numpy()
    yy        = df.iloc[:, ycols].to_numpy()
    errorbars = df.iloc[:, yerr_cols].to_numpy() if yerr_cols else None
    xname  = df.index.name if not opt.seqnum else 'seqnum'
    ynames = df.columns.to_list()
    #print(f'XXX df.index.name={df.index.name} opt.seqnum={opt.seqnum} xname={xname} df.columns={df.columns} ynames={ynames}')
    ny = len(ycols)
    if opt.filtery: # like y-y**3
        for iy in range(ny):
            yy[:, iy] = eval(re.sub(r'\by\b', 'yy[:, iy]', filtery))
    ystats = None
    if opt.stats:
        ystats = ComputeStats(opt.stats, yy)
        if 'P' in opt.stats:
            PrintStats(ystats, ynames)
            ystats = None
    if opt.cumsum:
        assert not errorbars, 'cumsum errorbars not supported'
        for iy in range(ny):
            if True: # new (20210429) to stop after y turns into all nans
                y = yy[:, iy]
                if np.isnan(y).any():
                    non_nan_idx = np.where(~np.isnan(y))
                    beg_of_all_nans = non_nan_idx[0][-1] + 1
                    yy[0:beg_of_all_nans, iy] = np.nancumsum(y[0:beg_of_all_nans])
                else:
                    yy[:, iy] = np.cumsum(y[:])    
            else: # this effectively replaced all trailing nans (for shorter dataset) with zeros
                yy[:, iy] = np.nancumsum(yy[:, iy])
    elif opt.diff:
        assert not errorbars, 'diff errorbars not supported'
        for iy in range(ny):
            yy[:, iy] = np.diff(yy[:, iy], prepend=yy[0, iy])
    if opt.ols_fit:
        icepts, betas = ComputeOlsFit(x, yy)
        yy_fit = np.zeros(yy.shape)
        #print(f'XXX x={x}')
        for iy in range(ny):
            #print(f'XXX iy={iy} icept={icepts[iy]} beta={betas[iy]} for yy={yy[:, iy]}')
            yy_fit[:, iy] = icepts[iy] + betas[iy]*x[:]
            ynames.append(f'fit({ynames[iy]})={FormatNum(betas[iy])}*{xname}+{FormatNum(icepts[iy])}')
            hows.append('l')
        yy = np.hstack((yy, yy_fit))
    title = opt.title
    if not opt.title:
        ew = 'EW 'if opt.equal_wt else ''
        if   opt.hgram: title = f'{int(len(df)/len(df.columns))}-bin {ew}hgram'
        elif opt.rgram: title = f'{len(df)}-bin {ew}rgram'
        else:           title = 'plot'
    if not opt.hgram and opt.yscale != 1:      title += f' yscale={opt.yscale}'
    if (opt.hgram or opt.rgram) and opt.xtrim: title += f' xtrim={opt.xtrim}'
    if opt.wts_col is not None:                title += f' wts={opt.wts_col}'
    if opt.filtery:                            title += f' ({filtery})'
    if opt.cumsum:                             title += ' cumsum'
    if opt.diff:                               title += ' diff'
    if opt.smooth:                             title += f' {opt.smooth}'
    if opt.rgram and opt.esm:                  title += ' (esm)'
    if opt.llr:
        SmoothY(x, yy, opt)
        title += f' LLR({opt.llr})'
    #print(f'XXX cols={cols} xcol={xcol} ycols={ycols} names={names}') # names[0] == None ??
    if   opt.driver in ('gnuplot', 'gnu'): ExecuteGnuplot(x, yy, xname, ynames, cms, hows, title, ystats, errorbars, opt)
    elif opt.driver in ('pyplot', 'plt'):  ExecutePyplot(x, yy, xname, ynames, cms, hows, title, ystats, errorbars, opt)
    else:                                  assert False, f'unsupported plotting driver "{driver}"'
    global global_plot_idx
    global_plot_idx += 1

def AllDataframeColumns(df, x_in_index):
    if x_in_index:
        return [df.index.name] + list(df.columns)
    else:
        return list(df.columns)
    
def GnuplotLinetype(how):
    for code, gnuplot_type, pyplot_type in _linetype_table:
        if how == code:
            return gnuplot_type
    return how

def PyplotLinetype(how):
    for code, gnuplot_type, pyplot_type in _linetype_table:
        if how == code:
            return pyplot_type
    return how

def PrintLinetypes():
    print('Columns in input are specified as:  <col_group>,<col_group2>,...')
    print('  A column group is either <col>, <col>:<linetype>, or <beg_col>-<end_col>[:<linetype>]')
    print('  Columns are zero-based numeric indices.  Ranges are inclusive.')
    print('Supported linetypes depend on the driver used:')
    print('  {:10s} {:18s} {:s}'.format('linetype', 'gnuplot', 'pyplot'))
    for code, gnuplot_type, pyplot_type in _linetype_table:
        print('  {:10s} {:18s} {}'.format(code, gnuplot_type, pyplot_type))

def LogRange(x):
    min, max = np.nanmin(x), np.nanmax(x)
    return max/min if min > 0 else 0

def GetYear(date):
    if type(date) is str:  # ex: '20210330' or '20210318:1600'
        return int(date[0:4])
    else: # ex: 20210319.0
        date = int(date)
        return int(date/10000)

def GnuplotBestOutputDateFormat(beg, end):
    years = GetYear(end) - GetYear(beg)
    #print(f'XXX {beg}-{end}: {years} years')
    if   years > 10: return '%Y'
    elif years > 1:  return '%y/%m'
    else:            return '%y/%m/%d'

#
# Gnuplot driver
#

def ComputeInputTimeFmt(x):
    import re
    if   re.search(r'\d{8}$',                x): return '%Y%m%d'
    elif re.search(r'\d{8}\.0$',             x): return '%Y%m%d.0'
    elif re.search(r'\d{4}-\d{2}-\d{2}$',    x): return '%Y-%m-%d'
    elif re.search(r'\d{4}/\d{2}/\d{2}$',    x): return '%Y-%m-%d'
    elif re.search(r'\d{4}-\d{2}-\d{2}T00:', x): return '%Y-%m-%dT00:00:00.000000000' # 2021-10-04T00:00:00.000000000
    else:                                        return '%Y%m%d' # XXX

def ExecuteGnuplot(x, yy, xname, ynames, cms, yhows, title, ystats, errorbars, opt):
    #print(f'XXX xname={xname} ynames={ynames}')
    xdate = not opt.seqnum and IsDate(xname)
    script = []
    script.append(f"set datafile missing 'nan'")
    if not opt.hgram:
        script.append(f"set xlabel '{xname}'")
        if xdate:
            script.append(f'set xdata time')
            #in_fmt = '%Y%m%d'
            in_fmt = ComputeInputTimeFmt(str(x[0]))
            script.append(f'set timefmt "{in_fmt}"')
            #print(f'XXX x={x}')
            out_fmt = GnuplotBestOutputDateFormat(x[0], x[-1])
            script.append(f'set format x "{out_fmt}"')
    if   opt.nolegend:  script.append('unset key')
    elif not opt.hgram: script.append('set key autotitle columnhead')
    ny = yy.shape[1]
    if opt.logscale:
        if 'X' in opt.logscale or 'Y' in opt.logscale:
            range_threshold = 1000
            auto_logscale = ''
            if 'X' in opt.logscale:
                r = LogRange(x)
                #print(f'XXX LogRange({xname})={r}')
                if r > range_threshold:
                    auto_logscale += 'x'
            if 'Y' in opt.logscale:
                for iy in range(ny):
                    if LogRange(yy[iy]) > range_threshold:
                        auto_logscale += 'y'
                        break
            if auto_logscale:
                script.append(f'set logscale {auto_logscale}')
        else:
            script.append(f'set logscale {opt.logscale}')
    if opt.xzeroaxis:
        script.append('set xzeroaxis')
    extra = f'smooth {opt.smooth}' if opt.smooth else ''
    if opt.output is not None:
        script.append('set terminal pdfcairo noenhanced')
        output = opt.output if opt.nplots == 1 else opt.output.replace('.pdf', f'.{global_plot_idx}.pdf')
        #print(f'# {{{output}}}')
        script.append(f"set output '{output}'")
        script.append(f"set key title '{title}'")
    else:
        script.append(f"set term qt noenhanced title '{title}'")
    if opt.noaxes:
        script.append('unset xtics')
        script.append('unset ytics')
        script.append('unset xlabel')
        script.append('unset border')
        script.append('unset title')
    with tempfile.NamedTemporaryFile(delete=not opt.keeptmp) as tmp:
        line = ' '.join([xname] + ynames) + '\n'
        tmp.write(line.encode())
        tokens = ['' for _ in range(1 + ny if errorbars is None else 1 + 2*ny)]
        for row, xval in enumerate(x): # XXX np.savetxt or df.write_csv would be faster
            tokens[0] = str(xval)
            for iy, yval in enumerate(yy[row]):
                tokens[1 + iy] = str(yval)
            if errorbars is not None:
                for iy, errval in enumerate(errorbars[row]):
                    tokens[1 + ny + iy] = str(errval)
            line = ' '.join(tokens) + '\n'
            tmp.write(line.encode())
        tmp.flush()
        items= []
        for iy in range(ny):
            xcol, ycol, yname = 0, 1 + iy, ynames[iy]
            fname = tmp.name if 0 == iy else ''
            how = yhows[iy]
            items.append(f"'{fname}' using {xcol + 1}:(${ycol + 1}) {extra} {GnuplotLinetype(how)}") # 1:($2) will make line breaks on missing data
            if   opt.hgram: items[-1] += f' title "{yname} {FormatHgramStats(cms, iy)}"'
            elif opt.rgram: items[-1] += f' title "{yname} {FormatRgramStats(cms, iy)}"' # cms actually contains mis, neff, fmean, fstd
            elif ystats:    items[-1] += f' title "{yname} {FormatYStats(ystats, iy)}"'
            #print(f'XXX {yname}({xname}) title={items[-1]}')
            if errorbars is not None:
                items.append(f"'{fname}' using {xcol + 1}:{ycol + 1}:{ycol + 1 + ny} with errorbars notitle")
        script.append(f"plot {', '.join(items)}")
        cmd_line = [opt.gnuplot, '-persist', '-e', '; '.join(script)]
        if opt.verbose:
            print(f'# cmd_line: {cmd_line}')
        subprocess.run(cmd_line, shell=False, stderr=None if opt.verbose else subprocess.DEVNULL)

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

def SmoothY(x, yy, opt): # works with nans
    bandwidth = opt.llr if opt.llr > 0 else -opt.llr*(x[-1] - x[0])
    for iy in range(yy.shape[1]):
        #print(f'XXX10 y[{iy}] before:\n{yy[:, iy]}')
        yy[:, iy] = LLR_values(x, yy[:, iy], opt.llr)
        #print(f'XXX10 y[{iy}] after:\n{yy[:, iy]}')
        if opt.hgram:
            np.clip(yy, 0, 1e20, out=yy)

def ScatterPointSize(n):
    if n < 100:   return None # pyplot default
    if n < 1000:  return 2
    if n < 10000: return 1
    return 0.5

def _ExecutePyplot(x, yy, xname, ynames, cms, hows, title, ystats, errorbars, opt):
    import matplotlib.pyplot as plt
    ny = yy.shape[1];
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if opt.logscale is not None:
        if 'x' in opt.logscale: plt.xscale('log')
        if 'y' in opt.logscale: plt.yscale('log')
    for iy in range(ny):
        y = yy[:, iy]
        how = hows[iy]
        lw = 1 if how == 'D' else 1.5 # default (2?) seems thick
        if opt.smooth is not None:
            from scipy.interpolate import interp1d
            #assert False, 'smooth option not supported by pyplot driver'
            y = interp1d(x, y, kind='cubic')(x)
        if opt.nolegend:
            label = None
        elif ystats is not None:
            label = f'{ynames[iy]} {FormatYStats(ystats, iy)}'
        else:
            label = ynames[iy]
        if errorbars is not None:
            plt.errorbar(x, y, errorbars[:, iy], label=label, linewidth=lw, linestyle='solid')
        else:
            if how == 'p':
                plt.scatter(x, y, label=label, s=ScatterPointSize(len(x)))
            else:
                plt.plot(x, y, label=label, linewidth=lw, linestyle=f'{PyplotLinetype(how)}')
    if opt.xzeroaxis:
        plt.plot(x, np.zeros(len(x)), '--', label='', color='grey', linewidth=1)
    xdate = False and nx and IsDate(names[0]) # XXX later
    if xdate:
        import datetime as dt
        x = [dt.datetime.strptime(str(int(d)),'%Y%m%d').date() for d in x]
    if opt.noaxes:
        plt.axis('off')
    else:
        if not opt.hgram:    ax.set_xlabel(xname)
        if not opt.nolegend: plt.legend()
    if opt.output is not None:
        output = opt.output if opt.nplots == 1 else opt.output.replace('.pdf', f'.{global_plot_idx}.pdf')
        print(f'# {output}')
        fig.savefig(output, bbox_inches='tight')
    else:
        plt.show(block=True) # must kill the pyplot window with 'q' or mouse, or with Ctr-C in the shell

def ExecutePyplot(x, yy, xname, ynames, cms, hows, title, ystats, errorbars, opt):
    if opt.output:
        _ExecutePyplot(x, yy, xname, ynames, cms, hows, title, ystats, errorbars, opt)
    elif True: # try fork
        pid = os.fork()
        if pid: # parent
            #import time
            #print(f'parent: child window as pid {pid}')
            #time.sleep(10)
            sys.exit(0)
        else: # child
            #print(f'enter child')
            #sys.stderr = open(os.devnull, 'w') 
            _ExecutePyplot(x, yy, xname, ynames, cms, hows, title, ystats, errorbars, opt)
            
#
# 2D data plotter
#   
def RunDataplot(fnames, columns, opt):
    assert opt.wts_col is None
    cols, hows = ParseColumns(columns, opt)
    do_merge_by_x = True
    no_wts_col = None
    df, cols = ReadDataframe(fnames, cols, no_wts_col, do_merge_by_x, opt) # as cols can be extended (open range)
    while len(hows) < len(cols):
        hows.append(hows[-1])
    # x in in df.index, df.columns are yy
    ny = len(cols) if opt.seqnum else len(cols) - 1
    assert len(df.columns) == ny*len(fnames), f'df.columns={df.columns} vs ny={ny} and {len(fnames)} files'
    if len(fnames) > 1:
        orig_hows = hows
        for _ in range(1, len(fnames)):
            hows.extend(orig_hows[1:])
    if not opt.title:
        opt.title = ' '.join(fnames)
    no_yerr = None
    no_cms = None
    if opt.yscale != 1:
        for iy in range(ny):
            df.iloc[:, iy] *= opt.yscale
    if not opt.separate:
        PlotDataframe(df, list(range(ny*len(fnames))), no_yerr, no_cms, hows, opt)
    else:
        opt.nplots = ny
        for iy in range(ny):
            PlotDataframe(df, [iy], no_yerr, no_cms, hows, opt)

def PrintExamples():
    print('# optionally: export PLOT_DRIVER=plt')
    print('  ls -l /usr/bin|grep -v total|plot -cnQw 4 # cumsum of file sizes (input without header)\n'
          + '  plot -t | plot 0-4                        # plot columns 1-4 vs column 0\n'
          + '  plot -t | plot 0-4 -cS msaD -s 500 -e 900 # cumulative plots for x in [500,900) with statistics in legend\n'
          + '  for i in 3 4 5; do plot -t1 $i > tmp$i.csv; done; plot 0,2-3 tmp{3,4,5}.csv -zc # cumulative plots from multiple files\n'
          + '  plot -t | plot 0-4 -zO bezier             # bezier-smoothed plots (gnuplot driver only)\n'
          + '  plot -t | plot 0,1,bar,baz,fee -zL 10     # same plots, mixed column spec, LLR smoothed\n'
          + '  plot -t | plot 1-4 -p                     # scatter plots\n'
          + '  plot -t | plot 1   -H                     # histograms of column 1 data\n'
          + '  plot -t | plot 2-4 -HL 0.5                # histograms of multiple columns, LLR-smoothed\n'
          + '  plot -t | plot 1-4 -Rz                    # regressograms of columns 2-4 (y) vs column 1 (x)\n'
          + '  plot -t | plot 1-4 -RzL 0.4               # regressograms smoothed with bandwidth 0.4\n'
          + '  plot -t | plot 1-4 -PRzW 5                # weighted regressograms in separate windows\n'
          + '  plot -t | plot 1-4 -rERW 5 -B 60 -L 0.4   # weighted regressograms with 60 equal-weight bins and errorbars\n'
          + '  plot -tN 500 -B 100|plot 1-500 -Rq        # spaghetti art\n'
          + '  plot -t3N 200|plot 1-199 -qxO bezier      # alpha art\n'
          + '  plot -X | head -16 | bash                 # run all of the above\n'
          + '  unplot                                    # kill all active gnuplot windows\n'
)
    
def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='Plot columnar data in files or stdin, normally separated by comma or whitespace and with a header line\n'
                                     + 'The columns to plot can be specified by name, int index, inclusive range of ints, or a mix thereof\n'
                                     + 'Column examples: foo | 1,0,5 | foo,bar,4-8,2 | 1,3- (trailing dash means through the end of available columns)\n'
                                     + 'Default action: Print columns in the input available for plotting')
    parser.add_argument('-X', '--examples',   action='store_true', help='Print use examples')
    parser.add_argument('-Y', '--show_linetypes', action='store_true', help='Print help on supported line types')
    parser.add_argument('-T', '--title',                           help='Use this window title [filename]')
    parser.add_argument('-q', '--nolegend',   action='store_true', help='Do not mark plotted line in legend')
    parser.add_argument('-z', '--xzeroaxis',  action='store_true', help='Add horizontal zero line')
    parser.add_argument('-w', '--whitespace', action='store_true', help='Split input by whitespace [comma]')
    parser.add_argument('-Q', '--seqnum',     action='store_true', help='Use sequence number for x data [first column]')
    parser.add_argument('-P', '--separate',   action='store_true', help='Make plot for each data column separate [together]')
    parser.add_argument('-n', '--noheader',   action='store_true', help='Indicate that data has no header.  Generate header names F0, F1, etc')
    parser.add_argument('-c', '--cumsum',     action='store_true', help='Plot cumulative sums')
    parser.add_argument('-O', '--smooth',     type=str,            help='Apply gnuplot smoothing, supported types: s?bezier|(c|mc|ac)splines|(f|c)normal|kdensity <bwidth>')
    parser.add_argument('-S', '--stats',      type=str,            help='Add stats (before any cumsum or diff) to plots: c=count, m=mean, s=sdev, a=sharpe, d=Drawdown, D=DrawdownToSdev, A=all, P=statPrintStdout')
    parser.add_argument('-s', '--start',      type=str,            help='Skip data for x < start')
    parser.add_argument('-e', '--end',        type=str,            help='Skip data for x >= end')
    parser.add_argument('-x', '--noaxes',     action='store_true', help='Generate plot without axes or border')
    parser.add_argument('-L', '--llr', metavar='bandwidth', type=float, help='Smooth data using local linear regression (LLR) over this bandwidth. Weird (buggy?) with --equal_wt option')
    parser.add_argument('-d', '--diff',       action='store_true', help='Plot differences')
    parser.add_argument('-f', '--filtery',    type=str,            help='Filter (transform) y data through an expression to be applied via eval().  Examples: "y-y**3", "np.sin(y)"')
    parser.add_argument('-p', '--points',     action='store_true', help='Plot with points (scatter plot) [with lines]')
    parser.add_argument('-D', '--dots',       action='store_true', help='Plot with small dots (scatter plot) [with lines]')
    parser.add_argument('-l', '--logscale',   type=str,            help='Use log scale for x|y|xy.  Uppercase means select automatically based on data.  Useful with --separate option')
    parser.add_argument('-F', '--ols_fit',    action='store_true', help='Add an OLS linear fit to each plotted y(x) line')
    parser.add_argument('-H', '--hgram',      action='store_true', help='Plot histogram of indicated data column(s)')
    parser.add_argument('-R', '--rgram',      action='store_true', help='Plot regressogram of data columns: first treated as x, rest as y(s)')
    parser.add_argument('-M', '--esm',        action='store_true', help='For regressogram with yerr: plot error bars of sample mean (x 1/sqrt(Neff)) rather than yerr')
    parser.add_argument('-a', '--swapxy',     action='store_true', help='For rgram: plot first column (x turned y) vs other indicated columns (yy turned xx) separately')
    parser.add_argument('-y', '--yscale',     type=float,          help='Scale y values for rgram and dataplot by this factor', default=1)
    parser.add_argument('-m', '--xtrim',      type=float,          help='For hgrm or rgram: trim this percentage of the x distribution tails on both sides [0.001]', default=0.001)
    parser.add_argument('-B', '--nbins',                           help=f'For histogram: use this many bins: sqrt: sqrt(size)/2, qbrt: size**(1/3), or <int> [{HIST_DEFAULT_NBINS}]', default=HIST_DEFAULT_NBINS)
    parser.add_argument('-W', '--wts_col',                         help='For hgram or rgram: use this column for weights')
    parser.add_argument('-r', '--yerr',       action='store_true', help='For regressogram: plot with yerrorbars')
    parser.add_argument('-E', '--equal_wt',   action='store_true', help='Use equal-weight (histo|regresso)gram bins. Implies a density plot for histogram [equal-size].  Can look weird for highly peaked histogram')
   #parser.add_argument('-g', '--grep_key',                        help='Skip input lines without this word (cna be done with a grep/csv pipe')
   #parser.add_argument('-M', '--multikey_col', nargs=1,           help='Plot separate lines for for each value in this column (TODO or can be done with csv -W pipe)')
    parser.add_argument('-o', '--output',                          help='Plot to this pdf file [qt window].  Seems to work better with -D plt option')
    parser.add_argument('-v', '--verbose',    action='store_true', help='Print gnuplot command line')
    parser.add_argument('-K', '--keeptmp',    action='store_true', help='Keep tmp file for gnuplot (helpful with -v option) [delete]')
    parser.add_argument('-t', '--test_csv',   action='store_true', help='Print a csv stream for testing')
    parser.add_argument('-3', '--test_parametric', action='store_true', help='For test_csv: Print data for parametric plots')
    parser.add_argument('-C', '--test_clean', action='store_true', help='For test_csv: Print simple clean data for plots')
    parser.add_argument('-1', '--seed',       type=int,            help='Use this seed for random test_csv data', default=0)
    parser.add_argument('-2', '--driver',                          help=f'Use this backend graphics driver: gnuplot (gnu) or pyplot (plt) [$PLOT_DRIVER (if set) or {_DEFAULT_DRIVER}]', default=DEFAULT_DRIVER)
    parser.add_argument('-g', '--gnuplot',                         help=f'Use this gnuplot binary [gnuplot]', default='gnuplot')
    parser.add_argument('-N', '--ntests',     type=int,            help=f'test_csv: print this many columns [{DEFAULT_NTESTS}]', default=DEFAULT_NTESTS)
    parser.add_argument('files_and_columns',  nargs='*')
    opt = parser.parse_args()

    if opt.rgram and opt.equal_wt:
        pass # assert not opt.llr, 'llr option conflicts with equal_wt (a bug to fix)'
    if opt.examples:
        PrintExamples()
    elif opt.show_linetypes:
        PrintLinetypes()
    elif opt.test_csv:
        PrintTestDataframe(opt)
    else:
        opt.how = 'p' if opt.points else 'D' if opt.dots else 'l' # plot with lines by default
        fnames = []
        columns = None
        for arg in opt.files_and_columns:
            if os.path.isfile(arg) or arg == '-':
                fnames.append(arg)
            else: # assume this is columns spec (numbers, commas, dashes, with optional linetypes)
                assert columns is None, f'multiple columns specs ({arg} after {columns}) not supported'
                columns = arg
        if not fnames:
            fnames.append('-')
        opt.nplots = 1 # initial/default
        if columns is None:          PrintColumns(fnames[0], opt)
        elif opt.hgram or opt.rgram: RunHRgram(fnames, columns, opt)
        else:                        RunDataplot(fnames, columns, opt)
        if opt.output is not None and opt.nplots > 1:
            from PyPDF2 import PdfFileMerger
            pdfs = [opt.output.replace('.pdf', f'.{idx}.pdf') for idx in range(opt.nplots)]
            merger = PdfFileMerger()
            count = 0
            for pdf in pdfs:
                if os.path.isfile(pdf) and os.stat(pdf).st_size > 0:
                    merger.append(pdf)
                    count += 1
            merger.write(opt.output)
            print(f'# {opt.output} from {count} files')
            merger.close()
            for pdf in pdfs:
                if os.path.isfile(pdf):
                    os.unlink(pdf)

if __name__ == '__main__':
    try:
        main()
    except BrokenPipeError:
        pass
