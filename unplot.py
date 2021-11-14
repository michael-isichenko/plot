#!/usr/local/miniconda3/bin/python

import os, psutil, argparse, re

_DEF_GREP_KEY = 'gnuplot_qt'
DEF_GREP_KEY = 'plot' if 'PLOT_DRIVER' in os.environ and os.environ['PLOT_DRIVER'] in ('pyplot', 'plt') else 'gnuplot_qt'

def KillProcs(opt):
    for proc in psutil.process_iter():
        try:
            # Get process name & pid from process object.
            procName = proc.name()
            if opt.grep_key != procName: # more restrictive than 'not in procName'
                continue
            if opt.verbose or opt.list:
                print(f'{proc.pid} :: {procName} :: {proc.cmdline()}')
            if not opt.list:
                proc.kill()
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='Kill or list your gnuplot processes',
                                     epilog='Examples:\n'
                                     + '  unplot -l # list processes\n'
                                     + '  unplot    # kill them all\n'
    )
    parser.add_argument('-l', '--list',    action='store_true', help="List gnuplot processes, don't kill any")
    parser.add_argument('-v', '--verbose', action='store_true', help="List gnuplot processes being killed")
    parser.add_argument('-g', '--grep_key',                     help=f'Filter processes by this substring [depends on $PLOT_DRIVER or {_DEF_GREP_KEY}]', default=DEF_GREP_KEY)
    args = parser.parse_args()
    KillProcs(args)

if __name__ == '__main__':
    main()
