#!/usr/local/miniconda3/bin/python

import psutil, argparse

def KillProcs(opt):
    for proc in psutil.process_iter():
        try:
            # Get process name & pid from process object.
            procName = proc.name()
            if opt.grep_key not in procName:
                continue
            if opt.verbose or opt.list:
                print(f'{proc.pid} ::: {proc.cmdline()}')
            if not opt.list:
                proc.kill()
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass

def main():
    def_grep_key = 'gnuplot_qt'
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='Kill or list your gnuplot processes',
                                     epilog='Examples:\n'
                                     + '  unplot -l # list processes\n'
                                     + '  unplot    # kill them all\n'
    )
    parser.add_argument('-l', '--list',    action='store_true', help="List gnuplot processes, don't kill any")
    parser.add_argument('-v', '--verbose', action='store_true', help="List gnuplot processes being killed")
    parser.add_argument('-g', '--grep_key',                     help=f'Filter processes by this substring [{def_grep_key}]', default=def_grep_key)
    args = parser.parse_args()
    KillProcs(args)

if __name__ == '__main__':
    main()
