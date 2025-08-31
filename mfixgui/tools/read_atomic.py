#!/usr/bin/env python

import os
import sys
import pickle
import time

def main():
    if len(sys.argv) != 3:
        print("Usage: %s infile outfile" % sys.argv[0])
        sys.exit(-1)

    infile_name = sys.argv[1]
    outfile_name = sys.argv[2]
    if os.path.exists(outfile_name):
        print("Not clobbering %s" % outfile_name)
        sys.exit(-1)

    infile = open(infile_name, 'r', encoding='ascii')
    outfile = open(outfile_name, 'wb')

    data = {}
    for line in infile:
        line = line.strip()
        if not (line and line[0].isdigit()):
            continue
        print(line)
        atno, sym, name, wt = line.split()
        if wt.startswith('['): # It's radioactive, baby
            wt = wt[1:-1]
        data[sym.strip()] = (int(atno), name.strip(), float(wt))

    pickle.dump(data, outfile)


if __name__ == '__main__':
    main()
