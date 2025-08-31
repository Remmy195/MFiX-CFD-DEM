#!/usr/bin/env python
# run this to generate 'atomic.inc':  ./read_atomic.py > atomic.inc

for line in open('ATOMIC.TXT'):
    if not line[0].isdigit():
        continue
    atno, sym, name, wt_str = line.split()
    atno = int(atno)
    if wt_str.startswith('['):
        rad = True
        wt = float(wt_str[1:-1])
    else:
        rad = False
        wt = float(wt_str)
    print('ATOMIC(%d,"%s","%s", "%s", %f, %s)%s &'%
          (atno, sym.upper(), name, wt_str, wt,
           '.True.' if rad else '.False.',
          ',' if atno<118 else ''))
