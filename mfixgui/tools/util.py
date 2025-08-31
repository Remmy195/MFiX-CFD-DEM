# -*- coding: utf-8 -*-

"""utility functions which do not depend on Qt"""

import copy
import glob
import math
import operator
import os
import re
import subprocess
from mfixgui.constants import *
from pdb import set_trace
from functools import reduce

def find_project_file(dirname):
    """ returns *.mfx or mfix.dat in dirname, or None if no project file exists in dirname"""
    mfx_files = glob.glob(os.path.join(dirname, '*.mfx'))
    if mfx_files:
        return mfx_files[0]
    dat_files = glob.glob(os.path.join(dirname, 'mfix.dat'))
    if dat_files:
        return dat_files[0]
    return None

def solver_name(solver):
    return ('Single-phase' if solver==SINGLE
            else 'TFM' if solver==TFM
            else 'DEM' if solver==DEM
            else 'CGP' if solver==CGP
            else 'SQP' if solver==SQP
            else 'GSP' if solver==GSP
            else 'PIC' if solver==PIC
            else 'Hybrid' if solver==HYBRID
            else '???')


def get_run_name_from_file(filename):
    """Get run name from file without loading it.
    For a loaded project, use get_value('run_name')"""

    run_name = ""
    if not filename:
        return run_name
    try:
        f = open(filename, encoding="utf-8", errors="replace")
    except:
        return run_name
    for line in f.readlines():
        line = line.strip()
        if line.lower().startswith("run_name") and "=" in line:
            tok = line.split("=", 1)[1].strip()
            # Remove quotes if present.
            # NB: Due to a bug, in some files, run_name may lack quotes
            if tok.startswith('"') or tok.startswith("'"):
                tok = tok[1:]
            if tok.endswith('"') or tok.endswith("'"):
                tok = tok[:-1]
            run_name = tok.strip().replace(" ", "_")
            break
    return run_name


def format_key_with_args(key, args=None):
    if args is not None and args != [] and args != ():
        if isinstance(args, int):
            args = [args]
        return "%s(%s)" % (key, ','.join(str(a) for a in args))
    else:
        return str(key)


def parse_key_with_args(string):
    # inverse of format_key_with_args
    if string.endswith(')'):
        key, args = string[:-1].split('(')
        args = [int(arg) for arg in args.split(',')]
    else:
        key = string
        args = []
    return key, args


def plural(n, word):
    fmt = "%d %s" if n == 1 else "%d %ss"
    return fmt % (n, word)


def drop_row_column_triangular(a, n, r):
    # Inputs:
    #  a :  upper-triangular n by n matrix represented as a list (n)(n+1)/2
    #  n :  size of input
    #  r :  index of row/column to drop, 1-based
    # Return:
    #  list representing matrix "a" with row/column "r" removed
    ret = []
    i = j = 1
    for x in a:
        if i != r and j != r:
            ret.append(x)
        j += 1
        if j > n:
            i += 1
            j = i
    return ret


def append_row_column_triangular(a, n, fill_value=None):
    # Append a row and column to an upper-triangular rank-n matrix
    ret = []
    i = j = 1
    for x in a:
        ret.append(x)
        if j == n:
            ret.append(fill_value)
        j += 1
        if j > n:
            i += 1
            j = i
    ret.append(fill_value)
    return ret

def get_unique_string(base, listofstrings):
    "uniquify a string"
    if base in listofstrings:
        # look for number at end
        nums = re.findall(r'[\d]+$', base)
        if nums:
            number = int(nums[-1]) + 1
            base = base.replace(nums[-1], '')
        else:
            number = 1
        base = get_unique_string(base + str(number), listofstrings)

    return base


def num_to_time(time, unit='s', outunit='time', digits=0):
    """Convert time with a unit to another unit."""
    unit = unit.lower()
    time = float(time)
    if not (isinstance(time, (float, int)) and math.isfinite(time)):
        return '∞'
    if time < 0:
        time = 0
    # convert time to seconds
    if unit in ['d', 'days']:
        time *= 60 * 60 * 24
    elif unit in ['h', 'hr']:
        time *= 60 * 60
    elif unit in ['m', 'min']:
        time *= 60

    if outunit == 'time':
        d,h,m,s = reduce(
            lambda ll, b: divmod(ll[0], b) + ll[1:], [(time, ), 60, 60, 24])

        t = '{:02d}:{:02d}:{:02d}:{:02d}'.format(*(int(t) for t in (d,h,m,s)))
        if digits > 0:
            fraction = round(time - int(time), digits)
            t += ('%.*f' % (digits, fraction))[1:] # strip leading 0
        return t

    elif outunit in ['d', 'days']:
        return time / (60.0 * 60.0 * 24.0)

    elif outunit in ['h', 'hr', 'hrs']:
        return time / (60.0 * 60.0)

    elif outunit in ['m', 'min', 'mins']:
        return time / (60.0)

    else:
        return time


def topological_sort(dependency_dict):
    '''
    Sort the dependency tree.
    Inspired by: http://code.activestate.com/recipes/578272-topological-sort/
    '''

    data = copy.deepcopy(dependency_dict)

    # Ignore self dependencies.
    for k, v in data.items():
        v.discard(k)
    # Find all items that don't depend on anything.
    extra_items_in_deps = reduce(set.union, data.values()) - set(data.keys())
    # Add empty dependences where needed
    data.update({item: set() for item in extra_items_in_deps})
    while True:
        ordered = set(item for item, dep in data.items() if not dep)
        if not ordered:
            break
        yield ordered
        data = {item: (dep - ordered)
                for item, dep in data.items()
                if item not in ordered}
    assert not data, "Cyclic dependencies exist among these items:\n%s" % '\n'.join(repr(x) for x in data.items())


def sort_dict(dict_, key, start=0):
    """given an dict of dicts and a key, sort the outside dict based on the
    value of one of the the internal dict's keys and return the sorted dict"""
    sorted_enum = enumerate(
        sorted(
            [(k, v[key]) for k, v in dict_.items()],
            key=operator.itemgetter(1)), start)
    return dict([(k, dict_[old_k]) for k, (old_k, v) in sorted_enum])


def extract_config(path):
    '''Given a path to a file, extract the section that starts with ## CONFIG
    and ends with ## END CONFIG'''
    config = []
    script = []
    in_config = False
    with open(path, encoding='utf-8', errors='replace') as f:
        for line in f:
            l = line.rstrip()
            if '## CONFIG' in l:
                in_config = True
            elif '## END CONFIG' in l:
                in_config = False
            elif in_config:
                config.append(l)
            else:
                script.append(l)
    return '\n'.join(config), '\n'.join(script + [''])


def replace_with_dict(string, replacements):
    '''given a string and a dict, replace all dict.key found in the string
    with the corresponding dict.value'''

    for key, value in replacements.items():
        string = string.replace('${' + key + '}', str(value))
    return string


def is_vnc():
    """determine if the gui is running in vnc"""
    if os.name == 'nt':
        return False
    try:
        xdpyinfo = subprocess.Popen(
            'xdpyinfo', stdout=subprocess.PIPE).communicate()[0]
        return 'vnc' in str(xdpyinfo)
    except:
        return False


def convert_string_to_python(string):
    """Attempt to convert strings to python types"""

    if not string:
        return ''

    # remove all quotes
    string = string.replace("'", '').replace('"', '')
    # remove any leading or trailing space, after removing quotes
    string = string.strip()
    # Remove comma separators if present
    if string.endswith(','):
        string = string[:-1]

    # lower-case version of string
    s_low = string.lower()

    if s_low in ('.t.', '.true.'):
        return True
    elif s_low in ('.f.', '.false.'):
        return False

    # Maybe it's a number (would a regex be better?)
    if any(char.isdigit() for char in string):
        try:
            return int(string)
        except ValueError:
            try:
                return float(string)
            except ValueError:
                pass

    # default - return string unchanged
    return string


def safe_float(value, default=0.0):
    """try to convert the value to a float, if ValueError or None, return default"""
    if value is None:
        return default
    try:
        return float(value)
    except ValueError:
        return default
    except TypeError:
        return default


def safe_int(value, default=0):
    """try to convert the value to a int, if ValueError or None, return default"""
    if value is None:
        return default
    try:
        return int(value)
    except ValueError:
        return default
    except TypeError:
        return default

def case_insensitive(pattern):
    def either(c):
        return '[%s%s]'%(c.lower(),c.upper()) if c.isalpha() else c
    return ''.join(map(either, pattern))


def break_string(s, w):
    if len(s) <= w:
        return (s, '')
    else:
        while w > 0 and s[w] != ' ':
            w -= 1
        if w == 0: # no space found, can't break
            return  (s, '')
        s, rest = s[:w], s[w:]
        if len(rest) <= 6: # Prefer a slightly long line to a 'dangler'
            return (s+rest, '')
        else:
            return (s, rest)

def dequote(s):
    s = s.strip()
    while s.startswith(('"', "'")):
        s = s[1:].strip()
    while s.endswith(('"', "'")):
        s = s[:-1].strip()
    return s

def unicode_greek(text):
    ### epsilon
    for t in 'epsilon', 'Epsilon': # 'EPSILON'
        text = text.replace(t, 'ε')
    ### sigma
    for t in 'sigma', 'Sigma':
        text = text.replace(t, 'σ')
    ### theta
    for t in 'theta', 'THETA', 'Theta':
        text = text.replace(t, 'θ')
    ### lambda
    for t in 'lambda', 'LAMBDA', 'Lambda':
        text = text.replace(t, 'λ')

    return text

def format_description(desc, highlight_term = None):
    desc = desc.strip()

    ### Note, order of replacements is important here!
    #We do not want to do replacements on our html tags, such as background-color

    # '>' and '<' will mess up HTML
    desc = desc.replace('<', '&lt;')
    desc = desc.replace('>', '&gt;')

    # Only SI is supported in the GUI
    desc = desc.replace(' in SI]', ']')

    # Remove default value (doesn't make sense once key is set)
    if desc.endswith(']') or desc.endswith('].'):
        dot = desc.endswith('.')
        desc = desc.rsplit('[', 1)[0]
        desc = desc.strip()
        if dot and not desc.endswith('.'): # Put it back
            desc += '.'

    # Bullets
    desc = re.sub(r'^\s*-', '<br/>&bull;', desc, flags=re.MULTILINE)
    desc = re.sub(' -- ', '<br/>&bull;', desc, flags=re.MULTILINE)

    # Non-breaking hyphen - do this before using 'background-color'
    desc = desc.replace('-', '&#8209;')

    # Don't split diff. eq's over multiple lines
    pat = re.compile('condition: *(\n.+[\n,])')
    match = pat.search(desc)
    if match or desc.startswith('Partial slip'):
        if match:
            text = match.group(1)
            repl = text
            if repl.endswith(','):
                repl = repl[:-1]

        else:
            text = re.search('dV[^,]*,', desc).group()
            repl = text
            if repl.endswith(','):
                repl = repl[:-1]

        repl = repl.replace('d(', 'd', 1)
        repl = repl.replace(')/', '/', 1)
        repl = repl.replace(' (', '(') # Squeeze out spaces before parens
        # Subscripts
        for sub in 'g', 's', 'M', 'ref':
            repl = repl.replace('_%s'%sub, '<sub>%s</sub>'%sub) # Make _g into a subscript
        repl = repl.replace('ScalarW', 'Sw')
        repl = repl.replace('Scalar', 'S') # Scalar -> S in equation (?)
        repl = repl.replace('d', '</i>d<i>') # de-italicize differential op
        repl = '<br><nobr><i>'+repl+'</i></nobr><br>'   # put it on its own line
        desc = desc.replace(text, repl)

        # Variable name before the BC
        pat = re.compile(', (.+), in', re.MULTILINE)
        match = pat.search(desc)
        if match:
            var = match.group(1)
            text = ', %s, in' % var
            repl = var
            repl = repl.replace('ScalarW', 'Sw')
            repl = repl.replace('Scalar', 'S')
            # Subscripts
            for sub in 'g', 's', 'M', 'ref':
                repl = repl.replace('_%s'%sub, '<sub>%s</sub>'%sub) # Make _g into a subscript
            repl = ' <i>%s</i> in' % repl
            desc = desc.replace(text, repl, 1)

        # Italicize normal vector to match equation display
        desc = desc.replace('where n is', 'where <i>n</i> is')
        if desc.startswith("Partial slip"):
            for s in 'Hw', 'vw':
                desc = desc.replace('%s'%s,
                                    '<i>%s</i>'%s)
    # remove cylindrical refs
    for x in ' (theta) direction', ' (r) direction':
        desc = desc.replace(x, '-direction')

    # Highlight search term
    if highlight_term:
        pat = re.compile('(%s)'%highlight_term, re.I)
        desc = pat.sub(r'<span style="background-color: #FFFF00">\1</span>', desc)

    # Show Greek characters
    desc = unicode_greek(desc)
    # T_ref
    desc = desc.replace("T_ref", "T<sub>ref</sub>")

    # Literature citation
    for x in ('\nSee ', ' See '):
        if x in desc:
            desc = desc.replace(x, '<br>See ')

    # Replace blank lines with break
    lines = [l if l.strip() else '<br>' for l in  desc.splitlines()]
    desc = '\n'.join(lines)

    # Issues/590 (workaround)
    desc = re.sub(r'\.\s*-', '<br/>&bull;', desc, flags=re.MULTILINE)

    # Lists
    if 'list:' in desc.lower():
        for n in range(1, 20):
            s = ' %d: ' % n
            desc = desc.replace(s, '<br/>&bull;%d: '%n)

    return desc



# Debugging hooks
def debug_trace():
    """Set a tracepoint in the Python debugger that works with Qt"""
    set_trace()


def main():
    def test_drop_add_triangular():
        a = [ 1, 2, 3, 4,
                 5, 6, 7,
                    8, 9,
                       10]
        b = drop_row_column_triangular(a, 4, 2)
        assert b == [1, 3, 4,
                        8, 9,
                           10]

        c = drop_row_column_triangular(b, 3, 3)
        assert c == [1, 3,
                        8]

        d = append_row_column_triangular(c, 2, 999)
        assert d == [1, 3, 999,
                        8, 999,
                           999]

    test_drop_add_triangular()

if __name__ == '__main__':
    main()
