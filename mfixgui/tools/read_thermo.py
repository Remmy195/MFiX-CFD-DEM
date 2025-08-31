#!/usr/bin/env python

# This reads THERMO DATA in Burcat or CHEMKIN formats


import os
import sys
import pickle
import time
import re
import warnings

re_float_exp = re.compile(r"""
    [-+]?           # optional leading sign
    (\d*[.])?       # optional zero or more digits in front of an optional decimal point
    \d+             # non-optional digits
    ([eE][-+]?\d+)? # optional exponent
""", re.VERBOSE)

name = None

def is_data_block(lines):
    if len(lines) != 4:
        return False
    if not all(75 <= len(line) <= 90 for line in lines):
        return False
    if (lines[0].endswith('1')
        and lines[1].endswith('2')
        and lines[2].endswith('3')
        and lines[3].endswith('4')):
        return True
    return False

def parse_section_burcat(section):
    global name
    section = [expand_tabs(l).strip() for l in section]
    n_lines = len(section)
    start_line = 0
    comment = ""
    comment_block = []
    while start_line < n_lines:
        tcom = 1000.0
        h_f = None
        if not is_data_block(section[start_line:start_line+4]):
            comment_block.append(section[start_line])
            start_line += 1
            continue
        #print("COMMENT:", '\n'.join(comment_block))
        #print("DATA:", '\n'. join(section[start_line:start_line+4]))
        l1, l2, l3, l4 = (l[:-1].strip() for l in section[start_line:start_line+4])
        comment = '\n'.join(comment_block)
        start_line += 4
        if not is_data_block(section[start_line:start_line+4]):
            # NB, single comment applies to multiple data blocks if not
            # separated by blank line
            comment_block = []

        name = l1[:18].strip()
        com = l1[18:44]
        lower = com.lower()
        if 'tcom=' in lower:
            idx = lower.index('tcom=')
            tok = lower[idx+5:].split()
            if tok:
                tok = tok[0].strip()
                tcom = float(tok)
        #elif len(l1) > 80: # I do not like this - see issues/1631 cgw 2022-05-10
        #    tok = l1.split()[-1]
        #    tcom = float(tok)
        phase = l1[44]
        if phase == ' ':
            phase = l1[43]  # C2HCl5
        words = l1[45:].split()
        tmin = float(words[0])
        tmax = float(words[1])
        mol_weight = words[-1]
        # There might be some non-numeric characters abutting the last field
        while mol_weight and not mol_weight[0].isdigit():
            mol_weight = mol_weight[1:]
        mol_weight = float(mol_weight)
        coeffs = list(extract_coeffs(' '.join((l2, l3, l4))))
        if len(coeffs) == 15:
            h_f = coeffs[14]
        coeffs = coeffs[7:14] + coeffs[:7]

        composition = {} # No composition data for Burcat

        yield (name, phase, [tmin, tcom, tmax], coeffs, composition, mol_weight, h_f, comment)

    if comment_block and any(l.strip() for l in comment_block):
        warnings.warn("Ignored extra lines in thermodynamic data:"+'\n'.join(comment_block))

    # TODO handle cross references

def parse_section_chemkin(section, atomic=None):
    global name
    section = [expand_tabs(l).strip() for l in section]
    n_lines = len(section)
    start_line = 0
    comment = "" # Do we have comments for Chemkin data?
    while start_line < n_lines:
        tcom = 1000.0 # ? default for Chemkin?
        composition = {}
        mol_weight = None
        h_f = None
        if not section[start_line]: # Blank line
            start_line += 1
            continue
        # Simple Chemkin format
        if is_data_block(section[start_line:start_line+4]):
            l1, l2, l3, l4 = (l[:-1].strip() for l in section[start_line:start_line+4])
            name = l1[:18].strip()
            phase = l1[44]
            words = l1[45:].split()
            tmin = float(words[0])
            tmax = float(words[1])
            tcom = float(words[2])
            temps = [tmin, tcom, tmax]
            date = l1[18:24] # Not used
            coeffs = list(extract_coeffs(' '.join((l2, l3, l4))))
            formula = l1[24:44]
            while formula:
                part, formula = formula[:5], formula[5:]
                sym = part[:2]
                num = part[2:]
                if not sym.strip():
                    break
                sym = sym.strip().title()
                if atomic:
                    if sym in atomic:
                        if float(num): # Filter out 0 counts
                            composition[sym] = float(num)
                    else:
                        raise ValueError("Unrecognized element %s in data for %s" %
                                         (sym, name))
            if atomic:
                mol_weight = sum(v*atomic[k][2] for (k,v) in composition.items())
                mol_weight = round(mol_weight, 8)
            if len(coeffs) == 15:
                h_f = coeffs[14]
            start_line += 4
            coeffs = coeffs[7:14] + coeffs[:7] # Low temp range first
            temps = [tmin, tcom, tmax]

        elif section[start_line].endswith('&'):
            l1 = section[start_line]
            name = l1[:18].strip()
            phase = l1[44]
            words = l1[45:].split()
            tmin = float(words[0])
            tmax = float(words[1])
            tcom = float(words[2])
            start_line += 1
            comp_lines = []
            while start_line < n_lines and section[start_line].endswith("&"):
                comp_lines.append(section[start_line][:-1])
                start_line += 1
            comp_lines.append(section[start_line]) # Last line does not have "&"
            if start_line < n_lines:
                start_line += 1
            else:
                raise ValueError("EOF reached reading thermo data for %s" % name)
            comp_text = ' '.join(comp_lines)
            words = comp_text.split()
            if len(words) % 2:
                raise ValueError(comp_text,
                                 "Odd number of terms in composition data for %s" % name)
            while words:
                sym, num, words = words[0], words[1], words[2:]
                sym = sym.strip().title()
                if atomic:
                    if sym in atomic:
                        num = float(num)
                        if num:
                            composition[sym] = float(num)
                    else:
                        raise ValueError("Unrecognized element %s in data for %s" %
                                         (sym, name))
            # Multiple temperature ranges
            temps = [tmin, tcom, tmax]
            if section[start_line].startswith("TEMP"):
                words = section[start_line].split()
                temps = [float(w) for w in words[1:]]
                temps.sort()
                if start_line < n_lines:
                    start_line += 1
                else:
                    raise ValueError("EOF reached reading thermo data for %s" % name)
                coeffs = []
                for n in range(2*(len(temps)-1)):
                    coeffs.extend(list(extract_coeffs(section[start_line])))
                    if start_line < n_lines:
                        start_line += 1
                    else:
                        raise ValueError("EOF reached reading thermo data for %s" % name)
                # Reverse coeffs in blocks of 7
                groups = [coeffs[i:i+7] for i in range(0, len(coeffs), 7)]
                coeffs = [coeff for group in reversed(groups) for coeff in group]

            else: # Single range
                lines = section[start_line:start_line+3]
                if not(len(lines)==3
                       and all (75 <= len(line) <= 90 for line in lines)
                       and lines[0].endswith('2')
                       and lines[1].endswith('3')
                       and lines[2].endswith('4')):
                    raise ValueError("Invalid thermo data: %s" % section[start_line])
                coeffs = list(extract_coeffs(' '.join((l2, l3, l4))))
                if len(coeffs) == 15:
                    h_f = coeffs[14]
                    coeffs = coeffs[7:14]+coeffs[:7]


        else:
            if start_line < n_lines:
                line = section[start_line].lower().strip()
                # Ignore global "THERMO" values
                if line == 'thermo':
                    start_line += 2
                    continue
                if '!' in line:
                    line = line.split('!')[0].strip()
                if line and line != 'end' and not line.startswith('thermo data'):
                    warnings.warn("Ignored extra line in thermodynamic data:"+'\n'+section[start_line])
            start_line += 1
            continue

        if atomic:
            mol_weight = sum(v*atomic[k][2] for (k,v) in composition.items())
        # Compute h_f
        T = 298.15
        a = coeffs
        h_f = a[5] + sum((a[i-1]*T**i)/i for i in range(1,6))


        yield (name, phase, temps, coeffs, composition, mol_weight, h_f, comment)


def extract_coeffs(text):
    '''extract E15.0 formatted floats from line of text.'''
    # Note that the values may lack separating whitespace!
    text = text.upper().replace('D', 'E')  # could be lower case
    text = re.sub(" +E", "E", text) # Really?  0.3332728 E+05 ??
    text = re.sub("E +", "E", text) # Really? .24753989E 01 ??
    while text:
        text = text.lstrip()
        if not text:
            break
        if text.startswith('N/A'):
            text = text[3:]
            yield None
            continue
        match = re_float_exp.match(text)
        if not match:
            msg = "Cannot parse entry for %s: %s" % (name, text)
            raise ValueError(msg)
        l = match.end()
        word = text[:l]
        val = float(word)
        text = text[l:]
        yield val

def expand_tabs(line):
    r = ''
    for c in line:
        if c == '\t':
            l = len(r)
            n = 8 * (1 + int(l / 8))
            r += ' ' * (n - l)
        else:
            r += c
    return r


def main():
    global data

    if len(sys.argv) != 3:
        print("Usage: %s infile outfile" % sys.argv[0])

    infile_name = sys.argv[1]
    outfile_name = sys.argv[2]
    if os.path.exists(outfile_name):
        print("Not clobbering %s" % outfile_name)
        sys.exit(-1)

    infile = open(infile_name, 'r', encoding='ascii')
    outfile = open(outfile_name, 'wb')

    prelude = True
    data = {}
    section = []
    for line in infile:
        line = expand_tabs(line).strip()
        if prelude:
            if line.startswith('THERMO DATA'):
                prelude = False
            continue

        if "REST IN PEACE" in line:
            time.sleep(1)

        if line:
            section.append(line.rstrip())
        else:  # sections separated by blanks
            if not section:  # extra blank lines
                continue
            for (name, phase, temps, coeffs, composition, mol_weight, h_f, comment) in parse_section_burcat(section):
                tmin, tcom, tmax = temps
                key = (name, phase, tmin, tmax)
                if key in data:
                    # uniquify names? not now, just skip
                    print("duplicate key %s, skipping" % str(key))
                else:
                    data[key] = (coeffs, mol_weight, h_f, tcom, comment)
            section = []

    pickle.dump(data, outfile)


if __name__ == '__main__':
    main()
