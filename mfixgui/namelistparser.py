#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Parse the MFiX init_namelist files to extract documentation of the keywords.
"""

import os
import re
import xml.etree.ElementTree

from glob import glob

import warnings

from mfixgui.tools import get_mfix_src
from mfixgui.arrows import LEFT_ARROW, RIGHT_ARROW, LEFT_RIGHT_ARROW

TO_DTYPE = {"REAL": "DP", "LOGICAL": "L", "INTEGER": "I", "CHARACTER": "C"}

FROM_DTYPE = {
    "DP": "double precision",
    "C": "character",
    "L": "logical",
    "I": "integer",
}

# Values from param_mod.f
# TODO get these from source if modified
DIMS = {
    # DES_MMAX
    'DIMENSION_BC': 500,
    # DIMENSION_CTRL
    'DIMENSION_IC': 500,
    'DIMENSION_IS': 500,
    'DIMENSION_PS': 5000,
    'DIMENSION_USR': 5,
    'DIMENSION_VTK': 100,
    'DIM_EQS': 10,
    'DIM_I': 5000,
    'DIM_J': 5000,
    'DIM_K': 5000,
    # DIM_LM
    'DIM_M': 10,
    'DIM_N_G': 100,
    'DIM_N_S': 100,
    'DIM_QUADRIC': 500,
    'DIM_SCALAR': 100,
    'N_SPX': 11,
}

_cached_keyword_doc = None

def getKeywordDoc():
    """ build keyword dict from mfix source"""

    global _cached_keyword_doc
    if _cached_keyword_doc:
        return _cached_keyword_doc

    searchPathList = [namelist_f for namelist_f in
                      glob(os.path.join(get_mfix_src(), 'model', '**', '*.f'), recursive=True)
                      if 'namelist' in os.path.basename(namelist_f).lower()]

    sortedPathListTemp = searchPathList

    for _ in range(2):
        sortedPathList = []
        for fname in sortedPathListTemp:
            base = os.path.basename(fname)
            base_lower = base.lower()
            if base_lower == "init_namelist.f":
                sortedPathList.insert(0, fname)
            elif base_lower == "cartesian_grid_init_namelist.f":
                sortedPathList.insert(1, fname)
            elif base_lower == "des_init_namelist.f":
                sortedPathList.insert(2, fname)
            else:
                sortedPathList.append(fname)
        sortedPathListTemp = sortedPathList

    mfixKeywordDict = {}

    for fname in sortedPathList:
        with open(fname, encoding="utf-8", errors="replace") as f:
            txt = f.read()
            with warnings.catch_warnings() as ws:
                warnings.simplefilter("error") ### turn warnings into exceptions
                mfixKeywordDict.update(parse(txt, fname))

    _cached_keyword_doc = mfixKeywordDict
    return mfixKeywordDict


def parse(text, fname=None):
    r"""
    Read mfix namelists to generate documentation.

    returns dictionary

    !<keyword category="category name" required="true/false"
    !    tfm="true/false" dem="true/false" pic="true/false"
    !                                    legacy="true/false">
    !  <description></description>
    !  <arg index="" id="" max="" min=""/>
    !  <dependent keyword="" value="DEFINED"/>
    !  <conflict keyword="" value="DEFINED"/>
    !  <valid value="" note="" alias=""/>
    !  <range min="" max="" />
      MFIX_KEYWORD = INIT_VALUE
    !</keyword>


    UNDEFINED <- double
    UNDEFINED_I <-integer
    UNDEFINED_C <- character
    ZERO <- double
    .true./.false. <- logical
    \d* <- integer
    \d*\.\d*[Dd]?\d? <- double

    """

    ret = {}
    lines = text.splitlines()
    offset = 0
    in_block = False
    block = []
    for (lineno, line) in enumerate(lines):
        if line.startswith('!') and '<keyword' in line:
            in_block = True
            offset = lineno
        if in_block:
            for (pat, repl) in ((' & ', ' &amp; '),
                                ('<-->', LEFT_RIGHT_ARROW),
                                ('-->', RIGHT_ARROW),
                                ('<--', LEFT_ARROW)):
                line = line.replace(pat,repl)
            # Non-breaking spaces around arrows and '+', keep chem EQ's on one line
            for s in (LEFT_RIGHT_ARROW, RIGHT_ARROW, LEFT_ARROW, '+'):
                line = line.replace(' %s '%s, '&#160;%s&#160;'%s)

            block.append(line)
            if line.startswith('!') and '</keyword>' in line:

                try:
                    key, kwdict = parse_block(block)
                    if key:
                        ret[key] = kwdict

                except Exception as e:
                    if isinstance(e, xml.etree.ElementTree.ParseError):
                        lno, col = e.position
                        txt = block[lno-1] + '\n' + ' '*(1+col) + '^'
                        msg = str(e)
                        # Extract parenthesized term to simplify message
                        if '(' in msg:
                            msg = msg.split('(',1)[-1].split(')',1)[0]
                        # Fix line and column number in message because block did not
                        # start on line 1
                        if ': line' in msg:
                            msg = msg[:msg.index(': line')]

                        warnings.warn(msg + ' in %s at line %s, column %s'%(fname, lno+offset, 1+col)
                                      + '\n' +txt)
                    else:
                        warnings.warn(str(e) + ' in %s near line %s'%(fname, offset))

                finally:
                    in_block = False
                    block = []

    return ret


def parse_block(block):
    """ keyword_block is a list of lines  '<keyword ...> ... </keyword>"""

    # strip Fortran comments to find keyword and initval

    keyword, initval = None, None

    for line in block:
        line = line.split("!")[0].strip()
        if line:
            tok = line.split('=')
            if len(tok) >= 2:
                keyword, initval = tok[0], tok[1]
                break
    if keyword is None or initval is None:
        return None, {}

    init = initval.strip(" -")
    is_negative = initval.strip(" ").startswith("-")
    keyword = keyword.strip()
    if "(" in keyword:
        keyword, shape = keyword.split("(", 1)
        shape = shape[:-1]  # Remove trailing paren
    else:
        shape = ""

    nargs = 1 + shape.count(",") if shape else 0
    #  Strip ! character from the first column,
    #  remove trailing ! comment after keyword=initval,
    #  and parse the rest as XML.

    xml_block = [line.lstrip("!").split("!")[0] for line in block]
    tree = xml.etree.ElementTree.fromstring('\n'.join(xml_block))

    validrange = {}
    for range_element in tree.findall("range"):
        for bound in ("min", "max"):
            mval = range_element.get(bound)
            if mval:
                validrange[bound] = mval
                if mval.startswith('+') and bound=='min':
                    validrange['exclude_min'] = True
                    # Note, no special syntax for exclude_max
                try:
                    validrange[bound] = float(mval)
                    if int(mval) == float(mval):
                        validrange[bound] = int(mval)
                except ValueError:
                    pass

    dtype = tree.get('dtype').upper()
    initpython = find_initpython(dtype, init, is_negative)

    def get_bool_attr(elem, attr):
        val = elem.get(attr)
        return bool(val and val.lower() != "false")

    def to_int(x):
        if x.isnumeric():
            return int(x)
        return DIMS.get(x, x)

    expected_args = dict({
        to_int(arg.get("index")): {
            "id": to_int(arg.get("id")),
            "min": to_int(arg.get("min")),
            "max": to_int(arg.get("max")),
        } for arg in tree.findall("arg")})

    expected_nargs = len(expected_args)
    if expected_nargs != nargs:
        print(block)
        raise TypeError("%s: expected %d args, found %s" % (
            keyword.upper(), expected_nargs, nargs))

    return keyword.lower(), {
        'init': init,
        'initpython': initpython,
        'dtype': TO_DTYPE[dtype.upper()],
        'description': unindented_description(keyword, tree),
        'category': tree.get('category').lower(),
        'required': get_bool_attr(tree, 'required'),
        'tfm': get_bool_attr(tree, 'tfm'),
        'dem': get_bool_attr(tree, 'dem'),
        'pic': get_bool_attr(tree, 'pic'),
        'legacy': get_bool_attr(tree, 'legacy'),
        'locked': get_bool_attr(tree, 'locked'),
        'args': expected_args,
        'dependents': dict({
            dep.get('keyword'):
            {
                'keyword': dep.get('keyword'),
                'value': dep.get('value')
            } for dep in tree.findall('dependent')}),
        'conflicts': dict({
            conflict.get('keyword'): {
                'keyword': conflict.get('keyword'),
                'value': conflict.get('value')
            } for conflict in tree.findall('conflict')}),
        'valids': dict({
            valid.get('value'): {
                'value': valid.get('value'),
                'note': space_to_newline(valid.get('note')),
                'alias': valid.get('alias') or None,
            } for valid in tree.findall('valid')}),
        'validrange': validrange,
    }

def space_to_newline(text):
    t0 = text
    text = ' '.join(text.split())
    # Double-dash used for lists inside valid values, since XML
    #  does not preserve newlines
    text = text.replace(' -- ','\n\n        - ')
    return text

def unindented_description(keyword, tree):
    desc = tree.find("description").text
    lines = desc.splitlines()
    if len(lines) == 1:
        return lines[0].lstrip()
    # Deal with XML in init_namelist where first line
    # starts right after <description> tag and is therefore
    # not indented, but the rest of the block is
    if lines[0].startswith(' '):
        indent = len(lines[0]) - len(lines[0].lstrip())
        l0 = lines[0].lstrip()
    else:
        indent = len(lines[1]) - len(lines[1].lstrip())
        l0 = lines[0]
    ret = '\n'.join([l0]+[l[indent:] for l in lines[1:]])
    if ret.startswith('\n'):
        ret = ret[1:]
    return ret

def find_initpython(dtype, init, is_negative=False):
    """ return a python value from the init """

    init_lower = init.lower()

    if dtype == "CHARACTER":
        if "undefined_c" in init_lower:
            initpython = None
        else:
            initpython = init.strip("'\"")  # remove extra quotes

    elif dtype == "INTEGER":
        initpython = (None if init_lower == "undefined_i"
            else (-int(init) if is_negative else int(init)))

    elif dtype == "LOGICAL":
        if init_lower in (".true.", ".false.", ".t.", ".f."):
            initpython = init_lower in (".t.", ".true.")
        else:
            raise ValueError(init)

    elif dtype == "REAL":
        initpython = {
            "zero": 0,
            "one": 1,
            "large_number": 1.0e32,
            "undefined": None,
        }.get(init_lower, False)
        if initpython is False:
            initpython = float(init_lower.replace("d", "e"))
        if is_negative and isinstance(initpython, float):
            initpython = -initpython
    else:
        raise TypeError(dtype)

    return initpython


def build_keyword_rst(directory=None):
    kw_doc = getKeywordDoc()
    keys_by_cat = {}
    for k in kw_doc:
        cat = kw_doc[k]['category']
        if cat not in keys_by_cat:
            keys_by_cat[cat] = []
        keys_by_cat[cat].append(k)
    if directory is None:
        directory = ""

    for cat, keys in keys_by_cat.items():
        cat = cat.replace(" ", "_")
        cat = cat.replace('(','')
        cat = cat.replace(')','')
        fname = os.path.join(directory, cat + "_gen.rst")
        with open(fname, "w", encoding="utf-8", errors="ignore") as rst:
            rst.write(
                """
.. role:: kwname
.. role:: kwtype
.. role:: required

"""
            )
            for key in sorted(keys):
                kw_data = kw_doc[key]
                write_keyword_heading(rst, key, kw_data)
                if kw_data.get('valids'):
                    write_keyword_valids(rst, kw_data)

# There are only a few diff. eqs so just convert them to TeX explicitly
eq_repl = {

    r"d(Scalar)/dn + Hw (Scalar - ScalarW) = C":
    r"{\mathrm{d}S}/{\mathrm{d}n} + Hw (S - Sw) = C",

    r"d(T\_g)/dn + Hw (T\_g - T\_ref) = C":
    r"{\mathrm{d}T_g}/{\mathrm{d}n} + Hw (T_g - T_{ref}) = C",

    r"d(T\_s)/dn + Hw (T\_s - T\_ref) = C":
    r"{\mathrm{d}T_s}/{\mathrm{d}n} + Hw (T_s - T_{ref}) = C",

    r"d(THETA\_M)/dn + Hw (THETA\_M - THETAw\_M) = C":
    r"{\mathrm{d}{\Theta}_m}/{\mathrm{d}n} + Hw ({\Theta}_m - {\Theta}w_m) = C",

    r"d(X\_g)/dn + Hw (X\_g - Xw\_g) = C":
    r"{\mathrm{d}X_g}/{\mathrm{d}n} + Hw (X_g - Xw_g) = C",

    r"d(X\_s)/dn + Hw (X\_s - Xw\_s) = C":
    r"{\mathrm{d}X_s}/{\mathrm{d}n} + Hw (X_s - Xw_s) = C",

    r"dV/dn + Hw (V - Vw) = 0":
    r"{\mathrm{d}V}/{\mathrm{d}n} + Hw (V - Vw) = 0"}


term_repl = {
    r"Scalar": r"S",
    r"ScalarW": r"Sw",
    r"T\_g": r"T_g",
    r"T\_ref": r"T_{ref}",
    r"T\_s": r"T_s",
    r"THETA_M": r"{\Theta}_m",
    r"THETAw\_M": r"{\Theta}w_m",
    r"X\_s": r"X_s",
    r"X\_g": r"X_g",
    r"Xw\_s": r"Xw_s" }


def fix_math(desc):
    pat = re.compile('condition: *([^,]*),')
    match = pat.search(desc)
    desc = desc.replace("Vw = 0",r" :math:`Vw=0`")
    desc = desc.replace(" Vw=0", r" :math:`Vw=0`")
    desc = desc.replace(" Hw=0", r" :math:`Hw=0`")
    desc = desc.replace("Hw = 0", r" :math:`Hw=0`")
    desc = desc.replace("Hw = +inf", r" :math:`Hw=\infty`")
    desc = desc.replace("Hw=+inf", r":math:`Hw=\infty`")
    if not match:
        for x in ('Hw', 'Uw', 'Vw', 'Ww'):
            desc = desc.replace(' '+x+' ', " :math:`"+x+"` ")
        return desc

    pre = desc[:match.start(1)]
    post = desc[match.end(1):]
    eq = match.group(1)
    parts = pre.rsplit(',', maxsplit=1)

    if len(parts) == 2:
        pre, term = parts[0].rsplit(maxsplit=1)
        pre = pre + " :math:`" + term_repl.get(term,term) + "`, " + parts[1]
    else:
        for x in ('Hw', 'Vw'):
            post = post.replace(' '+x+' ', " :math:`"+x+"` ")
    post = post.replace(" n ", " :math:`n` ")
    eq = eq_repl[eq.strip()]
    return pre + " :math:`" + eq + "`" + post
    return desc

def write_keyword_heading(rst, keyword, data):
    args = data['args']
    datatype = FROM_DTYPE[data['dtype']].upper()

    if args:
        index_range = r"- :math:`%s \le %s \le %s`"
        indices = "\n".join(
            index_range % (v["min"], v["id"], str(v["max"]).replace("_", r"{\_}"))
            for v in args.values()
        )
        kwargs = "(" + ", ".join(a["id"] for a in args.values()) + ")"
    else:
        indices = ""
        kwargs = ""
    required = ":required:`Required`" if data["required"] else ""
    kw_heading = f":kwname:`{keyword}{kwargs}`".upper()
    description = data["description"].replace("*", r"\*").replace("_", r"\_")
    description = fix_math(description)
    flags = (
        (
            "Applies to solids model(s):  "
            + ", ".join(
                [f"**{flag.upper()}**" for flag in ("tfm", "dem") if data[flag]]
            )
        )
        if data["tfm"] or data["dem"]
        else ""
    )
    rst.write(
        f"""

.. _{keyword.upper()}:

{kw_heading}
{"~" * len(kw_heading)}

Data Type: :kwtype:`{datatype.upper()}`

{required}

{indices}

{flags}

{description}

"""
    )


def write_keyword_valids(rst, keyword_data):
    valids = keyword_data['valids']
    rst.write(
        """

.. list-table:: Valid values
   :widths: 10, 5, 40
   :header-rows: 1

   * - Name
     - Default?
     - Description """
    )
    unquoted_init = keyword_data["init"].strip("'").strip('"')
    for valid, valid_data in valids.items():
        default = "|fisheye|" if unquoted_init.lower() == valid.lower() else ""
        validdata = fix_math(valid_data["note"])
        alias = valid_data["alias"]
        if alias:
            valid += " ("+alias+")"

        rst.write(
            f"""
   * - ``{valid}``
     - {default}
     - {validdata} """
        )


def main():
    build_keyword_rst()


if __name__ == '__main__':
    main()
