"""Functions that return paths for MFiX project templates, documentation,
solver source code, etc."""

import platform
import sys, os


# Base dir of Python 'mfixgui' package is parent directory of 'tools' where this file lives
SCRIPT_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
CONDA_PREFIX = os.getenv("CONDA_PREFIX")

def get_mfix_templates():
    """ return directory storing MFiX tutorials """
    if CONDA_PREFIX:
        x = CONDA_PREFIX+"/share/mfix/templates"
        if os.path.exists(x):
            return x
    parent = os.path.dirname(SCRIPT_DIR)
    if os.path.exists(parent+'/tutorials'):
        return parent
    raise RuntimeError("Unable to find MFiX templates")


mfix_doc_html = None
def get_mfix_doc_html():
    """ returns the directory containing html documentation """
    global mfix_doc_html
    if CONDA_PREFIX:
        x = CONDA_PREFIX+'/share/mfix/doc/html'
        if os.path.isdir(x):
            mfix_doc_html = x
            return x
    parent = os.path.dirname(SCRIPT_DIR)
    for d in ("doc/html", # pip/tarball
              "doc/user_manual/build/html"): # source
        x = parent + '/'+d
        if os.path.isdir(x):
            mfix_doc_html = x
            return x

    return None

mfix_src = None
def get_mfix_src():
    """ returns the directory containing solver source (model, postmfix, build_aux) """
    global mfix_src
    if mfix_src:
        return mfix_src
    parent = os.path.dirname(SCRIPT_DIR)
    if os.path.isdir(parent+"/model"):
        mfix_src = os.path.normpath(parent)
        return mfix_src
    if CONDA_PREFIX:
        x = CONDA_PREFIX+"/share/mfix/src"
        if os.path.isdir(x+"/model"):
            mfix_src = os.path.normpath(x)
            return x

    raise RuntimeError("Unable to find MFiX source")


def get_full_path(project_dir, path):
    # This matches some screwy code in error_manager_mod.f
    # and should probably not be here
    # Files in the project dir start with the final component of the project
    # dir name

    if os.path.isabs(path):
        fullpath = path
    else:
        A = path.split(os.path.sep)
        B = project_dir.split(os.path.sep)

        if A[0] == B[-1]:
            fullpath = os.path.join(*[project_dir]+A[1:])
        else:
            fullpath = os.path.join(get_mfix_src(), "model", path) #?

    fullpath = os.path.realpath(fullpath)
    return fullpath
