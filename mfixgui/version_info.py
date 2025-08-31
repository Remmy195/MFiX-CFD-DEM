#!/usr/bin/env python

import os
import platform
import sys
import subprocess

from mfixgui.version import get_version
from mfixgui.tools import get_mfix_src

MFIX_VERSION=get_version()

def get_version_info():
    # pylint is so dumb sometimes.
    try:
        from qtpy import QT_VERSION # pylint: disable=import-error
    except ImportError:
        QT_VERSION = 'Unavailable'

    try:
        from qtpy import __version__ as QTPY_VERSION
    except ImportError:
        QTPY_VERSION = 'Unavailable'

    try:
        import numpy as np # pylint: disable=import-error
        NUMPY_VERSION = np.__version__
    except ImportError:
        NUMPY_VERSION = 'Unavailable'

    # try to get the vtk version
    VTK_VERSION = 'Unavailable'
    try:
        from vtk import vtkVersion # pylint: disable=import-error
        VTK_VERSION = vtkVersion.GetVTKVersion()
    except ImportError:
        pass

    # VTK 9+
    try:
        from vtkmodules.vtkCommonCore import vtkVersion # pylint: disable=import-error
        VTK_VERSION = vtkVersion.GetVTKVersion()
    except:
        pass

    OPENGL_VERSION = None
    try:
        from vtk import vtkRenderingOpenGL # pylint: disable=import-error
        OPENGL_VERSION = 1
    except ImportError:
        pass
    try:
        from vtk import vtkRenderingOpenGL2 # pylint: disable=import-error
        OPENGL_VERSION = 2
    except ImportError:
        pass

    try:
        from nodeworks import __version__ as NODEWORKS_VERSION # pylint: disable=import-error
    except ImportError:
        NODEWORKS_VERSION = 'Unavailable'

    SYSTEM_INFO = platform.platform()

    GFORTRAN_VERSION = 'Unavailable'
    try:
        status, output = subprocess.getstatusoutput(
            'gfortran --version')
        if status==0:
            GFORTRAN_VERSION = output.splitlines()[0]
    except:
        pass

    OPENMPI_VERSION = 'Unavailable'
    try:
        status, output = subprocess.getstatusoutput(
            'mpirun --version')
        if status==0:
            OPENMPI_VERSION = output.splitlines()[0].replace("mpirun (Open MPI) ","")
    except:
        pass

    CONDA_VERSION = 'Unavailable'
    try:
        status, output = subprocess.getstatusoutput(
            'conda --version')
        if status==0:
            CONDA_VERSION = output.splitlines()[0]
    except:
        pass

    VERSIONS = [
        ('MFiX', MFIX_VERSION),
        ("Conda", CONDA_VERSION),
        ('Python', sys.version.replace(os.linesep, ' ')),
        ('GFortran', GFORTRAN_VERSION),
        ('OpenMPI', OPENMPI_VERSION),
        ('Qt', QT_VERSION),
        ('Qtpy', QTPY_VERSION),
        ('Numpy', NUMPY_VERSION),
        ('Nodeworks', NODEWORKS_VERSION),
        ('VTK', VTK_VERSION),
        ('OpenGL backend', OPENGL_VERSION),
        ('Solver source', get_mfix_src()),
        ('System info', SYSTEM_INFO)]

    for label, version in VERSIONS:
        spaces = ' '*(20-len(label))
        yield "%s:%s%s" % (label, spaces, str(version).strip())

def main():
    for line in get_version_info():
        print(line)

if __name__ == '__main__':
    main()
