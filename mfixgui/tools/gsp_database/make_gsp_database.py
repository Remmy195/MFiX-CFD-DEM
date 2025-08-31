import mfixgui.widgets.gs_designer as gsd

import sys, os
import errno
import shutil
import signal
from glob import glob
# Note, we should honor the global 'novtk' flag.  But
#   this module is only loaded on-demand when the user clicks
#   'gs_designer', so if they try to run this without VTK it
#   will trigger an error
import vtk
import pyqtgraph as pg
import numpy as np
from math import floor, log10, hypot

from qtpy.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog, QFileDialog,
                            QFrame, QGridLayout, QGroupBox, QHBoxLayout, QLabel,
                            QProgressBar, QPushButton, QRadioButton, QSizePolicy,
                            QSlider, QSpacerItem, QSplitter, QVBoxLayout, QWidget)

from qtpy.QtGui import QColor, QBrush, QMouseEvent

from mfixgui.widgets.base import LineEdit

from mfixgui.tools.qt import get_icon, sub_icon_size, SETTINGS, find_gui, make_toolbutton

#import these as needed so as not to slow down MFiX startup
#from mfixgui.tools.gsp_generator import generate_gsp, set_progress_callback

from qtpy.QtCore import Signal, Qt, QPoint

from vtkmodules.vtkRenderingCore import (vtkActor,
                                         vtkPolyDataMapper,
                                         vtkRenderWindow,
                                         vtkRenderWindowInteractor,
                                         vtkRenderer)
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import vtkCommand
from vtkmodules.vtkCommonMath import vtkMatrix4x4
from vtkmodules.vtkFiltersSources import vtkSuperquadricSource, vtkSphereSource
from vtkmodules.vtkIOGeometry import vtkSTLReader
from vtkmodules.vtkRenderingAnnotation import vtkCubeAxesActor

import pandas as pd

def main():

# Read csv file containing all gsp settings and save in a dataframe
    df = pd.read_csv('gsp_data.csv')
    print(df)

    args = sys.argv
    qapp = QApplication(args)
    qapp.setStyle("fusion")
    p = gsd.GSPopup()
    q = [1,0,0,0]
    for index, row in df.iterrows():
        p.setup(None,row['label'])
        p.a = row['a']
        p.b = row['b']
        p.c = row['c']
        p.m = row['m']
        p.n = row['n']
        p.nx = row['nx']
        p.ny = row['ny']
        p.nz = row['nz']
        p.decimation = True # Always use decimation, seems faster and more robust
        p.show()
        p.generate_gsp()
        p.set_view(row['view'])
        p.save()
    # exit with Ctrl-C at the terminal
    print('Done. Press Ctrl-C to exit.')
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    qapp.exec_()
    qapp.deleteLater()
    sys.exit()

if __name__ == '__main__':
    main()
