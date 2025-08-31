#!/usr/bin/env python

#  Integrated version of gs_designer

import sys, os
import errno
import shutil
import signal
import traceback
from glob import glob
# Note, we should honor the global 'novtk' flag.  But
#   this module is only loaded on-demand when the user clicks
#   'gs_designer', so if they try to run this without VTK it
#   will trigger an error

# 'import vtk' should not be needed, but if we don't import it,
#  SetCurrentStyleToTrackballCamera fails - cgw 2024-04-04
import vtk

import numpy as np
from math import floor, log10

from qtpy.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog, QFileDialog,
                            QFrame, QGridLayout, QGroupBox, QHBoxLayout, QLabel,
                            QProgressBar, QPushButton, QRadioButton, QSizePolicy,
                            QSlider, QSpacerItem, QSplitter, QVBoxLayout, QWidget)

from qtpy.QtGui import QColor, QBrush, QMouseEvent

from mfixgui.widgets.base import LineEdit

from mfixgui.tools.qt import get_icon, sub_icon_size, SETTINGS, find_gui, make_toolbutton

#import these as needed so as not to slow down MFiX startup
#from mfixgui.tools.gsp_generator import generate_gsp, set_progress_callback

from qtpy.QtCore import Qt, QPoint

from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderer
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import vtkCommand
from vtkmodules.vtkCommonMath import vtkMatrix4x4
from vtkmodules.vtkFiltersSources import vtkSuperquadricSource, vtkSphereSource
from vtkmodules.vtkIOGeometry import vtkSTLReader
from vtkmodules.vtkRenderingAnnotation import vtkCubeAxesActor

blue = QBrush(QColor('blue'))
green = QBrush(QColor('green'))
small = 10
big = 20

SQ_TYPE, STL_TYPE = 0, 1

class GSPopup(QDialog):
    def __init__(self, parent=None):
        super(GSPopup, self).__init__(parent)
        self.parent = parent
        self.gsp_dir = 'GSP'

        # Defaults
        self.type = SQ_TYPE
        self.a, self.b, self.c = 0.01, 0.01, 0.01
        self.m = self.n = 3
        self.nx, self.ny, self.nz = 2, 2, 2
        self.wall_extend_ratio = 0.001
        self.r_ratio = 0.9
        self.decimation = False
        self.remove_small_spheres = False
        self.stl_scale = 1.0
        self.ss_ratio = 0.001
        self.ss_fraction = 0.001
        self.stl_file = None
        self.sq_resolution = 30

        # Used to show summary info in solids pane
        self.nspheres = 0
        self.equiv_diam = self.bounding_diam = self.area = self.vol = self.com = None
        self.surface_deviation = self.vol_deviation = None

        # VTK stuff
        self.sphere_actors = [] # sphere actors, TODO replace with glyph rendering
        self.base_opacity = 0.5
        self.gsp_opacity = 1.0
        self.exp = {'a':fexp(self.a),
                    'b':fexp(self.b),
                    'c':fexp(self.c),
                    'stl_scale':fexp(self.stl_scale)}

        self.q = np.array((1,0,0,0), dtype=float)
        self.setWindowTitle("Glued-sphere designer")
        self.layout = layout = QHBoxLayout()
        layout.setContentsMargins(5,5,5,5)
        self.setLayout(layout)
        splitter = self.splitter = QSplitter()
        layout.addWidget(splitter)
        self.widgets = {}
        self.P = None # Active solids phase
        self.callback = None
        self.running = False
        self.cancelled = False
        self.set_stop_flag = None
        self.view_flip = [False]*3
        self.invalidate()

        # Leftmost part of window:  the VTK viewer
        w = QWidget()
        left = QVBoxLayout()
        left.setContentsMargins(0,0,5,0)
        w.setLayout(left)
        self.splitter.addWidget(w)
        # Button bar
        bbar = QWidget(self)
        bbar_layout = QHBoxLayout(bbar)
        bbar_layout.setContentsMargins(0, 0, 0, 0)
        bbar_layout.setSpacing(0)
        bbar.setLayout(bbar_layout)
        size = sub_icon_size()
        if size.width() == 16:
            size *= 2
        for  (name,callback,tooltip) in (
                ('overscan', self.reset_view, 'Reset view'),
                ('xy', lambda: self.set_view('xy'), 'XY view'),
                ('yz', lambda: self.set_view('yz'), 'YZ view'),
                ('xz', lambda: self.set_view('xz'), 'XZ view'),
                ('rotate_left', lambda: self.rotate(90), 'Rotate counter-clockwise'),
                ('rotate_right', lambda: self.rotate(-90), 'Rotate clockwise'),
                ('perspective', lambda ignore: self.perspective(), 'Perspective'),
                ('axes', self.toggle_axes, 'Hide axes')):
            btn = make_toolbutton(name+'.svg', callback, tooltip, size)
            self.widgets['tb_'+name] = btn
            bbar_layout.addWidget(btn)
            btn.setAutoRaise(True)
            btn.setFocusPolicy(Qt.ClickFocus)
        self.widgets['tb_axes'].setCheckable(True)
        self.widgets['tb_axes'].setChecked(True)
        bbar_layout.addItem(QSpacerItem(20,10,
                                        QSizePolicy.Expanding,
                                        QSizePolicy.Minimum))

        left.addWidget(bbar)

        self.vtk = QVTKRenderWindowInteractor(self)
        self.vtk.setMinimumHeight(200)
        left.addWidget(self.vtk)

        def mk_label(text, tooltip=None):
            l = QLabel(text)
            if tooltip:
                l.setToolTip(tooltip)
            return l

        hbox = QHBoxLayout()
        b = QRadioButton('Rotate object')
        b.setToolTip("Apply rotation to base object")
        b.setChecked(False)
        b.clicked.connect(lambda val: self.handle_button(1))
        hbox.addWidget(b)
        b = QRadioButton('Rotate axes')
        b.setToolTip("Apply rotation to view")
        b.setChecked(True)
        b.clicked.connect(lambda val: self.handle_button(2))
        hbox.addWidget(b)
        hbox2 = QHBoxLayout()
        hbox2.addWidget(mk_label("Base opacity", "Opacity of base object"))
        s = QSlider(1)
        self.widgets['s_base_opacity'] = s
        s.setValue(int(100*self.base_opacity))
        s.setMaximumWidth(50)
        s.valueChanged.connect(self.set_base_opacity)
        hbox2.addWidget(s)
        hbox2.addStretch()
        hbox.addLayout(hbox2)
        hbox2 = QHBoxLayout()
        hbox2.addWidget(mk_label("GSP opacity", "Opacity of glued-sphere representation"))
        s = QSlider(1)
        self.widgets['s_gsp_opacity'] = s
        s.setValue(int(100*self.gsp_opacity))
        s.setMaximumWidth(50)
        s.valueChanged.connect(self.set_gsp_opacity)
        hbox2.addWidget(s)
        hbox2.addStretch()
        hbox.addLayout(hbox2)
        for x in range(4):
            hbox.setStretch(x,2)
        left.addLayout(hbox)

        ren = self.ren = vtkRenderer()
        rw = self.rw = self.vtk.GetRenderWindow()
        rw.AddRenderer(ren)
        iren = self.iren = rw.GetInteractor()
        style = self.style = iren.GetInteractorStyle()
        style.SetCurrentStyleToTrackballCamera()
        self.rotate_camera = True
        self.rotate_object = False
        iren._MouseMoveEvent = iren.MouseMoveEvent
        iren.MouseMoveEvent = lambda: iren._MouseMoveEvent() if not (self.rotate_object and iren.GetShiftKey()) else None

        # Disable keyboard commands
        iren.RemoveObservers('CharEvent')

        # Colors
        colors = vtkNamedColors()
        #bg_color = colors.GetColor3d("DarkSlateGray")
        object_color = colors.GetColor3d("DarkGreen")
        sphere_color = self.sphere_color = colors.GetColor3d("LightSkyBlue")
        axis_colors = [colors.GetColor3d("Black"),
                       colors.GetColor3d("Black"),
                       colors.GetColor3d("Black")]

        # Create a superquadric
        sq_source = self.sq_source = vtkSuperquadricSource()
        sq_source.SetAxisOfSymmetry(2)
        sq_source.SetSize(1)
        sq_source.SetThetaResolution(500)
        sq_source.SetPhiResolution(500)

        sq_source.SetScale(self.a, self.b, self.c)
        sq_source.SetThetaRoundness(2/self.m)
        sq_source.SetPhiRoundness(2/self.n)
        sq_source.SetCenter(0,0,0)
        sq_source.Update()  # needed to GetBounds later

        mapper = self.sq_mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(sq_source.GetOutputPort())
        actor = self.sq_actor = vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetOpacity(self.base_opacity)
        actor.GetProperty().SetDiffuseColor(object_color)
        actor.GetProperty().SetDiffuse(.7)
        actor.GetProperty().SetSpecular(.7)
        actor.GetProperty().SetSpecularPower(50.0)

        stl_source = self.stl_source = vtkSTLReader()
        mapper = self.stl_mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(None)
        actor = self.stl_actor = vtkActor()
        actor.SetScale(self.stl_scale)
        actor.SetMapper(mapper)
        actor.GetProperty().SetOpacity(self.base_opacity)
        actor.GetProperty().SetDiffuseColor(object_color)
        actor.GetProperty().SetDiffuse(.7)
        actor.GetProperty().SetSpecular(.7)
        actor.GetProperty().SetSpecularPower(50.0)

        # Setup axes
        actor = self.axes_actor = vtkCubeAxesActor()
        #actor.SetUseTextActor3D(1) less legible
        actor.SetCamera(self.ren.GetActiveCamera())
        for i,col in enumerate(axis_colors):
            actor.GetTitleTextProperty(i).SetFontSize(48)
            actor.GetTitleTextProperty(i).SetColor(col)
            actor.GetLabelTextProperty(i).SetColor(col)
        actor.DrawXGridlinesOn()
        actor.DrawYGridlinesOn()
        actor.DrawZGridlinesOn()
        actor.SetGridLineLocation(actor.VTK_GRID_LINES_FURTHEST)
        actor.XAxisMinorTickVisibilityOff()
        actor.YAxisMinorTickVisibilityOff()
        actor.ZAxisMinorTickVisibilityOff()

        #https://vtk.org/doc/nightly/html/classvtkCubeAxesActor.html#details
        #actor.SetFlyModeToStaticEdges()
        #actor.SetFlyModeToOuterEdges()
        actor.SetFlyModeToClosestTriad()
        #actor.SetFlyModeToFurthestTriad()

        # Put the VTK scene together
        ren.AddActor(self.axes_actor)
        ren.AddActor(self.sq_actor)
        ren.AddActor(self.stl_actor)

        #ren.GetActiveCamera().Azimuth(30)
        ren.ResetCamera()
        c1 = QColor(SETTINGS.value('vtk_bck_color1', '#aaffff'))
        c2 = QColor(SETTINGS.value('vtk_bck_color2', '#ffffff'))

        #ren.SetBackground(bg_color)
        ren.GradientBackgroundOn()
        ren.SetBackground(c1.getRgbF()[:3])
        ren.SetBackground2(c2.getRgbF()[:3])

        #ren.GetActiveCamera().Zoom(0.8)
        iren.Initialize()
        self.suppress_update = False
        self.reset_view()

        # Right part:  GUI controls
        row = 0
        grid = QGridLayout()
        grid.setVerticalSpacing(2)
        grid.setContentsMargins(5,5,5,5)
        w = QWidget()
        w.setLayout(grid)
        self.splitter.addWidget(w)

        self.sq_actor.AddObserver(vtkCommand.ModifiedEvent, self.handle_orientation)
        self.stl_actor.AddObserver(vtkCommand.ModifiedEvent, self.handle_orientation)

        self.type = SQ_TYPE

        row += 1
        hbox = QHBoxLayout()
        label = QLabel("Base shape ")
        label.setStyleSheet("font-weight: bold;");
        hbox.addWidget(label)
        cb = QComboBox()
        cb.addItem("Superquadric")
        cb.addItem("STL")
        self.widgets['type'] = cb
        cb.currentIndexChanged.connect(self.handle_type)
        hbox.addWidget(cb)
        b = self.button_choose_stl = QPushButton("Choose STL file")
        b.setAutoDefault(False)
        b.clicked.connect(lambda arg:self.choose_stl_file())
        hbox.addWidget(b)
        l = self.label_stl_file = QLabel(' No file selected')
        l.setStyleSheet("font-style: italic;")
        hbox.addWidget(l)
        hbox.addItem(QSpacerItem(20,10,
                                 QSizePolicy.Expanding,
                                 QSizePolicy.Minimum))

        grid.addLayout(hbox, row, 0, 1, 3)

        row += 1
        gb = self.groupbox_sq = QGroupBox()
        grid.addWidget(gb, row, 0, 1, 4)
        subgrid = QGridLayout()
        subgrid.setVerticalSpacing(2)
        subgrid.setContentsMargins(5,0,5,0)
        gb.setLayout(subgrid)

        def mk_le(key, val, min=None, max=None, tooltip=None, exclude_min=False, dtype=float):
            le = LineEdit()
            if key != 'd':
                le.required = True
            le.key = key
            le.dtype = dtype
            le.saved_value = le.default_value = val
            le.setText(str(val))
            if min is not None:
                le.min = min
            if min == 0:
                le.exclude_min = True
            if exclude_min is not None:
                le.exclude_min = exclude_min
            if max is not None:
                le.max = max
            #le.setMaximumWidth(50)
            le.setMaximumWidth(le.fontMetrics().boundingRect('9.999999').width())
            le.value_updated.connect(self.handle_le)
            if tooltip:
                le.setToolTip(tooltip)
            self.widgets['le_'+key] = le
            return le

        def mk_s(key, val, min=None, max=None, tooltip=None):
            s = QSlider(1) #horiz
            width = 400
            s.setMaximum(width)
            s.min = min
            s.max = max
            if min is not None and max is not None:
                s.setValue(int(width * (val-min)/(max-min)))
            elif key in ('a','b','c','stl_scale'):
                s.setValue(int(40*fman(val)))
                self.exp[key] = fexp(val)

            s.valueChanged.connect(lambda val, key=key: self.handle_s(key, val))
            if tooltip:
                s.setToolTip(tooltip)
            self.widgets['s_'+key] = s
            return s

        subrow = 0
        for (c,v) in (('X','a'),('Y','b'),('Z','c')):
            subgrid.addWidget(QLabel("%s semiaxis"%c), subrow, 0)
            subgrid.addWidget(mk_le(v, self.a, 0), subrow, 1)
            subgrid.addWidget(mk_s(v, self.a), subrow, 2)
            subgrid.addWidget(QLabel("m"), subrow, 3)
            subrow += 1

        subgrid.addWidget(QLabel("Shape exponent m"), subrow, 0)
        subgrid.addWidget(mk_le('m', 2, 2, 24), subrow, 1)
        subgrid.addWidget(mk_s('m', 2, 2, 24), subrow, 2, 1, 2)

        subrow +=1
        subgrid.addWidget(QLabel("Shape exponent n"), subrow, 0)
        subgrid.addWidget(mk_le('n', 2, 2, 24), subrow, 1)
        subgrid.addWidget(mk_s('n', 2, 2, 24), subrow, 2, 1, 2)

        subrow +=1
        subgrid.addWidget(mk_label("STL scale factor"), subrow, 0)
        subgrid.addWidget(mk_le('stl_scale', self.stl_scale), subrow, 1)
        subgrid.addWidget(mk_s('stl_scale', self.stl_scale), subrow, 2, 1, 2)
        self.row_scale = subrow

        row += 1
        gb = QGroupBox("Orientation (quaternion)")
        gb.setStyleSheet("QGroupBox{font-weight: bold;}");
        gb.setFlat(True)
        grid.addWidget(gb, row, 0, 1, 4)
        subgrid = QGridLayout()
        subgrid.setVerticalSpacing(2)
        subgrid.setContentsMargins(5,0,5,0)
        gb.setLayout(subgrid)
        for i in range(4):
            key = 'q%d'%(i+1)
            le = mk_le(key, 1 if i==0 else 0, -1, 1)
            f = le.font()
            le.setMaximumWidth(100)
            le.value_updated.disconnect()
            le.arg = i
            le.value_updated.connect(self.handle_q)
            label = QLabel("Q%d"%(i+1))
            subgrid.addWidget(label, i//2, i%2 * 3) # gsp_q1..q4, 1-based)
            subgrid.addWidget(le, i//2, i%2 * 3 + 1)
            subgrid.addItem(QSpacerItem(20,10,
                                        QSizePolicy.Expanding,
                                        QSizePolicy.Minimum),
                            i//2, i%2 * 3 + 2)
            box = QHBoxLayout()
            subgrid.addLayout(box,3,0,1,6)
            b = QPushButton("Reset")
            b.setToolTip("Reset orientation quaternion to default.")
            b.setAutoDefault(False)
            b.clicked.connect(self.reset_q)
            box.addWidget(b)
            b = QPushButton("Apply quaternion")
            b.setToolTip("Normalize quaternion to unit and apply rotation to object.")
            b.setAutoDefault(False)
            b.clicked.connect(self.apply_q)
            box.addWidget(b)
            box.addItem(QSpacerItem(20,10,
                                    QSizePolicy.Expanding,
                                    QSizePolicy.Minimum))
        row += 1
        gb = self.groupbox_disc = QGroupBox("Discretization")
        gb.setStyleSheet("QGroupBox{font-weight: bold;}");
        grid.addWidget(gb, row, 0, 1, 4)
        subgrid = QGridLayout()
        subgrid.setVerticalSpacing(2)
        subgrid.setContentsMargins(5,0,5,0)
        gb.setLayout(subgrid)

        for subrow, c in enumerate('xyz'):
            subgrid.addWidget(QLabel("Spheres in %s direction"%c.upper()), subrow, 0)
            subgrid.addWidget(mk_le('n'+c, 2, 1, 20, dtype=int), subrow, 1)
            subgrid.addWidget(mk_s('n'+c, 2, 1, 20), subrow, 2)

        subrow += 1
        cb = QCheckBox("Mesh decimation")
        cb.setToolTip("Enable decimation of STL mesh before discretization.")
        cb = self.widgets['decimation']  = cb
        cb.clicked.connect(self.set_decimation)
        subgrid.addWidget(cb, subrow, 0)

        subrow += 1
        subgrid.addWidget(mk_label("Decimation ratio"), subrow, 0)
        subgrid.addWidget(mk_le('r_ratio', self.r_ratio, 0.01, 0.99), subrow, 1)
        subgrid.addWidget(mk_s('r_ratio', self.r_ratio, 0.01, 0.99), subrow, 2)
        self.row_decimation = subrow

        subrow += 1
        # Clear as mud...
        tip = 'Increase the bounding box size to ensure the entire geometry is enclosed.\nThe size is increased by this ratio times the inter-sphere spacing.'
        subgrid.addWidget(mk_label("Wall extension ratio",tip), subrow, 0)
        subgrid.addWidget(mk_le('wall_extend_ratio', self.wall_extend_ratio, 0, 0.1, tip, exclude_min=False), subrow, 1)
        subgrid.addWidget(mk_s('wall_extend_ratio', self.wall_extend_ratio, 0, 0.1, tip), subrow, 2)

        subrow +=1
        tip = "Resolution for superquadric meshing.\nHigher resolution will yield better fit but will be slower"
        subgrid.addWidget(mk_label("Mesh resolution", tip), subrow, 0)
        subgrid.addWidget(mk_le('sq_resolution', self.sq_resolution, 20, 100, tip, exclude_min=False,), subrow, 1)
        subgrid.addWidget(mk_s('sq_resolution', self.sq_resolution, 20, 100, tip), subrow, 2)
        self.row_resolution = subrow

        subrow += 1
        cb = QCheckBox("Remove small spheres")
        self.widgets['remove_small_spheres'] = cb
        cb.setToolTip("Enable small sphere removal during glued-sphere generation.")
        cb.clicked.connect(self.set_remove_small_spheres)
        subgrid.addWidget(cb, subrow, 0)

        subrow += 1
        tip = "Volume ratio below which a sphere is considered small and is a candidate for removal."
        subgrid.addWidget(mk_label("Small sphere ratio", tip), subrow, 0)
        subgrid.addWidget(mk_le('ss_ratio', self.ss_ratio, 0, 0.1, tip),subrow, 1)
        subgrid.addWidget(mk_s('ss_ratio', self.ss_ratio, 0, 0.1, tip), subrow, 2)
        self.row_small_spheres = subrow

        subrow += 1
        tip = "Maximum fraction of total volume loss to allow during small sphere removal."
        subgrid.addWidget(mk_label("Maximum removal fraction", tip), subrow, 0)
        subgrid.addWidget(mk_le('ss_fraction', self.ss_fraction, 0, 0.1, tip), subrow, 1)
        subgrid.addWidget(mk_s('ss_fraction', self.ss_fraction, 0, 0.1, tip), subrow, 2)

        subrow += 1
        b = self.pushbutton_generate = QPushButton("Generate GSP")
        b.setAutoDefault(False)
        b.clicked.connect(self.generate_gsp)
        subgrid.addWidget(b, subrow, 0, 1, 3)

        row += 1
        grid.addItem(QSpacerItem(20,10,QSizePolicy.Minimum,QSizePolicy.Expanding),row,0)

        row += 1
        pb = self.progress_bar = QProgressBar()
        grid.addWidget(pb, row, 0, 1, 4)

        row +=1
        hbox = QHBoxLayout()
        hbox.addItem(QSpacerItem(20,10,
                                 QSizePolicy.Expanding,
                                 QSizePolicy.Minimum))

        b = QPushButton("Cancel")
        b.clicked.connect(self.cancel)
        b.setAutoDefault(False)
        hbox.addWidget(b)
        b = QPushButton("OK")
        b.clicked.connect(self.save)
        b.setAutoDefault(False)
        hbox.addWidget(b)

        grid.addLayout(hbox, row, 0, 1, 4)

    def invalidate(self):
        self.generated = False
        try:
            #os.remove(self.gsp_dir+'/gsp_config.dat')
            pass
        except Exception as e:
            if e.errno != errno.ENOENT:
                raise

    def cancel(self):
        if self.running and self.set_stop_flag:
            self.running = False
            self.set_stop_flag()
            self.cancelled = True
            self.progress_bar.setFormat("Cancelled")
            dir = 'GSP-%s'%self.P if self.P else 'GSP'
            try:
                shutil.rmtree(dir)
            except Exception as e:
                self.error(str(e))
            self.progress_bar.setValue(0)
        else:
            self.close()

    # Event handlers
    def handle_type(self, arg, force=False):
        if self.type == arg and not force:
            return
        self.invalidate()
        self.clear_gsp_view()
        self.progress_bar.reset()
        self.type = arg
        self.widgets['type'].setCurrentIndex(arg)
        gb = self.groupbox_sq
        for i in (self.row_resolution,):
            for j in (0,1,2):
                w = self.groupbox_disc.layout().itemAtPosition(i,j).widget()
                w.setVisible(arg==0)

        if self.type == SQ_TYPE:
            # Scale control is in groupbox_sq, should it be?
            #self.groupbox_sq.setEnabled(True)
            for i in range(gb.layout().rowCount()):
                for j in (0,1,2):
                    w = gb.layout().itemAtPosition(i,j).widget()
                    w.setEnabled(i != self.row_scale)

            self.button_choose_stl.hide()
            self.label_stl_file.hide()
            self.sq_actor.SetVisibility(True)
            self.stl_actor.SetVisibility(False)
            self.axes_actor.SetBounds(self.sq_actor.GetBounds())

        elif self.type == STL_TYPE:
            # Scale control is in groupbox_sq, should it be?
            #self.groupbox_sq.setEnabled(False)
            for i in range(gb.layout().rowCount()):
                for j in (0,1,2):
                    w = gb.layout().itemAtPosition(i,j).widget()
                    w.setEnabled(i == self.row_scale)

            self.button_choose_stl.show()
            self.label_stl_file.show()

            self.stl_actor.SetVisibility(True)
            self.stl_actor.SetScale(self.stl_scale)
            self.sq_actor.SetVisibility(False)
            self.axes_actor.SetBounds(self.stl_actor.GetBounds())
        self.reset_q()
        self.reset_view()

    def choose_stl_file(self, fname=None):
        if glob('*.stl') or glob('*.STL'):
            start_dir = '.'
        else:
            start_dir =  os.path.dirname(os.path.dirname(__file__))+'/gsp_sample_data'

        if fname is None:
            f = QFileDialog.getOpenFileName(self,
                                            caption='Choose STL file for particle shape',
                                            dir=start_dir,
                                            filter='STL files (*.STL *.stl)')

            fname = f[0]
            if not fname: # User hit Cancel
                return

            self.stl_scale = 1.0 # Reset scaling
            w = self.widgets['le_stl_scale']
            w.updateValue(None, self.stl_scale)
            self.stl_actor.SetScale(1.0)

            basename= os.path.basename(fname)

            if not shutil._samefile(os.path.dirname(fname), os.getcwd()):
                try:
                    shutil.copyfile(fname, basename)
                    fname = basename
                    self.msg("Copying %s to current directory." % basename)
                except Exception as e:
                    self.warn("Error copying %s to current directory: %s" % (basename, e))
            else:
                fname = basename # It's in the current directory

        elif not os.path.exists(fname):
            self.warn("File %s not found" % fname)
            fname = None

        self.invalidate()
        self.clear_gsp_view()
        self.progress_bar.reset()

        self.stl_file = fname
        if fname:
            self.label_stl_file.setText(' '+os.path.basename(fname))
            self.label_stl_file.setStyleSheet("");
            self.stl_actor.SetScale(1.0)
            le = self.widgets['le_stl_scale']
            le.updateValue(None,1.0) # Set lineedit
            self.handle_le(le) # Propagate to sliders
            self.reset_q()
            self.stl_source.SetFileName(fname)
            self.stl_source.Update()
            self.stl_mapper.SetInputConnection(self.stl_source.GetOutputPort())
            # this SetPosition is needed if one STL object is rotated
            # and another STL is loaded
            self.stl_actor.SetPosition(0,0,0)
            self.axes_actor.SetBounds(self.stl_actor.GetBounds())
            self.stl_actor.SetVisibility(True)
            # Should we move center of mass to origin?
            #vcom = vtk.vtkCenterOfMass()
            #vcom.SetInputConnection(self.stl_source.GetOutputPort())
            #vcom.Update()
            #com = vcom.GetCenter()
        else:
            self.label_stl_file.setText(' No file selected')
            self.label_stl_file.setStyleSheet("font-style: italic;")
            self.stl_actor.SetVisibility(False)

        self.reset_view()
        if fname:
            self.axes_actor.SetBounds(self.stl_actor.GetBounds())
        self.rw.Render()

    def handle_button(self, arg):
        iren = self.iren
        if arg == 1:
            self.rotate_object = True
            self.rotate_camera = False
            self.style.SetCurrentStyleToTrackballActor()
            #iren.RemoveObservers('MiddleButtonPressEvent')
            #iren.RemoveObservers('RightButtonPressEvent')
            iren.RemoveObservers('CharEvent')
        else:
            self.rotate_object = False
            self.rotate_camera = True
            self.style.SetCurrentStyleToTrackballCamera()
            #iren.RemoveObservers('MiddleButtonPressEvent')
            iren.RemoveObservers('CharEvent')
            # allow right button zoom

    def set_mn(self, m, n):
        self.m = m
        self.n = n
        self.widgets['le_m'].updateValue(None, m)
        self.widgets['le_n'].updateValue(None, n)
        self.suppress_update = True
        self.widgets['s_m'].setValue(int(16*(m-2)))
        self.widgets['s_n'].setValue(int(16*(n-2)))
        self.suppress_update = False
        self.sq_source.SetThetaRoundness(2/self.m)
        self.sq_source.SetPhiRoundness(2/self.n)
        self.sq_source.Update()

    def set_base_opacity(self, val):
        val = (val+1)/100
        val = min(val, 1.0)
        val = max(val, 0.0)
        self.base_opacity = val
        if self.suppress_update:
            return
        self.sq_actor.GetProperty().SetOpacity(val)
        self.stl_actor.GetProperty().SetOpacity(val)
        self.rw.Render()

    def set_gsp_opacity(self, val):
        val = (val+1)/100
        val = min(val, 1.0)
        val = max(val, 0.0)
        if self.suppress_update:
            return
        for sa in self.sphere_actors:
            sa.GetProperty().SetOpacity(val)
            sa.GetProperty().SetOpacity(val)
        self.gsp_opacity = val
        self.rw.Render()

    def set_decimation(self, arg):
        arg = bool(arg)
        self.decimation = arg
        self.widgets['decimation'].setChecked(arg)
        for i in (self.row_decimation,):
            for j in (0,1,2):
                w = self.groupbox_disc.layout().itemAtPosition(i,j).widget()
                w.setEnabled(arg)
        self.invalidate()
        self.clear_gsp_view()

    def set_remove_small_spheres(self,arg):
        arg = bool(arg)

        self.remove_small_spheres = arg
        self.widgets['remove_small_spheres'].setChecked(arg)
        for i in (self.row_small_spheres, self.row_small_spheres+1):
            for j in (0,1,2):
                w = self.groupbox_disc.layout().itemAtPosition(i,j).widget()
                w.setEnabled(arg)
        self.invalidate()
        self.clear_gsp_view()

    def msg(self, text):
        if self.parent:
            self.parent.print_internal(text, color='blue')
        else:
            print(text)

    def warn(self, text):
        if self.parent:
            self.parent.warn(text, popup=True)
        else:
            print("Warning: " + text)

    def error(self, text):
        if self.parent:
            self.parent.error(text, popup=True)
        else:
            print("Error: " + text)

    def setup(self, callback, P):
        self.callback = callback
        self.P = P
        if P is not None:
            gsp_dir = 'GSP-%s' % P
        else:
            gsp_dir = 'GSP'
        if not os.path.exists(gsp_dir):
            self.msg("Creating "+gsp_dir)
        os.makedirs(gsp_dir, exist_ok=True)
        self.gsp_dir = gsp_dir
        fname = gsp_dir+'/gsp_properties.dat'
        if os.path.exists(fname):
            self.msg("Reading "+fname)
            with open(fname) as f:
                for line in f:
                    tok = line.split('=')
                    if len(tok) == 2:
                        k,v = tok
                        k = k.strip()
                        v = v.strip()
                        try:
                            if k in ('q', 'com'):
                                v = np.array(list(map(float, v.split())))
                            else:
                                v = (None if v=="None" else
                                     bool(v=="True") if v in ("True","False")
                                     else v if k=='stl_file'
                                     else float(v) if k=='n'
                                     else int(v) if k.startswith(('n','type')) else float(v))

                            setattr(self, k, v)
                        except Exception as e:
                            print(e)

        self.handle_type(self.type, force=True)
        if self.type == STL_TYPE:
            # choose_stl_file resets stl_scale to 1.0
            stl_scale = self.stl_scale
            self.choose_stl_file(fname=self.stl_file)
            self.stl_scale = stl_scale

        self.set_mn(self.m, self.n)
        self.set_decimation(self.decimation)
        self.set_remove_small_spheres(self.remove_small_spheres)

        for k in ('a','b','c', 'stl_scale', 'nx','ny','nz',
                  'r_ratio', 'wall_extend_ratio', 'ss_ratio', 'ss_fraction'):
            v = getattr(self, k)
            le = self.widgets['le_'+k]
            le.updateValue(None, v) # Set lineedit value
            self.handle_le(le)  # Propagate event to sliders and model

        q = self.q
        for i in range(4):
            self.widgets['le_q%d'%(i+1)].updateValue(None, q[i])
        self.apply_q()
        self.update_gsp_view(scale=self.bounding_diam, vol=self.vol, com=self.com)
        self.reset_q()
        self.reset_view()
        self.invalidate()
        if self.type == STL_TYPE: # Repaint axes, not sure why this is needed
            self.axes_actor.SetBounds(self.stl_actor.GetBounds())
        elif self.type == SQ_TYPE:
            self.axes_actor.SetBounds(self.sq_actor.GetBounds())

    def save(self):
        self.thumbnail()
        if self.generated:
            self.save_gsp_properties()
        if self.callback:
            self.callback(self.P, self.bounding_diam)
            # We get other params from the gsp_properties file
        self.close()


    def reset(self):
        self.progress_bar.reset()
        self.invalidate()
        self.clear_gsp_view()
        self.pushbutton_generate.setDown(False)


    def handle_le(self, w):
        self.progress_bar.reset()
        self.invalidate()
        self.clear_gsp_view()
        k,v = w.key, w.value
        #print(k,v)
        setattr(self, k, v)
        self.set_slider_pos(k, v)
        if k in ('m', 'n'):
            self.sq_source.SetThetaRoundness(2/self.m)
            self.sq_source.SetPhiRoundness(2/self.n)
            self.sq_source.Update()

        if k in ('a','b','c'):
            self.sq_source.SetScale(self.a, self.b, self.c)
            self.sq_source.Update()

        if k == 'stl_scale':
            self.stl_actor.SetScale(self.stl_scale)
            self.axes_actor.SetBounds(self.stl_actor.GetBounds())
            self.reset_view()

        if k in ('a','b','c','m','n'):
            self.axes_actor.SetBounds(self.sq_actor.GetBounds())

        self.rw.Render()

        #self.reset_view()

    def set_slider_pos(self, key, val):
        self.suppress_update = True
        s = self.widgets['s_'+key]
        mx = int(s.maximum())

        if key in ('m', 'n'):
            x = mx if val==24 else 0 if val==2 else (val-1)*16
            s.setValue(int(x))
        elif key in ('nx', 'ny', 'nz'):
            x = (val-1)*(1+mx/s.max)
            s.setValue(int(x))
        elif key in ('a','b','c','stl_scale'):
            man, exp = fman(val), fexp(val)
            self.exp[key] = exp
            s.setValue(int(40*man))
        elif 'opacity' in key:
            s.setValue(int(mx * val))
        else:
            v = int((1 + mx)*(val-s.min)/s.max)
            v = min(v,mx)
            s.setValue(v)
        self.suppress_update = False

    def handle_s(self, key, val):
        if self.suppress_update:
            return
        self.progress_bar.reset()
        self.invalidate()
        self.clear_gsp_view()
        s = self.widgets['s_'+key]
        mx = s.maximum()
        if key in ('m', 'n'):
            x = min(2 + val/16, 24)
            if key == 'm':
                self.m = x
                self.sq_source.SetThetaRoundness(2/x)
            else:
                self.n = x
                self.sq_source.SetPhiRoundness(2/x)
            self.widgets['le_'+key].updateValue(None,x)
            self.sq_source.Update()
            self.rw.Render()

        elif key in ('a', 'b', 'c', 'stl_scale'):
            if val == mx:
                # Hit right end of slider.  Fake a mouse release event
                # so we don't move right forever
                e = QMouseEvent(3, QPoint(0,0), Qt.MouseButton(Qt.LeftButton),
                                Qt.MouseButtons(Qt.NoButton),
                                Qt.KeyboardModifiers(Qt.NoModifier))
                QApplication.instance().sendEvent(s, e)
                val = int(mx/10)
                s.setSliderPosition(val)

                self.exp[key] += 1

            elif val == 0:
                # Hit left end of slider.  Fake a mouse release event
                # so we don't move left forever
                e = QMouseEvent(3, QPoint(0,0), Qt.MouseButton(Qt.LeftButton),
                                Qt.MouseButtons(Qt.NoButton),
                                Qt.KeyboardModifiers(Qt.NoModifier))
                QApplication.instance().sendEvent(s, e)
                s.setSliderPosition(40)
                val = 40
                self.exp[key] -= 1
            x = (val/40) * 10**self.exp[key]
            setattr(self, key, x)
            self.widgets['le_'+key].updateValue(None, x)
            if key in ('a','b','c'):
                self.sq_source.SetScale(self.a, self.b, self.c)
                self.sq_source.Update()
            if key == 'stl_scale':
                self.stl_actor.SetScale(self.stl_scale)
                self.axes_actor.SetBounds(self.stl_actor.GetBounds())
                self.reset_view()
            else:
                self.sq_source.Update()
            self.rw.Render()
        elif key in ('nx', 'ny', 'nz'):
            x = int (1+ val/21)
            setattr(self, key, x)
            self.widgets['le_'+key].updateValue(None, x)
        elif key.endswith(('ratio', 'fraction', 'resolution')):
            x = s.min + (s.max-s.min)*val/400
            if key.endswith('resolution'): x = int(x)
            setattr(self, key, x)
            self.widgets['le_'+key].updateValue(None, x)
        if key in ('a','b','c','m','n'):
            self.axes_actor.SetBounds(self.sq_actor.GetBounds())


    def toggle_axes(self, arg):
        self.widgets['tb_axes'].setToolTip("Hide axes" if arg else "Show axes")
        self.axes_actor.SetVisibility(arg)
        self.rw.Render()


    def handle_orientation(self, actor, event):
        # Interactive rotation handler
        if not self.rotate_object:
            return
        actor.ComputeMatrix()
        m = actor.GetMatrix()
        M = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                M[i,j] = m.GetElement(i,j)
        d = np.linalg.det(M)
        M /= d**(1/3)
        self.q = q = m_to_q(M)
        for i in range(4):
            self.widgets['le_q%d'%(i+1)].updateValue(None, q[i])
        actor = self.stl_actor if self.type == STL_TYPE else self.sq_actor
        self.progress_bar.reset()
        self.invalidate()
        self.clear_gsp_view()
        self.axes_actor.SetBounds(actor.GetBounds())
        return True

    def handle_q(self, w, val_dict, args):
        key = w.key
        val = val_dict[key]
        self.q[w.arg] = val
        if np.sqrt(np.dot(self.q, self.q)) < 0.001:
            self.q = np.array((1,0,0,0), dtype=float)
            for i in range(4):
                self.widgets['le_q%d'%i].updateValue(None, self.q[i],8)

    def perspective(self, parallel=None):
        """change the perspective of the vtk scene"""
        cam = self.ren.GetActiveCamera()

        if parallel is None:
            parallel = not cam.GetParallelProjection()

        tb = self.widgets['tb_perspective']
        if parallel:
            cam.ParallelProjectionOn()
            tb.setIcon(get_icon('parallel.svg'))
        else:
            cam.ParallelProjectionOff()
            tb.setIcon(get_icon('perspective.svg'))

        self.rw.Render()

    def set_view(self, view='xy'):
        #self.perspective(parallel=True)
        cam = self.ren.GetActiveCamera()
        if view == 'xy':
            cam.SetPosition(0, 0, 10000000)
            cam.SetViewUp(0, 1, 0)
            #if self.view_flip[1]:
            #    cam.Azimuth(180)
            #self.view_flip[1] = not self.view_flip[1]
        elif view == 'xz':
            cam.SetPosition(0, 10000000, 0)
            cam.SetViewUp(1, 0, 0)
            #if self.view_flip[2]:
            #    cam.Azimuth(180)
            #self.view_flip[2] = not self.view_flip[2]
        elif view == 'yz':
            cam.SetPosition(10000000, 0, 0)
            cam.SetViewUp(0, 1, 0)
            #if self.view_flip[2]:
            #    cam.Azimuth(180)
            #self.view_flip[2] = not self.view_flip[2]

        self.reset_view()

    def rotate(self, degrees):
        cam = self.ren.GetActiveCamera()
        cam.Roll(degrees)
        self.rw.Render()

    def reset_view(self):
        self.ren.ResetCamera()
        self.ren.GetActiveCamera().Zoom(0.8)
        self.rw.Render()

    def reset_q(self):
        vt = vtk.vtkTransform()
        vm = vtkMatrix4x4()
        vm.Identity()
        vt.SetMatrix(vm)
        self.sq_actor.SetOrientation(vt.GetOrientation())
        self.stl_actor.SetOrientation(vt.GetOrientation())
        self.stl_actor.SetPosition(0,0,0)
        self.q = np.array([1,0,0,0], dtype=float)
        self.widgets['le_q1'].setText('1')
        self.widgets['le_q2'].setText('0')
        self.widgets['le_q3'].setText('0')
        self.widgets['le_q4'].setText('0')
        self.rw.Render()
        self.reset_view()

    def apply_q(self):
        m = q_to_m(self.q)
        vm = vtkMatrix4x4()
        vm.Identity()
        for i in range(3):
            for j in range(3):
                vm.SetElement(i,j,m[i,j])
        vt = vtk.vtkTransform()
        vt.SetMatrix(vm)
        self.sq_actor.SetOrientation(vt.GetOrientation())
        #self.sq_source.Update()
        #self.axes_actor.SetBounds(self.sphere_actor.GetBounds())
        self.rw.Render()

    def generate_gsp(self):
        if self.running:
            return
        self.running = True
        self.cancelled = False
        self.progress_bar.reset()
        self.progress_bar.setFormat("Generating... please wait")
        self.progress_bar.setValue(0)
        self.pushbutton_generate.setDown(True)
        self.pushbutton_generate.setEnabled(False)
        self.clear_gsp_view()
        ret = None
        for ext in ('dat', 'csv'):
            try:
                os.remove(self.gsp_dir+'/gsp_config.'+ext)
            except Exception as e:
                self.pushbutton_generate.setDown(False)
                self.pushbutton_generate.setEnabled(True)
                if e.errno != errno.ENOENT:
                    raise
        try:
            from mfixgui.tools.gsp_generator import generate_gsp, set_progress_callback, set_stop_flag
            self.set_stop_flag = set_stop_flag
            set_stop_flag(False)
            set_progress_callback(lambda p: (self.progress_bar.setValue(int(p)),
                                             self.progress_bar.setFormat("Generating... %p%"),
                                             QApplication.instance().processEvents()))

            args = dict(
                file_type = 'sq' if self.type==SQ_TYPE else 'stl',
                usr_input_file_names = [self.stl_file],
                units = [self.stl_scale],
                sq_inputs = [[self.a, self.b, self.c, self.m, self.n]],
                qs = [self.q],
                centers = [self.stl_actor.GetCenter() if self.type == STL_TYPE else
                           self.sq_actor.GetCenter()],
                discretizations = [[self.nx, self.ny, self.nz]],
                wall_extend_ratios = [self.wall_extend_ratio] ,
                mesh_decimations=[self.decimation],
                r_ratios=[self.r_ratio],
                isRemoveTinyVolume=self.remove_small_spheres,
                critical_vol_fraction=self.ss_ratio,
                critical_totalvol_fraction=1-self.ss_fraction,
                isUsrVels = [False],
                isTp = False,
                isRandomQ = True,
                startPhase = self.P,
                sq_resolution = int(self.sq_resolution),
                save_dir = self.gsp_dir)
            ret = generate_gsp(**args)
            self.generated = True
        except Exception as e:
            self.error(''.join(traceback.format_exception(e)[1:]))
        finally:
            self.running = False
            self.pushbutton_generate.setDown(False)
            self.pushbutton_generate.setEnabled(True)
        if ret:
            nspheres, equiv_diam, bounding_diam, area, vol, com, surface_deviation, vol_deviation = ret
            self.progress_bar.setFormat("Done")
        else:
            self.progress_bar.setFormat("Cancelled" if self.cancelled else "Failed")
            nspheres = equiv_diam = bounding_diam = area = vol = com = None
            surface_deviation = vol_deviation = None
        # surface/volume_deviation is a percentage, check if they are larger than 1%
        if (surface_deviation and surface_deviation >= 1
            or vol_deviation and vol_deviation >= 1):
            self.warn("Surface/volume deviation is larger than 1%. Adjust discretization or enable mesh decimation.")
        self.nspheres = nspheres
        self.equiv_diam = equiv_diam
        self.bounding_diam = bounding_diam
        self.area = area
        self.vol = vol
        self.com = com
        self.surface_deviation = surface_deviation
        self.vol_deviation = vol_deviation
        self.running = self.cancelled = False
        self.update_gsp_view(scale=bounding_diam, vol=vol, com=com)

    def clear_gsp_view(self):
        for sa in self.sphere_actors:
            self.ren.RemoveActor(sa)
            del(sa)
        self.sphere_actors = []
        if self.base_opacity <= 0.01:
            self.set_slider_pos('base_opacity', 50)
            self.set_base_opacity(50)
        self.rw.Render()

    def update_gsp_view(self, scale=1, vol=None, com=None):
        fname = self.gsp_dir + '/' + 'gsp_config.csv'
        if scale is None:
            scale = 1
        try:
            data = []
            with open(fname, 'r') as f:
                header_seen = False
                while True:
                    line = f.readline()
                    if not line:
                        break
                    if not line.strip(): # skip blank lines
                        continue
                    if line.startswith('#'):# skip comments
                        continue
                    if not header_seen:
                        header_seen = True
                        continue
                    x,y,z,d,rest = line.split(',', maxsplit=4)
                    data.append([float(x),float(y),float(z),float(d)])
        except IOError as e:
            if e.errno == errno.ENOENT:
                self.clear_gsp_view()
                return
            raise

        data = np.array(data)

        self.clear_gsp_view()
        sphere_actors = []
        for row in data:
            x,y,z,d = row*scale
            #print(x,y,z,d)
            m = vtkPolyDataMapper()
            s = vtkSphereSource()
            s.SetThetaResolution(100)
            s.SetPhiResolution(100)
            if com is not None:
                x += com[0]
                y += com[1]
                z += com[2]
            s.SetCenter(x, y, z)
            s.SetRadius(d/2)
            m.SetInputConnection(s.GetOutputPort())
            sa = vtkActor()
            sphere_actors.append(sa)
            sa.SetVisibility(True)
            sa.SetPickable(False)
            sa.SetMapper(m)
            sa.GetProperty().SetOpacity(self.gsp_opacity)
            sa.GetProperty().SetDiffuseColor(self.sphere_color)
            sa.GetProperty().SetDiffuse(.7)
            sa.GetProperty().SetSpecular(.7)
            sa.GetProperty().SetSpecularPower(50.0)
            #self.rw.Render() # Pop up spheres one at a time
            self.ren.AddActor(sa)

        self.rw.Render() # Show all spheres at once
        self.sphere_actors = sphere_actors

    def save_gsp_properties(self):
        fname = self.gsp_dir+'/gsp_properties.dat'
        self.msg("Saving "+fname)
        try:
            with open(fname, 'w') as f:
                f.write("version = 1\n")
                for k in ('type', 'a', 'b', 'c', 'm', 'n',
                          'stl_file', 'stl_scale',
                          'nx', 'ny', 'nz', 'wall_extend_ratio', 'r_ratio',
                          'ss_ratio', 'ss_fraction', 'decimation',
                          'remove_small_spheres',
                          'nspheres', 'equiv_diam', 'area', 'vol',
                          'bounding_diam'):
                    v = getattr(self,k)
                    f.write('%s = %s\n' % (k, getattr(self,k)))
                if self.com is None:
                    f.write('com = None\n')
                else:
                    f.write('com = %s %s %s\n' % tuple(self.com))
                f.write('q = %s %s %s %s\n' % tuple(self.q))
                f.write('surface_deviation = %s\n' % self.surface_deviation)
                f.write('vol_deviation = %s\n' % self.vol_deviation)
        except Exception as e:
            self.warn(str(e))
            raise


    def closeEvent(self, ev):
        if self.parent:
            self.parent.ui.solids.groupbox_gsp_shape.setDisabled(False)
        #self.close()
        self.hide()
        ev.ignore()

    def thumbnail(self):
        #make a pixmap of the particle
        # I do not know why setting these directly results in a blank thumbnail - cgw
        #self.axes_actor.SetVisibility(False)
        #self.stl_actor.SetVisibility(False)
        #self.sq_actor.SetVisibility(False)
        if not self.sphere_actors: # Don't create blank thumbnail
            return

        # Save state
        base_opacity = self.base_opacity
        gsp_opacity = self.gsp_opacity
        axes = self.widgets['tb_axes'].isChecked()
        # Set up for the picture
        self.set_base_opacity(0)
        self.set_gsp_opacity(100)
        self.ren.GetActiveCamera().Zoom(1.0)
        self.toggle_axes(False)
        self.ren.GradientBackgroundOff()
        c1 = QColor('#ffffffff')
        self.ren.SetBackground(c1.getRgbF()[:3])
        try:
            window_image = vtk.vtkWindowToImageFilter()
            window_image.SetInput(self.rw)
            window_image.ReadFrontBufferOff() # cargo
            if hasattr(window_image, 'SetScale'):
                window_image.SetScale(1)
            else:
                window_image.SetMagnification(1)
            writer = vtk.vtkPNGWriter()
            fname = self.gsp_dir + "/thumbnail.png"
            writer.SetFileName(fname)
            writer.SetInputConnection(window_image.GetOutputPort())
            window_image.Update()
            writer.Write()
        except Exception as e:
            self.error(str(e))

        finally: # Put everything back
            c1 = QColor(SETTINGS.value('vtk_bck_color1', '#aaffff'))
            c2 = QColor(SETTINGS.value('vtk_bck_color2', '#ffffff'))
            self.ren.GradientBackgroundOn()
            self.ren.SetBackground(c1.getRgbF()[:3])
            self.ren.SetBackground2(c2.getRgbF()[:3])
            self.set_gsp_opacity(100*gsp_opacity)
            self.set_base_opacity(100*base_opacity)
            self.toggle_axes(axes)

def q_to_m(q):
    h = np.sqrt(np.dot(q,q))
    q = np.array(q) / h # normalize to unit quaternion
    M = np.zeros((3,3))
    M[0,0] = 2*(q[0]**2 + q[1]**2) - 1
    M[1,1] = 2*(q[0]**2 + q[2]**2) - 1
    M[2,2] = 2*(q[0]**2 + q[3]**2) - 1
    M[1,0] = 2*(q[1]*q[2] + q[0]*q[3])
    M[0,1] = 2*(q[1]*q[2] - q[0]*q[3])
    M[0,2] = 2*(q[1]*q[3] + q[0]*q[2])
    M[2,0] = 2*(q[1]*q[3] - q[0]*q[2])
    M[2,1] = 2*(q[2]*q[3] + q[0]*q[1])
    M[1,2] = 2*(q[2]*q[3] - q[0]*q[1])
    return M

def m_to_q(M):
    w = np.sqrt(max(0, 1 + M[0,0] + M[1,1] + M[2,2])) / 2
    x = np.sqrt(max(0, 1 + M[0,0] - M[1,1] - M[2,2])) / 2;
    y = np.sqrt(max(0, 1 - M[0,0] + M[1,1] - M[2,2])) / 2;
    z = np.sqrt(max(0, 1 - M[0,0] - M[1,1] + M[2,2])) / 2;
    x = np.copysign(x, M[2,1] - M[1,2])
    y = np.copysign(y, M[0,2] - M[2,0])
    z = np.copysign(z, M[1,0] - M[0,1])
    q = np.array((w,x,y,z))
    h = np.sqrt(np.dot(q,q))
    q /= h
    # Less stable,  see https://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    #r = np.math.sqrt(1+M[0,0]+M[1,1]+M[2,2])*0.5
    #i = (M[2,1]-M[1,2])/(4*r)
    #j = (M[0,2]-M[2,0])/(4*r)
    #k = (M[1,0]-M[0,1])/(4*r)
    #q = np.array((r,i,j,k))

    for i in range(4):
        if abs(q[i]) < 1e-6:
            q[i] = 0
    return q

# Exponent and mantissa
def fexp(f):
    return int(floor(log10(abs(f)))) if f != 0 else 0

def fman(f):
    return f/10**fexp(f)

def main():
    args = sys.argv
    qapp = QApplication(args)
    qapp.setStyle("fusion")
    p = GSPopup()
    p.setup(None,None)
    #q = [0.34934191, 0.50357956, 0.5967449, 0.51794149]
    q = [1,0,0,0]
    p.show()
    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    qapp.exec_()
    qapp.deleteLater()
    sys.exit()

if __name__ == '__main__':
    main()
