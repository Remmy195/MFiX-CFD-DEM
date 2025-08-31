#!/usr/bin/env python

#  Integrated version of sq_designer

# TODO put "parameters" and "properties" in group box

import sys
import signal
# Note, we should honor the global 'novtk' flag.  But
#   this module is only loaded on-demand when the user clicks
#   'sq_designer', so if they try to run this without VTK it
#   will trigger an error

# 'import vtk' should not be needed, but if we don't import it,
#  SetCurrentStyleToTrackballCamera fails - cgw 2024-04-04
import vtk

import pyqtgraph as pg
import numpy as np
from math import floor, log10, hypot

from qtpy.QtWidgets import (QApplication, QCheckBox, QDialog, QFrame,
                            QGridLayout, QGroupBox, QHBoxLayout, QLabel,
                            QPushButton, QRadioButton, QSizePolicy,
                            QSlider, QSpacerItem, QSplitter, QVBoxLayout, QWidget)

from qtpy.QtGui import QColor, QBrush, QMouseEvent

from mfixgui.tools import sqp_d0, sqp_volume
from mfixgui.tools.qt import get_icon, sub_icon_size, SETTINGS, find_gui, make_toolbutton

from mfixgui.widgets.base import LineEdit
from qtpy.QtCore import  Qt, QPoint

from vtkmodules.vtkRenderingCore import vtkActor, vtkPolyDataMapper, vtkRenderer
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkCommonCore import vtkCommand
from vtkmodules.vtkCommonMath import vtkMatrix4x4
from vtkmodules.vtkFiltersSources import vtkSuperquadricSource, vtkSphereSource
from vtkmodules.vtkRenderingAnnotation import vtkCubeAxesActor

blue = QBrush(QColor('blue'))
green =QBrush(QColor('green'))
small = 10
big = 20

class SQPopup(QDialog):
    def __init__(self, parent=None):
        super(SQPopup, self).__init__(parent)
        self.parent = parent
        self.setWindowTitle("Superquadric designer")
        self.layout = layout = QHBoxLayout()
        layout.setContentsMargins(5,5,5,5)
        self.setLayout(layout)
        splitter = self.splitter = QSplitter()
        layout.addWidget(splitter)
        self.widgets = {}
        self.P = None # Active solids phase
        self.callback = None
        self.rotate_object = False
        self.rotate_camera = True
        self.view_flip = [False]*3
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

        hbox = QHBoxLayout()
        b = QRadioButton('Rotate object')
        b.setChecked(False)
        b.clicked.connect(lambda val: self.handle_button(1))
        hbox.addWidget(b)
        b = QRadioButton('Rotate axes')
        b.setChecked(True)
        b.clicked.connect(lambda val: self.handle_button(2))
        hbox.addWidget(b)
        b = QCheckBox('Show bounding sphere')
        b.clicked.connect(self.show_bounding_sphere)
        b.setChecked(True)
        hbox.addWidget(b)
        left.addLayout(hbox)

        ren = self.ren = vtkRenderer()
        rw = self.rw = self.vtk.GetRenderWindow()
        rw.AddRenderer(ren)
        iren = self.iren = rw.GetInteractor()
        style = self.style = iren.GetInteractorStyle()
        iren._MouseMoveEvent = iren.MouseMoveEvent
        iren.MouseMoveEvent = lambda: iren._MouseMoveEvent() if not (self.rotate_object and iren.GetShiftKey()) else None

        style.SetCurrentStyleToTrackballCamera()

        # Disable keyboard commands
        iren.RemoveObservers('CharEvent')

        # Colors
        colors = vtkNamedColors()
        #bg_color = colors.GetColor3d("DarkSlateGray")
        object_color = colors.GetColor3d("DarkGreen")
        sphere_color = colors.GetColor3d("LightSkyBlue")
        axis_colors = [colors.GetColor3d("Black"),
                       colors.GetColor3d("Black"),
                       colors.GetColor3d("Black")]

        # Create a superquadric
        sq_source = self.sq_source = vtkSuperquadricSource()
        sq_source.SetAxisOfSymmetry(2)
        sq_source.SetSize(1)
        sq_source.SetThetaResolution(500)
        sq_source.SetPhiResolution(500)

        self.m = self.n = 2
        self.a = self.b = self.c  = 0.01 # 1 cm

        self.d = sqp_d0(self.a, self.b, self.c, self.m, self.n)
        self.exp = {'a':fexp(self.a),
                    'b':fexp(self.b),
                    'c':fexp(self.c),
                    'd':fexp(self.d)}

        self.q = np.array((1,0,0,0), dtype=float)
        sq_source.SetScale(self.a, self.b, self.c)
        sq_source.SetThetaRoundness(2/self.m)
        sq_source.SetPhiRoundness(2/self.n)
        sq_source.SetCenter(0,0,0)
        sq_source.Update()  # needed to GetBounds later

        sphere_source = self.sphere_source = vtkSphereSource()
        sphere_source.SetThetaResolution(500)
        sphere_source.SetPhiResolution(500)
        sphere_source.SetRadius(self.d/2)
        sphere_source.SetCenter(0,0,0)
        sphere_source.Update()  # needed to GetBounds later

        mapper = self.mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(sq_source.GetOutputPort())
        self.sq_actor = sq_actor = vtkActor()
        sq_actor.SetMapper(mapper)
        sq_actor.GetProperty().SetDiffuseColor(object_color)
        sq_actor.GetProperty().SetDiffuse(.7)
        sq_actor.GetProperty().SetSpecular(.7)
        sq_actor.GetProperty().SetSpecularPower(50.0)

        mapper2 = self.mapper2 = vtkPolyDataMapper()
        mapper2.SetInputConnection(sphere_source.GetOutputPort())
        self.sphere_actor = sphere_actor = vtkActor()
        sphere_actor.SetPickable(False)
        sphere_actor.SetMapper(mapper2)
        sphere_actor.GetProperty().SetOpacity(0.3)
        sphere_actor.GetProperty().SetDiffuseColor(sphere_color)
        sphere_actor.GetProperty().SetDiffuse(.7)
        sphere_actor.GetProperty().SetSpecular(.7)
        sphere_actor.GetProperty().SetSpecularPower(50.0)

        # Setup axes
        self.axes_actor = axes_actor = vtkCubeAxesActor()
        # axes_actor.SetUseTextActor3D(1) less legible
        axes_actor.SetCamera(self.ren.GetActiveCamera())
        for i,col in enumerate(axis_colors):
            axes_actor.GetTitleTextProperty(i).SetFontSize(48)
            axes_actor.GetTitleTextProperty(i).SetColor(col)
            axes_actor.GetLabelTextProperty(i).SetColor(col)
        axes_actor.SetBounds(sphere_actor.GetBounds())
        axes_actor.DrawXGridlinesOn()
        axes_actor.DrawYGridlinesOn()
        axes_actor.DrawZGridlinesOn()
        axes_actor.SetGridLineLocation(axes_actor.VTK_GRID_LINES_FURTHEST)
        axes_actor.XAxisMinorTickVisibilityOff()
        axes_actor.YAxisMinorTickVisibilityOff()
        axes_actor.ZAxisMinorTickVisibilityOff()

        #https://vtk.org/doc/nightly/html/classvtkCubeAxesActor.html#details
        #axes_actor.SetFlyModeToStaticEdges()
        #axes_actor.SetFlyModeToOuterEdges()
        axes_actor.SetFlyModeToClosestTriad()
        #axes_actor.SetFlyModeToFurthestTriad()

        # Put the VTK scene together
        ren.AddActor(axes_actor)
        ren.AddActor(sq_actor)
        ren.AddActor(sphere_actor)
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
        self.suppress_diameter = True
        self.reset_view()

        # Right part:  GUI controls
        row = 0
        grid = QGridLayout()
        grid.setVerticalSpacing(2)
        grid.setContentsMargins(5,5,5,5)
        w = QWidget()
        w.setLayout(grid)
        self.splitter.addWidget(w)

        def mk_le(key, val, min=None, max=None):
            le = LineEdit()
            if key != 'd':
                le.required = True
            le.key = key
            le.dtype = float
            le.saved_value = le.default_value = val
            le.setText(str(val))
            if min is not None:
                le.min = min
            if min == 0:
                le.exclude_min = True
            if max is not None:
                le.max = max
            le.setMaximumWidth(le.fontMetrics().boundingRect('9.999999').width())
            #le.setMaximumWidth(50)
            le.value_updated.connect(self.handle_le)
            self.widgets['le_'+key] = le
            return le

        def mk_s(key, val, min=None, max=None):
            s = QSlider(1) #horiz
            width = 400
            s.setMaximum(width)
            s.width = width
            if min is not None and max is not None:
                s.setValue(int(width * (val-min)/(max-min)))
            elif key in ('a','b','c','d'):
                s.setValue(int(40*fman(val)))
                self.exp[key] = fexp(val)

            s.valueChanged.connect(lambda val, key=key: self.handle_s(key, val))

            if key in ('m','n'):
                s.sliderReleased.connect(self.finish_s)
            self.widgets['s_'+key] = s
            return s

        row += 1
        grid.addWidget(QLabel("X semiaxis"), row, 0)
        grid.addWidget(mk_le('a', self.a, 0), row, 1)
        grid.addWidget(mk_s('a', self.a), row, 2)
        grid.addWidget(QLabel("m"), row, 3)

        row += 1
        grid.addWidget(QLabel("Y semiaxis"), row, 0)
        grid.addWidget(mk_le('b', self.b, 0), row, 1)
        grid.addWidget(mk_s('b', self.b), row, 2)
        grid.addWidget(QLabel("m"), row, 3)

        row += 1
        grid.addWidget(QLabel("Z semiaxis"), row, 0)
        grid.addWidget(mk_le('c', self.c, 0), row, 1)
        grid.addWidget(mk_s('c', self.c), row, 2)
        grid.addWidget(QLabel("m"), row, 3)

        row += 1
        grid.addWidget(QLabel("Bounding diameter"), row, 0)
        grid.addWidget(mk_le('d', self.d, 0), row, 1)
        grid.addWidget(mk_s('d', self.d), row, 2)
        grid.addWidget(QLabel("m"), row, 3)

        row += 1
        grid.addWidget(QLabel("Shape exponent m"), row, 0)
        grid.addWidget(mk_le('m', 2, 2, 24), row, 1)
        grid.addWidget(mk_s('m', 2, 2, 24), row, 2, 1, 2)

        row += 1
        grid.addWidget(QLabel("Shape exponent n"), row, 0)
        grid.addWidget(mk_le('n', 2, 2, 24), row, 1)
        grid.addWidget(mk_s('n', 2, 2, 24), row, 2, 1, 2)

        row += 1
        cb = QCheckBox("Constrain exponents to favorable values")
        cb.setToolTip("Simulations will run faster when certain exponent pairs are used.")
        self.constrain = True
        cb.setChecked(self.constrain)
        cb.toggled.connect(self.handle_constrain)
        grid.addWidget(cb, row, 0, 1, 3)

        row += 1
        b = QPushButton("Show exponent map")
        b.setAutoDefault(False)
        self.plot = self.mk_plot()
        self.plot_shown = False
        self.plot_button = b
        grid.addWidget(b, row, 0, 1, 4)
        b.clicked.connect(self.show_hide_plot)

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
            le.setMaximumWidth(100)
            le.value_updated.disconnect()
            le.arg = i
            le.value_updated.connect(self.handle_q)
            subgrid.addWidget(QLabel("Q%d"%(i+1)), i//2, i%2 * 3) # sqp_q1..q4, 1-based
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
        sq_actor.AddObserver(vtkCommand.ModifiedEvent, self.handle_orientation)

        row +=1
        grid.addWidget(QLabel("Volume"), row, 0)
        self.le_volume = LineEdit()
        self.le_volume.setReadOnly(True)
        #self.le_volume.setMaximumWidth(50)
        grid.addWidget(self.le_volume, row, 1,1,2)
        grid.addWidget(QLabel('m³'), row, 3)

        row += 1
        grid.addItem(QSpacerItem(20,10,QSizePolicy.Minimum,QSizePolicy.Expanding),row,0)

        row +=1
        hbox = QHBoxLayout()
        hbox.addItem(QSpacerItem(20,10,
                                 QSizePolicy.Expanding,
                                 QSizePolicy.Minimum))
        b = QPushButton("Cancel")
        b.clicked.connect(lambda arg: (self.close(), self.show_hide_plot(None,False)))
        b.setAutoDefault(False)
        hbox.addWidget(b)
        b = QPushButton("OK")
        b.clicked.connect(self.save)
        b.setAutoDefault(False)
        hbox.addWidget(b)
        grid.addLayout(hbox, row, 0, 1, 4)


    def mk_plot(self):
        # Interactive exponent plot
        p = pg.PlotWidget()
        p.setBackground('#efefef') # match mfix background
        p.setLabel('bottom', 'm')
        p.setLabel('left', 'n')
        p.setAspectLocked(1)
        p.setMouseEnabled(x=False, y=False)
        p.showGrid(x=True, y=True, alpha=1)
        p.scene().sigMouseClicked.connect(self.handle_click)
        # Override drag events
        class MyGraphItem(pg.GraphItem):
            def __init__(self,parent):
                pg.GraphItem.__init__(self)
                self.parent = parent
            def mouseDragEvent(self, ev):
                if ev.button() != Qt.LeftButton:
                    return
                ev.accept()
                self.parent.handle_click(ev)
        it = MyGraphItem(self)
        self.scene = sc = it.scatter


        p.addItem(it)
        xs = []
        ys = []
        # want m and n/m to be dyadic
        for m4 in range(8,97):
            m = m4 / 4
            for n4 in range(8, 97):
                n = n4/4
                if int(m4/n) == m4/n:
                    xs.append(m)
                    ys.append(n)
        self.ms = np.array(xs)
        self.ns = np.array(ys)
        self.move_dot(self.m, self.n)

        #it = pg.GridItem()
        #p.addItem(it)

        p.setWindowTitle("Favorable shape exponents")
        def intercept_close(ev):
            self.show_hide_plot(False)
            ev.ignore()
        p.closeEvent = intercept_close
        return p


    # Event handlers
    def handle_button(self, arg):
        iren = self.iren
        if arg == 1:
            self.style.SetCurrentStyleToTrackballActor()
            self.rotate_camera = False
            self.rotate_object = True
            #iren.RemoveObservers('MiddleButtonPressEvent')
            #iren.RemoveObservers('RightButtonPressEvent')
            iren.RemoveObservers('CharEvent')
        else:
            self.rotate_camera = True
            self.rotate_object = False
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
        self.update_volume()

    def update_volume(self):
        try:
            text = "%.6g" % sqp_volume(self.a, self.b, self.c, self.m, self.n)
        except:
            text = ''
        self.le_volume.setText(text)

    def setup(self, callback, P, a, b, c, d, m, n, q):
        self.suppress_diameter = True
        self.callback = callback
        self.P = P
        self.set_mn(m, n)
        self.move_dot(m, n)
        for k,v in ('a',a), ('b',b), ('c',c):
            le = self.widgets['le_'+k]
            le.updateValue(None, v) # Set lineedit value
            self.handle_le(le)  # Propagate event to sliders and model
        if d is None:
            d = sqp_d0(a,b,c,m,n)
        self.suppress_diameter = False
        self.set_diameter(d)
        le = self.widgets['le_d']
        le.updateValue(None, d) # Set lineedit value
        self.q = np.array(q, dtype=float)
        for i in range(4):
            self.widgets['le_q%d'%(i+1)].updateValue(None, q[i])
        self.apply_q()

        self.reset_view()

        # This should not be necessary but it causes the
        # z-axis labels to show up correctly on first view
        r = d/2
        r += 0.001
        self.axes_actor.SetBounds(-r,r,-r,r,-r,r)


    def save(self):
        if self.callback:
            self.callback(self.P, self.a, self.b, self.c, self.d,
                          self.m, self.n, self.q)
        self.show_hide_plot(None,False)
        self.close()

    def handle_le(self, w):
        k,v = w.key, w.value
        if k == 'd' and w.value in ('', None):
            v = sqp_d0(self.a, self.b, self.c, self.m, self.n)
            w.setText(str(v))
        setattr(self, k, v)
        self.set_slider_pos(k, v)
        if k in ('m', 'n'):
            if self.constrain:
                self.m, self.n = self.find_closest(self.m, self.n)
                self.widgets['le_m'].updateValue(None,self.m)
                self.widgets['le_n'].updateValue(None,self.n)
            self.sq_source.SetThetaRoundness(2/self.m)
            self.sq_source.SetPhiRoundness(2/self.n)
            self.sq_source.Update()
            self.move_dot(self.m, self.n)

        if k in ('a','b','c'):
            self.sq_source.SetScale(self.a, self.b, self.c)

        if k != 'd' and not self.suppress_diameter:
            le = self.widgets['le_d']
            self.d = sqp_d0(self.a, self.b, self.c, self.m, self.n)
            le.saved_value = self.d
            le.updateValue(None, self.d)
            self.set_slider_pos('d', self.d)
        self.set_diameter(self.d)
        self.sq_source.Update()
        self.rw.Render()
        self.update_volume()
        #self.reset_view()

    def set_diameter(self, d):
        if self.suppress_diameter:
            return
        self.d = d
        r = d/2
        self.sphere_source.SetRadius(r)
        self.sphere_source.Update()
        r += 0.001
        self.axes_actor.SetBounds(-r,r,-r,r,-r,r)

    def set_slider_pos(self, key, val):
        self.suppress_update = True
        s = self.widgets['s_'+key]
        if key in ('m', 'n'):
            x = (val-2)*16
            s.setValue(int(x))
        else:
            man, exp = fman(val), fexp(val)
            self.exp[key] = exp
            s.setValue(int(40*man))
        self.suppress_update = False

    def handle_s(self, key, val):
        if self.suppress_update:
            return

        if key in ('m', 'n'):
            x = min(2 + val/16, 24)
            if key == 'm':
                self.m = x
                self.sq_source.SetThetaRoundness(2/x)
            else:
                self.n = x
                self.sq_source.SetPhiRoundness(2/x)
            self.widgets['le_'+key].updateValue(None,x)
            self.move_dot(self.m, self.n)
        elif key in ('a', 'b', 'c', 'd'):
            s = self.widgets['s_'+key]
            if val == 400:
                # Hit right end of slider.  Fake a mouse release event
                # so we don't move right forever
                e = QMouseEvent(3, QPoint(0,0), Qt.MouseButton(Qt.LeftButton),
                                Qt.MouseButtons(Qt.NoButton),
                                Qt.KeyboardModifiers(Qt.NoModifier))
                QApplication.instance().sendEvent(s, e)
                s.setSliderPosition(40)
                val = 40
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
            if key == 'd':
                self.d = x
        if key != 'd' and not self.suppress_diameter:
            self.d = sqp_d0(self.a, self.b, self.c, self.m, self.n)
            self.set_slider_pos('d', self.d)
            le = self.widgets['le_d']
            le.saved_value = self.d
            le.updateValue(None, self.d)
        self.sq_source.Update()
        self.set_diameter(self.d)
        self.rw.Render()
        self.update_volume()

    def finish_s(self):
        if not self.constrain:
            return
        self.m, self.n = self.find_closest(self.m, self.n)
        self.widgets['le_m'].updateValue(None, self.m)
        self.widgets['le_n'].updateValue(None, self.n)
        self.handle_le(self.widgets['le_m'])
        self.handle_le(self.widgets['le_n'])

    def toggle_axes(self, arg):
        self.widgets['tb_axes'].setToolTip("Hide axes" if arg else "Show axes")
        self.axes_actor.SetVisibility(arg)
        self.rw.Render()

    def handle_constrain(self, arg):
        self.constrain = bool(arg)
        if self.constrain:
            self.finish_s()

    def find_closest(self, m, n):
        a = np.argmin(np.hypot(self.ms-m, self.ns-n))
        m,n = self.ms[a], self.ns[a]
        return (m,n)

    def show_hide_plot(self, ev, state=None):
        if state is not None:
            self.plot_shown = state
        else:
            self.plot_shown = not self.plot_shown
        if self.plot_shown:
            # Show plot to right of our window
            self.plot_button.setText("Hide exponent map")
            g = self.geometry()
            self.plot.setGeometry(g.left()+g.width(), g.top(), g.height(), g.height())
            self.plot.show()
            self.plot.update()
        else:
            self.plot_button.setText("Show exponent map")
            self.plot.hide()
        return False

    def show_bounding_sphere(self, val):
        self.sphere_actor.SetVisibility(val)
        self.rw.Render()

    def handle_orientation(self, actor, event):
        # Interactive rotation handler
        actor.ComputeMatrix()
        m = actor.GetMatrix()
        M = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                M[i,j] = m.GetElement(i,j)
        self.q = q = m_to_q(M)
        for i in range(4):
            self.widgets['le_q%d'%(i+1)].updateValue(None, q[i])
        #self.axes_actor.SetBounds(self.sphere_actor.GetBounds())
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
        self.rw.Render()

    def reset_q(self):
        vt = vtk.vtkTransform()
        vm = vtkMatrix4x4()
        vm.Identity()
        vt.SetMatrix(vm)
        self.sq_actor.SetOrientation(vt.GetOrientation())
        self.q = np.array([1,0,0,0], dtype=float)
        self.widgets['le_q1'].setText('1')
        self.widgets['le_q2'].setText('0')
        self.widgets['le_q3'].setText('0')
        self.widgets['le_q4'].setText('0')
        self.rw.Render()

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

    def handle_click(self, ev):
        if ev.button() != Qt.LeftButton:
            return
        sc = self.scene
        pt = sc.mapFromScene(ev.scenePos())
        m, n = pt.x(), pt.y()
        def clip(lo,x,hi):
            return lo if x<lo else hi if x>hi else x
        m = clip(2,m,24)
        n = clip(2,n,24)
        if self.constrain:
            m,n = self.find_closest(m, n)
        self.set_mn(m,n) # Update GUI and model
        self.d = sqp_d0(self.a, self.b, self.c, self.m, self.n)
        self.set_diameter(self.d)
        self.set_slider_pos('d', self.d)
        le = self.widgets['le_d']
        le.saved_value = self.d
        le.updateValue(None, self.d)

        self.move_dot(m,n)
        self.rw.Render()

    def move_dot(self, m,n):
        sc = self.scene
        sc.setData(list(self.ms)+[m],
                   list(self.ns)+[n],
                   tip=lambda x,y,data: '(%s,%s)'%(round(x,6),round(y,6)),
                   hoverable=True)
        sc.setSize([small]*len(self.ms) + [big])
        sc.setBrush([blue]*len(self.ms) + [green])



    def closeEvent(self, ev):
        if self.parent:
            self.parent.ui.solids.groupbox_sqp_shape.setDisabled(False)
        self.show_hide_plot(None,False)
        #self.close()
        self.hide()
        ev.ignore()

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
    sq_popup = SQPopup()
    q = [1,0,0,0]
    sq_popup.setup(None, None, 1,2,3,None,3,3,q)

    sq_popup.show()
    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    qapp.exec_()
    qapp.deleteLater()
    sys.exit()

if __name__ == '__main__':
    main()
