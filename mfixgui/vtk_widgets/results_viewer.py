# -*- coding: utf-8 -*-
import copy
import glob
import os
from collections.abc import Mapping
from qtpy import QtCore, QtGui, QtWidgets
import numpy as np

QIcon = QtGui.QIcon

from mfixgui.colormaps.color_maps import build_qicons
from mfixgui.tools import safe_int, safe_float
from mfixgui.tools.qt import get_icon, get_image_path, sub_icon_size, get_ui, SETTINGS
from mfixgui.vtk_widgets.pvd_reader import get_pvd_reader, PvdFile
from mfixgui.vtk_widgets.base import BaseVtkWidget
from mfixgui.vtk_widgets.constants import *
from mfixgui.vtk_widgets.tools import  safe_combo, update_combo, remove_vtk_objects
from mfixgui.vtk_widgets.colormap_dialog import ColorMapPopUp
from mfixgui.vtk_widgets.dialogs import CreateMoviePopUp

# graphics libraries
try:
    import vtk
    from mfixgui.colormaps.color_maps import build_vtk_lookup_tables
    VTK_AVAILABLE = True
    LOOKUP_TABLES = build_vtk_lookup_tables()

    # check for point gaussian, vtk 7+
    POINT_GAUSSIAN = hasattr(vtk, 'vtkPointGaussianMapper')
except ImportError:
    vtk = None
    VTK_AVAILABLE = False
    LOOKUP_TABLES = {}
    POINT_GAUSSIAN = False

DEFAULT_MAXIMUM_POINTS = 100000
MAX_PLAYBACK_DELAY = 1000  # ms
DEFAULT_PLAYBACK_DELAY = 100  # ms

GEOMETRY = ['cells', 'points', 'geometry', 'color_bar', 'time_label', 'axes']
INDEX_DICT = {'x': 0, 'y': 1, 'z': 2, 'mag': -1}
DEFAULT_OBJECT_DATA = {
    'type': '', # stl, vtu, vtp
    'visible': False,
    'opacity': 1.0,
    'variable_color': {},
    'color': DEFAULT_GEO_COLOR,
}
DEFAULT_COLOR_MAP = 'viridis'
DEFAULT_VARIABLE_DATA = {
    'component': 'mag',
    'single_color': False,
    'color': DEFAULT_GEO_COLOR,
    'color_map': DEFAULT_COLOR_MAP,
    'from': None,
    'to': None,
    'range': [[0, 1]],
    'magnitude': [[0, 1]],
    'reversed': False
}

def find_nearest(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()


def clean_dict(dirty_dict):
    """remove qcolor objects save the hex values"""
    cd = {}
    for k, v in dirty_dict.items():
        if isinstance(v, dict):
            cd[k] = clean_dict(v)
        elif isinstance(v, QtGui.QColor):
            cd[k] = v.name()
        elif 'vtk' in str(v):
            continue
        elif isinstance(v, PvdFile):
            continue
        else:
            cd[k] = v
    return cd


def drop_none(d):
    """remove keys that have None values"""
    clean_d = {}
    for k, v in d.items():
        if isinstance(v, dict):
            clean_d[k] = drop_none(v)
        elif v is None:
            continue
        else:
            clean_d[k] = v
    return clean_d


def minmax(list1, list2):
    r = []
    for a, b in zip(list1, list2):
        r.append([min(a[0], b[0]), max(a[1], b[1])])
    return r


def qcolor_dict(d):
    """the reverse of clean_dict, change hex back to qcolor"""
    qd = {}
    for k, v in d.items():
        if isinstance(v, dict):
            qd[k] = qcolor_dict(v)
        elif isinstance(v, str) and v.startswith('#'):
            qd[k] = QtGui.QColor(v)
        else:
            qd[k] = v
    return qd


def update(d, u):
    for k, v in u.items():
        if isinstance(v, Mapping):
            r = update(d.get(k, {}), v)
            d[k] = r
        else:
            d[k] = u[k]
    return d


class GraphicsVtkWidget(BaseVtkWidget):
    """vtk widget for showing results"""
    def __init__(self, parent=None, load=False):
        BaseVtkWidget.__init__(self, parent)

        self.loading = load
        self.frame_index = -1
        self.time = 0.0
        self.time_map = []
        self.time_format = '%.2f s'
        self.time_label_color = DEFAULT_TEXT_COLOR
        self.pipeline = {}
        self.active_objects = []
        self.color_maps = {'vtp':{}, 'vtu':{}}

        # PVD reader
        self.pvd_reader = get_pvd_reader()
        self.pvd_reader.signal_new_files.connect(self.refresh_pvd_list)
        self.pvd_reader.signal_removed_files.connect(self.remove_pvd_list)

        # dialogs
        self.color_dialog = ColorMapPopUp(self)
        self.color_dialog.applyEvent.connect(self.change_variable_color)
        self.movie_dialog = CreateMoviePopUp(self)

        self.init_sidebar()
        self.init_top_toolbar()
        self.init_vtk()
        self.refresh_pvd_list()

        # start the timers
        self.play_timer = QtCore.QTimer()
        self.play_timer.timeout.connect(self.handle_next)


    def init_sidebar(self):
        # build to sidebar
        sb = self.sidebar = get_ui('results_viewer_sidebar.ui')
        self.grid_layout.addWidget(self.sidebar, 1, 1)
        self.hide_show_sidebar(False)

        format_help_text = "printf-style format string such as %.2f, %.3g, or %.1e"
        # tweak scrollarea
        # monkey patching
        def customResize(event):
            sc = sb.scrollAreaWidgetContents
            sa = sb.scrollArea
            sbar = sa.verticalScrollBar()
            size_hint = sc.minimumSizeHint()
            width = size_hint.width()
            if size_hint.height() > sa.height():
                width += sbar.width()
            sa.setMinimumWidth(width)

        sb.scrollAreaWidgetContents.resizeEvent = customResize

        # hide/Show
        tb = sb.toolButton_hide_show
        tb.setIcon(get_icon('left.svg'))
        tb.setToolTip('Open sidebar')
        tb.clicked.connect(lambda ignore: self.hide_show_sidebar())
        tb.setIconSize(sub_icon_size())

        def hide_show_options_widget(tb, opt_widget, visible=None):
            if visible is None:
                visible = not opt_widget.isVisible()
            opt_widget.setVisible(visible)
            tb.setArrowType(QtCore.Qt.UpArrow if visible else QtCore.Qt.DownArrow)

        # toggles
        for tb, opt in [
                (sb.toolButton_objects_toggle, sb.widget_object_opts),
                (sb.toolButton_color_bar_toggle, sb.widget_color_bar_opts),
                (sb.toolButton_time_label_toggle, sb.widget_time_label_opts),
                (sb.toolButton_playback_toggle, sb.widget_playback_opts),
                (sb.toolButton_image_stack_toggle, sb.widget_image_stack_opts),
                ]:
            hide_show_options_widget(tb, opt, False)
            tb.clicked.connect(
                lambda down, tb=tb, w=opt: hide_show_options_widget(tb, w))

        # visible buttons
        tb_list = [sb.toolButton_time_label,
                   sb.toolButton_axes_vis, sb.toolButton_image_stack,
                   ]
        # icons
        icon_size = sub_icon_size()
        for tb, icon in zip(tb_list, ['time', 'axes', 'camera_stack']):
            tb.setIcon(get_icon(icon+'.svg'))
            tb.setIconSize(icon_size)
        # Use disabled toolbuttons for icon holders.  They are not clickable.
        # We do this so they scale with font size changes.
        for tb, icon_name in zip([sb.toolButton_label_objects,
                                  sb.toolButton_label_colorbar,
                                  sb.toolButton_label_playback],
                                 ['geometry', 'gradient', 'speed']):
            icon = get_icon(icon_name)
            pixmap = icon.pixmap(icon_size)
            icon.addPixmap(pixmap, QIcon.Disabled, QIcon.Off) # Avoid greyed-out icon
            tb.setIcon(icon)
            tb.setIconSize(icon_size)
            tb.setAutoRaise(True)
            tb.setEnabled(False) # It's not clickable, just an icon holder


        # -- pipeline --
        sb.treeWidget_pipeline.setStyleSheet(
            "QTreeView::indicator:unchecked {image: url(%s);}"
            "QTreeView::indicator:checked {image: url(%s);}"
            % (get_image_path('visibilityofflight.svg'),
               get_image_path('visibility.svg'))
            )

        sb.treeWidget_pipeline.itemSelectionChanged.connect(
            self.object_selection_changed)
        sb.treeWidget_pipeline.itemClicked.connect(self.object_clicked)

        # variables
        sb.comboBox_points_variable.activated.connect(self.variable_selection_changed)
        sb.comboBox_cells_variable.activated.connect(self.variable_selection_changed)

        # component
        sb.comboBox_points_component.activated.connect(self.component_selection_changed)
        sb.comboBox_cells_component.activated.connect(self.component_selection_changed)

        # color
        sb.toolButton_geometry_color.clicked.connect(self.geometry_color_pressed)

        # opacity
        sb.doubleSpinBox_points_opacity.valueChanged.connect(self.opacity_changed)
        sb.doubleSpinBox_cells_opacity.valueChanged.connect(self.opacity_changed)
        sb.doubleSpinBox_geometry_opacity.valueChanged.connect(self.opacity_changed)

        # particle mapper
        sb.comboBox_particle_mapper.activated.connect(lambda _: self.change_particle_mapper())
        sb.comboBox_particle_shader.addItems(SPLAT_SHADERS.keys())
        sb.comboBox_particle_shader.activated.connect(self.change_shader)
        sb.comboBox_glyph.addItems(PRIMITIVE_DICT.keys())
        sb.comboBox_glyph.activated.connect(self.change_glyph)

        sb.lineedit_maximum_glyph_points.textChanged.connect(self.change_max_points)

        # force chain
        sb.doubleSpinBox_forcechain_radius.valueChanged.connect(self.change_forcechain_radius)
        sb.doubleSpinBox_forcechain_factor.valueChanged.connect(self.change_forcechain_factor)

        # VTU: cell style
        sb.comboBox_cells_style.activated.connect(self.cell_style_changed)

        # geometry
        sb.comboBox_geometry_style.activated.connect(self.geometry_style_changed)

        # -- color bar --
        sb.treeWidget_colorbars.setStyleSheet(
            "QTreeView::indicator:unchecked {image: url(%s);}"
            "QTreeView::indicator:checked {image: url(%s);}"
            % (get_image_path('visibilityofflight.svg'),
               get_image_path('visibility.svg'))
            )
        sb.treeWidget_colorbars.itemSelectionChanged.connect(
            self.colorbar_selection_changed)
        sb.treeWidget_colorbars.itemClicked.connect(self.colorbar_clicked)

        sb.toolButton_color_scale.clicked.connect(self.change_color_pressed)
        self.set_color_button(sb.toolButton_color_scale, {'single':False})
        sb.lineedit_color_bar_title.textChanged.connect(self.change_color_bar_title)
        sb.spinBox_color_bar_title_size.valueChanged.connect(self.change_color_bar_title_size)
        sb.spinBox_color_bar_n_labels.valueChanged.connect(self.change_color_bar_n_labels)
        sb.lineedit_color_bar_format.textChanged.connect(self.change_color_bar_format)
        sb.lineedit_color_bar_format.setToolTip(format_help_text + " for color bar label")
        sb.spinBox_color_bar_label_size.valueChanged.connect(self.change_color_bar_label_size)
        sb.comboBox_color_bar_position.activated.connect(self.change_color_bar_loc)
        sb.toolButton_color_bar_color.clicked.connect(self.change_color_bar_color)
        sb.toolButton_color_bar_color.setStyleSheet(
            "QToolButton{{ background: {};}}".format(DEFAULT_TEXT_COLOR.name()))
        sb.doubleSpinBox_color_bar_opacity.valueChanged.connect(self.change_color_bar_opacity)

        # -- time label --
        le = sb.lineedit_time_label_format
        le.setText(self.time_format)
        le.textChanged.connect(self.handle_time_label_format)
        le.setToolTip(format_help_text + " for time label")

        sb.comboBox_time_label_loc.currentIndexChanged.connect(self.change_time_label_loc)
        sb.toolButton_time_label.toggled.connect(lambda vis: self.change_visibility(self.time_label, vis))
        sb.toolButton_time_label_color.clicked.connect(self.change_time_label_color)
        sb.toolButton_time_label_color.setStyleSheet(
            "QToolButton{{ background: {};}}".format(DEFAULT_TEXT_COLOR.name()))
        sb.spinBox_time_label_size.valueChanged.connect(self.change_time_label_text_size)

        # --- ---
        sb.toolButton_axes_vis.toggled.connect(lambda vis: self.change_visibility(self.axes, vis))


        # playback
        slider = self.speed_slider = sb.horizontalSlider_playback_speed
        slider.setRange(0, MAX_PLAYBACK_DELAY-10) # delay = MAX - speed
        slider.setValue(MAX_PLAYBACK_DELAY - DEFAULT_PLAYBACK_DELAY)
        slider.sliderReleased.connect(self.handle_speed_changed)

        # image stack
        project_path = os.path.dirname(self.gui.get_project_file())
        sb.lineedit_image_stack_dir.setText(project_path)
        sb.toolButton_image_stack_browse.setIcon(get_icon('folder.svg'))
        sb.toolButton_image_stack_browse.clicked.connect(self.browse_image_stack)

        # hide resolution widgets
        if not safe_int(SETTINGS.value('enable_screenshot_res', 0),0):
            for wid in [sb.label_image_stack_width, sb.lineedit_image_stack_width,
                        sb.label_image_stack_height, sb.lineedit_image_stack_height]:
                wid.setVisible(False)

        # config default view
        sb.stackedWidget_object_options.setCurrentIndex(0)

    def init_top_toolbar(self):

        # create and add buttons to the toolbar
        self.init_base_toolbar()

        # --- video capture ---
        self.toolbutton_create_movie = QtWidgets.QToolButton()

        # --- play/stop/forward/backward controls ---
        self.toolbutton_first = QtWidgets.QToolButton()
        self.toolbutton_back = QtWidgets.QToolButton()
        self.toolbutton_play = QtWidgets.QToolButton()
        self.toolbutton_stop = QtWidgets.QToolButton()
        self.toolbutton_stop.setEnabled(False)
        self.toolbutton_next = QtWidgets.QToolButton()
        self.toolbutton_last = QtWidgets.QToolButton()
        self.toolbutton_repeat = QtWidgets.QToolButton()
        self.toolbutton_repeat.setCheckable(True)

        icon_size = sub_icon_size()
        for tb, callback, icon, tooltip in [
                (self.toolbutton_create_movie, self.create_movie, 'videocam', 'Create a video'),
                (self.toolbutton_first, self.handle_first, 'first', 'First'),
                (self.toolbutton_back, self.handle_back, 'back', 'Previous'),
                (self.toolbutton_play, self.handle_play, 'play', 'Play'),
                (self.toolbutton_stop, self.handle_stop, 'stop', 'Stop'),
                (self.toolbutton_next, self.handle_next, 'next', 'Next'),
                (self.toolbutton_last, self.handle_last, 'last', 'Last'),
                (self.toolbutton_repeat, None, 'autorenew', 'Repeat from beginning')]:

            if callback is not None:
                tb.clicked.connect(callback)
            tb.setIcon(get_icon(icon + '.svg'))
            tb.setIconSize(icon_size)
            tb.setToolTip(tooltip)
            tb.setAutoRaise(True)
            tb.setFocusPolicy(QtCore.Qt.ClickFocus)

        hspacer = QtWidgets.QSpacerItem(999, 10, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Maximum,)

        self.time_lineedit = QtWidgets.QLineEdit()
        self.time_lineedit.editingFinished.connect(self.handle_time_edited)
        self.time_lineedit.setFixedWidth(100)
        t_lbl = QtWidgets.QLabel('s')

        for btn in [self.toolbutton_create_movie,
                    self.toolbutton_first, self.toolbutton_back,
                    self.toolbutton_play, self.toolbutton_stop,
                    self.toolbutton_next, self.toolbutton_last,
                    self.toolbutton_repeat,
                    hspacer, self.time_lineedit, t_lbl,
                    ]:
            if btn == hspacer:
                self.button_bar_layout.addSpacerItem(btn)
            else:
                self.button_bar_layout.addWidget(btn)

        self.button_bar_layout.addStretch()

    @property
    def project_dir(self):
        return self.gui.get_project_dir()

    def refresh_pvd_list(self):
        sb = self.sidebar

        fnames = self.pvd_reader.pvd_filenames()
        stl_files = [os.path.basename(p) for p in glob.glob(os.path.join(self.project_dir, '*.stl'))]

        new_names = list(fnames) + stl_files

        n_list = []
        for new_name in new_names:
            if new_name not in self.active_objects:
                self.add_object(new_name)
                n_list.append(new_name)
                self.change_frame(0)

        # show something when vtk output files exist
        if not self.loading and n_list and not self.active_objects:
            self.default_view()

        self.active_objects += n_list

        self.make_time_series()
        self.update_sq_glyph_params()

    def remove_pvd_list(self, removed):
        for name in removed:
            self.remove_object(name)
            if name in self.active_objects:
                self.active_objects.remove(name)

    def hide_show_sidebar(self, visible=None):
        """hide/show the sidebar widgets"""
        sb = self.sidebar
        if visible is None:
            visible = not sb.widget_color_bar.isVisible()

        tb = sb.toolButton_hide_show
        tb.setIcon(get_icon('right.svg' if visible else 'left.svg'))
        tb.setToolTip('Collapse sidebar' if visible else 'Open sidebar')
        tb.setIconSize(sub_icon_size())
        tb.setDown(False)

        for wid in [sb.widget_objects, sb.widget_color_bar,
                    sb.widget_time_label, sb.label_axes,
                    sb.widget_playback, sb.widget_image_stack]:
            wid.setVisible(visible)
        size_hint = sb.scrollAreaWidgetContents.minimumSizeHint()
        width = size_hint.width()
        sb.setMaximumWidth(width if not visible else 16777215)

    def set_state(self, state):
        '''load a saved vtk state'''
        self.defer_render = True # defer rendering vtk until load complete
        self.loading = True
        sb = self.sidebar
        if not state:
            self.default_view()
        else:
            state = qcolor_dict(state)

            # look for result viewer 2 data
            if 'pipeline' in state:
                # load objects
                for k, v in drop_none(state.get('pipeline')).items():
                    if k in self.pipeline:
                        self.pipeline[k].update(v)
                        if v.get('visible', False):
                            item = sb.treeWidget_pipeline.findItems(k, QtCore.Qt.MatchExactly|QtCore.Qt.MatchRecursive)[0]
                            item.setCheckState(0, QtCore.Qt.Checked)
                            self.object_clicked(item)
                        obj_data = self.pipeline[k]
                        self.change_opacity(obj_data)
                        var = obj_data.get('variable')
                        if obj_data.get('type', '') != 'stl':
                            self.change_color(k, obj_data.get('variable_color', {}.get(var, {})))
                # load color maps
                for k, v in drop_none(state.get('color_maps')).items():
                    if k in self.color_maps:
                        for k1, v1 in v.items():
                            vdict = self.color_maps[k].setdefault(k1, {})
                            vdict.update(v1)
                            self.change_variable_color([k, k1], v1)
                            if v1.get('visible', False):
                                item = sb.treeWidget_colorbars.findItems(k1, QtCore.Qt.MatchExactly|QtCore.Qt.MatchRecursive)[0]
                                item.setCheckState(0, QtCore.Qt.Checked)
                                self.colorbar_clicked(item)

            # try to load old viz data
            else:
                visible = state.get('visible', {})
                names = [state.get(p) for p, v in [['vtu_pattern', 'cells'], ['vtp_pattern', 'points']] if visible.get(v)]
                if visible.get('geometry'):
                    names.append('geometry.stl')
                root = sb.treeWidget_pipeline.invisibleRootItem()
                for i in range(root.childCount()):
                    item = root.child(i)
                    if not item:
                        continue
                    text = item.text(0)
                    if text.lower() in names:
                        item.setCheckState(0, QtCore.Qt.Checked)
                        self.object_clicked(item)

            # time label
            time_label_format = state.get('time_label_format', '%.2f s')
            # Translate from str.format to printf-style
            if not '%' in time_label_format:
                time_label_format = time_label_format.replace('{:', '%')
                time_label_format = time_label_format.replace('}', ' ')
            sb.lineedit_time_label_format.setText(time_label_format)
            sb.comboBox_time_label_loc.setCurrentText(state.get('time_label_pos', 'top right'))
            sb.spinBox_time_label_size.setValue(state.get('time_label_text_size', 24))
            self.change_time_label_color(color=state.get('time_label_color',  DEFAULT_TEXT_COLOR))

            # camera
            camera_state = state.get('camera', None)
            if camera_state is not None:
                self.set_camera_state(camera_state)

            # image stack
            project_path = os.path.dirname(self.gui.get_project_file())
            sb.lineedit_image_stack_dir.setText(state.get('image_stack_dir', project_path))
            sb.lineedit_image_stack_prefix.setText(state.get('image_stack_prefix', 'frame_'))
            sb.lineedit_image_stack_width.setText(state.get('image_stack_width', '1920'))
            sb.lineedit_image_stack_height.setText(state.get('image_stack_height', '1080'))
            safe_combo(sb.comboBox_image_stack_type, state.get('image_stack_type', 'png'))
            sb.checkBox_trans_back.setChecked(state.get('image_stack_trans', False))

        # set the current frame
        time = state.get('time', None)
        if time:
            self.change_frame(time)
        elif self.time_series:
            frame = state.get('frame', 0)
            if frame < len(self.time_series):
                time = self.time_series[frame]
            self.change_frame(time)

        # update widgets
        self.object_selection_changed()
        self.colorbar_selection_changed()

        self.loading = False
        self.render(defer_render=False)  # render

    def default_view(self):
        # new vis, look for some items to show
        sb = self.sidebar
        root = sb.treeWidget_pipeline.invisibleRootItem()
        picked_pvd = False
        for i in range(root.childCount()):
            item = root.child(i)
            if not item:
                continue
            text = item.text(0).lower()
            if 'geometry.stl' in text:
                item.setCheckState(0, QtCore.Qt.Checked)
                self.object_clicked(item)
            elif not picked_pvd and not text.endswith('.stl'):
                item.setCheckState(0, QtCore.Qt.Checked)
                self.object_clicked(item)
                picked_pvd = True
        # show a color bar
        item = sb.treeWidget_colorbars.currentItem()
        if item:
            item.setCheckState(0, QtCore.Qt.Checked)
            self.colorbar_clicked(item)
        self.reset_view()
        self.render()

    def get_state(self):
        '''collect a dictionary of values to save'''
        sb = self.sidebar
        state = {
            'time':   self.time,
            'camera': self.get_camera_state(),

            'pipeline': self.pipeline,
            'color_maps': self.color_maps,

            # time label
            'time_label_format': sb.lineedit_time_label_format.text(),
            'time_label_pos': sb.comboBox_time_label_loc.currentText(),
            'time_label_color': self.time_label_color,
            'time_label_text_size': sb.spinBox_time_label_size.value(),

            # image stack
            'image_stack_dir': sb.lineedit_image_stack_dir.text(),
            'image_stack_prefix': sb.lineedit_image_stack_prefix.text(),
            'image_stack_width': sb.lineedit_image_stack_width.text(),
            'image_stack_height': sb.lineedit_image_stack_height.text(),
            'image_stack_type': sb.comboBox_image_stack_type.currentText(),
            'image_stack_trans': sb.checkBox_trans_back.isChecked(),
        }

        return clean_dict(state)

    def reset(self):
        self.play_timer.stop()
        if not self.deleting:
            self.vtkRenderWindow.Finalize()

    def init_vtk(self):

        self.actors = {'time_label': self.time_label,
                       'axes':       self.axes}
        self.mappers = {}
        self.lookuptables = {}

        self.time_label.SetVisibility(True)
        self.time_label.SetInput(self.time_format % self.time)

    def create_movie(self):
        self.toolbutton_repeat.setChecked(False)
        self.movie_dialog.popup()

    def showEvent(self, event):
        # has to be called after the widget is visible
        self.vtkiren.Initialize()

    def hideEvent(self, event):
        self.stop()

    def close(self):
        BaseVtkWidget.close(self)
        self.reset()
        # clean up timer(s)
        self.play_timer.stop()

    def browse_image_stack(self):
        sb = self.sidebar
        path = sb.lineedit_image_stack_dir.text()

        filename = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Select a image stack directory", path)
        if isinstance(filename, (tuple, list)):
            filename = filename[0]
        if not filename:
            return
        sb.lineedit_image_stack_dir.setText(filename)

    def handle_speed_changed(self):
        if self.play_timer.isActive():
            self.play_timer.stop()
            self.handle_play()

    def stop(self):
        self.play_timer.stop()
        self.toolbutton_play.setEnabled(True)
        self.toolbutton_stop.setEnabled(False)

    def handle_stop(self):
        self.stop()

    def handle_play(self):
        delay_ms = max(10, MAX_PLAYBACK_DELAY-self.speed_slider.value())
        self.play_timer.start(delay_ms)
        self.toolbutton_play.setEnabled(False)
        self.toolbutton_stop.setEnabled(True)

    def make_time_series(self):
        # TODO: handle changes in dt
        self.time_series = [0]
        dt = np.inf
        # look for the pvd file with the shortest dt
        for pvd in self.pvd_reader:
            if pvd.dt < dt:
                dt = pvd.dt
                self.time_series = pvd.time_series

    def handle_first(self):
        if self.time_series:
            self.change_frame(self.time_series[0])

    def handle_back(self):
        if self.time_series and self.time is not None:
            idx = find_nearest(self.time_series, self.time)
            idx = max(idx, 0)
            self.change_frame(self.time_series[idx - 1])

    def handle_next(self):
        if self.time_series and self.time is not None:
            idx = find_nearest(self.time_series, self.time)
            if idx >= len(self.time_series)-1:
                if self.toolbutton_repeat.isChecked():
                    idx = -1
                else:
                    return
            self.change_frame(self.time_series[idx + 1])

    def handle_last(self):
        if self.time_series:
            self.change_frame(self.time_series[-1])

    def change_frame(self, time, force=False):
        self.time = time
        for pvd in self.pvd_reader:
            pvd.change_time(time)
        if time is not None:
            self.time_lineedit.setText('%.3f' % time)
            self.time_label.SetInput(self.time_format % time)

        obj_data = self.current_object_data()
        if obj_data is not None and obj_data['type'] == 'vtp':
            pvd = obj_data.get('pvd')
            if pvd is not None:
                array_info = pvd.get_arrayinfo()
                var = list(array_info.keys())[0]
                if hasattr(array_info[var], 'number_of_tuples'):
                    self.sidebar.lineedit_point_count.setText('{}'.format(
                        array_info[var].number_of_tuples))
        self.render()

        sb = self.sidebar
        if sb.toolButton_image_stack.isChecked():
            path = sb.lineedit_image_stack_dir.text()
            prefix = sb.lineedit_image_stack_prefix.text()
            width = safe_int(sb.lineedit_image_stack_width.text(), 1920)
            height = safe_int(sb.lineedit_image_stack_height.text(), 1080)
            type_ = sb.comboBox_image_stack_type.currentText()
            trans = sb.checkBox_trans_back.isChecked()
            index = self.time_series.index(time)
            if os.path.exists(path):
                self.screenshot(
                    True,
                    fname=os.path.join(
                        path, f'{prefix}{index:06d}.{type_}'),
                    size=(width, height),
                    transparent=trans)

    def handle_time_edited(self):
        try:
            time = float(self.time_lineedit.text())
        except ValueError:
            self.time_lineedit.setText(f'{self.time:.3f}')
            return
        idx = find_nearest(self.time_series, time)
        time = self.time_series[idx]
        self.change_frame(time)

    # --- pipeline functions ---
    def add_object(self, obj_name=None):
        sb = self.sidebar

        pvd = None
        if obj_name.endswith('.stl'):
            icon = 'geometry.svg'
            obj_type = 'stl'
        else:
            pvd = self.pvd_reader.pvd(obj_name)
            if pvd is None:
                # how did we get here? See issue/1573
                self.gui.print_internal(f"Error: Could not get PVD file: {obj_name}.",
                                        color='red')
                return
            obj_type = pvd.type_
            icon = 'dem.svg' if obj_type == 'vtp' else 'grid.svg'

        # get existing object data - saved from a project reset
        if obj_name in self.pipeline:
            obj_data = self.pipeline.get(obj_name)
            obj_data['pvd'] = pvd
            self.build_pipeline(obj_data)
        # create new object data
        else:
            obj_data = self.pipeline.setdefault(obj_name, copy.deepcopy(DEFAULT_OBJECT_DATA))
            obj_data['type'] = obj_type
            obj_data['obj_name'] = obj_name
            obj_data['pvd'] = pvd
            if obj_type == 'stl':
                obj_data['opacity'] = 0.4

        # add tree widget item
        # find existing tree item and enable
        item = sb.treeWidget_pipeline.findItems(obj_name, QtCore.Qt.MatchExactly|QtCore.Qt.MatchRecursive)
        if item:
            item[0].setDisabled(False)
        # add new tree item
        else:
            toplevel = QtWidgets.QTreeWidgetItem([obj_name])
            toplevel.setIcon(0, get_icon(icon))
            toplevel.setFlags(toplevel.flags() | QtCore.Qt.ItemIsUserCheckable)
            toplevel.setCheckState(0, QtCore.Qt.Unchecked)

            sb.treeWidget_pipeline.addTopLevelItem(toplevel)
            sb.treeWidget_pipeline.setCurrentItem(toplevel)

    def remove_object(self, obj_name):
        obj_data = self.pipeline.get(obj_name)
        actor = obj_data.get('actor', None)
        if actor is not None:
            self.vtkrenderer.RemoveActor(actor)
        self.pipeline[obj_name] = remove_vtk_objects(obj_data)
        self.pipeline[obj_name].pop('pvd', None)
        item = self.sidebar.treeWidget_pipeline.findItems(obj_name, QtCore.Qt.MatchExactly|QtCore.Qt.MatchRecursive)[0]
        item.setDisabled(True)

        self.render()

    def build_pipeline(self, obj_data):
        obj_type = obj_data.get('type')
        obj_name = obj_data.get('obj_name')
        pvd = obj_data.get('pvd')
        if pvd is None and obj_type != 'stl':
            return

        # make sure the object is selected
        # item = self.sidebar.treeWidget_pipeline.findItems(obj_name, QtCore.Qt.MatchExactly|QtCore.Qt.MatchRecursive)
        # self.sidebar.treeWidget_pipeline.setCurrentItem(item[0])

        if obj_type == 'vtp':
            self.build_vtp_pipeline(pvd, obj_data)
        elif obj_type == 'vtu':
            self.build_vtu_pipeline(pvd, obj_data)
        elif obj_type == 'stl':
            self.build_stl_pipeline(obj_data)

        if obj_type != 'stl':
            # set the initial color
            array_info = pvd.get_arrayinfo()

            # variable
            vars = list(array_info.keys())
            if vars:
                var = obj_data.setdefault('variable', vars[0])
                var_data = obj_data['variable_color'].setdefault(var, copy.deepcopy(DEFAULT_VARIABLE_DATA))

                self.change_color(obj_name, var_data)

        # issues/1835
        if 'opacity' in obj_data:
            self.change_opacity(obj_data)

    def make_actor(self, shiny=True):
        actor = vtk.vtkActor()
        if shiny:
            self.set_actor_shiny_props(actor)
        return actor

    def set_actor_shiny_props(self, actor):
        # actor.SetForceOpaque(True)
        actor.GetProperty().SetDiffuse(1.0)
        actor.GetProperty().SetSpecular(.7)
        actor.GetProperty().SetSpecularPower(50.0)
        # actor.GetProperty().SetInterpolationToGouraud()
        return actor

    def build_vtu_pipeline(self, pvd, obj_data):
        reader = pvd.reader

        # filter out unused arrays
        pass_array = obj_data['entry'] = vtk.vtkPassArrays()
        pass_array.SetInputConnection(reader.GetOutputPort())

        # extract surface
        surface_extract = vtk.vtkDataSetSurfaceFilter()
        surface_extract.SetInputConnection(pass_array.GetOutputPort())

        ugrid_cell_mapper = obj_data['mapper_cells'] = vtk.vtkDataSetMapper()
        ugrid_cell_mapper.SetInputConnection(surface_extract.GetOutputPort())
        ugrid_cell_mapper.SetScalarVisibility(True)
        ugrid_cell_mapper.SetScalarModeToUseCellFieldData()

        actor = obj_data['actor'] = obj_data['actor_cells'] = self.make_actor(shiny=False)
        actor.SetMapper(ugrid_cell_mapper)

        self.vtkrenderer.AddActor(actor)
        cell_style = obj_data.get('cell_style', 'cells') == 'cells'
        actor.SetVisibility(cell_style)  # hide cells by default

        # nodes
        ugrid_cell_to_points = vtk.vtkCellDataToPointData()
        ugrid_cell_to_points.SetInputConnection(surface_extract.GetOutputPort())

        ugrid_node_mapper = obj_data['mapper_nodes'] = vtk.vtkDataSetMapper()
        ugrid_node_mapper.SetInputConnection(ugrid_cell_to_points.GetOutputPort())
        ugrid_node_mapper.SetScalarVisibility(True)
        ugrid_node_mapper.SetScalarModeToUsePointFieldData()
        #ugrid_node_mapper.UseLookupTableScalarRangeOn()
        #ugrid_node_mapper.InterpolateScalarsBeforeMappingOn()

        actor = obj_data['actor_nodes'] = self.make_actor(shiny=False)
        actor.SetMapper(ugrid_node_mapper)
        actor.SetVisibility(not cell_style)  # hide nodes by default

        self.vtkrenderer.AddActor(actor)

    def build_vtp_pipeline(self, pvd, obj_data):
        reader = pvd.reader
        array_info = pvd.get_arrayinfo()
        mtext = obj_data.get('mtext', None)
        if mtext is None:
            if 'FORCE_CHAIN_LENGTH' in array_info:
                mtext = 'force chain'
            elif array_info.get('field_data', {}):
                mtext = 'glyphs'
            else:
                mtext = 'point gaussian'
        actor = obj_data['actor'] = self.make_actor()

        self.change_particle_mapper(mtext, obj_data)
        self.vtkrenderer.AddActor(actor)
        safe_combo(self.sidebar.comboBox_particle_mapper, mtext)
        self.show_particle_mapper_options(mtext, obj_data)

    def change_particle_mapper(self, mtext=None, obj_data=None):
        sb = self.sidebar
        if obj_data is None:
            obj_name = self.current_object()
            obj_data = self.pipeline.get(obj_name)
        else:
            obj_name = self.pipeline.get('obj_name')
        pvd = obj_data.get('pvd')
        if pvd is None:
            return
        reader = pvd.reader
        array_info = pvd.get_arrayinfo()

        # Make a point mask to limit points
        mask = obj_data.get('point_mask', None)
        if mask is None:
            mask = vtk.vtkMaskPoints()
            mask.SetInputConnection(reader.GetOutputPort())
            mask.RandomModeOn()
            mask.SetRandomModeType(2)
            mask.SetMaximumNumberOfPoints(DEFAULT_MAXIMUM_POINTS)
            obj_data['point_mask'] = mask

        if mtext is None:
            mtext = sb.comboBox_particle_mapper.currentText()

        if mtext == 'force chain':
            pg_mapper, tube_filter = self.make_force_chain_mapper(mask)
            obj_data['mapper'] = pg_mapper
            obj_data['tube_filter'] = tube_filter
            mtext = 'force chain'
        elif mtext == 'glyphs':
            # look for user selected glyph
            glyph_name = obj_data.get('glyph_name', None)
            # guess glyph
            if glyph_name is None:
                glyph_name = 'superquadric' if array_info.get('field_data', {}) else 'sphere'
                obj_data['glyph_name'] = glyph_name
            pg_mapper = obj_data['mapper'] = self.make_GlyphMapper(mask, array_info, pvd, glyph_name=glyph_name)
            mtext = 'glyphs'
        else:
            # setup the point gaussian mapper by default
            pg_mapper = obj_data['mapper'] = self.make_PointGaussianMapper(mask, array_info)
            mtext = 'point gaussian'
        actor = obj_data.get('actor')
        if not actor:
            return
        actor.SetMapper(pg_mapper)
        obj_data['mtext'] = mtext

        # update color map
        var = obj_data.get('variable')
        var_data = obj_data.get('variable_color', {}).get(var, {})
        self.change_color(obj_name, var_data)

        self.show_particle_mapper_options(mtext, obj_data)
        self.render()

    def show_particle_mapper_options(self, mtext=None, obj_data=None):
        sb = self.sidebar
        if obj_data is None:
            obj_data = self.current_object_data()
        if mtext is None:
            mtext = sb.comboBox_particle_mapper.currentText()
        show = [False]*10  # shader, glyph, force chain, force chain
        if mtext == 'glyphs':
            show[1] = True
        elif mtext == 'force chain':
            show[2] = True
            show[3] = True
        else:
            show[0] = True

        wlist = [
            (sb.label_shader, sb.comboBox_particle_shader),
            (sb.label_glyph, sb.comboBox_glyph),
            (sb.label_forcechain_radius, sb.doubleSpinBox_forcechain_radius),
            (sb.label_forcechain_radius_factor, sb.doubleSpinBox_forcechain_factor),
        ]
        for show, ws in zip(show, wlist):
            for w in ws:
                w.setVisible(show)

    def change_shader(self):
        sb = self.sidebar
        shader = sb.comboBox_particle_shader.currentText()
        obj_data = self.current_object_data()
        mapper = obj_data['mapper']
        mapper.SetSplatShaderCode(SPLAT_SHADERS.get(shader))
        self.render()

    def change_glyph(self):
        sb = self.sidebar
        obj_data = self.current_object_data()
        mapper = obj_data['mapper']
        glyph_name = sb.comboBox_glyph.currentText()

        pvd = obj_data.get('pvd')
        array_info = pvd.get_arrayinfo()

        # Update glyph
        self.update_glyph(mapper, glyph_name, array_info, pvd)

        obj_data['glyph_name'] = glyph_name
        self.render()

    def update_sq_glyph_params(self):
        for k, obj_data in self.pipeline.items():
            if obj_data.get('mtext', '') == 'glyphs':
                mapper = obj_data.get('mapper', None)
                glyph_name = obj_data.get('glyph_name', 'sphere')
                pvd = obj_data.get('pvd', None)
                if any(v is None for v in [mapper, pvd]):
                    continue
                array_info = pvd.get_arrayinfo()
                self.update_glyph(mapper, glyph_name, array_info, pvd)

    def change_forcechain_radius(self):
        sb = self.sidebar
        obj_data = self.current_object_data()
        v = sb.doubleSpinBox_forcechain_radius.value()

        tb = obj_data.get('tube_filter')
        if tb is None:
            return
        tb.SetRadius(v)
        self.render()

    def change_forcechain_factor(self):
        sb = self.sidebar
        obj_data = self.current_object_data()
        v = sb.doubleSpinBox_forcechain_factor.value()

        tb = obj_data.get('tube_filter')
        if tb is None:
            return
        tb.SetRadiusFactor(v)
        self.render()

    def change_max_points(self):
        sb = self.sidebar
        obj_data = self.current_object_data()
        if not obj_data:
            return
        v = safe_int(sb.lineedit_maximum_glyph_points.text(), DEFAULT_MAXIMUM_POINTS)

        mask = obj_data.get('point_mask')
        if mask is None:
            return

        mask.SetMaximumNumberOfPoints(v)
        self.render()

    def make_PointGaussianMapper(self, reader, array_info, shader='sphere'):
        pg = vtk.vtkPointGaussianMapper()
        pg.EmissiveOff()

        # scale with particle diameter info
        diameter = self.gui.project.get_value('d_p0', args=[1], default=0.5)
        if 'Diameter' in array_info:
            pg.SetScaleArray('Diameter')
            pg.SetScaleFactor(.5)
        else:
            pg.SetScaleFactor(diameter/2.0)

        pg.SetScalarVisibility(True)
        pg.SetScalarModeToUsePointFieldData()
        pg.SetInputConnection(reader.GetOutputPort())
        pg.SetSplatShaderCode(SPLAT_SHADERS.get(shader))
        return pg

    def make_GlyphMapper(self, reader, array_info, pvd, glyph_name='sphere'):
        glyph = vtk.vtkOpenGLGlyph3DMapper()
        glyph.SetScalarVisibility(True)
        glyph.SetScalarModeToUsePointFieldData()
        glyph.SetInputConnection(reader.GetOutputPort())
        if 'QUATERNION' in array_info:
            glyph.SetOrientationArray('QUATERNION')
            glyph.SetOrientationModeToQuaternion()
        elif 'Orientation' in array_info:
            glyph.SetOrientationArray('Orientation')
            glyph.SetOrientationModeToDirection()
        elif 'Velocity' in array_info:
            glyph.SetOrientationArray('Velocity')
            glyph.SetOrientationModeToDirection()

        # Update glyph
        self.update_glyph(glyph, glyph_name, array_info, pvd)

        sb = self.sidebar
        safe_combo(sb.comboBox_glyph, glyph_name)
        return glyph

    def update_glyph(self, glyph, glyph_name, array_info, pvd):
        # scale with particle diameter info
        diameter = self.gui.project.get_value('d_p0', args=[1], default=0.5)
        # superquadric
        if glyph_name == 'superquadric':
            if 'SQP_POLY_SCALE' in array_info:
                glyph.SetScaleFactor(1)
                glyph.SetScaleModeToScaleByMagnitude()
                glyph.SetScaleArray('SQP_POLY_SCALE')
                glyph.ScalingOn()
            else:
                glyph.ScalingOff()
        else:
            if 'Diameter' in array_info:
                glyph.SetScaleFactor(1)
                glyph.SetScaleModeToScaleByMagnitude()
                glyph.SetScaleArray('Diameter')
            else:
                glyph.SetScaleFactor(diameter)
            glyph.ScalingOn()

        if glyph_name == 'superquadric':
            if 'Particle_Phase_ID' in array_info:
                glyph.SetSourceIndexArray('Particle_Phase_ID')
                dummy = PRIMITIVE_DICT['superquadric']()
                dummy.Update()
                glyph.SetSourceConnection(0, dummy.GetOutputPort())
                idx_range = array_info['Particle_Phase_ID'].range_[0]
                for i in range(safe_int(idx_range[0],1), safe_int(idx_range[1],1)+1):
                    gs = PRIMITIVE_DICT.get(glyph_name, PRIMITIVE_DICT['superquadric'])()
                    gs = self.update_sq_glyph(gs, array_info, pvd, idx=i)
                    gs.Update()
                    glyph.SetSourceConnection(i, gs.GetOutputPort())
                glyph.SourceIndexingOn()
            else:
                gs = PRIMITIVE_DICT.get(glyph_name, PRIMITIVE_DICT['superquadric'])()
                gs = self.update_sq_glyph(gs, array_info, pvd)
                glyph.SetSourceConnection(gs.GetOutputPort())
                glyph.SourceIndexingOff()
        else:
            gs = PRIMITIVE_DICT.get(glyph_name, PRIMITIVE_DICT['sphere'])()
            gs.Update()
            glyph.SetSourceConnection(gs.GetOutputPort())
            glyph.SourceIndexingOff()

    def update_sq_glyph(self, g, array_info, pvd, idx=None):
        field_data = array_info.get('field_data')
        if idx is None:
            idx = 1
            # infer phase from file name
            fname = pvd.pvdfilename.stem
            if 'PHASE' in fname:
                idx = safe_int(fname.split('PHASE')[-1].split('.')[0], 1) # seems hacky

        a, b, c, m, n = field_data.get(
            idx, [0.3E-02, 0.3E-02, 0.3E-02, 0.4E+01, 0.4E+01])

        g.SetCenter(0.0, 0.0, 0.0)
        g.SetAxisOfSymmetry(2)
        g.SetSize(1)
        g.SetThetaResolution(10)
        g.SetPhiResolution(10)
        g.SetScale(a,b,c)
        g.SetPhiRoundness(2/n)      # 2/n
        g.SetThetaRoundness(2/m)    # 2/m
        return g

    def make_OpenGLSphereMapper(self, reader):
        pg = vtk.vtkOpenGLSphereMapper()
        pg.SetScaleArray('Diameter')
        pg.SetScalarVisibility(True)
        pg.SetScalarModeToUsePointFieldData()
        pg.SetInputConnection(reader.GetOutputPort())
        return pg

    def make_OpenGLStickMapper(self, reader):
        pg = vtk.vtkOpenGLStickMapper()
        # pg.SetScaleArray('FORCE_CHAIN_LENGTH')
        # pg.SetScalarVisibility(True)
        # pg.SetScalarModeToUsePointFieldData()
        pg.SetInputConnection(reader.GetOutputPort())
        return pg

    def make_force_chain_mapper(self, reader):

        gm = vtk.vtkGlyph3D()
        gm.SetInputConnection(reader.GetOutputPort())
        gm.SetInputArrayToProcess(
            0, 0, 0, 0,
            'FORCE_CHAIN_LENGTH')
        gm.SetInputArrayToProcess(
            1, 0, 0, 0,
            'FORCE_CHAIN_ORIENTATION')
        gm.ScalingOn()
        gm.SetScaleFactor(1)
        gm.SetScaleModeToScaleByScalar()
        gm.OrientOn()
        gm.SetVectorModeToUseVector()

        # line source
        ls = vtk.vtkLineSource()
        ls.SetResolution(1)
        ls.Update()
        gm.SetSourceConnection(ls.GetOutputPort())

        # post tube
        tb = vtk.vtkTubeFilter()
        tb.SetInputConnection(gm.GetOutputPort())
        tb.SetRadius(0.0005)
        tb.SetRadiusFactor(5)
        tb.SetNumberOfSides(3)
        tb.CappingOff()
        tb.SetVaryRadiusToVaryRadiusByScalar()
        tb.SetInputArrayToProcess(
            vtk.vtkDataSetAttributes.SCALARS, 0, 0,
            vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,
            'FORCE_CHAIN_FN')

        mp = vtk.vtkPolyDataMapper()
        mp.SetInputConnection(tb.GetOutputPort())
        # mp.SelectColorArray('FORCE_CHAIN_FN_LENGTH')
        return mp, tb


    def build_stl_pipeline(self, obj_data):
        reader = vtk.vtkSTLReader()
        filename = obj_data.get('obj_name', '')
        reader.SetFileName(filename)
        reader.MergingOn()

        mapper = obj_data['mapper'] = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(reader.GetOutputPort())

        actor = obj_data['actor'] = self.make_actor()
        actor.SetMapper(mapper)

        obj_data.setdefault('opacity', 0.4)
        obj_data.setdefault('color', copy.deepcopy(DEFAULT_GEO_COLOR))
        rep = obj_data.setdefault('style', 'solid')

        self.vtkrenderer.AddActor(actor)

        self.change_opacity(obj_data)
        self.set_representation(actor, rep)
        self.change_geometry_color(obj_data)


    def object_clicked(self, item):
        """tree object was clicked"""
        tw = self.sidebar.treeWidget_pipeline
        if not item:
            return
        obj_name = item.text(0)
        # if obj_name not in self.active_objects:
        #     return
        tw.setCurrentItem(item)
        self.gui.app.processEvents()
        obj_data = self.pipeline.get(obj_name)
        actor = obj_data.get('actor', None)
        visible = item.checkState(0) == QtCore.Qt.Checked
        obj_data['visible'] = visible
        if visible and actor is None:
            self.build_pipeline(obj_data)
        self.change_visibility(actor, visible)

    def current_object(self):
        item = self.sidebar.treeWidget_pipeline.currentItem()
        if item:
            return item.text(0)
        else:
            return None

    def current_object_data(self):
        # get the current object name/data
        return self.pipeline.get(self.current_object())

    def object_selection_changed(self):
        """New tree object selected"""
        obj_name = self.current_object()
        self.update_object_widgets(obj_name)

    def update_object_widgets(self, obj_name):
        sb = self.sidebar
        sw = sb.stackedWidget_object_options
        # get the current object name/data
        obj_data = self.pipeline.get(obj_name)
        if obj_data is None:
            return

        obj_type = obj_data.get('type', '')

        # show the correct widgets based on type
        new_index = 0
        for i in range(sw.count()):
            widget = sw.widget(i)
            if str(widget.objectName()).endswith(obj_type):
                new_index = i
                break
        sw.setCurrentIndex(new_index)

        # populate widget data
        if obj_type == 'stl':
            color = obj_data.get('color', DEFAULT_GEO_COLOR)
            btn = sb.toolButton_geometry_color
            btn.setStyleSheet("QToolButton{{ background: {};}}".format(color.name()))
            safe_combo(sb.comboBox_geometry_style, obj_data.get('style', 'solid'))
            sb.doubleSpinBox_geometry_opacity.setValue(obj_data.get('opacity', 0.4))
        else:
            pvdfile = self.pvd_reader.pvd(obj_data.get('obj_name'))
            if pvdfile is None:
                return
            array_info = pvdfile.get_arrayinfo()

            # variable
            vars = list(array_info.keys())
            if vars:
                var = obj_data.setdefault('variable', vars[0])
                var = update_combo(
                    sb.comboBox_points_variable if obj_type == 'vtp' else sb.comboBox_cells_variable,
                    vars, var)
                var_data = obj_data['variable_color'].setdefault(var, copy.deepcopy(DEFAULT_VARIABLE_DATA))

                self.update_variable_widgets(var_data, obj_type, array_info, var)
                obj_data['variable'] = var

                if obj_type == 'vtp' and hasattr(array_info[var], 'number_of_tuples'):
                    sb.lineedit_point_count.setText('{}'.format(array_info[var].number_of_tuples))

            # opacity
            spin = sb.doubleSpinBox_points_opacity if obj_type == 'vtp' else sb.doubleSpinBox_cells_opacity
            spin.setValue(obj_data.get('opacity', 1.0))

            if obj_type == 'vtu':
                safe_combo(sb.comboBox_cells_style, obj_data.get('cell_style', 'cells'))

            # particle mapper
            mtext = obj_data.get('mtext', 'point gaussian')
            safe_combo(sb.comboBox_particle_mapper, mtext)
            self.show_particle_mapper_options(mtext)


    def variable_selection_changed(self):
        sb = self.sidebar

        # get the current object name/data
        unique_obj_name = self.current_object()
        obj_data = self.current_object_data()
        obj_type = obj_data.get('type', '')

        pvdfile = self.pvd_reader.pvd(obj_data.get('obj_name'))
        if pvdfile is None:
            return
        array_info = pvdfile.get_arrayinfo()

        cb = sb.comboBox_points_variable if obj_type == 'vtp' else sb.comboBox_cells_variable
        var = obj_data['variable'] = cb.currentText()
        var_data = obj_data['variable_color'].setdefault(var, copy.deepcopy(DEFAULT_VARIABLE_DATA))
        self.update_variable_widgets(var_data, obj_type, array_info, var)

        self.change_color(unique_obj_name, var_data)

    def update_variable_widgets(self, var_data, obj_type, array_info, var):
        sb = self.sidebar

        # component
        comp = var_data.get('component', 'mag')
        comp_cb = sb.comboBox_points_component if obj_type == 'vtp' else sb.comboBox_cells_component
        if hasattr(array_info[var], 'components'):
            num_comp = array_info[var].components
            comp_cb.setEnabled(num_comp > 1)
        safe_combo(comp_cb, comp)

    def component_selection_changed(self):
        sb = self.sidebar

        # get the current object name/data
        unique_obj_name = self.current_object()
        obj_data = self.current_object_data()
        obj_type = obj_data.get('type', '')

        var = obj_data.get('variable')
        if var is None:
            return

        comp_cb = sb.comboBox_points_component if obj_type == 'vtp' else sb.comboBox_cells_component
        comp = comp_cb.currentText()
        obj_data.get('variable_color', {}).get(var, {})['component'] = comp

        var_data = obj_data.get('variable_color', {}).get(var, {})
        self.change_color(unique_obj_name, var_data)

        # update the component in the other onjects that use the same colormap
        for name in self.color_maps.get(obj_type, {}).get(var, {}).get('obj_list', []):
            if name == unique_obj_name:
                continue
            obj_data = self.pipeline.get(name)
            obj_data.get('variable_color', {}).get(var, {})['component'] = comp

    def change_color_pressed(self):
        sb = self.sidebar

        # get the current colorbar name
        sb = self.sidebar
        item = sb.treeWidget_colorbars.currentItem()
        if not item:
            return
        var = item.text(0)
        type_ = item.text(1)
        cbar_data = self.color_maps.get(type_, {}).get(var)

        # loop over objects, collecting range info
        cbar_var_data = {}
        # ParticleArray(
        #     number_of_tuples=2400,
        #     components=3,
        #     range_=[(-1.0186058282852173, 0.8821541666984558),
        #             (-0.5045252442359924, 2.5463428497314453),
        #             (0.0, 0.0)],
        #     magnitude=(0.002775479783295483, 2.5483403396637114))

        for obj_name in cbar_data.get('obj_list', []):
            obj_data = self.pipeline.get(obj_name)
            pvd = obj_data.get('pvd')
            array_info = pvd.get_arrayinfo()
            var_data = array_info.get(var, {})
            if not cbar_var_data:
                cbar_var_data = {
                    'range': var_data.range_ if hasattr(var_data, 'range_') else [(0, 1)],
                    'magnitude': var_data.magnitude if hasattr(var_data, 'magnitude') else [0, 1],
                }
            else:
                cbar_var_data['range'] = minmax(cbar_var_data['range'], var_data.range_)

        # add in the color info
        for k in ['color', 'single_color', 'color_map', 'reversed', 'from', 'to', 'component']:
            if cbar_data.get(k) is not None:
                cbar_var_data[k] = cbar_data.get(k)

        comp = cbar_var_data.get('component', 'mag')
        comp_num = INDEX_DICT[comp]
        self.color_dialog.set_((type_, var), cbar_var_data, comp_num)
        self.color_dialog.popup(var)

    def change_variable_color(self, type_var, new_var_data):
        cbar_data = self.color_maps.get(type_var[0], {}).get(type_var[1])
        cbar_data['from'] = new_var_data.get('from')
        cbar_data['to'] = new_var_data.get('to')

        self.make_lut(new_var_data.get('color_map'), cbar_data)

        for unique_obj_name in copy.deepcopy(cbar_data.get('obj_list', [])):
            self.change_color(unique_obj_name, new_var_data)

        # save info
        for k in ['color', 'single_color', 'color_map', 'reversed', 'from', 'to']:
            cbar_data[k] = new_var_data.get(k)

        # color the button
        btn = self.sidebar.toolButton_color_scale
        self.set_color_button(btn, new_var_data)

    def make_lut(self, color_map, cmap_data):
        new_lut = vtk.vtkLookupTable()
        new_lut.DeepCopy(LOOKUP_TABLES.get(color_map, LOOKUP_TABLES[DEFAULT_COLOR_MAP]))
        # Fix for bug introduced in VTK7, bug fixed in VTK8
        # https://gitlab.kitware.com/vtk/vtk/issues/16966
        #if hasattr(new_lut, 'BuildSpecialColors'):
        #    new_lut.BuildSpecialColors()
        #new_lut.Build()

        cmap_data['lut'] = lut = new_lut
        cmap_data['cmap'] = color_map

        # update the color bar with the new lut
        cbar_actor = cmap_data.get('actor', None)
        if cbar_actor is not None:
            cbar_actor.SetLookupTable(lut)
        self.update_colorbar_tree()
        return lut


    def change_color(self, unique_obj_name, var_data):
        sb = self.sidebar
        obj_data = self.pipeline.get(unique_obj_name)
        if obj_data is None or obj_data.get('pvd') is None:
            return
        obj_type = obj_data.get('type', '')
        var = obj_data.get('variable')

        # check the other color map object lists and remove
        for vname, vmap in self.color_maps.get(obj_type, {}).items():
            if vname == var:
                continue
            obj_list = vmap.get('obj_list', [])
            if unique_obj_name in obj_list: obj_list.remove(unique_obj_name)

        # add to the correct color map object list
        cmap_data = self.color_maps.get(obj_type, {}).setdefault(var, {})
        obj_list = cmap_data.setdefault('obj_list', [])
        if unique_obj_name not in obj_list:
            obj_list.append(unique_obj_name)

        if var_data.get('from', None) is None:
            pvdfile = obj_data.get('pvd')
            array_info = pvdfile.get_arrayinfo()
            var_info = array_info.get(var, None)
            if var_info is None:
                return
            mag = var_data['magnitude'] = var_info.magnitude if hasattr(var_info, 'magnitude') else [0, 1]
            var_data['from'] = min(mag)
            var_data['to'] = max(mag)
            var_data['range'] = var_info.range_ if hasattr(var_info, 'range_') else [(0, 1)]

        if obj_type == 'vtu':
            # array filter
            pass_array = obj_data.get('entry')
            if pass_array is not None:
                pass_array.ClearArrays()
                pass_array.AddCellDataArray(var)

        # build color map
        single_color = var_data.get('single_color', False)

        # get existing colormap
        lut = cmap_data.get('lut', None)
        if lut is None:
            color_map = var_data.get('color_map', DEFAULT_COLOR_MAP)
            lut = self.make_lut(color_map, cmap_data)

        # set component
        cmap_data['component'] = comp = var_data.get('component', 'mag')
        if comp != 'mag':
            cmap_data['comp'] = f'({comp})'
        else:
            cmap_data['comp'] = f''
        cbar = cmap_data.get('actor')
        if cbar is not None:
            self.update_colorbar(cbar, **cmap_data)
        comp_num = INDEX_DICT[comp]
        if comp_num >= 0:
            lut.SetVectorModeToComponent()
            lut.SetVectorComponent(comp_num)
        else:
            lut.SetVectorModeToMagnitude()

        # set the mappers
        vmin = cmap_data.get('from', None)
        if vmin is None:
            vmin = safe_float(var_data.get('from', 0.0), 0.0)
        vmax = cmap_data.get('to', None)
        if vmax is None:
            vmax = safe_float(var_data.get('to', 1.0), 1.0)
        for k in ['mapper', 'mapper_cells', 'mapper_nodes']:
            mapper = obj_data.get(k)
            if mapper is not None:
                mlist = self.color_maps[obj_type][var].setdefault('mappers', [])
                if mapper not in mlist:
                    mlist.append(mapper)
                mapper.SelectColorArray(var)
                mapper.SetScalarRange(vmin, vmax)
                if single_color:
                    mapper.ScalarVisibilityOff()
                else:
                    mapper.ScalarVisibilityOn()
                    mapper.SetLookupTable(lut)

        # set the actors
        color = var_data.get('color', DEFAULT_GEO_COLOR)
        for k in ['actor', 'actor_cells', 'actor_nodes']:
            actor = obj_data.get(k)
            if actor is not None:
                if single_color:
                    # actor.GetProperty().SetDiffuseColor(color.getRgbF()[:3])
                    # actor.GetProperty().SetSpecularColor(color.getRgbF()[:3])
                    actor.GetProperty().SetColor(color.getRgbF()[:3])
                # self.set_actor_shiny_props(actor)

        self.render()

    def set_color_button(self, btn, var_data):
        single_color = var_data.get('single_color', False)

        if single_color:
            color = var_data.get('color', QtCore.Qt.white)
            btn.setStyleSheet("QToolButton{{ background: {};}}".format(color.name()))
            btn.setIcon(QtGui.QIcon())
        else:
            color_map = var_data.get('color_map', DEFAULT_COLOR_MAP)
            btn.setIcon(build_qicons().get(color_map, {}).get('icon', QtGui.QIcon()))
            btn.setStyleSheet("QToolButton{{ background: {};}}".format(None))

    def geometry_color_pressed(self):
        obj_data = self.current_object_data()

        color = obj_data.get('color', DEFAULT_GEO_COLOR)
        color = QtWidgets.QColorDialog.getColor(color)
        if not color.isValid():
            return

        obj_data['color'] = color

        btn = self.sidebar.toolButton_geometry_color
        btn.setStyleSheet("QToolButton{{ background: {};}}".format(color.name()))

        self.change_geometry_color(obj_data)

    def change_geometry_color(self, obj_data):
        color = obj_data.get('color', DEFAULT_GEO_COLOR)
        actor = obj_data.get('actor', None)
        if actor is None:
            return
        actor.GetProperty().SetColor(color.getRgbF()[:3])
        self.render()

    def opacity_changed(self):
        sb = self.sidebar

        # get the current object name/data
        obj_data = self.current_object_data()
        obj_type = obj_data.get('type', '')

        if obj_type == 'vtp':
            spin = sb.doubleSpinBox_points_opacity
        elif obj_type == 'vtu':
            spin = sb.doubleSpinBox_cells_opacity
        elif obj_type == 'stl':
            spin = sb.doubleSpinBox_geometry_opacity
        new_op = spin.value()
        if new_op == obj_data['opacity']:
            return

        obj_data['opacity'] = new_op
        self.change_opacity(obj_data)

    def cell_style_changed(self):
        sb = self.sidebar

        # get the current object name/data
        obj_data = self.current_object_data()
        obj_data['cell_style'] = cell_style = sb.comboBox_cells_style.currentText()

        cells = cell_style == 'cells'
        vis = obj_data.get('visible', True)
        self.change_visibility(obj_data.get('actor_cells'), cells and vis)
        self.change_visibility(obj_data.get('actor_nodes'), not cells and vis)
        obj_data['actor'] = obj_data.get('actor_cells') if cells else obj_data.get('actor_nodes')
        self.render()

    def geometry_style_changed(self):
        sb = self.sidebar
        rep = sb.comboBox_geometry_style.currentText()
        obj_data = self.current_object_data()
        obj_data['style'] = rep
        actor = obj_data.get('actor', None)
        if actor is None:
            return
        self.set_representation(actor, rep)
        self.render()

    # --- vtk functions ---
    def change_visibility(self, actor, visible):
        if actor is None:
            return
        actor.SetVisibility(visible)
        self.render()

    def change_opacity(self, obj_data):
        opacity = obj_data.get('opacity', 1.0)
        actor = obj_data.get('actor', None)
        if actor is None:
            return
        actor.GetProperty().SetOpacity(opacity)
        self.render()

    # --- Color Bar ----
    def colorbar_selection_changed(self):
        sb = self.sidebar
        item = sb.treeWidget_colorbars.currentItem()
        if not item:
            return
        var = item.text(0)
        type_ = item.text(1)
        cbar_data = self.color_maps.get(type_, {}).get(var)
        if cbar_data is None:
            return
        cbar_data = DEFAULT_COLOR_BAR | cbar_data # Merge defaults
        sb = self.sidebar
        lbl = cbar_data.get('label')
        if lbl is not None:
            sb.lineedit_color_bar_title.setText(lbl)
        sb.spinBox_color_bar_title_size.setValue(cbar_data.get('title_size', 12))
        sb.spinBox_color_bar_n_labels.setValue(cbar_data.get('n_labels', 11))
        sb.lineedit_color_bar_format.setText(cbar_data.get('label_fmt', '%.2f'))
        sb.spinBox_color_bar_label_size.setValue(cbar_data.get('label_size', 12))
        safe_combo(sb.comboBox_color_bar_position, cbar_data.get('position', 'right'))
        # color the button
        btn = self.sidebar.toolButton_color_scale
        self.set_color_button(btn, cbar_data)
        sb.doubleSpinBox_color_bar_opacity.setValue(cbar_data['opacity'])
        text_color = cbar_data.get('color', DEFAULT_TEXT_COLOR)
        if text_color is None:
            text_color = DEFAULT_TEXT_COLOR
        sb.toolButton_color_bar_color.setStyleSheet(
            "QToolButton{{ background: {};}}".format(text_color.name()))

    def colorbar_value_changed(self, key, value):
        sb = self.sidebar
        item = sb.treeWidget_colorbars.currentItem()
        if not item:
            return
        var = item.text(0)
        type_ = item.text(1)
        cbar_data = self.color_maps.get(type_, {}).get(var)
        cbar_data[key] = value
        cbar = cbar_data.get('actor')
        if cbar is not None:
            self.update_colorbar(cbar, **cbar_data)

    def colorbar_clicked(self, item):
        if not item:
            return
        var = item.text(0)
        type_ = item.text(1)
        tw = self.sidebar.treeWidget_colorbars
        tw.setCurrentItem(item)
        cbar_data = self.color_maps.get(type_, {}).get(var)
        visible = item.checkState(0) == QtCore.Qt.Checked
        cbar_data['visible'] = visible

        if visible and cbar_data.get('actor') is None:
            lut = cbar_data.get('lut')
            if lut:
                cbar_data['actor'] = self.make_colorbar(lut)
                for k, v in DEFAULT_COLOR_BAR.items():
                    if k not in cbar_data:
                        cbar_data[k] = copy.deepcopy(v)
                cbar_data['label'] = cbar_data.get('label', var)

        actor = cbar_data.get('actor', None)
        if actor is not None:
            actor.SetVisibility(visible)
            self.update_colorbar(actor, **cbar_data)

        self.render()

    def make_colorbar(self, lut):
        cbar_actor = vtk.vtkScalarBarActor()
        cbar_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        cbar_actor.AnnotationTextScalingOff()
        cbar_actor.SetLookupTable(lut)
        cbar_actor.UnconstrainedFontSizeOn()
        self.vtkrenderer.AddActor(cbar_actor)
        return cbar_actor

    def update_colorbar(self, cbar, position=None, label=None, comp=None, color=None,
                     shadow=None, italic=None, n_labels=None, label_fmt=None,
                     title_size=None, label_size=None,
                     **kwargs):

        if position is not None:
            geo = [0.08, 0.7]
            if position == 'left':
                pos = [0.02, 0.15]
                cbar.SetOrientationToVertical()
            elif position == 'bottom':
                pos = [0.15, 0.02]
                geo = list(reversed(geo))
                cbar.SetOrientationToHorizontal()
            elif position == 'top':
                pos = [0.15, 0.9]
                geo = list(reversed(geo))
                cbar.SetOrientationToHorizontal()
            else: # default to right
                pos = [0.9, 0.15]
                cbar.SetOrientationToVertical()

            cbar.GetPositionCoordinate().SetValue(*pos)
            cbar.SetWidth(geo[0])
            cbar.SetHeight(geo[1])

        # title
        if label is not None:
            if comp is not None:
                cbar.SetComponentTitle(comp)
            else:
                cbar.SetComponentTitle('')
            cbar.SetTitle(label)
            cbar.SetVerticalTitleSeparation(5)

        # title
        if title_size is not None:
            p_title = cbar.GetTitleTextProperty()
            p_title.SetFontSize(title_size)

        # number labels
        if label_size is not None:
            p_label = cbar.GetLabelTextProperty()
            p_label.SetFontSize(label_size)

        if n_labels is not None:
            cbar.SetNumberOfLabels(n_labels)

        if color is not None:
            if isinstance(color, QtGui.QColor):
                color = color.getRgbF()[:3]
            for label in [cbar.GetLabelTextProperty(), cbar.GetTitleTextProperty()]:
                label.SetColor(color)

        if shadow is not None:
            for label in [cbar.GetLabelTextProperty(), cbar.GetTitleTextProperty()]:
                label.SetShadow(shadow)

        if italic is not None:
            for label in [cbar.GetLabelTextProperty(), cbar.GetTitleTextProperty()]:
                label.SetItalic(italic)

        if label_fmt is not None:
            cbar.SetLabelFormat(label_fmt)

        self.render()

    def change_color_bar_title(self):
        title = self.sidebar.lineedit_color_bar_title.text()
        self.colorbar_value_changed('label', title)

    def change_color_bar_title_size(self):
        size = self.sidebar.spinBox_color_bar_title_size.value()
        self.colorbar_value_changed('title_size', size)

    def change_color_bar_label_size(self):
        size = self.sidebar.spinBox_color_bar_label_size.value()
        self.colorbar_value_changed('label_size', size)

    def change_color_bar_n_labels(self):
        n_labels = self.sidebar.spinBox_color_bar_n_labels.value()
        self.colorbar_value_changed('n_labels', n_labels)

    def change_color_bar_format(self):
        le = self.sidebar.lineedit_color_bar_format
        fmt = le.text()

        try:
            fmt % 3.14
            color = 'black'
            self.colorbar_value_changed('label_fmt', fmt)
        except:
            color = 'red'
        le.setStyleSheet("color: " + color)

    def change_color_bar_loc(self):
        pos = self.sidebar.comboBox_color_bar_position.currentText()
        self.colorbar_value_changed('position', pos)

    def change_color_bar_color(self, checked=False, color=None):
        if color is None:
            sb = self.sidebar
            item = sb.treeWidget_colorbars.currentItem()
            if not item:
                return
            var = item.text(0)
            type_ = item.text(1)
            cbar_data = self.color_maps.get(type_, {}).get(var)
            cur_color = cbar_data.get('color',  DEFAULT_TEXT_COLOR)
            if cur_color is None:
                cur_color = DEFAULT_TEXT_COLOR
            color = QtWidgets.QColorDialog.getColor(cur_color, parent=self)
            if not color.isValid():
                return

        self.sidebar.toolButton_color_bar_color.setStyleSheet(
            "QToolButton{{ background: {};}}".format(color.name()))
        self.colorbar_value_changed('color', color)

    def change_color_bar_opacity(self):
        sb = self.sidebar
        item = sb.treeWidget_colorbars.currentItem()
        if not item:
            return
        var = item.text(0)
        type_ = item.text(1)
        cbar_data = self.color_maps.get(type_, {}).get(var)
        op = cbar_data['opacity'] = sb.doubleSpinBox_color_bar_opacity.value()
        self.change_opacity(cbar_data)

    def update_colorbar_tree(self):
        tw = self.sidebar.treeWidget_colorbars
        last_var = tw.currentItem()
        if last_var is not None:
            last_var = last_var.text(0)

        tw.clear()
        toplevel = None
        for type_, vars in self.color_maps.items():
            for var in vars.keys():
                # add tree widget item
                toplevel = QtWidgets.QTreeWidgetItem([var, type_])
                toplevel.setFlags(toplevel.flags() | QtCore.Qt.ItemIsUserCheckable)
                icon = 'dem.svg' if type_ == 'vtp' else 'grid.svg'
                toplevel.setIcon(0, get_icon(icon))
                vis = self.color_maps.get(type_, {}).get(var, {}).get('visible', False)
                checked = QtCore.Qt.Checked if vis else QtCore.Qt.Unchecked
                toplevel.setCheckState(0, checked)

                tw.addTopLevelItem(toplevel)

                if var == last_var:
                    tw.setCurrentItem(toplevel)

        # make sure something is selected
        if not self.loading and last_var is None and toplevel is not None:
            tw.setCurrentItem(toplevel)

        # self.colorbar_selection_changed()

    # --- time label ---
    def handle_time_label_format(self, text):
        try:
            text % 1.34
            self.time_format = text
            self.set_timelabel(text=self.time_format % self.time)
            self.render()
            color = 'black'
        except:
            color = 'red'
        self.sidebar.lineedit_time_label_format.setStyleSheet("color: " + color)

    def change_time_label_color(self, checked=False, color=None):
        """Change the color of the geometry actor"""
        if color is None:
            color = QtWidgets.QColorDialog.getColor(parent=self)
            if not color.isValid():
                return
        self.time_label_color = color

        self.sidebar.toolButton_time_label_color.setStyleSheet(
            "QToolButton{{ background: {};}}".format(color.name()))

        self.set_timelabel(color=color.getRgbF()[:3])
        self.render()

    def change_time_label_loc(self):
        cb = self.sidebar.comboBox_time_label_loc.currentText()
        self.set_timelabel(pos=cb)

        self.render()

    def change_time_label_text_size(self):
        size = self.sidebar.spinBox_time_label_size.value()
        self.set_timelabel(fontsize=size)

        self.render()
