from qtpy import QtCore, QtGui, QtWidgets

from mfixgui.tools import safe_float, safe_int
from mfixgui.tools.qt import get_icon, get_image_path, sub_icon_size, get_ui, SETTINGS
from mfixgui.colormaps.color_maps import build_qicons
from mfixgui.vtk_widgets.tools import  safe_combo, parse_pvd_file, update_combo

DEFAULT_GEO_COLOR = QtGui.QColor(224, 224, 224)
DEFAULT_COLOR_MAP = 'viridis'


class ColorMapPopUp(QtWidgets.QDialog):
    applyEvent = QtCore.Signal(tuple, dict)

    def __init__(self, parent=None):
        QtWidgets.QDialog.__init__(self, parent)

        self.color = None
        self.org_range = [0, 1]
        ui = self.ui = get_ui('color_map.ui', self)

        self.setWindowTitle('Color Map')

        ui.toolbutton_select_color.clicked.connect(self.select_color)
        ui.toolbutton_select_color.setStyleSheet(
            "QToolButton{{ background: {};}}".format(DEFAULT_GEO_COLOR.name()))
        ui.lineedit_from.dtype = float
        ui.lineedit_to.dtype = float

        for name, icons in build_qicons().items():
            if not name.endswith('_reversed'):
                ui.combobox_color_map.addItem(
                    icons.get('bar', QtGui.QIcon()), name)
        self.set_color_map(DEFAULT_COLOR_MAP)

        btn = ui.buttonBox.button(QtWidgets.QDialogButtonBox.Apply)
        btn.clicked.connect(self.emit_apply_event)
        btn = ui.buttonBox.button(QtWidgets.QDialogButtonBox.Ok)
        btn.clicked.connect(self.emit_apply_event)

        ui.toolbutton_auto_scale.clicked.connect(self.auto_scale)

    def emit_apply_event(self):
        self.collect_values()
        self.applyEvent.emit(self.type_var, self.var_data)

    def set_(self, type_var, var_data, comp):
        self.var_data = var_data
        self.comp = comp
        self.type_var = type_var
        self.color = self.var_data.get('color', DEFAULT_GEO_COLOR)
        self.set_color(self.color)
        self.set_color_map(self.var_data.get('color_map', DEFAULT_COLOR_MAP))
        self.ui.checkbox_reversed.setChecked(self.var_data.get('reversed', False))

        d_range = self.var_data.get('range', [[0, 1]])
        if len(d_range) == 3:
            if comp >= 0:
                d_range = d_range[comp]
            else:
                d_range = self.var_data.get('magnitude', [0, 1])
        else:
            d_range = d_range[0]
        self.ui.lineedit_from.updateValue(
            None, self.var_data.get('from', '{:.3g}'.format(d_range[0])))
        self.ui.lineedit_to.updateValue(
            None, self.var_data.get('to', '{:.3g}'.format(d_range[1])))
        self.org_range = d_range
        single_color = self.var_data.get('single_color', False)
        self.ui.checkbox_single_color.setChecked(single_color)
        self.ui.widget_color_map.setEnabled(not single_color)

    def collect_values(self):
        color_map = self.ui.combobox_color_map.currentText()
        reverse = self.ui.checkbox_reversed.isChecked()

        f = self.ui.lineedit_from.value
        f = f if f is not None else self.org_range[0]

        t = self.ui.lineedit_to.value
        t = t if t is not None else self.org_range[1]
        rng = [safe_float(f, self.org_range[0]),
               safe_float(t, self.org_range[1])]

        if reverse:
            color_map += '_reversed'
        new_dict = {
            'color':        self.color,
            'single_color': self.ui.checkbox_single_color.isChecked(),
            'color_map':    color_map,
            'reversed':     reverse,
            'from':         min(rng),
            'to':           max(rng),
            }

        self.var_data.update(new_dict)

    def select_color(self):
        col = QtWidgets.QColorDialog.getColor(self.color)
        if not col.isValid():
            return
        self.color = col
        self.set_color(col)

    def set_color(self, color):
        if isinstance(color, QtGui.QColor):
            self.ui.toolbutton_select_color.setStyleSheet("QToolButton{{ background: {};}}".format(
                color.name()))

    def set_color_map(self, color_map):
        color_map = color_map.replace('_reversed', '')
        self.ui.combobox_color_map.setCurrentText(color_map)

    def auto_scale(self):
        comp = self.comp
        d_range = self.var_data.get('range', [[0, 1]])
        if len(d_range) == 3:
            if comp >= 0:
                d_range = d_range[comp]
            else:
                d_range = self.var_data.get('magnitude', [0, 1])
        else:
            d_range = d_range[0]
        self.ui.lineedit_from.updateValue(None, '{:.3g}'.format(d_range[0]))
        self.ui.lineedit_to.updateValue(None, '{:.3g}'.format(d_range[1]))

    def popup(self, text):

        self.setWindowTitle('Change color map for '+text)

        self.show()
        self.raise_()
        self.activateWindow()
