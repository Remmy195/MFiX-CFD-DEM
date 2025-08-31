"""
Widget to show histograms.
"""
import logging
import sys

from qtpy.QtWidgets import (
    QVBoxLayout,
    QGridLayout,
    QWidget,
    QToolButton,
    QLabel,
    QComboBox,
    QSizePolicy,
    QSpacerItem,
    QSpinBox,
    QAbstractSpinBox,
)
from qtpy.QtCore import QTimer, Qt

from mfixgui.tools.qt import get_icon, sub_icon_size
from mfixgui.tools.util import safe_int
from mfixgui.vtk_widgets.tools import update_combo, safe_combo
from mfixgui.vtk_widgets.pvd_reader import get_pvd_reader

# graphics libraries
from mfixgui.tools.pyqtgraph import pg, PlotWidget, DEFAULT_PENS, DEFAULT_BRUSHS

LOG = logging.getLogger(__name__)


class HistogramViewer(QWidget):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)

        self.name = 'histogram'

        self._pvd_reader = get_pvd_reader()
        self._pvd_reader.signal_new_files.connect(self.update_pvd_combo)
        self._hgram = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(5, 5, 5, 5)

        # controls
        btn_bar = QWidget()
        self.button_bar_layout = QGridLayout(btn_bar)
        layout.addWidget(btn_bar)
        self.build_toolbar()

        # plot
        self.histogram_plot = PlotWidget()
        layout.addWidget(self.histogram_plot)

        # play timer
        self.play_timer = QTimer()
        self.play_timer.timeout.connect(self.handle_next)

        # check for files on init
        self.update_pvd_combo(self._pvd_reader.pvd_filenames())

    @property
    def frame_index(self):
        return self.frame_spinbox.value()

    @frame_index.setter
    def frame_index(self, index):
        index = max(0, min(index, self.frame_index_max()))
        self.frame_spinbox.setValue(index)
        self.frame_refresh()

    def frame_index_max(self):
        pvdfile = self.current_pvd()
        if pvdfile is None:
            LOG.warning("PVD not initialized; defaulting to maximum int")
            return sys.maxsize

        return pvdfile.max_index()

    def current_pvd(self):
        return self._pvd_reader.pvd(self.combox_pvd.currentText())

    def get_state(self):
        state = {
            'pvd': self.combox_pvd.currentText(),
            'var': self.combox_variable.currentText(),
            'comp': self.combobox_component.currentText(),
            'bins': self.spinbox_bins.value(),
            'style': self.combobox_style.currentText(),
        }

        return state

    def set_state(self, state):
        for combo, key, d in [
            (self.combox_pvd, 'pvd', ''),
            (self.combox_variable, 'var', ''),
            (self.combobox_component, 'comp', 'x'),
            (self.combobox_style, 'style', 'bar'),
        ]:
            safe_combo(combo, state.get(key, d))

        self.spinbox_bins.setValue(safe_int(state.get('bins', 100), 100))

    def build_toolbar(self):
        layout = self.button_bar_layout
        layout.setContentsMargins(5, 5, 5, 5)
        hspacer = QSpacerItem(99999, 10, QSizePolicy.Expanding, QSizePolicy.Maximum)

        # PVD select
        lbl = QLabel('File pattern')
        layout.addWidget(lbl, 0, 0)
        self.combox_pvd = QComboBox()
        self.combox_pvd.activated.connect(self.frame_refresh)
        layout.addWidget(self.combox_pvd, 0, 1)

        lbl = QLabel('Variable')
        layout.addWidget(lbl, 0, 2)
        self.combox_variable = QComboBox()
        self.combox_variable.activated.connect(self.frame_refresh)
        layout.addWidget(self.combox_variable, 0, 3)

        lbl = QLabel('Component')
        layout.addWidget(lbl, 0, 4)
        self.combobox_component = QComboBox()
        self.combobox_component.addItems(['x', 'y', 'z', 'mag'])
        self.combobox_component.activated.connect(self.frame_refresh)
        layout.addWidget(self.combobox_component, 0, 5)

        self.button_bar_layout.addItem(hspacer, 0, 6)

        lbl = QLabel('Bins')
        layout.addWidget(lbl, 1, 0)
        self.spinbox_bins = QSpinBox()
        self.spinbox_bins.setRange(2, 9999)
        self.spinbox_bins.setValue(10)
        self.spinbox_bins.valueChanged.connect(self.frame_refresh)
        layout.addWidget(self.spinbox_bins, 1, 1)

        lbl = QLabel('Style')
        layout.addWidget(lbl, 1, 2)
        self.combobox_style = QComboBox()
        self.combobox_style.addItems(['bar', 'line'])
        self.combobox_style.activated.connect(self.plot_refresh)
        layout.addWidget(self.combobox_style, 1, 3)

        # --- play/stop/next/back controls ---
        for i, (icon, callback, tooltip) in enumerate(
            [
                ('first.svg', self.handle_first, 'First'),
                ('back.svg', self.handle_back, 'Previous'),
                ('play.svg', self.handle_play_stop, 'Play'),
                ('next.svg', self.handle_next, 'Next'),
                ('last.svg', self.handle_last, 'Last'),
                ('autorenew.svg', None, 'Repeat from beginning'),
            ]
        ):
            btn = QToolButton()
            if tooltip == 'Play':
                self.toolbutton_play = btn
            elif tooltip.startswith('Repeat'):
                btn.setCheckable(True)
                self.toolbutton_repeat = btn
            if callback is not None:
                btn.clicked.connect(callback)
            btn.setIcon(get_icon(icon))
            btn.setIconSize(sub_icon_size())
            btn.setToolTip(tooltip)
            btn.setAutoRaise(True)
            btn.setFocusPolicy(Qt.ClickFocus)
            layout.addWidget(btn, 0, 7 + i)

        self.button_bar_layout.addItem(hspacer, 0, 15)

        self.frame_spinbox = QSpinBox()
        self.frame_spinbox.editingFinished.connect(self.frame_refresh)
        self.frame_spinbox.setMaximum(9999999)
        self.frame_spinbox.setButtonSymbols(QAbstractSpinBox.NoButtons)

        layout.addWidget(self.frame_spinbox, 0, 16)

    def handle_first(self):
        self.frame_index = 0

    def handle_last(self):
        self.frame_index = self.frame_index_max()

    def handle_play_stop(self):
        if self.play_timer.isActive():
            self.stop()
        else:
            self.toolbutton_play.setIcon(get_icon('stop.svg'))
            self.play_timer.start(100)

    def stop(self):
        self.toolbutton_play.setIcon(get_icon('play.svg'))
        self.play_timer.stop()

    def handle_next(self):
        index = self.frame_index + 1
        if (
            self.play_timer.isActive()
            and self.toolbutton_repeat.isChecked()
            and self.frame_index_max() > 0
        ):
            index = index % self.frame_index_max()
        self.frame_index = index

    def handle_back(self):
        self.frame_index -= 1

    def frame_refresh(self):
        pvdfile = self.current_pvd()
        if pvdfile is None:
            LOG.warning("PVD not initialized")
            return

        pvdfile.change_index(self.frame_index)
        array_info = pvdfile.get_arrayinfo()

        var = update_combo(self.combox_variable, array_info.keys())
        var_info = array_info.get(var)
        num_components = 1
        if var_info is not None:
            num_components = array_info[var].components
        self.combobox_component.setEnabled(num_components > 1)

        comp = self.combobox_component.currentText()
        self._hgram = pvdfile.histogram(var, comp, self.spinbox_bins.value())
        self.plot_refresh()

    def plot_refresh(self):
        hgram = self._hgram
        if hgram is None:
            return

        centers = (hgram.edges[1:] + hgram.edges[:-1]) / 2
        width = hgram.edges[1] - hgram.edges[0]

        style = self.combobox_style.currentText()

        pw = self.histogram_plot
        pw.clear()
        pw.showGrid(True, True, 0.5)
        if style == 'line':
            pw.plot(centers, hgram.counts, pen=DEFAULT_PENS[0])
        else:
            bg = pg.BarGraphItem(
                x=centers,
                width=width,
                height=hgram.counts,
                pen=pg.mkPen(color=[255, 255, 255], width=1),
                brush=DEFAULT_BRUSHS[0],
            )
            pw.addItem(bg)

        # pw.setXRange(min_, max_)
        pw.setLabel('bottom', self.combox_variable.currentText())
        pw.setLabel('left', 'Count')

    @property
    def project_dir(self):
        return self.gui.get_project_dir()

    def update_pvd_combo(self, pvd_filenames):
        first = bool(self.combox_pvd.currentIndex())
        var = update_combo(self.combox_pvd, pvd_filenames)
        if first:
            self.frame_index = 0

    def reset(self):
        self.histogram_plot.clear()
        self.play_timer.stop()
        self.frame_index = 0

    def close(self):
        # clean up timer
        self.play_timer.stop()

    def hideEvent(self, _event):
        self.stop()
