from qtpy import QtGui, QtCore
from .regexes import re_valid_run_name_qt
from .constants import *

## TODO implement limits for dt_min/dt_max

class Run:
    def init_run(self):
        ui = self.ui.run
        # Set input validator for run_name
        le = ui.lineedit_keyword_run_name
        le.setMaxLength(60)
        le.setValidator(QtGui.QRegExpValidator(
            QtCore.QRegExp(re_valid_run_name_qt)))
        le.required = True
        ui.radiobutton_steady_state.clicked.connect(self.set_run_steady_state)
        ui.radiobutton_steady_state.keys = ['dt_min', 'dt_max', 'dt_fac']
        ui.radiobutton_unsteady_state.clicked.connect(self.set_run_unsteady_state)
        ui.radiobutton_unsteady_state.keys = ['dt_min', 'dt_max', 'dt_fac']
        ui.radiobutton_constant_dt.clicked.connect(self.set_run_constant_dt)

        ui.radiobutton_constant_ratio.clicked.connect(self.set_run_constant_ratio)
        ui.radiobutton_constant_ratio.keys = ['dt_min', 'dt_max', 'dt_fac']
        ui.radiobutton_nice_dt.clicked.connect(self.set_run_nice_dt)
        self.add_tooltip(ui.radiobutton_nice_dt, 'nice_dt')

        ui.lineedit_keyword_dt.exclude_min = True
        ui.lineedit_keyword_dt.required = True
        ui.lineedit_keyword_dt.saved_value = 1e-6

        ui.lineedit_keyword_dt_fac.exclude_max = True

        ui.lineedit_keyword_dt_fac.required = True
        ui.lineedit_keyword_dt_fac.saved_value = 0.9
        ui.lineedit_keyword_dt_min.required = True
        ui.lineedit_keyword_dt_min.saved_value = 1.0e-6
        ui.lineedit_keyword_dt_max.required = True
        ui.lineedit_keyword_dt_max.saved_value = 1.0

        ui.lineedit_keyword_dt.post_update = self.check_dt_range

    def check_dt_range(self):
        ui = self.ui.run
        if self.project.get_value('dt_fac') == 1.0:
            dt = self.project.get_value('dt')
            if dt is not None:
                self.update_keyword('dt_min', dt)
                self.update_keyword('dt_max', dt)
                txt = ui.lineedit_keyword_dt.text()
                ui.lineedit_keyword_dt_min.setText(txt)
                ui.lineedit_keyword_dt_max.setText(txt)


    def set_run_steady_state(self):
        ui = self.ui.run
        ui.groupbox_adaptive_time_step.setEnabled(False)
        dt_fac = self.project.get_value('dt_fac')
        for k in ('dt', 'tstop', 'dt_min', 'dt_max', 'dt_fac'):
            if not(dt_fac == 1.0 and k in ('dt_min', 'dt_max')):
                self.retain_keyword(k)
            self.unset_keyword(k)
            getattr(ui, 'lineedit_keyword_'+k).setText('')
        self.retain_keyword('nice_dt')
        self.unset_keyword('nice_dt')

        for w in (ui.radiobutton_constant_dt,
                  ui.radiobutton_constant_ratio,
                  ui.radiobutton_nice_dt):
            w.setAutoExclusive(False)
            w.setChecked(False)
            w.setAutoExclusive(True)
        for w in (ui.label_tstop, ui.lineedit_keyword_tstop):
            w.setEnabled(False)

    def set_run_unsteady_state(self):
        ui = self.ui.run
        self.enable_widget(ui.groupbox_adaptive_time_step, ignore_lock=True)
        nice_dt = self.project.get_value('nice_dt')
        if nice_dt is None:
            nice_dt = self.get_retained_keyword('nice_dt')
            if nice_dt:
                self.update_keyword('nice_dt', True)
        keys = ['dt_min', 'dt_max', 'tstop']
        defaults = {'dt_fac':0.9,
                    'dt_min':1e-6,
                    'dt_max':1.0}
        if not nice_dt:
            keys += ['dt', 'dt_fac']
        for k in keys:
            val = self.project.get_value(k)
            if val is None:
                val = self.get_retained_keyword(k, default=defaults.get(k))
                if val is not None:
                    self.update_keyword(k, val)
            w = getattr(ui, 'lineedit_keyword_' + k)
            w.setText('' if val is None else str(val))
        for w in (ui.label_tstop, ui.lineedit_keyword_tstop):
            self.enable_widget(w)
        # Set radiobutton state.  Cannot call setup_run here
        # FIXME this is messy
        nice_dt = self.project.get_value('nice_dt')
        dt_fac = self.project.get_value('dt_fac')
        nice_dt = self.project.get_value('nice_dt')
        if nice_dt: #Nice DT
            ui.radiobutton_nice_dt.setChecked(True)
            ui.lineedit_keyword_dt_fac.setText('')
            ui.lineedit_keyword_dt_fac.setEnabled(False)
            self.set_run_nice_dt()
        elif dt_fac == 1.0: #Constant DT
            ui.radiobutton_constant_dt.setChecked(True)
            ui.lineedit_keyword_dt_fac.setEnabled(False)
            ui.lineedit_keyword_dt_fac.setText('1.0')
            self.set_run_constant_dt()
        else: #Constant ratio
            ui.radiobutton_constant_ratio.setChecked(True)
            self.enable_widget(ui.lineedit_keyword_dt_fac)
            self.set_run_constant_ratio()


    def set_run_constant_dt(self):
        ui = self.ui.run
        self.unset_keyword('nice_dt')
        self.clear_retained_keyword('nice_dt')
        dt_fac = self.project.get_value('dt_fac')
        if dt_fac != 1.0:
            self.retain_keyword('dt_fac')
        self.update_keyword('dt_fac',1.0)
        ui.lineedit_keyword_dt_fac.setText('1.0')
        val = self.project.get_value('dt')
        if dt_fac != 1.0:
            self.retain_keyword('dt_min')
        self.update_keyword('dt_min',val)
        if dt_fac != 1.0:
            self.retain_keyword('dt_max')
        self.update_keyword('dt_max',val)
        for w in (ui.label_dt_fac,
                  ui.lineedit_keyword_dt_fac,
                  ui.label_dt_min,
                  ui.lineedit_keyword_dt_min,
                  ui.label_dt_max,
                  ui.lineedit_keyword_dt_max):
            w.setEnabled(False)

    def set_run_constant_ratio(self):
        ui = self.ui.run
        self.unset_keyword('nice_dt')
        self.clear_retained_keyword('nice_dt')
        # restore dt_fac
        val = dt_fac = self.project.get_value('dt_fac')
        if val is None or val==1.0:
            val = self.get_retained_keyword('dt_fac', default=0.9)
            if val == 1.0:
                val = 0.9
            self.update_keyword('dt_fac', val)
        ui.lineedit_keyword_dt_fac.setText(str(val))

        # restore dt/dt_min/dt_max
        val = self.project.get_value('dt')
        if val is None:
            val = self.get_retained_keyword('dt_min', default=1e-6)
        if val is not None:
            self.update_keyword('dt', val)
            ui.lineedit_keyword_dt.setText(str(val))

        val = self.project.get_value('dt_min')
        if val is None or dt_fac == 1.0:
            val = self.get_retained_keyword('dt_min', default=1e-6)
        if val is not None:
            self.update_keyword('dt_min', val)
            ui.lineedit_keyword_dt_min.setText(str(val))

        val = self.project.get_value('dt_max')
        if val is None or dt_fac == 1.0:
            val = self.get_retained_keyword('dt_max', default=1.0)
        if val is not None:
            self.update_keyword('dt_max', val)
            ui.lineedit_keyword_dt_max.setText(str(val))

        for w in (ui.label_dt_fac,
                  ui.lineedit_keyword_dt_fac,
                  ui.label_dt_min,
                  ui.lineedit_keyword_dt_min,
                  ui.label_dt_max,
                  ui.lineedit_keyword_dt_max):
            self.enable_widget(w)


    def set_run_nice_dt(self):
        ui = self.ui.run
        self.update_keyword('nice_dt', True)
        if self.project.get_value('dt_fac') != 1.0:
            self.retain_keyword('dt_fac')
        self.unset_keyword('dt_fac')
        ui.lineedit_keyword_dt_fac.setEnabled(False)
        ui.label_dt_fac.setEnabled(False)
        #Restore dt/dt_min/dt_max
        dt = self.project.get_value('dt')
        if dt is None:
            dt = self.get_retained_keyword('dt', default=1e-6)
            self.update_keyword('dt', dt)
        ui.lineedit_keyword_dt.setText(str(dt) if dt is not None else '')

        dt_min = self.project.get_value('dt_min')
        if dt_min is None or (dt_min == dt):
            dt_min = self.get_retained_keyword('dt_min')
            if dt_min is None or (dt_min == dt):
                dt_min = 1e-6
        self.update_keyword('dt_min', dt_min)
        ui.lineedit_keyword_dt_min.setText(str(dt_min) if dt_min is not None else '')

        dt_max = self.project.get_value('dt_max')
        if dt_max is None or (dt_max == dt):
            dt_max = self.get_retained_keyword('dt_max')
            if dt_max is None or (dt_max == dt):
                dt_max = 1e-6
        self.update_keyword('dt_max', dt_max)
        ui.lineedit_keyword_dt_max.setText(str(dt_max) if dt_max is not None else '')

        for w in (ui.label_dt_min,
                  ui.lineedit_keyword_dt_min,
                  ui.label_dt_max,
                  ui.lineedit_keyword_dt_max):
            self.enable_widget(w)
        for w in (ui.label_dt_fac,
                  ui.lineedit_keyword_dt_fac):
            w.setEnabled(False)
        ui.lineedit_keyword_dt_fac.setText('')


    def setup_run(self, allow_disabled_tab=False):
        ui = self.ui.run
        # "Time marching" groupbox
        rb = ui.radiobutton_steady_state
        enable = self.project.solver in [SINGLE, TFM]
        self.set_widget_enabled(rb, enable, reason="Only available for TFM and single-phase models.")

        if not enable and rb.isChecked():
            rb.setChecked(False)
            ui.radiobutton_unsteady_state.setChecked(True)
            self.set_run_unsteady_state()
        if self.project.get_value('dt') is not None or self.project.get_value('tstop') is not None:
            ui.radiobutton_unsteady_state.setChecked(True)
            self.set_run_unsteady_state()
        elif enable:
            rb.setChecked(True)
            self.set_run_steady_state()


    def reset_run(self):
        ui = self.ui.run
        ui.lineedit_keyword_dt_fac.saved_value = 0.9
        ui.lineedit_keyword_dt_min.saved_value = 1.0e-6
        ui.lineedit_keyword_dt_max.saved_value = 1.0
