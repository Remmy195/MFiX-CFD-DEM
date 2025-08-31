# Methods to deal with solids DEM tab, split off from solids_handler.py
# This also handles SQP,CGP, and GSP

import os, shutil

from qtpy.QtWidgets import (QComboBox, QLabel, QLineEdit,
                            QWidget, QTableWidgetItem)
from qtpy.QtCore import Qt

from .widgets.base import LineEdit
from .project import Equation
from .tools.qt import (get_combobox_item, get_selected_row,
                       set_combobox_tooltip, set_item_enabled,
                       set_item_noedit, widget_iter)
from . import default_values

from .constants import *
# These should be in constants.py so project.py can normalize them

DES_INTG_METHODS = ['EULER', 'ADAMS_BASHFORTH']
DES_COLL_MODELS = ['LSD', 'HERTZIAN']
DES_INTERP_SCHEMES = ['NONE', 'SQUARE_DPVM', 'DPVM_SATELLITE', 'GARG_2012']

NONE, SQUARE_DPVM, DPVM_SATELLITE, GARG_2012 = range(4)


class SolidsDEM:
    def init_solids_dem(self):
        ui = self.ui.solids
        ui.groupbox_gener_part_config.clicked.connect(self.set_gener_part_config)
        ui.combobox_des_intg_method.activated.connect(self.set_des_intg_method)
        ui.combobox_des_intg_method.key = 'des_intg_method'
        ui.combobox_des_coll_model.activated.connect(self.set_des_coll_model)
        ui.combobox_des_coll_model.key = 'des_coll_model'
        ui.combobox_des_oneway_coupled.activated.connect(self.set_des_oneway_coupled)
        ui.checkbox_keyword_des_explicitly_coupled.post_update = self.setup_solids_dem_tab
        ui.combobox_des_interp.activated.connect(self.set_des_interp)
        ui.combobox_des_interp_scheme.activated.connect(self.set_des_interp_scheme)
        ui.checkbox_enable_des_diffuse_width.clicked.connect(self.enable_des_diffuse_width)
        ui.combobox_cohesion_model.activated.connect(self.set_cohesion_model)
        ui.combobox_cohesion_model.keys = ['use_cohesion', 'van_der_waals']
        ui.combobox_des_neighbor_search.activated.connect(self.set_des_neighbor_search)
        ui.lineedit_keyword_des_diffuse_width.setdtype('dp')
        ui.combobox_des_conv_corr.activated.connect(self.set_des_conv_corr)
        self.add_tooltip(ui.combobox_des_interp_scheme, key='des_interp_scheme')
        self.add_tooltip(ui.combobox_des_conv_corr, key='des_conv_corr')
        self.add_tooltip(ui.checkbox_enable_des_diffuse_width, key='des_diffuse_width')

        ui.groupbox_gener_part_config.key = 'gener_part_config'

        self.solids_dem_saved_solids_names = [] # keep track of when we need to rebuild per-phase gui items

        for w in widget_iter(ui.groupbox_filter_particles):
            if isinstance(w, LineEdit):
                w.hidden_ctrl = ui.groupbox_filter_particles
                w.post_update = [self.setup_solids_dem_tab, self.update_particle_preview]

        # GSP Explicit/Implicit
        cb = ui.combobox_gsp_model
        cb.keys = ['gsp_model']
        self.add_tooltip(get_combobox_item(cb,0), key='gsp_model', value='EXPLICIT')
        self.add_tooltip(get_combobox_item(cb,1), key='gsp_model', value='IMPLICIT')
        cb.currentIndexChanged.connect(self.solids_handle_combobox_gsp_model)
        self.gsp_hertzian_warning_printed = False

        # Tooltips for combobox items
        key = 'des_intg_method'
        cb = ui.combobox_des_intg_method
        for i, m in enumerate(DES_INTG_METHODS):
            item = get_combobox_item(cb, i)
            self.add_tooltip(item, key, value=m)

        key = 'des_coll_model'
        cb = ui.combobox_des_coll_model
        self.add_tooltip(cb, 'des_coll_model')
        for i, m in enumerate(DES_COLL_MODELS):
            item = get_combobox_item(cb, i)
            self.add_tooltip(item, key, value=m)

        key = 'des_oneway_coupled'
        cb = ui.combobox_des_oneway_coupled
        cb.key = key
        for i, m in enumerate((True, False)):
            item = get_combobox_item(cb, i)
            self.add_tooltip(item, key, value=m)
        label = ui.label_des_oneway_coupled
        self.add_tooltip(label, key)

        key = 'des_interp_scheme'
        cb = ui.combobox_des_interp_scheme
        for i, m in enumerate(DES_INTERP_SCHEMES):
            item = get_combobox_item(cb, i)
            self.add_tooltip(item, key, value=m)

        cb = ui.combobox_des_interp
        text = ''
        cb.keys = ['des_interp_on', 'des_interp_mean_fields']
        for key in cb.keys:
            doc = self.keyword_doc[key]
            if text:
                text += '<br>'
            text += '<b>%s</b>: %s' % (key, doc['description'])
        ui.label_des_interp_on.setToolTip(text)
        ui.label_des_interp_on.keys = cb.keys
        for i in range(4):
            key = 'des_interp_on'
            val = not bool(i>>1)
            vkey = '.TRUE.' if val else '.FALSE.'
            text = '<b>%s=%s:</b> %s' % (key, val,
                                         self.keyword_doc[key]['valids'][vkey]['note'])
            key = 'des_interp_mean_fields'
            val = not bool(i%2)
            vkey = '.TRUE.' if val else '.FALSE.'
            text += '<br><b>%s=%s:</b> %s' % (key, val,
                                              self.keyword_doc[key]['valids'][vkey]['note'])
            get_combobox_item(cb, i).setToolTip(text)

        key = 'des_conv_corr'
        cb = ui.combobox_des_conv_corr
        for i, m in enumerate(DES_CONV_CORR_VALUES):
            item = get_combobox_item(cb, i)
            self.add_tooltip(item, key, value=m)

        cb = ui.combobox_time_step_option
        cb.key = 'dtsolid'
        item = get_combobox_item(cb, 0) # Automatic
        item.setToolTip('<html>Solids time step is automatically computed based on particle properties. It will be set to the collision time divided by the DEM time step factor.</html>')
        item = get_combobox_item(cb, 1) # Constant
        item.setToolTip('<html>Solids time step does not depend on particle properties. It must be set to sufficiently small value (less than collision time) to avoid unphysical results.</html>')
        cb.currentIndexChanged.connect(self.solids_handle_combobox_time_step_option)

        self.fixup_solids_table(ui.tablewidget_psd)

        ui.groupbox_filter_particles.clicked.connect(
            self.solids_handle_groupbox_filter_particles)

        ui.lineedit_keyword_particles.post_update = self.solids_handle_particles
        ui.groupbox_filter_particles_checked = None # Tristate
        # Set up widgets for _exclude keys
        for w in widget_iter(ui.groupbox_filter_particles):
            w.hidden_ctrl = ui.groupbox_filter_particles
            if not isinstance(w, (QComboBox, QLabel)):
                continue
            name = w.objectName()
            if name.startswith('combobox_filter_'):
                w.key = 'part_in_' +  name.split('combobox_filter_', 1)[-1] + '_exclude'
                w.args = None
                self.add_tooltip(w, w.key)
                w.currentIndexChanged.connect(lambda idx, w=w:
                                              self.solids_handle_part_in_exclude(w, idx))
            if name.startswith('label_filter_'):
                if name.endswith('_units'):
                    name = name[:-6]
                w.key = 'part_in_' +  name.split('label_filter_', 1)[-1] + '_max'
                w.args = None
                self.add_tooltip(w, w.key)

        # Particle size distribution, this is common to DEM/CGP/SQP/GSP, and PIC #1356
        ui.groupbox_psd.setToolTip("") # why?
        tb = ui.toolbutton_add_psd
        tb.clicked.connect(self.solids_psd_add)
        tb = ui.toolbutton_delete_psd
        tb.clicked.connect(self.solids_psd_delete)
        tb = ui.toolbutton_edit_psd
        tb.clicked.connect(self.solids_psd_edit)


    def set_gener_part_config(self, val, use_retained=True):
        ui = self.ui.solids
        enabled = not val
        prev_val = self.project.get_value('gener_part_config', default=False)
        for item in (ui.label_particles,
                     ui.lineedit_keyword_particles,
                     ui.groupbox_filter_particles):
            self.set_widget_enabled(item, enabled)

        if val:
            self.update_keyword('gener_part_config', val)
        else:
            if self.project.get_value('gener_part_config', default=False) != False:
                self.unset_keyword('gener_part_config') # Don't unset if set to False
        if val: # Automatic particle generation
            if use_retained:
                self.retain_keyword("particles")
            self.unset_keyword("particles")
            ui.lineedit_keyword_particles.setText('')
            ui.groupbox_filter_particles.setChecked(False)
            self.solids_handle_groupbox_filter_particles(False, use_retained=use_retained)

        else:   # Use particle_input.dat
            if use_retained:
                ret = self.get_retained_keyword("particles")
                if ret:
                    self.update_keyword("particles", ret)
            if any(k[0].startswith('part_in_') for k in self.retained_keys):
                ui.groupbox_filter_particles.setChecked(True)
                self.solids_handle_groupbox_filter_particles(True, use_retained=use_retained)

        if val != prev_val:
            self.update_particle_preview()


    def set_des_intg_method(self, idx):
        key = 'des_intg_method'
        val = DES_INTG_METHODS[idx]
        self.update_keyword(key, val)
        self.setup_solids_dem_tab()

    def set_des_coll_model(self, idx):
        key = 'des_coll_model'
        val = DES_COLL_MODELS[idx]
        self.update_keyword(key, val)
        self.setup_solids_dem_tab()

    def set_des_oneway_coupled(self, val):
        self.update_keyword('des_oneway_coupled', [True,False][val])
        self.setup_solids_dem_tab()

    def set_des_interp(self, val):
        des_interp_on = not bool(val>>1)
        des_interp_mean_fields = not bool(val%2)
        self.update_keyword('des_interp_on', des_interp_on)
        self.update_keyword('des_interp_mean_fields', des_interp_mean_fields)
        self.setup_solids_dem_tab()

    def set_des_interp_scheme(self, val):
        ui = self.ui.solids
        des_interp_scheme = DES_INTERP_SCHEMES[val]

        self.update_keyword('des_interp_scheme', des_interp_scheme)
        key = 'des_interp_width' # Only valid for SQUARE_DPVM
        if des_interp_scheme != 'SQUARE_DPVM':
            self.retain_keyword(key)
            self.unset_keyword(key)
            ui.lineedit_keyword_des_interp_width.setText('')
        else:
            value = self.get_retained_keyword(key)
            if value is not None:
                self.update_keyword(key, value)
                ui.lineedit_keyword_des_interp_width.setText(str(value))
        self.setup_solids_dem_tab()


    def enable_des_diffuse_width(self, val):
        ui = self.ui.solids
        enabled = val
        for item in (ui.label_des_diffuse_width, ui.lineedit_keyword_des_diffuse_width,
                     ui.label_des_diffuse_width_units):
            self.set_widget_enabled(item, enabled)
        if not enabled:
            self.unset_keyword('des_diffuse_width')
        else: #Restore value
            self.update_keyword('des_diffuse_width',
                                ui.lineedit_keyword_des_diffuse_width.value)


    def set_cohesion_model(self, val):
        for kw in ('use_cohesion', 'van_der_waals'):
            self.update_keyword(kw, bool(val))
        self.setup_solids_dem_tab()

    def set_des_neighbor_search(self, val):
        self.update_keyword('des_neighbor_search', 4 if val==0 else 1)
        # No setup_solids_dem_tab needed

    def set_des_conv_corr(self, val):
        self.update_keyword('des_conv_corr', value=DES_CONV_CORR_VALUES[val])

    def setup_solids_dem_tab(self):
        # Ensures all constraints (items enabled/disabled) are set

        # After any item in this tab is edited, we call setup_solids_dem_tab again to make
        # sure all constraints (disabled/enabled inputs) are enforced
        # Therefore, 'setup_solids_dem_tab' must not call any of the setters

        ui = self.ui.solids
        ui.DEM.setEnabled(any(b.isEnabled() for b in
                              (ui.pushbutton_solids_DEM,
                               ui.pushbutton_solids_CGP,
                               ui.pushbutton_solids_SQP,
                               ui.pushbutton_solids_GSP))) # We may be just visiting due to locate_keyword

        #MFIX-UI_SRS
        #Enable automatic particle generation
        # Enabled sets keyword GENER_PART_CONFIG to true
        # Disabled enables the user to specify number of entries in particle input file
        # Default value is 0
        default = 0
        # Sets keyword PARTICLES
        gener_part_config = self.project.get_value('gener_part_config',
                                                   default=False)
        particles = self.project.get_value('particles')
        if gener_part_config:
            if particles: # Should not both be set
                self.warning("gener_part_config set, particles=%s" % particles)
                particles = default
        else:
            if particles is not None and particles < 0:
                self.warning("Invalid particles %s" % particles)
                particles = default
        self.update_keyword('particles', particles)

        enabled = not gener_part_config
        for item in (ui.label_particles, ui.lineedit_keyword_particles):
            self.set_widget_enabled(item, enabled)

        # Set groupbox to correct state
        gb = ui.groupbox_filter_particles

        if ui.groupbox_filter_particles_checked is None:
            ui.groupbox_filter_particles_checked = any(k.key.startswith('part_in_')
                                                       for k in self.project.keywordItems())

        gb.setChecked(enabled and ui.groupbox_filter_particles_checked)

        energy_eq = self.project.get_value('energy_eq', default=True)
        for w in (ui.combobox_filter_temp,
                  ui.label_filter_temp,
                  ui.label_filter_temp_units,
                  ui.lineedit_keyword_part_in_temp_min,
                  ui.lineedit_keyword_part_in_temp_max):
            self.set_widget_enabled(w, energy_eq, reason='Requires energy equations.')

        # PSD table is shared between DEM and PIC per #1356 (should it be in Materials?)
        gb = ui.groupbox_psd
        par = gb.parent()
        if par != ui.frame_dem: # Widget has moved to PIC pane
            gb.hide()
            par.layout().removeWidget(gb)
            layout = ui.frame_dem.layout()
            layout.insertWidget(2, gb)
            gb.show()

        # Filter particle_input.dat
        gb = ui.groupbox_filter_particles
        layout = gb.layout()
        last = ui.lineedit_keyword_part_in_temp_min
        # Dynamic widgets for phase, species and usr vars
        n_species = max((self.project.get_value('nmax_s', default=0, args=[p])
                         for p in range(1, 1+len(self.solids))),
                        default=0)
        key = 'part_in_x_s_exclude'
        if key not in ui.dynamic_widgets.keys():
            ui.dynamic_widgets[key] = []
        dws = ui.dynamic_widgets[key]
        while len(dws) < n_species:
            n = 1+len(dws)
            cb = QComboBox()
            cb.addItem("Include")
            cb.addItem("Exclude")
            cb.key = key
            cb.args = [n]
            cb.currentIndexChanged.connect(lambda idx, w=cb:
                                           self.solids_handle_part_in_exclude(w, idx))

            self.add_tooltip(cb, cb.key)
            layout.addWidget(cb)
            QWidget.setTabOrder(last, cb)
            last = cb
            l1 = QLabel("Species %s mass fraction" % n)
            l1.key = 'part_in_x_s_min'
            l1.args = [n]
            self.add_tooltip(l1, l1.key)
            layout.addWidget(l1)
            le1 = LineEdit()
            le1.dtype = float
            le1.key = 'part_in_x_s_min'
            le1.args = [n]
            le1.min = 0.0
            le1.max = 1.0
            le1.value_updated.connect(self.project.submit_change)
            le1.post_update = [self.setup_solids_dem_tab, self.update_particle_preview]
            self.add_tooltip(le1, le1.key)
            layout.addWidget(le1)
            QWidget.setTabOrder(last, le1)
            last = le1
            self.add_tooltip(le1, le1.key)
            le2 = LineEdit()
            le2.dtype = float
            le2.key = 'part_in_x_s_max'
            le2.args = [n]
            le2.min = 0.0
            le2.max = 1.0
            le2.value_updated.connect(self.project.submit_change)
            le2.post_update = [self.setup_solids_dem_tab, self.update_particle_preview]
            self.add_tooltip(le2, le2.key)
            layout.addWidget(le2)
            QWidget.setTabOrder(last, le2)
            last = le2
            l2 = QLabel('') # Unit placeholder
            layout.addWidget(l2)
            dws.append((cb, l1, le1, le2, l2))

        while len(dws) > n_species:
            for kw in ('part_in_x_s_min', 'part_in_x_s_max', 'part_in_x_s_exclude'):
                self.unset_keyword(kw, args=[len(dws)])
            row = dws.pop()
            for w in row:
                layout.removeWidget(w)
                w.deleteLater()

        def as_str(x):
            return '' if x is None else str(x)

        for i in range(1, n_species+1):
            val = self.project.get_value('part_in_x_s_min', args=[i])
            dws[i-1][2].setText(as_str(val))
            val = self.project.get_value('part_in_x_s_max', args=[i])
            dws[i-1][3].setText(as_str(val))

        # User vars
        n_vars = self.project.get_value('des_usr_var_size', default=0)
        key = 'part_in_usr_var_exclude'
        if key not in ui.dynamic_widgets.keys():
            ui.dynamic_widgets[key] = []
        dws = ui.dynamic_widgets[key]
        while len(dws) < n_vars:
            n = 1+len(dws)
            cb = QComboBox()
            cb.addItem("Include")
            cb.addItem("Exclude")
            cb.key = key
            cb.args = [n]
            cb.currentIndexChanged.connect(lambda idx, w=cb:
                                           self.solids_handle_part_in_exclude(w, idx))
            self.add_tooltip(cb, cb.key)
            layout.addWidget(cb)
            QWidget.setTabOrder(last, cb)
            last = cb
            l1 = QLabel("User scalar %s" % n)
            l1.key = 'part_in_usr_var_min'
            l1.args = [n]
            self.add_tooltip(l1, l1.key)
            layout.addWidget(l1)
            le1 = LineEdit()
            le1.dtype = float
            le1.key = 'part_in_usr_var_min'
            le1.args = [n]
            le1.value_updated.connect(self.project.submit_change)
            le1.post_update = [self.setup_solids_dem_tab, self.update_particle_preview]
            self.add_tooltip(le1, le1.key)
            layout.addWidget(le1)
            QWidget.setTabOrder(last, le1)
            last = le1
            le2 = LineEdit()
            le2.dtype = float
            le2.key = 'part_in_usr_var_max'
            le2.args = [n]
            le2.value_updated.connect(self.project.submit_change)
            le2.post_update = [self.setup_solids_dem_tab, self.update_particle_preview]
            self.add_tooltip(le2, le2.key)
            layout.addWidget(le2)
            QWidget.setTabOrder(last, le2)
            last = le2
            l2 = QLabel('') # Unit placeholder
            layout.addWidget(l2)
            dws.append((cb, l1, le1, le2, l2))

        while len(dws) > n_vars:
            for kw in ('part_in_usr_var_min', 'part_in_usr_var_max', 'part_in_usr_var_exclude'):
                self.unset_keyword(kw, args=[len(dws)])
            row = dws.pop()

            for w in row:
                layout.removeWidget(w)
                w.deleteLater()

        for i in range(1, n_vars+1):
            val = self.project.get_value('part_in_usr_var_min', args=[i])
            dws[i-1][2].setText(as_str(val))
            val = self.project.get_value('part_in_usr_var_max', args=[i])
            dws[i-1][3].setText(as_str(val))

        key = 'part_in_phase'
        n_phases = len(self.solids)
        if key not in ui.dynamic_widgets.keys():
            ui.dynamic_widgets[key] = []
        dws = ui.dynamic_widgets[key]
        while len(dws) < n_phases:
            n = 1+len(dws)
            cb = QComboBox()
            cb.addItem("Include")
            cb.addItem("Exclude")
            cb.key = key
            cb.args = [n]
            cb.currentIndexChanged.connect(lambda idx, w=cb:
                                           self.solids_handle_part_in_phase(w, idx))
            self.add_tooltip(cb, cb.key)
            row = layout.rowCount()
            layout.addWidget(cb, row, 0, 1, 1)
            QWidget.setTabOrder(last, cb)
            last = cb
            l = QLabel('')
            l.args = [n]
            self.add_tooltip(l, key)
            layout.addWidget(l, row, 1, 1, -1)
            dws.append((cb, l))

        while len(dws) > n_phases:
            self.unset_keyword('part_in_phase', args=[len(dws)])
            row = dws.pop()
            for w in row:
                layout.removeWidget(w)
                w.deleteLater()

        for (i, name) in enumerate(self.solids):
            dws[i][1].setText(name)

        # Set comboboxes to correct state
        for w in widget_iter(ui.groupbox_filter_particles):
            if not isinstance(w, QComboBox):
                continue
            key = getattr(w, 'key', None)
            if key is None:
                continue
            args = w.args
            default = not key.endswith('_exclude')
            val = self.project.get_value(key, args=args, default=default)
            # Default  Val     Idx
            # False    False   0
            # False    True    1
            # True     False   1
            # True     True    0
            w.setCurrentIndex(int(val != default))

        # Set groupbox to correct state, do this last to handle dynamic_widgets
        ui.groupbox_gener_part_config.setChecked(gener_part_config)
        self.set_widget_enabled(ui.groupbox_gener_part_config, True) # Lock during run
        self.set_gener_part_config(gener_part_config, use_retained=False) # UI update
        self.solids_handle_groupbox_filter_particles(gb.isChecked(), use_retained=False) # collapse/uncollapse
        ui.groupbox_filter_particles.setToolTip('Select which particles to read from particle_input.dat')
        # Particle size distribution, issues/1158
        self.update_solids_psd_table()
        self.fixup_solids_table(ui.tablewidget_psd)

        #Select numerical integration method
        # Selection always available
        # Available selections
        #Euler [DEFAULT]
        # Selection always available
        # Sets keyword DES_INTG_METHOD to 'EULER'
        #Adams-Bashforth
        # Selection always available
        # Sets keyword DES_INTG_METHOD to 'ADAMS_BASHFORTH'1
        key = 'des_intg_method'
        cb = ui.combobox_des_intg_method
        des_intg_method = self.project.get_value(key, default='EULER')
        if des_intg_method not in DES_INTG_METHODS:
            self.warn("Invalid des_intg_method %s" % des_intg_method)
            des_intg_method = 'EULER'
        cb.setCurrentIndex(DES_INTG_METHODS.index(des_intg_method))
        cb.setToolTip(get_combobox_item(cb, cb.currentIndex()).toolTip())

        #Selection collision model
        # Selection always available
        # Available selections
        #Linear Spring-Dashpot [DEFAULT]
        # Selection always available
        # Sets keyword DES_COLL_MODEL to 'LSD'
        #Hertzian
        # Selection always available
        # Sets keyword DES_COLL_MODEL to 'HERTZIAN'
        # Issues/1646:  Force DES_COLL_MODEL to 'HERTZIAN' for SQP cases
        key = 'des_coll_model'
        cb = ui.combobox_des_coll_model
        item = get_combobox_item(cb, DES_COLL_MODELS.index('LSD'))
        des_coll_model = self.project.get_value(key, default='LSD')
        des_explicitly_coupled = self.project.get_value('des_explicitly_coupled', default=False)
        if self.project.solver == SQP:
            self.disable_widget(item, reason='Not available for SQP model.')
            if des_coll_model != 'HERTZIAN':
                self.warn("Setting Hertzian collision model for SQP", popup=True)
                des_coll_model = 'HERTZIAN'
                self.update_keyword(key, des_coll_model)
        else:
            self.enable_widget(item)

        # Issues/1870: Enable DES_COLL_MODEL to 'HERTZIAN' for GSP cases
        # Pop up warning to notify the user
        item = get_combobox_item(cb, DES_COLL_MODELS.index('HERTZIAN'))
        if self.project.solver == GSP:
            self.enable_widget(item)
            item.setToolTip(item.tooltip0 + '<br>Not recommended for GSP.')
            if des_coll_model != 'LSD':
                if not self.gsp_hertzian_warning_printed:
                    self.gsp_hertzian_warning_printed = True
                    self.warn("LSD collision model is recommended for GSP", popup=True)
            else:
                self.gsp_hertzian_warning_printed = False
        else:
            item.setToolTip(item.tooltip0)
        set_combobox_tooltip(cb)

        # Show gsp_model for GSP only
        for w in (ui.label_gsp_model, ui.combobox_gsp_model):
            w.setVisible(self.project.solver == GSP)
        if self.project.solver == GSP:
            gsp_model = self.project.get_value('gsp_model', default='EXPLICIT')
            if gsp_model not in ("EXPLICIT", "IMPLICIT"):
                self.warning('Invalid value "%s" for gsp_model, setting to EXPLICIT'%
                             gsp_model,
                             popup=True)
                gsp_model = 'EXPLICIT'
                self.update_keyword('gsp_model', gsp_model)
            ui.combobox_gsp_model.setCurrentIndex(0 if gsp_model=="EXPLICIT"
                                                  else 1)
            set_combobox_tooltip(ui.combobox_gsp_model)

        if des_coll_model not in DES_COLL_MODELS:
            self.warn("Invalid des_coll_model %s" % des_coll_model)
            des_coll_model = 'LSD'
        cb.setCurrentIndex(DES_COLL_MODELS.index(des_coll_model))
        set_combobox_tooltip(cb)

        #Select gas-solids coupling scheme:
        # Selection unavailable if fluid model is disabled
        # Available selections:
        # One-way Coupled
        # Selection always available
        # Sets keyword DES_ONEWAY_COUPLED true
        # Fully Coupled
        # Selection always available
        # Sets keyword DES_ONEWAY_COUPLED false
        enabled = not self.fluid_solver_disabled
        for item in (ui.label_des_oneway_coupled, ui.combobox_des_oneway_coupled):
            self.set_widget_enabled(item, enabled)
        key = 'des_oneway_coupled'
        cb = ui.combobox_des_oneway_coupled
        des_oneway_coupled = self.project.get_value(key, default=False)
        if des_oneway_coupled not in (True, False):
            self.warn("Invalid des_oneway_coupled %s" % des_oneway_coupled)
            des_oneway_coupled = False
            self.update_keyword(key, des_oneway_coupled)
        cb.setCurrentIndex(0 if des_oneway_coupled else 1)
        cb.setToolTip(get_combobox_item(cb, cb.currentIndex()).toolTip())

        #Optional to enable explicitly coupled simulation
        # Unavailable for GARG_2012 interpolation
        des_interp_scheme = self.project.get_value('des_interp_scheme', default='NONE')
        enabled = (des_interp_scheme!='GARG_2012')
        self.set_widget_enabled(ui.checkbox_keyword_des_explicitly_coupled,
                                enabled)

        #Select interpolation framework:
        # Selection always available
        # Available selections:
        # Field-to-Particle and Particle-to-Field [DEFAULT]
        #  Sets keyword DES_INTERP_ON to true
        #  Sets keyword DES_INTERP_MEAN_FIELDS to true
        # Field-to-Particle only
        #  Sets keyword DES_INTERP_ON to true
        #  Sets keyword DES_INTERP_MEAN_FIELDS to false
        # Particle-to-Field only
        #  Sets keyword DES_INTERP_ON to false
        #  Sets keyword DES_INTERP_MEAN_FIELDS to true
        # No Interpolation
        #  Sets keyword DES_INTERP_ON to false
        #  Sets keyword DES_INTERP_MEAN_FIELDS to false
        #
        # issues/116 must also set DES_INTERP_SCHEME to None when no-interpolation
        #
        # issues/1646 this must be on for SQP
        cb = ui.combobox_des_interp
        sqp = self.project.solver == SQP
        des_interp_on = self.project.get_value('des_interp_on', default=False)
        des_interp_mean_fields = self.project.get_value('des_interp_mean_fields', default=False)

        item = get_combobox_item(cb, 3)
        item.setEnabled(not sqp)
        if sqp:
            item.setToolTip("Interpolation must be enabled for SQP cases.")
        else:
            key = 'des_interp_on'
            text = '<b>%s=%s:</b> %s' % (key, False,
                                         self.keyword_doc[key]['valids']['.FALSE.']['note'])
            key = 'des_interp_mean_fields'
            text += '<br><b>%s=%s:</b> %s' % (key, False,
                                              self.keyword_doc[key]['valids']['.FALSE.']['note'])
            item.setToolTip(text)


        if des_interp_on not in (True, False):
            self.warn("Invalid des_interp_on %s" % des_interp_on)
            des_interp_on = True if sqp else False
            self.update_keyword('des_interp_on', des_interp_on)

        if des_interp_mean_fields not in (True, False):
            self.warn("Invalid des_interp_mean_fields %s" % des_interp_mean_fields)
            des_interp_mean_fields = True if sqp else False
            self.update_keyword('des_interp_mean_fields', des_interp_mean_fields)

        interp_enabled = des_interp_on or des_interp_mean_fields # not no-interp
        if sqp and not interp_enabled:
            self.warn("Enabling interpolation for SQP", popup=True)
            des_interp_mean_fields = True
            des_interp_on = True
            self.update_keyword('des_interp_on', True)
            self.update_keyword('des_interp_mean_fields', True)
            interp_enabled = True

        index = 2*(1-des_interp_on) + (1-des_interp_mean_fields)
        cb.setCurrentIndex(index)
        cb.setToolTip(get_combobox_item(cb, cb.currentIndex()).toolTip())

        #Select interpolation scheme:
        # Selection available except when no-interpolation framework is selected
        # Available selections:
        #  None [locked default for no-interpolation framework]
        #  Selection always available
        #  Sets keyword DES_INTERP_SCHEME='NONE'
        # Square DPVM
        #  Selection always available
        #  Requires an interpolation width, DES_INTERP_WIDTH
        #  Sets keyword DES_INTERP_SCHEME='SQUARE_DPVM'
        # Garg 2012
        #  Selection not available with explicit coupling enabled
        #  issues/557  not available with energy equations or species equations
        #  Sets keyword DES_INTERP_SCHEME='GARG_2012'
        # Issues/1646:  Force DES_INTERP_SCHEME to 'DPVM_SATELLITE' for SQP cases
        cb = ui.combobox_des_interp_scheme
        label = ui.label_des_interp_scheme
        if self.project.solver == SQP and des_interp_scheme != 'DPVM_SATELLITE':
            des_interp_scheme = 'DPVM_SATELLITE'
            self.warn("Setting DES_INTERP_SCHEME to 'DPVM_SATELLITE' for SQP", popup=True)
            self.update_keyword('des_interp_scheme', des_interp_scheme)


        for item in (cb, label):
            self.set_widget_enabled(item, interp_enabled)
        if not interp_enabled or des_interp_scheme not in DES_INTERP_SCHEMES:
            des_interp_scheme = 'NONE'
            self.update_keyword('des_interp_scheme', des_interp_scheme)
        cb.setCurrentIndex(DES_INTERP_SCHEMES.index(des_interp_scheme))
        cb.setToolTip(get_combobox_item(cb, cb.currentIndex()).toolTip())

        # per-item enable flags
        energy_eq = self.project.get_value('energy_eq', default=True)
        mmax = self.project.get_value('mmax', default=1)
        species_eq = [self.project.get_value('species_eq', default=True, args=[i])
                      for i in range(0, 1+mmax)] # including index 0 = fluid

        key = 'des_interp_scheme'

        item = get_combobox_item(cb, NONE)
        if sqp:
            item.setEnabled(False)
            item.setToolTip("Not compatible with SQP")
        else:
            item.setEnabled(True)
            self.add_tooltip(item, key, value='NONE')

        item = get_combobox_item(cb, SQUARE_DPVM)
        if sqp:
            item.setEnabled(False)
            item.setToolTip("Not compatible with SQP")
        else:
            item.setEnabled(True)
            self.add_tooltip(item, key, value='SQUARE_DPVM')

        item = get_combobox_item(cb, DPVM_SATELLITE)
        if sqp:
            item.setEnabled(True)
            self.add_tooltip(item, key, value='DPVM_SATELLITE')
        else:
            item.setEnabled(False)
            item.setToolTip('SQP only.')

        item = get_combobox_item(cb,GARG_2012)
        enabled = not (sqp
                      or energy_eq
                      or des_explicitly_coupled or any(species_eq))
        item.setEnabled(enabled)
        if not enabled:
            self.add_tooltip(item,
                             key,
                             description="GARG_2012 interpolation is not compatible with "
                             "SQP, energy equation, species equation, or explicit coupling")
        else:
            self.add_tooltip(item,
                             key,
                             value='GARG_2012')
        cb.setToolTip(get_combobox_item(cb, cb.currentIndex()).toolTip())

        idx = DES_INTERP_SCHEMES.index(des_interp_scheme)
        # Don't stay on a disabled item
        if not get_combobox_item(cb, idx).isEnabled():
            des_interp_scheme = 'DPVM_SATELLITE' if sqp else 'NONE' if not interp_enabled else 'SQUARE_DPVM'
            self.warn("Setting 'DES_INTERP_SCHEME' to '%s'" % des_interp_scheme, popup=True)
            self.update_keyword('des_interp_scheme', des_interp_scheme)
            idx = DES_INTERP_SCHEMES.index(des_interp_scheme)

        cb.setCurrentIndex(idx)

        #Define interpolation width (DPVM only) (required)
        # Specification only available with SQUARE_DPVM interpolation scheme
        # Sets keyword DES_INTERP_WIDTH
        # TODO default?
        key = 'des_interp_width'
        enabled = interp_enabled and (des_interp_scheme=='SQUARE_DPVM') #?
        for item in (ui.label_des_interp_width, ui.lineedit_keyword_des_interp_width,
                     ui.label_des_interp_width_units):
            self.set_widget_enabled(item, enabled)

        #Option to enable diffusion of particle data
        # Selection unavailable with GARG_2012 interpolation scheme
        # No keyword is set by this option
        # Enables the user to specify a diffusion width
        # Sets keyword DES_DIFFUSE_WIDTH
        key = 'des_diffuse_width'
        enabled = (des_interp_scheme!='GARG_2012')
        self.set_widget_enabled(ui.checkbox_enable_des_diffuse_width,
                                enabled)
        if not enabled:
            ui.checkbox_enable_des_diffuse_width.setChecked(False)
            self.unset_keyword(key)
            ui.lineedit_keyword_des_diffuse_width.clear() # ??? FIXME
        enabled = ui.checkbox_enable_des_diffuse_width.isChecked()
        for item in (ui.label_des_diffuse_width, ui.lineedit_keyword_des_diffuse_width,
                     ui.label_des_diffuse_width_units):
            self.set_widget_enabled(item, enabled)
            if enabled:
                self.update_keyword(key, ui.lineedit_keyword_des_diffuse_width.value)


        #Specify friction coefficient
        # Specification always required
        # Sets keyword MEW (MEW_W)
        pass

        # Solids time step option (issues/1855)
        idx = 1 if self.project.get_value('dtsolid') is not None else 0
        ui.combobox_time_step_option.setCurrentIndex(idx)
        self.solids_handle_combobox_time_step_option(idx)

        #Specify normal spring constant
        # Only available for LSD collision model
        # Sets keyword KN (KN_W)
        #Specify tangential spring constant factor
        # Only available for LSD collision model
        # Sets keyword KT_FAC (KT_W_FAC)
        # Default values of 2.0/7.0
        #Specify tangential damping coefficient factor
        # Only available for LSD collision model
        # Sets keyword DES_ETAT_FAC (DES_ETAT_W_FAC)
        # Default values of 0.5
        enabled = (des_coll_model=='LSD')
        for item in (ui.label_kn,
                     ui.lineedit_keyword_kn,
                     ui.lineedit_keyword_kn_w,
                     ui.label_kn_units,
                     ui.label_kt_fac,
                     ui.lineedit_keyword_kt_fac,
                     ui.lineedit_keyword_kt_w_fac,
                     ui.label_des_etat_fac,
                     ui.lineedit_keyword_des_etat_fac,
                     ui.lineedit_keyword_des_etat_w_fac):
            self.set_widget_enabled(item, enabled)


        if enabled: # TODO set these defaults at load-time, not when this tab is shown
            for (key, default) in [('kt_fac', Equation('2/7')), ('kt_w_fac', Equation('2/7')),
                                   ('des_etat_fac', 0.5), ('des_etat_w_fac', 0.5)]:
                if self.project.get_value(key) is None:
                    self.update_keyword(key, default)

        # Unset keywords if not enabled?
        #Specify Young's modulus
        # Only available for Hertzian collision model
        # Sets keyword E_YOUNG (EW_YOUNG)
        #Specify Poisson ratio:
        # Only available for Hertzian collision model
        # Sets keyword V_POISSON (VW_POISSON)

        layout = ui.gridlayout_dem_parameters
        enabled = (des_coll_model=='HERTZIAN')
        for item in (ui.groupbox_young,
                     ui.groupbox_poisson):
            self.set_widget_enabled(item, enabled)

        key = 'des_conv_corr'
        enabled = (energy_eq)
        for item in (ui.label_des_conv_corr,
                     ui.combobox_des_conv_corr):
            self.set_widget_enabled(item, enabled, reason='Requires energy equations.')


        val = self.project.get_value(key, default=DEFAULT_DES_CONV_CORR)
        if val not in DES_CONV_CORR_VALUES:
            self.error('Invalid value for "%s", resetting to %s' %
                       (key, DES_CONV_CORR_VALUES),
                       popup=True)
            val = DEFAULT_DES_CONV_CORR
            self.update_keyword(key, val)
        ui.combobox_des_conv_corr.setCurrentIndex(DES_CONV_CORR_VALUES.index(val))

        solids_names = list(self.solids.keys())
        row = 0
        if self.solids_dem_saved_solids_names != solids_names:

            for (key, gb) in (('e_young', ui.groupbox_young),
                              ('v_poisson', ui.groupbox_poisson),
                              ('e_young_actual', ui.groupbox_young_actual),
                              ('v_poisson_actual', ui.groupbox_poisson_actual)):
                layout = gb.layout()
                # Delete all the old ones
                while layout.count() > 0:
                    item = layout.itemAt(0)
                    if not item:
                        continue
                    w = item.widget()
                    if not w:
                        continue
                    if isinstance(w, LineEdit):
                        self.project.unregister_widget(w)
                    w.hide()
                    layout.removeWidget(w)
                    w.setParent(None)
                    w.deleteLater()

                # ...and make new ones
                for (p, name) in enumerate(self.solids.keys(), 1):
                    row = p-1
                    label = QLabel(name)
                    #label.setObjectName('label_%s_args_%s' % (key,p))
                    #setattr(ui, label.objectName(), label)
                    label.args = [p]
                    self.add_tooltip(label, key)
                    layout.addWidget(label, row, 0, 1, 1)

                    le = LineEdit()
                    #le.setMaximumWidth(150) # matches ui file
                    le.key = key
                    le.args = [p]
                    le.dtype = float
                    if 'young' in key:
                        le.setValInfo(min=0.0, exclude_min=True)
                    #le.setObjectName('lineedit_keyword_%s_args_%s' % (key, p))
                    #setattr(ui, le.objectName(), le)
                    self.add_tooltip(le, key)
                    layout.addWidget(le, row, 1, 1, 1)
                    val = self.project.get_value(key, args=[p])
                    if val is not None:
                        le.updateValue(key, val)
                    self.project.register_widget(le, keys=[key], args=[p])

                    if 'young' in key:
                        label = QLabel('Pa')
                        #label.setObjectName('label_%s_units_args_%s' % (key, p))
                        #setattr(ui, label.objectName(), label)
                        layout.addWidget(label, row, 3, 1, 1)

                # put back 'Wall' row
                key = key[0] + 'w' + key[1:]  # e_young -> ew_young, etc
                row += 1
                name = 'Wall'
                label = QLabel(name)
                self.add_tooltip(label, key)
                layout.addWidget(label, row, 0, 1, 1)
                le = LineEdit()
                #le.setMaximumWidth(150) # matches ui file
                le.key = key
                le.args = None
                le.dtype = float
                if 'young' in key:
                    le.setValInfo(min=0.0, exclude_min=True)
                #le.setObjectName('lineedit_keyword_%s_args' % key)
                #setattr(ui, le.objectName(), le)
                self.add_tooltip(le, key)
                layout.addWidget(le, row, 1, 1, 1)
                val = self.project.get_value(key)
                if val is not None:
                    le.updateValue(key, val)
                self.project.register_widget(le, keys=[key], args=None)
                if 'young' in key:
                    label = QLabel('Pa')
                    #label.setObjectName('label_%s_units_args_%s' % (key, p))
                    #setattr(ui, label.objectName(), label)
                    layout.addWidget(label, row, 3, 1, 1)


        #Specify normal restitution coefficient
        # Specification always required
        # Sets keyword DES_EN_INPUT (DES_EN_WALL_INPUT)
        # Input given as an upper triangular matrix
        mmax = self.project.get_value('mmax', default=len(self.solids))
        tw = ui.tablewidget_des_en_input
        def make_item(str):
            item = QTableWidgetItem(str)
            set_item_noedit(item)
            set_item_enabled(item, False)
            return item

        if (self.solids_dem_saved_solids_names != solids_names
            or tw.rowCount() != mmax+1
            or tw.columnCount() != mmax):

            # Clear out old lineedit widgets
            for row in range(tw.rowCount()):
                for col in range(tw.columnCount()):
                    w = tw.cellWidget(row, col)
                    if w:
                        self.project.unregister_widget(w)
                        w.deleteLater()
            tw.clearContents()

            # Make a new batch
            tw.setRowCount(mmax+1) # extra row for "Wall"
            tw.setColumnCount(mmax)
            tw.setHorizontalHeaderLabels(solids_names)
            tw.setVerticalHeaderLabels(solids_names + ['Wall'])

            arg = 1 # One-based
            key = 'des_en_input'
            for row in range(mmax):
                for col in range(mmax):
                    if col < row:
                        tw.setItem(row, col, make_item('--'))
                    else:
                        le = LineEdit()
                        #le.setMaximumWidth(150)
                        le.key = key
                        le.args = [arg]
                        le.setdtype('dp')
                        self.add_tooltip(le, key)
                        tw.setCellWidget(row, col, le)
                        val = self.project.get_value(key, args=[arg])
                        if val is not None:
                            le.updateValue(key, val)
                        self.project.register_widget(le, keys=[key], args=[arg])
                        arg += 1
            arg = 1
            key = 'des_en_wall_input'
            row = mmax
            for col in range(mmax):
                le = LineEdit()
                #le.setMaximumWidth(150)
                le.key = key
                le.args = [arg]
                le.setdtype('dp')
                self.add_tooltip(le, key)
                tw.setCellWidget(row, col, le)
                val = self.project.get_value(key, args=[arg])
                if val is not None:
                    le.updateValue(key, val)
                self.project.register_widget(le, keys=[key], args=[arg])
                arg += 1

        self.fixup_solids_table(tw, stretch_column=mmax-1)
        # This makes the table look a little nicer
        tw.setShowGrid(False)
        # Move column headers to left so they line up with lineedits
        for i in range(tw.columnCount()):
            item = tw.horizontalHeaderItem(i)
            if item:
                item.setTextAlignment(Qt.AlignLeft)

        #Specify tangential restitution coefficient
        # Specification available for Hertzian collision model
        # Sets keyword DES_ET_INPUT (DES_ET_WALL_INPUT)
        # Input given as an upper triangular matrix
        enabled = (des_coll_model=='HERTZIAN')
        self.set_widget_enabled(ui.label_des_et_input, enabled, reason="Requires Hertzian collision model.")
        tw = ui.tablewidget_des_et_input
        # note - this is too much of a duplicate of des_en_input above
        if not enabled:
            # Clear out old lineedit widgets
            for row in range(tw.rowCount()):
                for col in range(tw.columnCount()):
                    w = tw.cellWidget(row, col)
                    if w:
                        self.project.unregister_widget(w)
                        w.deleteLater()
            tw.clearContents()
            tw.setRowCount(0)
            tw.setColumnCount(0)

        if enabled:
            if (self.solids_dem_saved_solids_names != solids_names
                or tw.rowCount() != mmax+1
                or tw.columnCount() != mmax):

                # Clear out old lineedit widgets
                for row in range(tw.rowCount()):
                    for col in range(tw.columnCount()):
                        w = tw.cellWidget(row, col)
                        if w:
                            self.project.unregister_widget(w)
                            w.deleteLater()
                tw.clearContents()
                # Make a new batch
                tw.setRowCount(mmax+1) # extra row for "Wall"
                tw.setColumnCount(mmax)
                tw.setHorizontalHeaderLabels(solids_names)
                tw.setVerticalHeaderLabels(solids_names + ['Wall'])

                arg = 1
                key = 'des_et_input'
                for row in range(mmax):
                    for col in range(mmax):
                        if col < row:
                            tw.setItem(row, col, make_item('--'))
                        else:
                            le = LineEdit()
                            #le.setMaximumWidth(150)
                            le.key = key
                            le.args = [arg]
                            le.setdtype('dp')
                            self.add_tooltip(le, key)
                            tw.setCellWidget(row, col, le)
                            val = self.project.get_value(key, args=[arg])
                            if val is not None:
                                le.updateValue(key, val)
                            self.project.register_widget(le, keys=[key], args=[arg])
                            arg += 1
                key = 'des_et_wall_input'
                row = mmax
                arg = 1
                for col in range(mmax):
                    le = LineEdit()
                    #le.setMaximumWidth(150)
                    le.key = key
                    le.args = [arg]
                    le.setdtype('dp')
                    tw.setCellWidget(row, col, le)
                    val = self.project.get_value(key, args=[arg])
                    if val is not None:
                        le.updateValue(key, val)
                    self.project.register_widget(le, keys=[key], args=[arg])
                    arg += 1
        self.fixup_solids_table(tw, stretch_column=mmax-1)
        # This makes the table look a little nicer
        tw.setShowGrid(False)
        # Move column headers to left so they line up with lineedits
        for i in range(tw.columnCount()):
            item = tw.horizontalHeaderItem(i)
            if item:
                item.setTextAlignment(Qt.AlignLeft)

        #Select cohesion model
        # Selection always available
        # Available selections
        #None [DEFAULT]
        #Selection always available
        #Sets keyword USE_COHESION to false
        #Sets keyword VAN_DER_WAALS to false
        #Van der Waals
        #Selection always available
        #Sets keyword USE_COHESION to true
        #Sets keyword VAN_DER_WAALS to true
        use_cohesion = self.project.get_value('use_cohesion')
        van_der_waals = self.project.get_value('van_der_waals')
        cb = ui.combobox_cohesion_model
        if use_cohesion:
            if not van_der_waals:
                self.warn('inconsistent value for keyword van_der_waals')
                self.unset_keyword('van_der_waals')
            cb.setCurrentIndex(1)
        else:
            if van_der_waals:
                self.warn('inconsistent value for keyword van_der_waals')
                self.update_keyword('van_der_waals', True)
            cb.setCurrentIndex(0)

        #Specify Hamaker constant
        # Specification only available for Van der Waals cohesion model
        # Sets keyword HAMAKER_CONSTANT (WALL_HAMAKER_CONSTANT)
        #Specify outer cutoff
        # Specification only available for Van der Waals cohesion model
        # Sets keyword VDW_OUTER_CUTOFF (WALL_OUTER_CUTOFF)
        #Specify inner cutoff
        # Specification only available for Van der Waals cohesion model
        # Sets keyword VDW_INNER_CUTOFF (WALL_INNER_CUTOFF)
        #Specify asperities
        # Specification only available for Van der Waals cohesion model
        # Sets keyword ASPERITIES
        enabled = bool(van_der_waals)
        self.set_widget_enabled(ui.groupbox_cohesion_parameters, enabled, reason="Requires Van der Walls cohesion model.")
        # (settings handled by keyword widgets.  TODO:
        #  decide if we want to unset keywords if not enabled

        #List the following options under an 'Advanced' section header.
        #Select Neighbor Search Method
        # Selection always available
        # Available selection
        #Grid-based [DEFAULT]
        #Selection always available
        #Sets keyword DES_NEIGHBOR_SEARCH 4
        #N-Square
        #Selection always available
        #Sets keyword DES_NEIGHBOR_SEARCH 1
        des_neighbor_search = self.project.get_value('des_neighbor_search', default=4)
        if des_neighbor_search not in (1, 4):
            self.warn("Invalid des_neighbor_search %s" % des_neighbor_search)
            des_neighbor_search = 4
            self.update_keyword('des_neighbor_search', des_neighbor_search)
        cb = ui.combobox_des_neighbor_search
        cb.setCurrentIndex(0 if des_neighbor_search==4 else 1)

        #Specify maximum steps between neighbor search
        #Specification always available
        # Sets keyword NEIGHBOR_SEARCH_N
        #Specify factor defining particle neighborhood
        #Specification always available
        # Sets keyword FACTOR_RLM
        #Specify neighborhood search radius ratio
        #Specification always available
        #Sets keyword NEIGHBOR_SEARCH_RAD_RATIO
        #Specify search grid partitions (optional)
        #Specification always available
        #Sets keyword DESGRIDSEARCH_IMAX
        #Sets keyword DESGRIDSEARCH_JMAX
        #Sets keyword DESGRIDSEARCH_KMAX
        pass # handled by keyword widgets


        #Define minimum distance for contact conduction (optional)
        #Unavailable if not solving energy equations
        #Define fluid lens proportionality constant (optional)
        # Unavailable if not solving energy equations
        enabled = self.project.get_value('energy_eq', default=True)
        self.set_widget_enabled(ui.groupbox_des_corrections, enabled, reason="Requires energy equations.")

        # Remember the names of solids phases, to track changes
        self.solids_dem_saved_solids_names = solids_names

        self.handle_solids_psd_selection() # button state
        # Fin!

    def solids_handle_groupbox_filter_particles(self, idx, use_retained=True):
        ui = self.ui.solids
        ui.groupbox_filter_particles_checked = idx
        for w in widget_iter(ui.groupbox_filter_particles):
            if isinstance(w, (QComboBox, QLabel, QLineEdit)):
                w.setVisible(idx)
                if isinstance(w, QLabel):
                    continue
                if not idx:
                    if use_retained:
                        val = self.project.get_value(w.key, args=w.args)
                        if val is not None:
                            self.retain_keyword(w.key, args=w.args)
                        else:
                            self.clear_retained_keyword(w.key, args=w.args)
                        self.unset_keyword(w.key, args=w.args)
                else:
                    if use_retained:
                        val = self.get_retained_keyword(w.key, args=w.args)
                        if val is not None:
                            self.update_keyword(w.key, val, args=w.args)
                    default=self.keyword_doc.get(w.key, {}).get('initpython')
                    val = self.project.get_value(w.key, args=w.args, default=default)
                    if isinstance(w, QComboBox):
                        w.setCurrentIndex(int(val != default))
                    elif isinstance(w, LineEdit):
                        w.updateValue(w.key, val, args=w.args)

    def solids_handle_particles(self):
        # Don't want to restore from retained keys if user cleared this intentionally
        val = self.project.get_value('particles')
        if val is None:
            self.clear_retained_keyword('particles')
        self.update_particle_preview()

    def solids_handle_part_in_exclude(self, w, idx):
        key = w.key
        args = getattr(w, 'args', None)
        if idx == 0: # Include
            self.unset_keyword(key, args=args) # Don't save default False
            #self.update_keyword(key, False, args=args)
        else:
            self.update_keyword(key, True, args=args)
        self.update_particle_preview()

    def solids_handle_part_in_phase(self, w, idx):
        key = w.key
        args = getattr(w, 'args', None)
        if idx: # Exclude
            self.update_keyword(key, False, args=args)
        else: # Don't save default True
            #self.update_keyword(key, True, args=args)
            self.clear_retained_keyword(key, args=args)
            self.unset_keyword(key, args=args)
        self.update_particle_preview()

    def solids_handle_combobox_time_step_option(self, idx):
        ui = self.ui.solids
        if idx == 0: # Automatic
            self.unset_keyword('dtsolid')
            ui.lineedit_keyword_dtsolid.required = False
            ui.lineedit_keyword_dtsolid.updateValue(None, None)
        else:
            ui.lineedit_keyword_dtsolid.required = True
        for w in (ui.lineedit_keyword_dtsolid,
                  ui.label_dtsolid_units):
            self.set_widget_enabled(w, bool(idx==1)) # Constant
        for w in (ui.label_dtsolid_factor,
                  ui.lineedit_keyword_dtsolid_factor,
                  ui.label_dtsolid_update_dt,
                  ui.lineedit_keyword_dtsolid_update_dt,
                  ui.label_dtsolid_update_dt_units):
            self.set_widget_enabled(w, bool(idx==0))
        set_combobox_tooltip(ui.combobox_time_step_option)


    def solids_handle_combobox_gsp_model(self, idx):
        ui = self.ui.solids
        self.update_keyword("gsp_model", value="EXPLICIT" if idx==0 else "IMPLICIT")
        set_combobox_tooltip(ui.combobox_gsp_model)

    def solids_psd_add(self):
        p = self.psd_popup
        p.clear_inputs()
        p.reset_signals()
        p.cancel.connect(self.solids_psd_revert)
        p.editing = False
        p.reserved_aliases.clear()
        p.reserved_aliases.update(self.psd.keys())
        p.save.connect(self.solids_psd_save)
        p.popup()


    def solids_psd_delete(self):
        # Find users of psd
        ui = self.ui.solids
        p = self.psd_popup
        tw = ui.tablewidget_psd
        row = get_selected_row(tw)
        if row is None:
            return
        name = tw.item(row, 0).text()
        nrefs = 0
        for ic_data in self.ics.values():
            # count() does not work on dict_values
            nrefs += list(ic_data.get('psd',{}).values()).count(name)
        for bc_data in self.bcs.values():
            nrefs += list(bc_data.get('psd',{}).values()).count(name)
        if nrefs:
            ret = self.message(text="%s has %d reference%s.\nAll references will be deleted\nContinue?" %
                               (name,
                                nrefs,
                                's' if nrefs!=1 else ''),
                               icon='question',
                              buttons=['ok', 'cancel'])
            if ret != 'ok':
                self.print_internal("Not deleting %s" % name)
                return
        self.psd.pop(name, None)
        self.update_solids_psd_table()
        for IC, ic_data in self.ics.items():
            for (P, psd) in list(ic_data.get('psd', {}).items()):
                if psd == name:
                    for k in ('ic_psd_type', 'ic_psd_mean_dp',
                              'ic_psd_stdev', 'ic_psd_min_dp',
                              'ic_psd_max_dp'):
                        self.unset_keyword(k, args=[IC,P])
                    ic_data['psd'].pop(P)
        for BC, bc_data in self.bcs.items():
            for (P, psd) in list(bc_data.get('psd', {}).items()):
                if psd == name:
                    for k in ('bc_psd_type', 'bc_psd_mean_dp',
                              'bc_psd_stdev', 'bc_psd_min_dp',
                              'bc_psd_max_dp'):
                        self.unset_keyword(k, args=[BC,P])
                    bc_data['psd'].pop(P)

        self.fixup_solids_table(ui.tablewidget_psd)


    def solids_psd_edit(self, col=None):
        ui = self.ui.solids
        tw = ui.tablewidget_psd
        row = get_selected_row(tw)
        if row is None:
            return
        name = tw.item(row, 0).text()
        p = self.psd_popup
        p.editing = True
        p.reserved_aliases.clear()
        p.reserved_aliases.update(self.psd.keys())
        p.reserved_aliases.discard(name)

        p.ui.lineedit_alias.setText(name)
        data = self.psd.get(name, {})
        t = data.get('type')
        if t == 'Normal':
            idx = 0
        elif t == 'Log-normal':
            idx = 1
        elif t == 'Custom':
            idx = 2
        else:
            self.warning("Undefined PSD type '%s'" % t, popup=True)
            idx = 0
        p.ui.combobox_type.setCurrentIndex(idx)

        if t == 'Custom':
            p.filename = data.get('filename')
            p.label_filename.setText(p.filename or '')
        else:
            for x in ('mean', 'sigma', 'min', 'max'):
                getattr(p.ui, 'lineedit_%s'%x).updateValue(None, data.get(x))
        plot_range =  data.get('plot_range', [0,0.1])
        p.ui.lineedit_plot_range_min.updateValue(None, plot_range[0])
        p.ui.lineedit_plot_range_max.updateValue(None, plot_range[1])

        # Is this needed?
        p.reset_signals()
        p.cancel.connect(self.solids_psd_revert)
        p.save.connect(self.solids_psd_save)
        p.popup()
        p.handle_type(idx)

        if col is not None:
            w = (p.ui.lineedit_alias, p.ui.combobox_type,
                 p.ui.lineedit_mean, p.ui.lineedit_sigma,
                 p.ui.lineedit_min, p.ui.lineedit_max)[col]
            w.setFocus(0)
            if col != 1:
                w.selectAll()


    def solids_psd_save(self):
        ui = self.ui.solids
        tw = ui.tablewidget_psd
        p = self.psd_popup
        data = {}
        name = p.ui.lineedit_alias.text().strip()
        if not name:
            return
        row = get_selected_row(tw)

        if p.editing and row is not None:
            old_name = tw.item(row,0).text()
            if name != old_name:
                # Name change
                del self.psd[old_name]
                for xcs in (self.ics.values(), self.bcs.values()):
                    for xc in xcs:
                        psd_data = xc.get('psd', {})
                        for (k,v) in psd_data.items():
                            if v == old_name:
                                psd_data[k] = name

        idx = p.ui.combobox_type.currentIndex()
        if idx == 0: #Normal
            data['type'] = 'Normal'
        elif idx == 1:
            data['type'] = 'Log-normal'
        else:
            data['type'] = 'Custom'
            if (os.path.isabs(p.filename)
                and not shutil._samefile(os.path.dirname(p.filename),os.getcwd())
                and p.ui.checkbox_copy_file.isChecked()):
                basename = os.path.basename(p.filename)
                shutil.copyfile(p.filename, basename)
                p.filename = basename
            data['filename'] = p.filename

        if idx in [0,1]:
            data['mean'] = p.ui.lineedit_mean.value
            data['sigma'] = p.ui.lineedit_sigma.value
            data['min'] = p.ui.lineedit_min.value
            data['max'] = p.ui.lineedit_max.value
            data['plot_range'] = [p.ui.lineedit_plot_range_min.value,
                                  p.ui.lineedit_plot_range_max.value]
        else:
            if p.custom_data:
                X, pdf, cdf, type_, µ, 𝜎 = p.custom_data
                data['mean'] = µ
                data['sigma'] = 𝜎
                data['min'] = X[0]
                data['max'] = X[-1]

        if name in self.psd and self.psd[name] == data:
            pass # No change
        else:
            self.psd[name] = data
            self.set_unsaved_flag()

            #update IC and BC keys
            for (x, xcs) in ('ic', self.ics), ('bc', self.bcs):
                for XC, xc_data in xcs.items():
                    psd_data = xc_data.get('psd')
                    if not psd_data:
                        continue
                    for (P, psd) in psd_data.items():
                        if psd == name:
                            self.update_keyword(x+'_psd_type',
                                                data['type'].replace('-','_').upper(), args=[XC,P])
                            if idx == 2:
                                data = {} # data has vals for custom PSD, but we don't want to set keys
                            self.update_keyword(x+'_psd_mean_dp', data.get('mean'), args=[XC,P])
                            self.update_keyword(x+'_psd_min_dp', data.get('min'), args=[XC,P])
                            self.update_keyword(x+'_psd_max_dp', data.get('max'), args=[XC,P])
                            self.update_keyword(x+'_psd_stdev', data.get('sigma'), args=[XC,P])

        self.update_solids_psd_table()
        self.fixup_solids_table(ui.tablewidget_psd)


    def solids_psd_revert(self):
        pass

    def update_solids_psd_table(self):
        ui = self.ui.solids
        tw = ui.tablewidget_psd
        tw.setRowCount(len(self.psd))
        def make_item(val):
            if isinstance(val, float):
                val = round(val, 6)
            item = QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item

        for (row,(name, data)) in enumerate(self.psd.items()):
            tw.setItem(row, 0, make_item(name))
            psd_type = data.get('type', '???')
            tw.setItem(row, 1, make_item(psd_type))
            tw.setItem(row, 2, make_item(data.get('mean')))
            tw.setItem(row, 3, make_item(data.get('sigma')))
            tw.setItem(row, 4, make_item(data.get('min')))
            tw.setItem(row, 5, make_item(data.get('max')))
