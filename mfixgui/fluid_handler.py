from copy import deepcopy

from qtpy.QtWidgets import QTableWidgetItem, QHeaderView, QLabel
from qtpy.QtCore import QTimer, Qt

from mfixgui import default_values
from mfixgui.constants import *
from mfixgui.species_handler import SpeciesHandler
from mfixgui.tools import format_key_with_args, keyword_args
from mfixgui.tools.qt import (clear_layout, get_combobox_item,
                              get_selected_row, set_combobox_tooltip,
                              set_item_noedit)
from mfixgui.tools.cantera_poly import create_cantera_polynomials
from mfixgui.widgets.base import LineEdit

cantera_keys = ('alpha_transport', 'config_transport', 'mu_transport', 'zrot_transport')
lennard_jones_keys = ('ljeps','ljsig')

class FluidHandler(SpeciesHandler):
    # Defaults
    def init_fluid_default_models(self):
        ui = self.ui.fluid
        self.fluid_density_model = IDEAL
        self.fluid_viscosity_model = SUTHERLAND # Sutherland
        self.fluid_mol_weight_model = MIXTURE
        self.fluid_specific_heat_model = MIXTURE
        self.fluid_kg_model = AIR
        self.fluid_diffusion_model = AIR
        for s in ('density', 'viscosity', 'mol_weight',
                  'specific_heat', 'kg', 'diffusion'):
            name = 'fluid_%s_model' % s
            getattr(ui, 'combobox_'+name).default(getattr(self, name))


    ## Fluid phase methods
    def handle_fluid_species_eq(self, enabled):
        ui = self.ui.fluid
        self.update_keyword('species_eq', enabled, args=[0])

        ui.frame_add_delete_copy_species.setVisible(enabled or bool(self.fluid_species))
        for item in (ui.combobox_fluid_diffusion_model,
                     ui.label_fluid_diffusion_model,
                     # more ?
                     ):
            if enabled:
                self.enable_widget(item)
            else:
                self.disable_widget(item, reason="Requires fluid species equations.")

        # dif_g0 == diffusion coeff model
        items = (ui.lineedit_keyword_dif_g0,
                 ui.label_dif_g0_units)
        for item in items:
            if enabled and (self.fluid_diffusion_model == CONSTANT):
                self.enable_widget(item)
            else:
                self.disable_widget(item, reason="Requires fluid species equations and constant diffusion model.")

        # avoid warnings from check_gas_phase.f
        item = get_combobox_item(ui.combobox_fluid_mol_weight_model, 0)  # CONSTANT
        if enabled: #species_eq
            #! MW_AVG is defined and the gas phase species equations are solved, then
            #! the user specified average molecular weight is ignored.
            self.set_fluid_mol_weight_model(MIXTURE)
            self.disable_widget(item, reason="Not available with fluid species equations.")
        else:
            ro_g0 = self.project.get_value('ro_g0')
            if ro_g0 is None:
                self.enable_widget(item)
            else:
                self.set_fluid_mol_weight_model(MIXTURE)
                self.disable_widget(item, reason="Not compatible with specified constant gas density.")


    def handle_fluid_momentum_eq(self, key, enabled):
        self.update_keyword(key, enabled, args=[0])


    def init_fluid_handler(self):
        self.fluid_species = {} # keyed by alias
        ui = self.ui.fluid
        cb = ui.checkbox_species_eq

        ui.lineedit_fluid_phase_name.default_value = self.fluid_phase_name = "Fluid"
        self.init_fluid_default_models()
        self.fluid_species_warnings_shown = set()

        # Handle a number of cases which are essentially the same
        #  (Note, fewer cases are the same now, see issues/1890)
        # see 'set_fluid_mol_weight_model' below to help understand this
        def make_fluid_model_setter(self, name, key):
            def setter(model, in_setup=False):
                ui = self.ui.fluid
                key_g0 = 'c_pg0' if key=='c_p' else key + '_g0'
                key_usr = 'usr_cpg' if key=='c_p' else 'usr_' + key + 'g'

                setattr(self, name, model) # self.fluid_<name>_model = model
                cb = getattr(ui, 'combobox_' + name)
                prev_model = cb.currentIndex()
                if model != prev_model:
                    cb.setCurrentIndex(model)
                # Make tooltip match item setting
                cb.setToolTip(get_combobox_item(cb, model).toolTip())
                #
                # Make combobox locatable
                cb.keys = [key_g0, key_usr]
                label = getattr(ui, 'label_' + name)
                label.keys = cb.keys
                # Enable lineedit for constant model
                lineedit = getattr(ui, 'lineedit_keyword_%s' % key_g0)
                label = getattr(ui, 'label_%s_units' % key_g0)

                for item in (lineedit, label):
                    if model==CONSTANT:
                        self.enable_widget(item)
                    else:
                        self.disable_widget(item, reason="Requires constant %s model."%name)

                # Workaround for disabled fluid solver
                if self.fluid_solver_disabled and key_g0 == 'ro_g0':
                    lineedit.minimum = 0
                    lineedit.saved_value = 0.0
                    return

                if model == CONSTANT:
                    if self.project.get_value(key_g0) is None:
                        value = self.get_retained_keyword(key_g0, default=getattr(default_values, key_g0, None))
                        if value != '' and value is not None:
                            self.update_keyword(key_g0, value) # Restore keyword value
                    self.unset_keyword(key_usr) # Issues/435
                elif model == UDF:
                    self.retain_keyword(key_g0)
                    self.unset_keyword(key_g0)
                    self.update_keyword(key_usr, True)
                    lineedit.updateValue(key_g0, None)
                else: # Ideal gas law, etc
                    self.retain_keyword(key_g0)
                    self.unset_keyword(key_g0)
                    self.unset_keyword(key_usr)
                    lineedit.updateValue(key_g0, None)

                if key_g0 == 'ro_g0':
                    # avoid warnings from check_gas_phase.f
                    const_mw_item = get_combobox_item(ui.combobox_fluid_mol_weight_model, 0)  # CONSTANT
                    species_eq = self.project.get_value('species_eq', args=[0], default=True)
                    if model == CONSTANT:
                        lineedit.required = True
                        lineedit.exclude_min = True
                        self.disable_widget(const_mw_item, reason="Not available with constant fluid density model.")
                        self.set_fluid_mol_weight_model(MIXTURE)
                    else:
                        lineedit.required = False
                        lineedit.exclude_min = False
                        if species_eq:
                            self.enable_widget(const_mw_item)
                        else:
                            self.disable_widget(const_mw_item, reason="Requires fluid species equations.")
                        self.set_fluid_mol_weight_model(CONSTANT) # ??XXXX
                    self.setup_fluid_diffusivity() # requires ideal gas law
                # anything else to do in this case? validation?
            return setter

        # Create setters for the cases which are similar (mol. wt. handled separately)
        # (This is getting less and less generic since the first implementation, see
        #  mol_weight, viscosity, and diffusion below)
        for (name, key) in (
                ('density', 'ro'),
                ('specific_heat', 'c_p')):
            model_name = 'fluid_%s_model' % name
            setattr(self, 'set_'+model_name, make_fluid_model_setter(self, model_name, key))

            # Set the combobox default value (?)
            cb = getattr(ui, 'combobox_'+model_name)
            cb.default_value = getattr(self, model_name)

            # Tooltips
            key_g0 = 'c_pg0' if key=='c_p' else key + '_g0'
            key_usr = 'usr_cpg' if key=='c_p' else 'usr_' + key + 'g'
            item =  get_combobox_item(cb, 0)
            self.add_tooltip(item, key=key_g0)
            item =  get_combobox_item(cb, 1)
            item.setToolTip('<b>%s</b>: Use MFiX default calculation.' % item.text()) #Full name of model
            item =  get_combobox_item(cb, 2)
            self.add_tooltip(item, key=key_usr)

        ui.combobox_fluid_density_model.keys = ['ro_g0', 'usr_rog']
        ui.label_fluid_density_model.keys =  ['ro_g0', 'usr_rog']
        ui.combobox_fluid_specific_heat_model.keys = ['c_pg0', 'usr_cpg']
        ui.label_fluid_specific_heat_model.keys = ['c_pg0', 'usr_cpg']

        # Viscosity
        cb = ui.combobox_fluid_viscosity_model
        cb.default_value = self.fluid_viscosity_model
        cb.keys = ['mu_g_model', 'mu_g0', 'sl_muref', 'sl_tref', 'sl_s',
                   'hb_tau0', 'hb_k0', 'hb_n', 'hb_gama_c'] # locatability
        key = 'mu_g_model'
        for (i,v) in enumerate(VISCOSITY_MODELS):
            self.add_tooltip(get_combobox_item(cb,i), key=key, value=v)
        set_combobox_tooltip(cb)

        # Mol weight
        cb = ui.combobox_fluid_mol_weight_model
        cb.default_value = self.fluid_mol_weight_model
        item = get_combobox_item(cb, 0)
        self.add_tooltip(item, key='mw_avg')
        item = get_combobox_item(cb, 1)
        item.setToolTip("<b>Mixture</b>: Use MFiX default calculation.")
        set_combobox_tooltip(cb)

        # Specific heat
        cb = ui.combobox_fluid_specific_heat_model
        set_combobox_tooltip(cb)

        # Thermal conductivity
        cb = ui.combobox_fluid_kg_model
        key = 'kg_model'
        self.add_tooltip(cb, key)
        for i,val in enumerate(KG_MODELS):
            self.add_tooltip(get_combobox_item(cb,i), key, value=val)
        set_combobox_tooltip(cb)

        # Diffusion
        cb = ui.combobox_fluid_diffusion_model
        cb.default_value = self.fluid_diffusion_model
        cb.keys = ['dif_g0', 'multi_component_diffusion', 'usr_difg'] # locatability
        self.add_tooltip(get_combobox_item(cb,0), key='dif_g0')
        item = get_combobox_item(cb,1)
        item.setToolTip('<b>%s</b>: Use MFiX default calculation.' % item.text()) #Full name of model
        self.add_tooltip(get_combobox_item(cb,2), key='multi_component_diffusion')
        self.add_tooltip(get_combobox_item(cb,3), key='usr_difg')
        set_combobox_tooltip(cb)

        cb = ui.combobox_dif_coeff_kt
        item = get_combobox_item(cb,0)
        self.add_tooltip(item, key='dif_coeff_kt', value=False)
        item = get_combobox_item(cb,1)
        self.add_tooltip(item, key='dif_coeff_kt', value=True)
        set_combobox_tooltip(cb)
        ui.combobox_dif_coeff_kt.currentIndexChanged.connect(self.handle_dif_coeff_kt)

        ui.lineedit_fluid_phase_name.value_updated.connect(
            self.handle_fluid_phase_name)

        cb = ui.checkbox_species_eq
        self.add_tooltip(cb, key='species_eq')
        cb.clicked.connect(self.handle_fluid_species_eq)

        for c in 'xyz':
            key = 'momentum_%s_eq' % c
            cb = getattr(ui, 'checkbox_'+key)
            self.add_tooltip(cb, key)
            cb.clicked.connect(lambda enabled, key=key: self.handle_fluid_momentum_eq(key, enabled))

        # Fluid phase models
        for name in ('density', 'viscosity', 'specific_heat', 'mol_weight',
                     'kg', 'diffusion'):
            model_name = 'fluid_%s_model' % name
            cb = getattr(ui, 'combobox_%s' % model_name)
            setter = getattr(self,'set_%s' % model_name)
            cb.currentIndexChanged.connect(setter)

        # Fluid species
        tb = ui.toolbutton_fluid_species_add
        tb.key = 'nmax_g' # locatability
        tb.clicked.connect(self.fluid_species_add)
        tb = ui.toolbutton_fluid_species_edit
        tb.clicked.connect(self.fluid_species_edit)
        self.disable_widget(tb)
        tb = ui.toolbutton_fluid_species_delete
        tb.key = 'nmax_g' # locatability
        self.disable_widget(tb)
        tb.clicked.connect(self.fluid_species_delete)
        tw = ui.tablewidget_fluid_species
        tw.itemSelectionChanged.connect(self.handle_fluid_species_selection)
        self.fixup_fluid_table()

        # Locatability
        # Note, with find-by-args this will have to select a table row, not the
        # groupbox
        # Should nmax_g be on +/- buttons or on the groupbox?
        ui.groupbox_species.keys = ['species_alias_g', 'species_g'] # Locatability
        #ui.lineedit_keyword_ro_g0.required = True # Only required when constant density
        #ui.widget_binary.keys = ['dabg']
        ui.widget_multicomponent_diffusion.keys = ['dabg',
                                                   'dif_thermal']
        ui.widget_multicomponent_diffusion.hidden_ctrl = ui.combobox_fluid_diffusion_model
        ui.widget_binary.hidden_ctrl = [ui.combobox_fluid_diffusion_model,
                                        ui.combobox_dif_coeff_kt]


    def fixup_fluid_table(self):
        hv = QHeaderView
        ui = self.ui.fluid
        tw = ui.tablewidget_fluid_species
        resize = tw.horizontalHeader().setSectionResizeMode
        for n in range(tw.columnCount()):
            resize(n, hv.ResizeToContents if n>0
                   else hv.Stretch)


    def set_fluid_kg_model(self, model, in_setup=False):
        ui = self.ui.fluid
        cb = ui.combobox_fluid_kg_model

        if in_setup:
            cb.setCurrentIndex(model)
        set_combobox_tooltip(cb)
        for w in ui.lineedit_keyword_k_g0, ui.label_k_g0_units:
            self.set_widget_enabled(w, model==CONSTANT)
        if in_setup:
            return
        self.update_keyword('kg_model', KG_MODELS[model])
        if model == CONSTANT:
            k_g0 = self.project.get_value('k_g0', default=None)
            if k_g0 is None:
                k_g0 = self.get_retained_keyword('k_g0')
                if k_g0 is not None:
                    self.update_keyword('k_g0', k_g0)
        else:
            self.retain_keyword('k_g0')
            self.unset_keyword('k_g0')
            ui.lineedit_keyword_k_g0.setText('')

        self.species_popup.enable_lennard_jones(self.lennard_jones_required())
        self.species_popup.enable_cantera(self.cantera_required())
        self.species_popup.check_data()
        self.update_cantera_polynomials()


    # molecular wt model only has 2 choices, and the key names don't
    # follow the same pattern, so create its setter specially
    def set_fluid_mol_weight_model(self, model):
        ui = self.ui.fluid
        self.fluid_mol_weight_model = model
        cb = ui.combobox_fluid_mol_weight_model
        # Make tooltip match setting (for longer names which are truncated)
        cb.setToolTip(get_combobox_item(cb, cb.currentIndex()).toolTip())
        prev_model = cb.currentIndex()
        if model != prev_model:
            cb.setCurrentIndex(model)
        # Enable lineedit for constant mol_weight model
        lineedit = ui.lineedit_keyword_mw_avg
        label = ui.label_mw_avg_units
        species_eq = self.project.get_value('species_eq', args=[0], default=True)
        # avoid warnings from check_gas_phase.f
        for item in (lineedit, label):
            item.setEnabled(model==CONSTANT and not species_eq)
        if model == CONSTANT and not species_eq:
            key = 'mw_avg'
            if self.project.get_value(key) is None:
                value = self.get_retained_keyword(key, default=default_values.mw_avg)
                self.update_keyword("mw_avg", value) # Restore keyword value
        else: # Mixture
            # TODO: validate, require mw for all component species
            self.retain_keyword("mw_avg")
            self.unset_keyword("mw_avg")

    # viscosity has some extra choices, and does not
    # follow the same pattern, so create its setter specially
    def set_fluid_viscosity_model(self, model, in_setup=False):
        ui = self.ui.fluid
        self.fluid_viscosity_model = model
        cb = ui.combobox_fluid_viscosity_model
        # Make tooltip match setting (for longer names which are truncated)
        cb.setToolTip(get_combobox_item(cb, cb.currentIndex()).toolTip())
        prev_model = cb.currentIndex()
        if model != prev_model:
            cb.setCurrentIndex(model)
        self.update_keyword('mu_g_model', VISCOSITY_MODELS[model])
        constant_items = [ui.lineedit_keyword_mu_g0, ui.label_mu_g0_units]
        for item in constant_items:
            if model==CONSTANT:
                self.enable_widget(item)
            else:
                self.disable_widget(item, reason="Requires constant fluid viscosity model.")
        key = 'mu_g0'
        if model != CONSTANT:
            self.retain_keyword(key)
            self.unset_keyword(key)
            ui.lineedit_keyword_mu_g0.setText('')
        else:
            val = self.get_retained_keyword(key)
            if val is not None:
                self.update_keyword(key, val)
                ui.lineedit_keyword_mu_g0.updateValue(key, val)

        sutherland_keys = ('sl_muref', 'sl_tref', 'sl_s')
        sutherland_items = [getattr(ui, pat%k)
                            for pat in ('label_%s', 'label_%s_units', 'lineedit_keyword_%s')
                            for k in sutherland_keys]
        for item in sutherland_items:
            item.setVisible(model==SUTHERLAND)
        hb_keys = ('hb_tau0', 'hb_k0', 'hb_n', 'hb_gama_c')
        hb_items = [getattr(ui, pat%k)
                    for pat in ('label_%s', 'label_%s_units', 'lineedit_keyword_%s')
                    for k in hb_keys]
        for item in hb_items:
            item.setVisible(model==HERSCHEL_BULKLEY)

        self.species_popup.enable_lennard_jones(self.lennard_jones_required())
        self.species_popup.enable_cantera(self.cantera_required())
        self.species_popup.check_data()
        self.update_cantera_polynomials()

    def lennard_jones_required(self):
        return bool(self.project.get_value('mu_g_model') in ('LENNARD_JONES', 'CANTERA_POLY')
                    or self.project.get_value('kg_model') in ('LENNARD_JONES', 'CANTERA_POLY')
                    or self.project.get_value('dif_coeff_kt'))

    def cantera_required(self):
        return bool(any(self.project.get_value(key)=='CANTERA_POLY'
                    for key in('mu_g_model', 'kg_model')))

    def check_lennard_jones_cantera(self, force=False):
        msg = ""
        if self.cantera_required():
            if any((self.project.get_value(k ,args=[i]) is None
                    for i in range(1,1+len(self.fluid_species))
                    for k in ('ljeps', 'ljsig', 'config_transport'))):
                msg = "Cantera and Lennard-Jones"
        elif self.lennard_jones_required():
            if any((self.project.get_value(k,args=[i]) is None
                   for i in range(1,1+len(self.fluid_species))
                   for k in ('ljeps', 'ljsig'))):
                msg = "Lennard-Jones"
        else:
                return True
        if not msg:
            return True
        msg += " parameters must be defined for all species. Use species editor to set values."
        self.species_popup.check_data()
        if force or (msg not in self.fluid_species_warnings_shown):
            self.warning(msg, popup=True)
            self.fluid_species_warnings_shown.add(msg)
        return False

    # diffusion has some extra choices, and does not
    # follow the same pattern, so create its setter specially
    def set_fluid_diffusion_model(self, model):
        ui = self.ui.fluid
        self.fluid_diffusion_model = model
        cb = ui.combobox_fluid_diffusion_model
        # Make tooltip match setting (for longer names which are truncated)
        cb.setToolTip(get_combobox_item(cb, cb.currentIndex()).toolTip())
        prev_model = cb.currentIndex()
        if model != prev_model:
            cb.setCurrentIndex(model)
        set_combobox_tooltip(cb)
        constant_items = [ui.lineedit_keyword_dif_g0, ui.label_dif_g0_units]
        for item in constant_items:
            if model==CONSTANT:
                self.enable_widget(item)
            else:
                self.disable_widget(item, reason="Requires constant fluid diffusion model.")
        key = 'dif_g0'
        if model != CONSTANT:
            self.retain_keyword(key)
            self.unset_keyword(key)
            ui.lineedit_keyword_dif_g0.setText('')
        else:
            val = self.get_retained_keyword(key)
            if val is not None:
                self.update_keyword(key, val)
                ui.lineedit_keyword_dif_g0.updateValue(key, val)

        if model==MULTICOMPONENT:
            self.update_keyword('multi_component_diffusion', True)
        else:
            self.unset_keyword('multi_component_diffusion')
        ui.widget_multicomponent_diffusion.setVisible(model==MULTICOMPONENT)

        # NB UDF is 2 in constants.py but we added MULTICOMPONENT
        if model == 3:
            self.update_keyword('usr_difg', True)
        else:
            self.unset_keyword('usr_difg')

        self.setup_fluid_diffusivity(model=model)

    def handle_dif_coeff_kt(self, idx):
        ui = self.ui.fluid
        if idx == 0: # constant binary coeffs
            self.unset_keyword('dif_coeff_kt')
            self.setup_widget_binary()
            ui.widget_binary.setVisible(True)
        else:
            self.update_keyword('dif_coeff_kt', True)
            ui.widget_binary.setVisible(False)
        self.species_popup.enable_lennard_jones(self.lennard_jones_required())
        set_combobox_tooltip(ui.combobox_dif_coeff_kt)
        self.species_popup.check_data()
        self.update_cantera_polynomials()

    def setup_widget_binary(self,force=False):
        ui = self.ui.fluid
        layout = ui.layout_binary
        if len(self.fluid_species) < 2:
            # Shouldn't be here
            return
        if force or layout.columnCount() != len(self.fluid_species)-1 + 2: # Upper tri, plus label and units
            clear_layout(layout)
            names = list(self.fluid_species.keys())
            # First row:  headers
            for col, name in enumerate(names[1:],1):
                label = QLabel(name)
                label.setToolTip(name)
                layout.addWidget(label, 0, col)
            for row, name in enumerate(names[:-1]):
                label = QLabel(name)
                label.setToolTip(name)
                layout.addWidget(label, row+1, 0)
                key = 'dabg'
                for col in range(row+1, len(names)):
                    le = LineEdit()
                    le.key = key
                    le.dtype = float
                    le.allow_parameters = True
                    le.args = [row+1, col+1]
                    le.min = 0
                    le.exclude_min = True
                    le.value_updated.connect(self.project.submit_change)
                    # Move le.updateValue  outside "if layout.columnCount() != " block ?
                    val = self.project.get_value(key, args=[row+1,col+1])
                    le.updateValue(key, val)
                    self.add_tooltip(le, key)
                    layout.addWidget(le, row+1, col)
                unit_label = QLabel("m²/s")
                layout.addWidget(unit_label, row+1, len(names))


    def handle_fluid_phase_name(self, widget, value_dict, args):
        ui = self.ui.fluid
        le = ui.lineedit_fluid_phase_name
        old_name = self.fluid_phase_name
        new_name = le.text()
        if new_name in self.solids: # Reject the input
            self.warning("%s: name is in use" % new_name, popup=True)
            le.setText(old_name)
        else:
            self.set_fluid_phase_name(new_name)


    def set_fluid_phase_name(self, value):
        if value != self.ui.fluid.lineedit_fluid_phase_name.text():
            self.ui.fluid.lineedit_fluid_phase_name.setText(value) # set GUI state
        if self.fluid_phase_name == value:
            return
        self.fluid_phase_name = value
        if self.project.mfix_gui_comments.get('fluid_phase_name') != value:
            self.project.mfix_gui_comments['fluid_phase_name'] = value
            self.set_unsaved_flag()

    def fluid_species_revert(self):
        pass

    def fluid_species_save(self):
        self.set_unsaved_flag()
        old_aliases = dict(enumerate(self.fluid_species.keys(), 1))
        self.set_unsaved_flag() # TODO check if anything really changed
        self.fluid_species = dict((alias, deepcopy(data))
            for (alias,data) in self.species_popup.defined_species.items())
        for k in lennard_jones_keys:
            for (i,name) in enumerate(self.fluid_species,1):
                val = self.fluid_species[name].pop(k, None)
                if val is not None:
                    self.update_keyword(k, val, args=[i])
                else:
                    self.unset_keyword(k, args=[i])
        for k in cantera_keys:
            for (i,name) in enumerate(self.fluid_species,1):
                val = self.fluid_species[name].pop(k, None)
                if val is not None:
                    if k != 'config_transport' and val == 0:
                        self.unset_keyword(k, args=[i])
                    else:
                        self.update_keyword(k, val, args=[i])
                else:
                    self.unset_keyword(k, args=[i])


        self.fluid_normalize_keys()
        self.update_fluid_species_table()
        self.setup_fluid_diffusivity(force=True)
        for (i,name) in enumerate(self.fluid_species.keys(),1):
            old_alias = old_aliases.get(i)
            new_alias = name
            if old_alias is None:
                continue
            if new_alias != old_alias:
                self.chemistry_rename_species(old_alias, new_alias)

        # Update IC_X_G for new species
        for i in self.ics:
            self.ics_set_default_keys(i)
        self.setup_fluid()

        self.bcs_check_wall_keys()
        for BC in self.bcs:
            self.bcs_set_default_keys(BC)
        for IC in self.ics:
            self.ics_set_default_keys(IC)
        for PS in self.pss:
            self.pss_set_default_keys(PS)
        for IS in self.iss:
            self.iss_set_default_keys(IS)

        self.update_cantera_polynomials()

        self.update_nav_tree() # Chemistry

    def update_cantera_polynomials(self):
        if not self.check_lennard_jones_cantera():
            return False

        k = self.project.get_value("kg_model")=="CANTERA_POLY"
        mu = self.project.get_value("mu_g_model")=="CANTERA_POLY"
        if k or mu:
            try:
                for (key, args, val) in create_cantera_polynomials(
                        self.fluid_species_copy_dict(),
                        k=k, mu=mu):

                    self.update_keyword(key, val, args=args)
            except Exception as e:
                self.error(str(e))
                return False

        # Clear unused poly coeffs
        for (var,key) in ((k,'poly_kg'), (mu,'poly_mu_g')):
            if not var:
                for args in self.project.get_key_indices(key):
                    self.unset_keyword(key, args=args)

        return True



    def update_fluid_species_table(self):
        """Update table in fluid pane.  Also set nmax_g, species_g and species_alias_g keywords,
        which are not tied to a single widget"""
        tw = self.ui.fluid.tablewidget_fluid_species
        tw.clearContents()
        if self.fluid_species is None:
            self.fixup_fluid_table()
            return
        nrows = len(self.fluid_species)
        tw.setRowCount(nrows)
        def make_item(val):
            if isinstance(val, float):
                val = round(val, 8)
            item = QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item
        old_nmax_g = self.project.get_value('nmax_g')
        nmax_g = len(self.fluid_species)
        if nmax_g > 0:
            self.update_keyword('nmax_g', nmax_g)
        else:
            self.unset_keyword('nmax_g')
        for (row, (alias,data)) in enumerate(self.fluid_species.items()):
            for (col, key) in enumerate(('alias', 'phase', 'mol_weight', 'h_f')):
                data['alias'] = alias # for make_item
                tw.setItem(row, col, make_item(data.get(key)))

        # Clear any keywords with indices above nmax_g NEEDED?
        if old_nmax_g is None:
            old_nmax_g = 0
        for i in range(nmax_g+1, old_nmax_g+1):
            self.unset_keyword('species_g', args=i)
            self.unset_keyword('species_alias_g', args=i)
        self.fixup_fluid_table()


    def fluid_normalize_keys(self):
        # issues/869
        for (sp, (alias,data)) in enumerate(self.fluid_species.items(), 1):
            species_g = self.species_burcat_name(alias, 0, sp)
            data['species'] = species_g # for thermo_data
            self.update_keyword('species_g',
                                species_g,
                                args=[sp])
            self.update_keyword('species_alias_g', alias, args=[sp])
            # We're avoiding mw_g in favor of the settings in THERMO DATA
            #self.update_keyword('mw_g', data['mol_weight'], args=row+1)#
        self.project.update_thermo_data(self.fluid_species)


    def handle_fluid_species_selection(self):
        ui = self.ui.fluid
        tw = ui.tablewidget_fluid_species
        row = get_selected_row(tw)
        enabled = (row is not None)
        self.set_widget_enabled(ui.toolbutton_fluid_species_delete, enabled)
        self.set_widget_enabled(ui.toolbutton_fluid_species_edit, enabled)
        if enabled:
            #tw.doubleClicked.connect(lambda: (time.sleep(0.25), self.fluid_species_edit()))
            tw.doubleClicked.connect(lambda: QTimer.singleShot(250, self.fluid_species_edit))
        else:
            try:
                tw.doubleClicked.disconnect() #self.fluid_species_edit)
            except:
                pass


    def fluid_species_add(self):
        sp = self.species_popup
        sp.set_phases('GL')
        sp.do_search('') # Init to full db
        # how to avoid this if dialog open already?
        sp.reset_signals()
        sp.cancel.connect(self.fluid_species_revert)
        sp.save.connect(self.fluid_species_save)
        sp.defined_species = self.fluid_species_copy_dict()
        sp.update_defined_species()
        sp.reserved_aliases = set(map(str.lower, self.species_all_aliases()))
        sp.setWindowTitle("Fluid species")
        sp.enable_density(False)
        sp.enable_lennard_jones(self.lennard_jones_required())
        sp.enable_cantera(self.cantera_required())
        sp.popup()


    def fluid_species_delete(self):
        ui = self.ui.fluid
        tw = ui.tablewidget_fluid_species
        row = get_selected_row(tw)
        if row is None: # No selection
            return

        alias = tw.item(row,0).text()

        # Warn if species in use
        refs = self.fluid_get_species_refs(alias)

        if refs:
            # Do we want to show the references?
            ret = self.message(text="%s has %d reference%s.\nAll references will be deleted\nContinue?" %
                               (alias,
                                len(refs),
                                's' if len(refs)!=1 else ''),
                               icon='question',
                              buttons=['ok', 'cancel'])
            if ret != 'ok':
                self.print_internal("Not deleting %s" % alias)
                return

        tw.clearSelection() #?
        species_index = 1 + list(self.fluid_species.keys()).index(alias)
        self.bcs_delete_fluid_species(species_index) # special handling for memoized eq_type
        self.fluid_species.pop(alias, None)
        self.fluid_delete_species_keys(species_index) # Must remove fluid species first (why?)
        self.update_fluid_species_table()
        # Sigh, we have to update the row in the popup too.
        # Should the popup just be modal, to avoid this?
        sp = self.species_popup
        sp.defined_species = deepcopy(self.fluid_species)
        sp.update_defined_species()
        self.chemistry_delete_reactions_of_species(alias)
        self.setup_fluid()
        self.update_nav_tree() # Chemistry


    def fluid_species_copy_dict(self):
        d = deepcopy(self.fluid_species)
        # Add LJ and Cantera keys
        for (i, name) in enumerate(self.fluid_species,1):
            for k in cantera_keys+lennard_jones_keys:
                val = self.project.get_value(k, args=[i], default=None)
                if val is not None:
                    d[name][k] = val
        return d

    def fluid_species_edit(self):
        ui = self.ui.fluid
        if self.lock_level == LOCKED: #prevent double-click from popping us up during run
            return
        tw = ui.tablewidget_fluid_species
        row = get_selected_row(tw)
        sp = self.species_popup
        sp.set_phases('GL')
        sp.reset_signals()
        sp.cancel.connect(self.fluid_species_revert)
        sp.save.connect(self.fluid_species_save)
        sp.defined_species = self.fluid_species_copy_dict()
        sp.update_defined_species()
        sp.reserved_aliases = set(map(str.lower, self.species_all_aliases()))
        if row is not None:
            sp.reserved_aliases.discard(tw.item(row,0).text().lower())
            sp.tablewidget_defined_species.setCurrentCell(row, 0)
            # Initialize search box to current species (?)
            sp.do_search(list(self.fluid_species.keys())[row])
        else:
            sp.do_search('')
            sp.tablewidget_defined_species.clearSelection()

        sp.setWindowTitle("Fluid species")
        sp.enable_density(False)
        sp.enable_lennard_jones(self.lennard_jones_required())
        sp.enable_cantera(self.cantera_required())
        sp.popup()
        col = tw.currentColumn()
        if col==0:
            sp.ui.lineedit_alias.setFocus()
        elif col==1:
            # Phase is not editable
            pass
        elif col==2:
            sp.ui.lineedit_mol_weight.setFocus()
        elif col==3:
            sp.ui.lineedit_h_f.setFocus()



    def fluid_get_species_refs(self, species):
        ret = []
        msg = self.chemistry_get_species_refs(species)
        if msg:
            ret.append("reaction %s" % msg)

        species_num = 1 + list(self.fluid_species.keys()).index(species) # :(

        for key in keyword_args.keys_by_type['species']:
            indices = self.project.get_key_indices(key)
            if not indices: #Keys not set
                continue
            arg_types = keyword_args.keyword_args.get(key,())
            if 'phase' in arg_types: # It's a solid species, not fluid
                continue
            if arg_types == ('species',): # This will be deleted
                continue
            species_pos = arg_types.index('species')
            for args in indices:
                if args[species_pos] != species_num:
                    continue
                if self.project.get_value(key, args=args): # Ignore settings of None, False, or 0
                    ret.append(format_key_with_args(key,args))

        return ret


    def fluid_delete_species_keys(self, species):
        """Delete all keywords associated with specified species,
        fixing up the resulting gap in sequence"""
        for key in keyword_args.keys_by_type['species']:
            indices = self.project.get_key_indices(key)
            if not indices:
                continue
            arg_types = keyword_args.keyword_args.get(key,())
            if 'phase' in arg_types: # solids species
                continue
            start = 0
            for n in range(arg_types.count('species')):
                species_pos = arg_types.index('species',start)
                for args in sorted(indices):
                    args_species = args[species_pos]
                    if args_species < species:
                        continue
                    args1 = list(args)
                    args1[species_pos] += 1
                    val = self.project.get_value(key, args=args1)
                    self.update_keyword(key, val, args=args)
                start = species_pos + 1


    def setup_fluid(self, allow_disabled_tab=False):
        # Called whenever we switch to fluid tab
        self.P = 0
        ui = self.ui.fluid

        #### Fluid phase
        self.fluid_solver_disabled = (self.project.get_value('ro_g0') == 0.0)
        self.disable_fluid_solver(self.fluid_solver_disabled)
        self.update_fluid_species_table() # Necessary?  Should be done when we show fluid pane.

        # fluid momentum and species eq. handled by _keyword_ widget

        # Moved from gui.py
        # handle a bunch of items which are essentially the same
        #  Note fewer of these are the same now due to changes
        for (setter, name) in ((self.set_fluid_density_model, 'ro'),
                               #(self.set_fluid_viscosity_model, 'mu'),
                               (self.set_fluid_specific_heat_model, 'c_p'), # inconsistent
                               ):
            name_g0 = 'c_pg0' if name=='c_p' else name+'_g0'
            name_usr = 'usr_cpg' if name=='c_p' else 'usr_'+name+'g'
            val_g0 = self.project.get_value(name_g0)
            val_usr = self.project.get_value(name_usr)

            if val_usr is not None and val_g0 is not None:
                self.print_internal('Warning: %s and %s are both set' % (name_g0, name_usr))
                # this is getting printed after error count ... should be included in # of errs

            setter(CONSTANT if val_g0 is not None
                   else UDF if val_usr is not None
                   else 1,
                   in_setup=True)

        # molecular weight model is the odd one (only 2 settings)
        if self.project.get_value('mw_avg') is not None:
            self.set_fluid_mol_weight_model(CONSTANT)
        else:
            # requires molecular weights for all species components, when should we validate?
            self.set_fluid_mol_weight_model(1)

        # Viscosity is special too
        if self.project.get_value('usr_mug'):
            self.unset_keyword('usr_mug')
            self.update_keyword('mu_g_model', 'USR')
        mu_g_model = self.project.get_value('mu_g_model')
        if mu_g_model is None:
            if self.project.get_value('mu_g0') is not None:
                mu_g_model = 'CONSTANT'
            elif self.project.get_value('usr_mug'):
                mu_g_model = 'USR'
        if mu_g_model is None:
            mu_g_model = 'SUTHERLAND' # default
        if mu_g_model.upper() not in VISCOSITY_MODELS:
            self.warning('Unknown viscosity model %s, setting to SUTHERLAND' % mu_g_model,
                         popup=True)
            mu_g_model = 'SUTHERLAND'
            self.update_keyword('mu_g_model', mu_g_model)
        for v in ('LENNARD_JONES', 'CANTERA_POLY'):
            item = get_combobox_item(ui.combobox_fluid_viscosity_model,
                                     VISCOSITY_MODELS.index(v))
            self.set_widget_enabled(item, bool(self.fluid_species),
                                    reason="Requires fluid species.")
        if mu_g_model in ('LENNARD_JONES', 'CANTERA_POLY'):
            if not self.fluid_species:
                self.warning('%s model requires fluid species, setting to SUTHERLAND' % mu_g_model,
                             popup=True)
                mu_g_model = 'SUTHERLAND'
                self.update_keyword('mu_g_model', mu_g_model)
        self.set_fluid_viscosity_model(VISCOSITY_MODELS.index(mu_g_model.upper()), in_setup=True)

        # Everybody is special now
        self.setup_fluid_diffusivity()

        # Thermal conductivity
        energy_eq = bool(self.project.get_value('energy_eq', default=True))
        species_eq = bool(self.project.get_value('species_eq', default=True, args=[0]))
        kg_model = self.project.get_value('kg_model', default="AIR")
        if kg_model and kg_model.upper() not in KG_MODELS:
            self.warning("Invalid kg_model '%s', setting to AIR"%kg_model,
                         popup=True)
            kg_model = "AIR"
            self.update_keyword('kg_model', kg_model)

        for v in ('LENNARD_JONES', 'CANTERA_POLY'):
            item = get_combobox_item(ui.combobox_fluid_kg_model,
                                     KG_MODELS.index(v))
            self.set_widget_enabled(item, bool(self.fluid_species),
                                    reason="Requires fluid species.")


        if kg_model in ('LENNARD_JONES', 'CANTERA_POLY'):
            if not self.fluid_species:
                self.warning('%s model requires fluid species, setting to AIR' % kg_model,
                             popup=True)
                kg_model = 'AIR'
                self.update_keyword('kg_model', kg_model)

        self.set_fluid_kg_model(KG_MODELS.index(kg_model.upper()), in_setup=True)

        enabled = energy_eq
        for item in (ui.label_fluid_specific_heat_model,
                     ui.combobox_fluid_specific_heat_model,
                     ui.label_fluid_kg_model,
                     ui.combobox_fluid_kg_model,
                     # more ?
                     ):
            if enabled:
                self.enable_widget(item)
            else:
                self.disable_widget(item, reason="Requires energy equations.")

        ui.checkbox_species_eq.setChecked(species_eq)
        self.handle_fluid_species_eq(species_eq)
        for c in 'xyz':
            key = 'momentum_%s_eq' % c
            cb = getattr(ui, 'checkbox_%s'%key)
            cb.setChecked(self.project.get_value(key, default=True, args=[0]))

        # c_pg0 == specific heat for fluid phase
        enabled = energy_eq
        cb = ui.combobox_fluid_specific_heat_model
        for w in (cb,
                  ui.label_fluid_specific_heat_model):
            if enabled:
                self.enable_widget(w)
            else:
                self.disable_widget(w, reason="Requires energy equations.")

        # Do not allow MIXTURE if no species defined
        if enabled:
            item = get_combobox_item(cb, MIXTURE)
            enable_item = bool(self.fluid_species)
            if enable_item:
                self.enable_widget(item)
            else:
                self.disable_widget(item, reason="Requires defined species, and specific heat of each.")
            if enable_item is False and cb.currentIndex() == MIXTURE:
                val = self.get_retained_keyword('c_pg0', default=default_values.c_pg0)
                self.update_keyword('c_pg0', val)
                self.set_fluid_specific_heat_model(CONSTANT)

        if not enabled:
            for item in (ui.lineedit_keyword_c_pg0,
                         ui.label_c_pg0_units):
                self.disable_widget(item, reason="Requires energy equations.")
            self.retain_keyword('c_pg0')
            self.unset_keyword('c_pg0')
            ui.lineedit_keyword_c_pg0.setText('')
            cb.setCurrentIndex(MIXTURE)

        if not enabled:
            self.retain_keyword('k_g0')
            self.unset_keyword('k_g0')
            self.retain_keyword('kg_model')
            self.unset_keyword('kg_model')

        # issues/533  do not allow species eq when energy eq not enabled and solver == DEM/CGP
        cb = ui.checkbox_species_eq
        if self.project.solver in (DEM,CGP,SQP,GSP) and not energy_eq:
            species_eq = False
            self.update_keyword('species_eq', False, args=[0])
            self.disable_widget(cb, reason="Requires energy equations when using %s model."%
                               ('DEM' if self.project.solver==DEM
                                else 'CGP' if self.project.solver==CGP
                                else 'SQP' if self.project.solver==SQP
                                else 'GSP' if self.project.solver==GSP
                                else '???'))
        else:
            self.enable_widget(cb)
        cb.setChecked(species_eq)

        # Autoselect if unique row
        tw = ui.tablewidget_fluid_species
        if get_selected_row(tw) is None and tw.rowCount() == 1:
            tw.setCurrentCell(0,0)

        ui.frame_add_delete_copy_species.setVisible(species_eq or bool(self.fluid_species))
        self.handle_fluid_species_selection() # buttons

    def setup_fluid_diffusivity(self, force=False, model=None):
        # Diffusivity is special too
        #TODO sanity check fluid diffusion model keys
        ui = self.ui.fluid
        cb = ui.combobox_fluid_diffusion_model
        if model is None:
            cb.setCurrentIndex(0 if self.project.get_value('dif_g0') is not None
                               else 2 if self.project.get_value('multi_component_diffusion')
                               else 3 if self.project.get_value('usr_difg')
                               else 1)
        set_combobox_tooltip(cb)

        item = get_combobox_item(cb, MULTICOMPONENT)
        self.set_widget_enabled(item, len(self.fluid_species)>1,
                                reason="Requires more than 1 fluid species.")
        # Don't stay on a disabled item
        if len(self.fluid_species) < 2 and self.project.get_value('multi_component_diffusion'):
            self.warning("Multi-component diffusion requires 2 or more fluid species, setting default diffusion model",
                         popup=True)
            self.set_fluid_diffusion_model(AIR)

        cb = ui.combobox_dif_coeff_kt
        item = get_combobox_item(cb, 1)
        dif_coeff_kt = bool(self.project.get_value('dif_coeff_kt', default=False))
        if not self.fluid_species:
            self.set_widget_enabled(item, False,
                                    reason="Requires fluid species.")
        else:
            self.set_widget_enabled(item, self.fluid_density_model==IDEAL,
                                    reason="Requires ideal gas law.")

        if dif_coeff_kt and not self.fluid_species:
            self.warning("Lennard-Jones model requires defined fluid species.",
                         popup=True)
            self.unset_keyword('dif_coeff_kt')
            dif_coeff_kt = False

        if dif_coeff_kt and self.fluid_density_model != IDEAL:
            self.warning("Lennard-Jones model requires ideal gas law.",
                         popup=True)
            self.unset_keyword('dif_coeff_kt')
            dif_coeff_kt = False

        if self.project.get_value('multi_component_diffusion'):
            ui.widget_multicomponent_diffusion.setVisible(True)
            cb.setCurrentIndex(1 if dif_coeff_kt else 0)
            ui.widget_binary.setVisible(not dif_coeff_kt)
            self.species_popup.enable_lennard_jones(self.lennard_jones_required())
            if not dif_coeff_kt:
                self.setup_widget_binary(force=force)
            #ui.checkbox_keyword_dif_thermal.setChecked(bool(self.project.get_value('dif_thermal')))
        else:
            ui.widget_multicomponent_diffusion.setVisible(False)


    def reset_fluid(self):
        # Set all fluid-related state back to default
        ui = self.ui.fluid
        self.fluid_phase_name = 'Fluid'
        self.fluid_species.clear()
        self.init_fluid_default_models()
        self.fluid_species_warnings_shown.clear()
        le = ui.lineedit_keyword_ro_g0
        le.required = False
        le.minimum = 0.0
        le.saved_value = default_values.ro_g0 # fallback
        le.updateValue('ro_g0', le.saved_value)
        # TODO remove dynamically created input widgets, although this should
        #  get handled next time we call 'setup'
