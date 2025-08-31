# -*- coding: utf-8 -*-

import os
from copy import deepcopy
from shutil import copyfile

from qtpy.QtWidgets import (QCheckBox, QComboBox, QFileDialog, QHeaderView,
                            QLabel, QMenu, QMessageBox, QTableWidgetItem,
                            QToolButton)
from qtpy.QtCore import QEvent
from qtpy.QtGui import QPalette, QValidator

from mfixgui.constants import *
from mfixgui.tools import dequote
from mfixgui.tools.read_chemkin import chemkin_to_mfix
from mfixgui.tools.qt import (get_combobox_item, get_selected_row,
                              get_selected_rows, get_icon, set_combobox_tooltip,
                              set_item_enabled, set_item_noedit,
                              sub_icon_height, sub_icon_size,
                              widget_iter)
from mfixgui.project import FloatExp
from mfixgui.widgets.base import LineEdit, ComboBox

from mfixgui.reaction_parser import ReactionParser
from mfixgui.arrows import RIGHT_ARROW

from json import JSONDecoder, JSONEncoder

MAX_RXNS = 2000 # param_mod.f

#Column numbers
COL_ENABLE, COL_RXN_NAME, COL_CHEM_EQ = 0,1,2
COL_PHASE, COL_SPECIES, COL_COEFF = 0,1,2

UNIMOL, BIMOL, PLOG = 0,1,2
LINDEMANN, TROE, SRI = 0,1,2
JAN, FIT1 = 0,1

# TODO remove .reaction_edited flag and just check working reaction

# Helper widgets for reactant and product tables
# Since these tables are populated with CellWidgets, normal selection
# mechanism doesn't work, so we capture focusInEvent
class ExtComboBox(ComboBox):
    def focusInEvent(self, ev):
        if self.side == 'reactants':
            self.parent.chemistry_handle_reactant_selection(row=self.row)
        else:
            self.parent.chemistry_handle_product_selection(row=self.row)
        ComboBox.focusInEvent(self, ev)

class ExtLineEdit(LineEdit):
    def focusInEvent(self, ev):
        if self.side == 'reactants':
            self.parent.chemistry_handle_reactant_selection(row=self.row)
        else:
            self.parent.chemistry_handle_product_selection(row=self.row)
        LineEdit.focusInEvent(self, ev)


def unicode_chem_eq(chem_eq):
    for x in '==', '-->':
        chem_eq = chem_eq.replace(x, RIGHT_ARROW)
    return chem_eq

def reaction_is_des(project, reaction):
    is_des = any(project.get_value('solids_model', default='TFM', args=[p]) != 'TFM'
                 for p in reaction.get('phases', []))
    return is_des


class Chemistry:
    #Chemistry Task Pane Window: This section allows a user to define chemical reaction input.

    def init_chemistry(self):
        ui = self.ui.chemistry
        ui.dynamic_widgets = {} # {keyword:  [(label, lineedit), ...] }
        tw = ui.tablewidget_reactions
        # data members
        self.current_reaction_name = None
        self.reaction_edited = False
        self.adding_reaction = False
        self.working_reaction = None
        self.reaction_mass_totals = [None, None]
        self.disabled_reactions = {}

        ui.checkbox_keyword_arrhenius_rrates_fluid.post_update = self.chemistry_chemkin_warning
        ui.checkbox_keyword_arrhenius_rrates_des.post_update = self.chemistry_chemkin_warning

        # Toolbuttons
        ui.toolbutton_add_reaction.clicked.connect(self.chemistry_add_reaction)
        ui.toolbutton_delete_reaction.clicked.connect(lambda checked: self.chemistry_delete_reaction())
        ui.toolbutton_add_reactant.clicked.connect(self.chemistry_add_reactant)
        ui.toolbutton_delete_reactant.clicked.connect(self.chemistry_delete_reactant)
        ui.toolbutton_add_product.clicked.connect(self.chemistry_add_product)
        ui.toolbutton_delete_product.clicked.connect(self.chemistry_delete_product)
        ui.toolbutton_ok.clicked.connect(self.chemistry_apply_changes)
        ui.toolbutton_cancel.clicked.connect(self.chemistry_revert_changes)

        #ui.toolbutton_load_reactions.clicked.connect(self.chemistry_load_reactions)
        ui.toolbutton_load_reactions.setPopupMode(QToolButton.InstantPopup)
        ui.toolbutton_load_reactions.setIcon(get_icon('load.svg'))
        ui.toolbutton_load_reactions.setIconSize(sub_icon_size()/2)
        menu = ui.load_menu = QMenu(ui)
        #https://stackoverflow.com/questions/44745459/remove-icon-space-from-qmenu
        #menu.setStyleSheet("QMenu::item {padding-left: 5px;}")
        menu.addAction("MFiX format", lambda: self.chemistry_load_reactions(chemkin=False))
        menu.addAction("CHEMKIN format", lambda: self.chemistry_load_reactions(chemkin=True))
        ui.toolbutton_load_reactions.setMenu(menu)

        ui.toolbutton_delete_reaction.setEnabled(False) # Need a selection
        ui.toolbutton_delete_reactant.setEnabled(False)
        ui.toolbutton_delete_product.setEnabled(False)
        ui.groupbox_dh.setEnabled(False)
        ui.toolbutton_ok.setEnabled(False)
        ui.toolbutton_cancel.setEnabled(False)


        # Tablewidgets
        tw = ui.tablewidget_reactions
        tw.itemSelectionChanged.connect(self.chemistry_handle_selection)
        tw.resizeEvent = (lambda old_method, tw=tw:
                          (lambda event:
                           (self.fixup_chemistry_table(tw, stretch_column=2),
                            old_method(event))[-1]))(tw.resizeEvent) # pylint: disable=undefined-variable


        for tw in ui.tablewidget_reactants, ui.tablewidget_products:
            tw.current_row = None

        ui.tablewidget_plog.horizontalHeader().setMinimumSectionSize(0)

        # Set reaction name
        ui.lineedit_reaction_name.editingFinished.connect(self.set_reaction_name)
        class RxnIdValidator(QValidator):
            #  Alphanumeric combinations (no special characters excluding underscores)
            #  Limited to 32 characters
            #  First character must be a letter
            #  No blank spaces
            def __init__(self, parent=None):
                super(RxnIdValidator, self).__init__()
                self.parent = parent

            def validate(self, text, pos):
                #self.parent.len = len(text)
                #if len(text) == 0:
                #    # How to reset the lineedit after user blanks out input?
                #    return (QValidator.Intermediate, text, pos)
                if len(text) == 0:
                    return (QValidator.Acceptable, text, pos)
                elif 1 <= len(text) <= 32 and text[0].isalpha() and all(c.isalnum() or c=='_' for c in text):
                    if text.lower() in self.parent.keyword_doc: # cannot use keywords as reaction names!
                        return (QValidator.Intermediate, text, pos)
                    else:
                        return (QValidator.Acceptable, text, pos)
                else:
                    return (QValidator.Invalid, text, pos)
        ui.lineedit_reaction_name.setValidator(RxnIdValidator(parent=self))

        ui.groupbox_dh.keys = ['dh', 'fracdh'] # Locatability (these are not really mfix keys)
        ui.checkbox_dh.keys = ui.groupbox_dh.keys
        ui.checkbox_dh.clicked.connect(self.handle_checkbox_dh)

        le =  ui.lineedit_dh
        le.key = 'dh'
        le.dtype = float
        le.value_updated.connect(self.handle_dh)
        self.add_tooltip(le, key='dh', description='Heat of reaction')
        label = ui.label_dh
        self.add_tooltip(label, key='dh', description='Heat of reaction')

        for le in (ui.lineedit_keyword_chem_min_species_fluid,
                   ui.lineedit_keyword_chem_min_species_solid,
                   ui.lineedit_keyword_stiff_chem_max_steps,
                   ui.lineedit_keyword_stiff_chem_abs_tol,
                   ui.lineedit_keyword_stiff_chem_rel_tol,
                   ui.lineedit_keyword_stiff_chem_min_rate,
                   ui.lineedit_keyword_stiff_chem_min_rate_dpm):
            le.required = True
        ui.checkbox_stiff_chemistry.clicked.connect(self.handle_checkbox_stiff_chemistry)
        self.add_tooltip(ui.checkbox_stiff_chemistry, key='stiff_chemistry')

        # Sub-panes
        for (i, btn) in enumerate((ui.pushbutton_reactions, ui.pushbutton_options)):
            btn.clicked.connect(lambda ignore, i=i:self.chemistry_change_tab(i))


### CHEMKIN/Arrhenius
        # Note some tooltips are set in chemistry.ui.  Tooltips here
        # are for items which cannot be set there (combobox items, etc)
        for w in (ui.lineedit_arrhenius_A,
                  ui.lineedit_pressure_arrhenius_A):
            w.key = 'A'
        for w in (ui.lineedit_arrhenius_B,
                  ui.lineedit_pressure_arrhenius_B):
            w.key = 'β'
        for w in (ui.lineedit_arrhenius_Ea,
                  ui.lineedit_arrhenius_Ea):
            w.key = 'Eₐ'
        for w in (ui.lineedit_arrhenius_A, ui.lineedit_arrhenius_B,
                  ui.lineedit_arrhenius_Ea):
            w.dtype = float
            w.value_updated.connect(self.chemistry_handle_lineedit_arrhenius)
            w.allow_parameters = False
        for w in (ui.lineedit_pressure_arrhenius_A, ui.lineedit_pressure_arrhenius_B,
                  ui.lineedit_pressure_arrhenius_Ea):
            w.dtype = float
            w.value_updated.connect(self.chemistry_update_press_rxn_param)
            w.allow_parameters = False

        # Checkboxes
        for x in ('reverse_reaction',
                  'reaction_order',
                  'third_body',
                  'pressure_dependent',
                  'landau_teller',
                  'rate_fit'):

            cb = getattr(ui, 'checkbox_' + x)
            f = getattr(self, 'chemistry_handle_checkbox_'+x)
            cb.clicked.connect(f)

        # Comboboxes
        for x in ('reverse_reaction',
                  'third_body',
                  'third_body_species',
                  'third_body_add_species',
                  'pressure_model',
                  'pressure_method',
                  'rate_fit_method'):

            cb = getattr(ui, 'combobox_' + x)
            f = getattr(self, 'chemistry_handle_combobox_' + x)
            cb.currentIndexChanged.connect(f)
        # lineedits
        for s,key in (('pressure_arrhenius_A','A'),
                      ('pressure_arrhenius_B','β'),
                      ('pressure_arrhenius_Ea','Eₐ'),
                      ('sri_a', 'a'),
                      ('sri_b', 'b'),
                      ('sri_c', 'c'),
                      ('sri_d', 'd'),
                      ('sri_e', 'e'),
                      ('troe_alpha', 'α'),
                      ('troe_tstar', 'T ⃰'),
                      ('troe_t2star', 'T ⃰  ⃰'),
                      ('troe_t3star', 'T ⃰  ⃰  ⃰')):
            w = getattr(ui, 'lineedit_' + s)
            w.dtype = float

            w.key = key
            w.value_updated.connect(self.chemistry_handle_pressure_coeff)
        ui.lineedit_sri_d.defaultValue = 1.0
        ui.lineedit_sri_e.defaultValue = 0.0
        ui.lineedit_troe_t2star.defaultValue = 0.0
        for s,key in (('landau_teller_B', 'B'),
                      ('landau_teller_C', 'C')):
            w = getattr(ui, 'lineedit_' + s)
            w.dtype = float
            w.key = key
            w.value_updated.connect(self.chemistry_handle_lineedit_landau_teller)
        subscripts = '₁₂₃₄₅₆₇₈₉'
        for i,s in enumerate(subscripts, 1):
            w = getattr(ui, 'lineedit_fit_b%s' % i)
            w.dtype = float
            w.key = 'b'+s
            w.value_updated.connect(self.chemistry_handle_lineedit_rate_fit)
        # Third-body
        le = ui.lineedit_third_body_effcy
        le.dtype = float
        le.min = 0
        ui.lineedit_third_body_effcy.value_updated.connect(self.chemistry_handle_lineedit_third_body)
        # PLOG table
        tw = ui.tablewidget_plog
        hv = QHeaderView
        resize = tw.horizontalHeader().setSectionResizeMode
        ncols = tw.columnCount()
        #for n in range(0, ncols):
        #    resize(n, hv.Stretch if n==0 else hv.ResizeToContents)
        tw.cellChanged.connect(self.chemistry_handle_plog)
        tw.horizontalHeaderItem(1).setToolTip("Pressure (Pa)")
        tw.horizontalHeaderItem(2).setToolTip(ui.lineedit_arrhenius_A.toolTip())
        tw.horizontalHeaderItem(3).setToolTip(ui.lineedit_arrhenius_B.toolTip())
        tw.horizontalHeaderItem(4).setToolTip(ui.lineedit_arrhenius_Ea.toolTip())
        cb = ui.combobox_pressure_model
        get_combobox_item(cb,2).setToolTip("<html>Specify Arrhenius coefficients at discrete temperatures</html>")

        cb = ui.combobox_reverse_reaction
        get_combobox_item(cb,0).setToolTip("<html>The rate constant of the reverse reaction is calculated based on the rate constant of the forward reaction (kf) and the equilibrium constant (Kc), kr=kf/Kc. kf is calculated based on the Arrhenius coefficients above.</html>")
        get_combobox_item(cb,1).setToolTip("<html>The rate constant of the reverse reaction is calculated directly from Arrhenius coefficients above.</html>")
        set_combobox_tooltip(cb)

        cb = ui.combobox_pressure_method
        get_combobox_item(cb,0).setToolTip("<html>See Lindemann, F.: Trans. Faraday Soc., 17 598 (1922)</html>")
        get_combobox_item(cb,1).setToolTip("<html>See Gilbert, R. G., Luther, K. and Troe, J.: Ber. Bunsenges. Phys. Chem., 87 169 (1983)</html>")
        get_combobox_item(cb,2).setToolTip("<html>See Stewart, P. H., Larson, C. W. and Golden, D.: Combustion and Flame, 75 25 (1989)</html>")
        set_combobox_tooltip(cb)

        ui.groupbox_landau_teller.setToolTip(ui.checkbox_landau_teller.toolTip())
        cb = ui.combobox_rate_fit_method
        get_combobox_item(cb,0).setToolTip("<html>Polynomial fit to the logarithm of temperature</html>")
        get_combobox_item(cb,1).setToolTip("<html>Power series with the exponential of a modified Arrhenius expression</html>")

    def chemistry_change_tab(self, tabnum):
        """switch stacked widget based on selected"""
        ui = self.ui.chemistry
        to_btn = [ui.pushbutton_reactions, ui.pushbutton_options][tabnum]
        self.animate_stacked_widget(
            ui.stackedwidget_chemistry,
            ui.stackedwidget_chemistry.currentIndex(),
            tabnum,
            line=ui.line_chemistry,
            to_btn=to_btn,
            btn_layout=ui.gridlayout_chemistry_tab_btns)

        for btn in [ui.pushbutton_reactions, ui.pushbutton_options]:
            btn.setChecked(btn == to_btn)
            font = btn.font()
            font.setBold(btn == to_btn)
            btn.setFont(font)

        if tabnum == 0: #Reactions
            self.chemistry_handle_selection()

    def chemistry_handle_selection(self):
        # selection callback for main table
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        rows = get_selected_rows(tw)
        enabled = (len(rows) == 1) # Enable widgets for single reaction selection only

        row = rows[0] if enabled else None
        self.reaction_mass_totals = [None, None]
        if enabled:
            self.current_reaction_name = tw.item(row, COL_RXN_NAME).text()
            self.working_reaction = deepcopy(self.project.reactions[self.current_reaction_name])
        else:
            self.current_reaction_name = None
            self.working_reaction = None
        self.clear_reaction_edited()
        reaction_disabled = self.working_reaction and self.working_reaction.get('chem_eq') == 'NONE'
        self.set_widget_enabled(ui.detail_pane, not reaction_disabled)

        # Warn user to specify arrhenius coeff!
        ## This is too aggressive, cannot cancel
        #if self.chemkin_enabled():
        #    arr = self.working_reaction.get("arrhenius_coeff")
        #    if not arr or 'None' in arr:
        #        self.set_reaction_edited()

        self.chemistry_update_detail_pane()
        ui.scrollarea_detail.ensureVisible(0, 0)

        # Restore selection since update_detail_pane does autoselect
        # This allows scrolling through reactions table
        #if row is not None:
        #    tw.setCurrentCell(row,0)
        #    tw.setFocus()
        # No longer needed, autoselect for reactants/products disabled (issues/1864)

    def chemistry_update_detail_pane(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        rows = get_selected_rows(tw)
        reaction = self.working_reaction
        enabled = (len(rows) == 1)
        row = rows[0] if enabled else None
        self.set_widget_enabled(ui.toolbutton_delete_reaction, bool(rows)) # Allow deleting multiple
        self.set_widget_enabled(ui.toolbutton_add_reactant,
                                bool(self.chemistry_find_available_species('reactants')) and enabled)
        self.set_widget_enabled(ui.toolbutton_add_product,
                                bool(self.chemistry_find_available_species('products')) and enabled)

        for w in (ui.label_reaction_name,
                  ui.lineedit_reaction_name,
                  ui.groupbox_reactants,
                  ui.groupbox_products,
                  ui.checkbox_dh,
                  ui.groupbox_arrhenius,
                  ui.groupbox_dh):
            self.set_widget_enabled(w, enabled)

        if reaction and 'PLOG' in reaction.get('press_rxn_param',''):
            ui.groupbox_arrhenius.setVisible(False)
        else:
            ui.groupbox_arrhenius.setVisible(True)

        if not enabled:
            self.chemistry_clear_tables()
            ui.lineedit_reaction_name.clear()
            ui.checkbox_dh.setChecked(False)
            ui.groupbox_dh.setVisible(False)
            for w in widget_iter(ui.groupbox_dh):
                if isinstance(w, LineEdit):
                    w.clear()

            for (cb, gb) in ((ui.checkbox_reverse_reaction,
                              ui.combobox_reverse_reaction),
                             (ui.checkbox_reaction_order,
                              ui.groupbox_reaction_order),
                             (ui.checkbox_third_body,
                              ui.groupbox_third_body),
                             (ui.checkbox_pressure_dependent,
                              ui.groupbox_pressure_dependent),
                             (ui.checkbox_landau_teller,
                              ui.groupbox_landau_teller),
                             (ui.checkbox_rate_fit,
                              ui.groupbox_rate_fit)):
                cb.setVisible(False)
                gb.setVisible(False)
            ui.label_reverse_reaction.setVisible(False)
            ui.checkbox_dh.setEnabled(False)
            ui.groupbox_dh.setVisible(False)
            ui.groupbox_arrhenius.setVisible(False)
            self.current_reaction_name = None
            return

        tw = ui.tablewidget_reactions
        name = tw.item(row, COL_RXN_NAME).text()
        ui.lineedit_reaction_name.setText(name)
        self.current_reaction_name = name

        def handle_phase(tw, cb, row, idx):
            ui = self.ui.chemistry
            # We have to replace the species widget
            old_species_cb = tw.cellWidget(row,  COL_SPECIES)
            if self.working_reaction is None:
                return
            reaction = self.working_reaction
            side = 'reactants' if tw == ui.tablewidget_reactants else 'products'
            species = self.chemistry_find_available_species(side, idx)
            reaction[side][row][0] = species
            item = make_species_item(tw, row, idx, species)
            tw.setCellWidget(row, COL_SPECIES, item)
            if old_species_cb: # necessary?
                try:
                    old_species_cb.activated.disconnect()
                except:
                    pass
                try:
                    old_species_cb.currentIndexChanged.disconnect()
                except:
                    pass
                old_species_cb.deleteLater()
            self.set_reaction_edited()
            self.chemistry_restrict_phases()
            self.chemistry_restrict_species()
            self.chemistry_check_fracdh(reaction) # sets 'phases'
            self.chemistry_update_detail_pane()

        def make_phase_item(tw, row, phase):
            ui = self.ui.chemistry
            side = 'reactants' if tw == ui.tablewidget_reactants else 'products'
            cb = ExtComboBox()
            cb.parent = self
            cb.row = row
            cb.side = side
            phases = [self.fluid_phase_name]
            for name in self.solids.keys():
                phases.append(name)
            for p in phases:
                cb.addItem(p)
            if phase is not None:
                cb.setCurrentIndex(phase)
            cb.currentIndexChanged.connect(lambda idx, tw=tw, cb=cb, row=row: # pylint: disable=undefined-variable
                                           handle_phase(tw, cb, row, idx)) # is this a pylint bug?
            return cb

        def handle_species(tw, cb, row, idx):
            ui = self.ui.chemistry
            #tw.cellWidget(row, COL_SPECIES).setText('1.0')
            self.set_reaction_edited()
            species = tw.cellWidget(row, COL_SPECIES).currentText()
            side = 'reactants' if tw == ui.tablewidget_reactants else 'products'
            #reaction = self.project.reactions[self.current_reaction_name]
            reaction = self.working_reaction
            reaction[side][row][0] = species
            self.chemistry_restrict_phases()
            self.chemistry_restrict_species()
            self.chemistry_update_totals()

        def make_species_item(tw, row, phase, species):
            cb = ExtComboBox()
            cb.parent = self
            idx = 0
            for s in self.species_of_phase(phase):
                cb.addItem(s)
                if s.lower() == species.lower():
                    cb.setCurrentIndex(idx)
                idx += 1
            cb.currentIndexChanged.connect(
                lambda idx, tw=tw, cb=cb, row=row: # pylint: disable=undefined-variable
                handle_species(tw, cb, row, idx))
            cb.side = 'reactants' if tw==ui.tablewidget_reactants else 'products'
            cb.row = row
            return cb

        def handle_coeff(widget, val, args):
            ui = self.ui.chemistry
            tw = ui.tablewidget_reactants if widget.side == 'reactants' else ui.tablewidget_products
            #reaction = self.project.reactions[self.current_reaction_name]
            reaction = self.working_reaction
            if not reaction:
                return
            val = widget.value
            row = widget.row
            tw.current_row = row
            if val in (None, ''):
                val = 1.0
                widget.setText('1.0')
            if val != reaction[widget.side][widget.row][1]:
                self.set_reaction_edited()
                reaction[widget.side][widget.row][1] = val
            self.chemistry_update_totals()

        def make_coeff_item(tw, row, val):
            le = ExtLineEdit()
            le.parent = self
            #le.setMaximumWidth(80) #?
            le.dtype = float
            le.min = 0
            le.setToolTip("Stoichometric coefficient")
            le.key = ''
            le.updateValue('', val)
            le.side = 'reactants' if tw == ui.tablewidget_reactants else 'products'
            le.tw = tw
            le.row = row
            le.value_updated.connect(handle_coeff)
            return le

        self.chemistry_clear_tables()

        for side in 'reactants', 'products':
            data = reaction.get(side,[])
            #data = self.project.reactions[name].get(side, [])
            # Add a "total" row, only if there is data
            tw = ui.tablewidget_reactants if side == 'reactants' else ui.tablewidget_products
            tw.setRowCount(len(data)+1 if data else 0)
            for (row, (species,coeff)) in enumerate(data):
                phase = self.find_species_phase(species)
                if phase is None:
                    self.error("Species %s not found in any phase" % species)
                tw.setCellWidget(row, COL_PHASE, make_phase_item(tw, row, phase))
                tw.setCellWidget(row, COL_SPECIES, make_species_item(tw, row, phase, species))
                tw.setCellWidget(row, COL_COEFF, make_coeff_item(tw, row, coeff))
            if data:
                item = QTableWidgetItem('Total mol. weight')
                set_item_noedit(item)
                tw.setItem(row+1, COL_SPECIES, item)
                item = QTableWidgetItem('0.0')
                set_item_noedit(item)
                font=item.font()
                font.setBold(True)
                item.setFont(font)
                set_item_noedit(item)
                tw.setItem(row+1, COL_COEFF, item)
                # Avoid selectable, editable cell in bottom corner
                item = QTableWidgetItem('')
                set_item_noedit(item)
                tw.setItem(row+1, COL_PHASE, item)

            self.fixup_chemistry_table(tw)
            row = tw.current_row
            # Avoid auto-select since it messes up table focus (issues/1864)
            #if row is None and len(data) == 1:
            #    row = 0 # Auto-select single row
            if side == 'reactants':
                self.chemistry_handle_reactant_selection(row)
            else:
                self.chemistry_handle_product_selection(row)

        self.fixup_chemistry_table(ui.tablewidget_reactions, stretch_column=COL_CHEM_EQ)
        self.chemistry_restrict_phases()
        self.chemistry_restrict_species()
        self.chemistry_update_totals()
        self.chemistry_update_dh()
        self.chemistry_update_arrhenius()
        self.chemistry_update_reverse_reaction()
        self.chemistry_update_reaction_order()
        self.chemistry_update_third_body()
        self.chemistry_update_pressure()
        self.chemistry_update_landau_teller()
        self.chemistry_update_rate_fit()

    def set_reaction_name(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        if row is None:
            return
        name = ui.lineedit_reaction_name.value
        if len(name) == 0: # Reset name if empty input
            name = tw.item(row, COL_RXN_NAME).text()
            ui.lineedit_reaction_name.setText(name)
            return
        tw.item(row, COL_RXN_NAME).setText(name)
        tw.resizeRowToContents(row)
        prev_name = self.current_reaction_name
        self.current_reaction_name = name # We can only edit the current reaction
        # update ordered dict, keeping order
        keys = list(self.project.reactions.keys())
        keys[row] = name
        self.project.reactions = dict(zip(keys, self.project.reactions.values()))
        # We might have renamed a disabled reaction
        if prev_name != name:
            if prev_name in self.disabled_reactions:
                self.disabled_reactions[name] = self.disabled_reactions.pop(prev_name)
            self.set_unsaved_flag()

    def chemistry_restrict_phases(self):
        # Set up comboboxes to only phases with defined species
        ui = self.ui.chemistry
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            side = 'reactants' if tw == ui.tablewidget_reactants else 'products'
            for row in range(tw.rowCount()-1):
                cb = tw.cellWidget(row, COL_PHASE)
                for i in range(cb.count()):
                    item = get_combobox_item(cb, i)
                    enabled = self.chemistry_find_available_species(side, match_phase=i)
                    set_item_enabled(item, enabled)
                    self.add_tooltip(item, key=None,
                                     description='No available species' if not enabled else None)



    def chemistry_restrict_species(self):
        # Set up comboboxes to ensure that no species is duplicated as a reactant or product
        ui = self.ui.chemistry
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            # Collect species info
            species = {}
            for row in range(tw.rowCount()-1): # Skip 'total'
                cb = tw.cellWidget(row, COL_SPECIES)
                species[row] = cb.currentText()
            # For each combobox item, determine whether setting the combobox to that value
            # results in duplicated species
            n_rows = tw.rowCount() - 1 # Skip 'total'
            for row in range(n_rows):
                cb = tw.cellWidget(row, COL_SPECIES)
                orig_species = cb.currentText()
                for i in range(cb.count()):
                    item = get_combobox_item(cb, i)
                    species[row] = item.text() # consider what setting this combobox would do...
                    enabled = (len(set(species.values())) == n_rows) # should be as many species as rows
                    set_item_enabled(item, enabled)
                species[row] = orig_species # ... and set it back


    def handle_checkbox_stiff_chemistry(self, enabled):
        ui = self.ui.chemistry
        key = 'stiff_chemistry'
        self.update_keyword(key, enabled)
        self.set_widget_enabled(ui.groupbox_stiff_chemistry, enabled,
                                ignore_lock=self.lock_level==PARTIAL_LOCK)
        ui.groupbox_stiff_chemistry.setVisible(enabled)


### Heat of reaction (dh/fracdh)
    def handle_checkbox_dh(self, enabled):
        ui = self.ui.chemistry
        key = 'fracdh'
        for w in ([v[1] for v in ui.dynamic_widgets.get(key, [])] +
                  [ui.label_dh, ui.lineedit_dh, ui.label_dh_units]):
            w.setEnabled(enabled)
        reaction = self.working_reaction
        if reaction is None:
            ui.groupbox_dh.setVisible(False)
            return
        if enabled:
            dh = reaction.get('dh')
            if dh is None:
                dh = ui.lineedit_dh.value #restore
                if dh == '':
                    dh = None
                if dh is not None:
                    reaction['dh'] = dh
            if dh is None:
                dh = 0.0 #Default
                reaction['dh'] = dh
            self.chemistry_check_fracdh(reaction)
        else:
            reaction.pop('dh', None)
            reaction.pop('fracdh', None)
        self.set_reaction_edited()
        self.chemistry_update_buttons()
        self.chemistry_update_dh()


    def handle_dh(self, widget, val, args):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if reaction is None:
            return
        reaction['dh'] = val['dh']
        self.set_reaction_edited()
        self.chemistry_update_buttons()


    def handle_fracdh(self, widget, vals, args):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if reaction is None:
            return
        if len(args) != 1:
            raise ValueError(args)
        arg = args[0]
        phases = self.chemistry_reaction_phases(reaction)
        num_phases = len(phases)
        key = 'fracdh'
        val = vals.get(key)
        if val in (None, ''):
            val = 0.0
        reaction[key][arg] = val
        if num_phases == 1:
            reaction[key][arg] = 1.0
        elif num_phases == 2:
            other = phases[1-phases.index(arg)]
            reaction[key][other] = 1.0 - val

        self.set_reaction_edited()
        self.chemistry_check_fracdh(reaction)
        self.chemistry_update_dh()
        self.chemistry_update_buttons()


    def chemistry_update_dh(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        key = 'dh'
        dh = reaction.get(key)
        ui.lineedit_dh.updateValue(key, dh)
        enabled = (dh is not None)
        ui.checkbox_dh.setChecked(enabled)
        key = 'fracdh'
        for w in ([v[1] for v in ui.dynamic_widgets.get(key, [])] + [
                ui.label_dh, ui.lineedit_dh, ui.label_dh_units]):
            self.set_widget_enabled(w, enabled)
        ui.groupbox_dh.setVisible(enabled)
        phases = self.chemistry_reaction_phases(reaction)
        phases = [p for p in phases if p is not None] # Filter out None's (only happens when species not in any phase)
        num_phases = len(phases)
        def phase_name(phase):
            if phase == 0:
                return self.fluid_phase_name
            else:
                return list(self.solids.keys())[phase-1]

        # Dynamically created inputs for fracdh
        key = 'fracdh'
        layout = ui.groupbox_dh.layout()
        ws = ui.dynamic_widgets.get(key, [])
        if len(ws) > num_phases:
            for (i, (label,le)) in enumerate(ws[num_phases:], num_phases):
                for w in (label, le):
                    layout.removeWidget(w)
                    layout.removeWidget(w)
                    w.setParent(None)
                    w.deleteLater()
            ui.dynamic_widgets[key] = ws = ws[:num_phases]

        while len(ws) < num_phases:
            i = len(ws)
            l = QLabel()
            layout.addWidget(l, i+1, 0)
            le = LineEdit()
            le.key = key
            le.dtype = float
            le.min, le.max = 0.0, 1.0
            le.value_updated.connect(self.handle_fracdh)
            layout.addWidget(le, i+1, 1, 1, 2)
            ws.append((l, le))

        fracdh = reaction.get(key, {})
        # remove any extra 'fracdh' settings (is this necessary?)
        for k, v in list(fracdh.items()):
            if k not in phases:
                fracdh.pop(k, None)
        reaction[key] = fracdh
        self.chemistry_check_fracdh(reaction)
        if ws:
            le = ws[0][1] # first fracdh lineedit
            enabled = len(ws) != 1
            le.setReadOnly(not enabled)
            #le.setEnabled(enabled) # better way to indicate readonly?
            if not enabled:
                reaction[key] = {phases[0]: 1.0}

        for ((label, le), phase) in zip(ws, phases):
            name = phase_name(phase)
            label.setText('Fraction assigned to %s' % name)
            descr = 'Fraction of heat of reaction for phase %s' % name
            label.args = [phase]
            le.args = [phase]
            self.add_tooltip(label, key=key, description=descr)
            self.add_tooltip(le, key=key, description=descr)
            le.updateValue(key, fracdh.get(phase))

        ui.dynamic_widgets[key] = ws


#### Chem eq
    def format_chem_eq(self, reactants, products):
        tmp = []
        for side in reactants, products:
            tmp.append(' + '.join(species if coeff == 1.0 else '%g*%s' % (coeff, species)
                                  for (species, coeff) in side))
        chem_eq = ' --> '.join(tmp)
        return chem_eq


    def chemistry_update_chem_eq(self, name):
        # Allow updating non-selected reaction, used when renaming species
        ui = self.ui.chemistry
        reaction = self.project.reactions.get(name)
        if reaction is None:
            return
        chem_eq = self.format_chem_eq(reaction.get('reactants',[]), reaction.get('products',[]))
        if reaction['chem_eq'].upper() == 'NONE': # Update disabled reaction
            self.disabled_reactions[name] = chem_eq
            display_text = 'Disabled'
        else: # Update active reaction
            reaction['chem_eq'] = chem_eq
            display_text = unicode_chem_eq(chem_eq)

        tw = ui.tablewidget_reactions
        for row in range(tw.rowCount()):
            if tw.item(row, COL_RXN_NAME).text() == name:
                ui.tablewidget_reactions.item(row, COL_CHEM_EQ).setText(display_text)
                break
        self.fixup_chemistry_table(tw, stretch_column=2) # chem eq


    def chemistry_update_totals(self):
        ui = self.ui.chemistry
        self.reaction_mass_totals = [None, None]
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            nrows = tw.rowCount()
            if nrows < 2: #
                continue
            tot = 0.0
            for row in range(nrows-1):
                species = tw.cellWidget(row, COL_SPECIES).currentText()
                m_w = self.species_mol_weight(species)
                if m_w is None: # Undefined mol. wt - species_mol_wt printed a warning
                    self.warning("Molecular weight for '%s' not found" % species)
                    continue
                coeff = tw.cellWidget(row, COL_COEFF).value
                if coeff in (None, ''):
                    continue # Empty input field
                tot += m_w * float(coeff)
            self.reaction_mass_totals[tw==ui.tablewidget_products] = tot
            tot = round(tot, 6)
            tw.item(nrows-1, COL_COEFF).setText(str(tot))
        self.chemistry_update_buttons()


    def chemistry_check_reaction_balance(self):
        ui = self.ui.chemistry
        if not self.working_reaction:
            return False, "No reaction"
        reaction = self.working_reaction
        reactants = reaction.get('reactants',[])
        products = reaction.get('products',[])
        phases = reaction.get('phases',[])
        if not reactants:
            return False, "No reactants defined"
        if not products:
            return False, "No products defined"
        if sorted(reactants) == sorted(products):
            return False, "Reaction cannot have the same species as both reactant and product."
        totals = self.reaction_mass_totals
        if any (t is None for t in totals):
            return False, "Total molecular weight unavailable"
        mass_reactants, mass_products = totals
        if mass_reactants == 0.0:
            return False, "Reactant mass total = 0"
        if mass_products == 0.0:
            return False, "Product mass total = 0"

        # Issues/559
        losers = []
        winners = []
        for p in phases:
            l = r = 0
            for (s,c) in reactants:
                if self.find_species_phase(s) != p:
                    continue
                l += c * self.species_mol_weight(s)
            for (s,c) in products:
                if self.find_species_phase(s) != p:
                    continue
                r += c * self.species_mol_weight(s)
            if r == 0:
                if l > 0.0001:
                    losers.append(p)
            elif l/r > 1.0001:
                losers.append(p)
            elif l/r < 0.9999:
                winners.append(p)
        if len(losers) > 1:
            return False, "More than one phase has net mass loss"
        if len(winners) > 1:
            return False, "More than one phase has net mass gain"

        if abs(mass_products/mass_reactants - 1.0) < 1e-4 :
            return True, "Apply changes?"
        else:
            return False, "Reaction is unbalanced"


    def chemistry_handle_reactant_selection(self, row=None):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactants
        tw.current_row = row
        enabled = (row is not None)
        self.set_widget_enabled(ui.toolbutton_delete_reactant, enabled)
        for r in range(tw.rowCount()):
            if r != row:
                for c in (COL_PHASE, COL_SPECIES, COL_COEFF):
                    w = tw.cellWidget(r, c)
                    if w:
                        if c == COL_COEFF:
                            w.deselect()
                        w.setStyleSheet("")
        if enabled:
            tw.setCurrentCell(row, COL_COEFF)
            for c in (COL_PHASE, COL_SPECIES, COL_COEFF):
                w = tw.cellWidget(row, c)
                if w:
                    if c == COL_COEFF:
                        #w.setStyleSheet('background-color: palette(highlight); color: palette(highlightedText)')
                        w.setStyleSheet('background-color: lightBlue')
                    else:
                        # Don't set color of dropdown, just the button
                        #w.setStyleSheet('QComboBox:!on{background-color: palette(highlight); color: palette(highlightedText)}')
                        w.setStyleSheet('QComboBox:!on{background-color: lightBlue}')

            #tw.cellWidget(row, COL_COEFF).setBackgroundRole(QPalette.HighlightedText)
            #tw.cellWidget(row, COL_COEFF).setForegroundRole(QPalette.HighlightedText)


    def chemistry_handle_product_selection(self, row=None):
        ui = self.ui.chemistry
        tw = ui.tablewidget_products
        tw.current_row = row
        enabled = (row is not None)
        self.set_widget_enabled(ui.toolbutton_delete_product, enabled)
        for r in range(tw.rowCount()):
            if r != row:
                for c in (COL_PHASE, COL_SPECIES, COL_COEFF):
                    w = tw.cellWidget(r, c)
                    if w:
                        if c == COL_COEFF:
                            w.deselect()
                        w.setStyleSheet("")
        if enabled:
            tw.setCurrentCell(row, COL_COEFF)
            for c in (COL_PHASE, COL_SPECIES, COL_COEFF):
                w = tw.cellWidget(row, c)
                if w:
                    if c == COL_COEFF:
                        #w.setStyleSheet('background-color: palette(highlight); color: palette(highlightedText)')
                        w.setStyleSheet('background-color: lightBlue')
                    else:
                        # Don't set color of dropdown, just the button
                        #w.setStyleSheet('QComboBox:!on{background-color: palette(highlight); color: palette(highlightedText)}')
                        w.setStyleSheet('QComboBox:!on{background-color: lightBlue}')

            #tw.cellWidget(row, COL_COEFF).setBackgroundRole(QPalette.HighlightedText)
            #tw.cellWidget(row, COL_COEFF).setForegroundRole(QPalette.HighlightedText)



    def chemistry_reaction_phases(self, reaction):
        """return sorted list of phase indices involved in specified reaction"""
        alias_list = [k[0] for k in reaction.get('reactants',[]) + reaction.get('products',[])]
        phases = list(set(map(self.find_species_phase, alias_list)))
        if None in phases: #Oops, couldn't find phase
            phases.remove(None)
        phases.sort()
        return phases


    def chemistry_update_enabled(self):
        disabled = False
        # Fewer than 2 species, no chemistry possible
        if len(self.fluid_species) + sum(map(len, self.solids_species.values())) < 2:
            disabled = True
        if self.project.reactions: # Don't ever disable pane if reactions are defined
            disabled = False
        item = self.find_navigation_tree_item("Chemistry")
        item.setDisabled(disabled)
        item.setToolTip(0, "Chemical reactions require at least 2 defined species" if disabled
                        else "Chemical reactions")

    def fixup_chemistry_table(self, tw, stretch_column=1): # species column, for reactant/product tables
        ui = self.ui.chemistry
        hv = QHeaderView
        resize = tw.horizontalHeader().setSectionResizeMode
        ncols = tw.columnCount()
        nrows = tw.rowCount()
        for n in range(0, ncols):
            resize(n, hv.Stretch if n == stretch_column else hv.ResizeToContents)
        for n in range(0, nrows):
            tw.resizeRowToContents(n)

        # trim excess vertical space - can't figure out how to do this in designer
        header_height = tw.horizontalHeader().height()

        # Note - scrollbar status can change outside of this function.
        # Do we need to call this every time window geometry changes?
        scrollbar_height = tw.horizontalScrollBar().isVisible() * (4+tw.horizontalScrollBar().height())
        nrows = tw.rowCount()
        if nrows == 0:
            row_height = 0
            height = header_height+scrollbar_height
        else:
            row_height = tw.rowHeight(0)
            height =  (header_height+scrollbar_height
                       + nrows*row_height)
        if tw == ui.tablewidget_reactions:
            table_min = row_height*min(nrows,3)
            icon_height = sub_icon_height() + 8
            ui.top_frame.setMaximumHeight(height+icon_height)
            ui.top_frame.setMinimumHeight(header_height+icon_height+table_min)
            tw.setMaximumHeight(height)
            tw.setMinimumHeight(header_height + table_min)
            tw.updateGeometry()
            ui.top_frame.updateGeometry()
            ui.bottom_frame.updateGeometry()
            #ui.splitter.setSizes((0,10000)) # https://stackoverflow.com/questions/31969465/how-to-resize-splitter-widgets-programmatically-in-qt

        else:
            tw.setMaximumHeight(height)
            tw.setMinimumHeight(height)
            tw.updateGeometry()


    def chemistry_add_reaction(self):
        ui = self.ui.chemistry
        aliases = list(self.species_all_aliases())
        if not aliases: # No species defined
            self.disable_widget(ui.toolbutton_add_reaction)
            return

        count = 1
        while 'Reaction_%s' % count in self.project.reactions:
            count += 1
        name = 'Reaction_%s' % count
        self.adding_reaction = True
        self.working_reaction =  {'reactants': [],
                                  'products': [],
                                  'chem_eq': ''}

        self.project.reactions[name] = self.working_reaction
        tw = ui.tablewidget_reactions
        self.setup_chemistry() # adds row to table
        row = tw.rowCount()-1
        #tw.cellWidget(row, COL_ENABLE).setChecked(True) #Enable new reaction
        cb = QCheckBox()
        cb.setChecked(True)

        tw.cellWidget(row, COL_ENABLE).setPixmap(cb.grab())
        tw.setCurrentCell(row, COL_RXN_NAME) # and select it
        self.chemistry_toggle_reaction(row, True) #Enable new reaction
        self.set_reaction_edited()
        self.chemistry_update_buttons()

    def chemistry_load_reactions(self, chemkin=False):
        if chemkin:
            filter = "Chemkin files (*.CKI);;All files (*)"
        else:
            filter = ""
        fname, filter = QFileDialog.getOpenFileName(self.ui.chemistry, "Load reactions", ".", filter)

        if not fname:
            return

        if os.path.abspath(os.path.dirname(fname)) != os.path.abspath(os.getcwd()):
            reply = QMessageBox.question(self, "Copy?", "Copy file to project directory?",
                                         QMessageBox.Yes|QMessageBox.No)

            if reply == QMessageBox.Yes:
                copy = True
                base = os.path.basename(fname)
                if os.path.exists(base):
                    reply = QMessageBox.question(self, "Replace?", "%s exists in project directory, replace it?" % base,
                                         QMessageBox.Yes|QMessageBox.No)
                    if reply != QMessageBox.Yes:
                        copy = False

                if copy:
                    copyfile(fname, base)
                    fname = base

        try:
            f = open(fname, encoding='utf8')
            lines = f.readlines()
        except Exception as e:
            self.error(str(e), popup=True)
            return

        if chemkin:
            try:
                lines = chemkin_to_mfix(lines)
            except Exception as e:
                self.error(str(e), popup=True)
                return

        rp = ReactionParser()
        for i,line in enumerate(lines):
            if line.startswith('@'):
                continue
            try:
                rp.parse(line)
            except Exception as e:
                self.error("Error at line %d of %s:\n %s" %
                           (i+1, fname, str(e)),
                           popup=True)
                return
        if len(rp.reactions) == 0:
            self.warning("No reactions found in file.",
                         popup=True)
        else:
            self.print_internal("Loaded %s reactions from %s." %
                                (len(rp.reactions), os.path.basename(fname)),
                                color='blue')
        count = 1+len(self.project.reactions)
        for name,reaction in rp.reactions.items():
            while "Reaction_%s"%count in self.project.reactions:
                count += 1
            new_name = "Reaction_%s"%count
            reaction['phases'] = self.chemistry_reaction_phases(reaction)
            self.project.reactions[new_name] = deepcopy(reaction)

        self.set_unsaved_flag()
        self.setup_chemistry()


    def chemistry_clear_tables(self):
        ui = self.ui.chemistry
        # Work from end, since removing down-shifts entries
        for tw in ui.tablewidget_reactants, ui.tablewidget_products:
            tw.current_row = None
            for row in range(tw.rowCount()-1, -1, -1):
                for col in (2, 1, 0):
                    w = tw.cellWidget(row, col)
                    if w:
                        try:
                            w.value_updated.disconnect()
                        except:
                            pass
                        try:
                            w.activated.disconnect()
                        except:
                            pass
                        try:
                            w.currentIndexChanged.disconnect()
                        except:
                            pass
                        tw.removeCellWidget(row, col)
                        w.deleteLater()
            tw.clearContents()
            tw.setRowCount(0)
            self.fixup_chemistry_table(tw)


    def chemistry_delete_reaction(self, row=None):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        if row is None:
            rows = [i.row() for i in tw.selectedIndexes()]
        else:
            rows = [row]
        if not rows:
            return
        self.chemistry_clear_tables()

        names = [tw.item(row, COL_RXN_NAME).text() for row in rows]
        for name in names:
            self.print_internal(name, font='strikeout')
            self.project.reactions.pop(name, None)
            self.disabled_reactions.pop(name, None)
            # TODO fix up (reindex) keys with 'rate' argument
        self.setup_chemistry()
        self.set_unsaved_flag()
        self.chemistry_update_buttons()


    def chemistry_update_buttons(self):
        # Update the apply/revert/delete buttons
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        label = ui.label_status
        ok = True

        if not self.reaction_edited:
            label.setText('') #"Reaction unmodified")
            label.setStyleSheet("")
        else:

            balanced, message = self.chemistry_check_reaction_balance()
            if not balanced:
                label.setText(message)
                ok = False
            elif self.chemkin_enabled():
                complete, message = self.chemistry_check_chemkin_params()
                if not complete:
                    label.setText(message)
                    ok = False
            if ok:
                label.setText("Apply changes?")

            label.setStyleSheet("color: %s" % ("green" if ok else "red"))

        ui.toolbutton_ok.setEnabled(ok and self.reaction_edited)
        ui.toolbutton_cancel.setEnabled(self.reaction_edited)
        self.set_widget_enabled(ui.toolbutton_add_reaction,
                                len(self.project.reactions)<MAX_RXNS and not self.reaction_edited)
        self.set_widget_enabled(ui.toolbutton_delete_reaction,
                                get_selected_row(tw) is not None and not self.reaction_edited)
        self.set_widget_enabled(tw, not self.reaction_edited, ignore_lock=True)


    def chemistry_apply_changes(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        self.enable_widget(ui.toolbutton_add_reaction)
        self.set_widget_enabled(ui.toolbutton_delete_reaction, row is not None)
        self.clear_reaction_edited()
        self.adding_reaction = False
        self.chemistry_update_buttons()
        if row is None or self.working_reaction is None:
            return
        self.chemistry_check_fracdh(self.working_reaction)
        # Make sure pressure-dependent reaction got updated,
        # because it is complex
        self.chemistry_update_press_rxn_param()
        self.project.reactions[self.current_reaction_name] = deepcopy(self.working_reaction)
        self.chemistry_update_chem_eq(self.current_reaction_name)
        self.chemistry_update_label_n_reactions()
        for line in self.project.format_reaction(self.current_reaction_name):
            if not line.endswith('\n'):
                line += '\n'
            self.print_internal(unicode_chem_eq(line),
                                font='Monospace')
        self.set_unsaved_flag()


    def chemistry_revert_changes(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        row = get_selected_row(tw)
        self.enable_widget(ui.toolbutton_add_reaction)
        if row is not None:
            chem_eq = tw.item(row, COL_CHEM_EQ)
            if (not chem_eq) or chem_eq.text().lower() in ('', 'none', 'disabled'):
                self.chemistry_delete_reaction()
        self.clear_reaction_edited()
        if self.adding_reaction:
            ui.tablewidget_reactions.removeRow(ui.tablewidget_reactions.rowCount()-1)
        self.adding_reaction = False
        self.chemistry_update_buttons()
        self.chemistry_handle_selection()
        self.print_internal('Reaction not saved', color='red')

    def chemistry_find_available_species(self, side, match_phase=None):
        # side is 'reactants' or 'products'
        # if match_phase is passed, species must belong to that phase
        if not self.current_reaction_name:
            return
        #reaction = self.project.reactions.get(self.current_reaction_name)
        reaction = self.working_reaction
        if reaction is None:
            return

        # Collect phase info
        if match_phase is None:
            phases = self.chemistry_reaction_phases(reaction)
        else:
            phases = []
        used = set(k[0] for k in reaction.get(side, []))

        for alias in self.species_all_aliases():
            # This is a bit inefficient - we should already know the phase
            # (but the species list should not be too long)
            species_phase = self.find_species_phase(alias)
            if match_phase is not None and species_phase != match_phase:
                continue
            if len(phases) > 1 and species_phase not in phases:
                continue
            if alias not in used:
                return alias

    def chemistry_add_reactant(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if reaction is None:
            return
        alias = self.chemistry_find_available_species('reactants')
        if not alias:
            return
        reaction['reactants'].append([alias, 1.0])
        self.set_reaction_edited()
        self.chemistry_update_detail_pane()
        # Select the new reaction
        self.chemistry_handle_reactant_selection(row=len(reaction['reactants'])-1)
        if not self.chemistry_find_available_species('reactants'):
            ui.toolbutton_add_reactant.setEnabled(False)


    def chemistry_delete_reactant(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactants
        row = tw.current_row
        if row is None:
            return
        reaction = self.working_reaction
        # reaction = self.project.reactions[self.current_reaction_name]
        del reaction['reactants'][row]
        self.set_reaction_edited()
        self.chemistry_update_detail_pane()
        # Move selection to last row
        n_reactants = len(reaction['reactants'])
        if n_reactants == 0:
            row = None
        elif row > n_reactants-1:
            row = n_reactants-1
        if row is not None:
            tw.setCurrentCell(row, COL_COEFF)
        tw.current_row = row
        self.chemistry_handle_reactant_selection(row=row)
        if self.chemistry_find_available_species('reactants'):
            self.enable_widget(ui.toolbutton_add_reactant)
        self.chemistry_update_label_n_reactions()
        self.chemistry_update_buttons()


    def chemistry_add_product(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if reaction is None:
            return
        alias = self.chemistry_find_available_species('products')
        if not alias:
            return
        reaction['products'].append([alias, 1.0])
        self.set_reaction_edited()
        self.chemistry_update_detail_pane()
        # Select the new product
        self.chemistry_handle_product_selection(row=len(reaction['products'])-1)
        if not self.chemistry_find_available_species('products'):
            ui.toolbutton_add_product.setEnabled(False)


    def chemistry_delete_product(self):
        ui = self.ui.chemistry
        tw = ui.tablewidget_products
        row = tw.current_row
        if row is None:
            return
        reaction = self.working_reaction
        del reaction['products'][row]
        self.set_reaction_edited()
        self.chemistry_update_detail_pane()
        # Move selection to last row
        n_products = len(reaction['products'])
        if n_products == 0:
            row = None
        elif row > n_products-1:
            row = n_products-1
        if row is not None:
            tw.setCurrentCell(row, COL_COEFF)
        tw.current_row = row
        self.chemistry_handle_product_selection(row=row)
        if self.chemistry_find_available_species('products'):
            self.enable_widget(ui.toolbutton_add_product)


    def chemistry_get_species_refs(self, alias):
        ret = []
        for (reaction_name, reaction) in self.project.reactions.items():
            for side in 'reactants', 'products':
                if any(v[0] == alias for v in reaction.get(side,[])):
                    ret.append(reaction_name)
        reaction = self.working_reaction
        if reaction:
            for side in 'reactants', 'products':
                if any(v[0] == alias for v in reaction.get(side,[])):
                    ret.append(self.current_reaction_name)
        return ret


    def chemistry_delete_reactions_of_species(self, alias):
        for name in self.chemistry_get_species_refs(alias):
            self.print_internal(name, font='strikeout')
            self.project.reactions.pop(name, None)
            self.disabled_reactions.pop(name, None)
        self.chemistry_update_buttons()

    def chemistry_rename_species(self, old_name, new_name):
        reaction = self.working_reaction
        if reaction is not None:
            for side in 'reactants', 'products':
                for v in reaction.get(side, []):
                    if v[0] == old_name:
                        v[0] = new_name

        for (reaction_name, reaction) in self.project.reactions.items():
            changed = False
            for side in 'reactants', 'products':
                for v in reaction.get(side, []):
                    if v[0] == old_name:
                        v[0] = new_name
                        changed = True
            if changed:
                if reaction_name is not None:
                    self.chemistry_update_chem_eq(reaction_name)


    def setup_chemistry(self, allow_disabled_tab=False):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions

        if self.reaction_edited:
            # We navigated away while editing.
            # Everything is already set up, leave it as-is
            return

        self.chemistry_update_label_n_reactions()
        # Note, because we clear and reconstruct this tab each time
        #  we lose the current selection
        old_selection = get_selected_row(tw)

        def make_item(sval):
            item = QTableWidgetItem(sval)
            set_item_noedit(item)
            return item

        tw.setRowCount(len(self.project.reactions))
        for row, (name, data) in enumerate(self.project.reactions.items()):
            chem_eq = data.get('chem_eq', '')
            enabled = bool(chem_eq and chem_eq.upper()!='NONE')
            #item = QCheckBox()
            #item.setCheckable(False)
            #item.setEnabled(False)
            cb = QCheckBox()
            cb.setChecked(enabled)
            item = QLabel()
            item.setPixmap(cb.grab())
            # Don't allow enable/disable reactions when paused
            #self.set_widget_enabled(item, self.lock_level==UNLOCKED)  # does not work, handle in event_filter
            item.eventFilter = lambda obj, event, row=row: self.chemistry_event_filter(row, event)
            item.installEventFilter(item)
            item.setToolTip('Enable/disable reaction')
            tw.setCellWidget(row, COL_ENABLE, item)

            item = make_item(name)
            tw.setItem(row, COL_RXN_NAME, item)

            if enabled:
                display_text = unicode_chem_eq(chem_eq)
            else:
                display_text = 'Disabled'
            item = make_item(display_text)
            tw.setItem(row, COL_CHEM_EQ, item)

        # Autoselect if only 1 row
        if tw.rowCount() == 1:
            tw.setCurrentCell(0, COL_RXN_NAME)
        elif old_selection is not None and old_selection < tw.rowCount():
            tw.setCurrentCell(old_selection, COL_RXN_NAME)

        self.chemistry_update_detail_pane()
        self.chemistry_update_buttons()
        self.fixup_chemistry_table(tw, stretch_column=2)#chem eq
        for tw in (ui.tablewidget_reactants, ui.tablewidget_products):
            self.fixup_chemistry_table(tw)

        # Issues/435 don't allow fracdh if not energy_eq
        if old_selection is not None:
            energy_eq = self.project.get_value('energy_eq', default=True)
            enabled = bool(energy_eq)
            ui.groupbox_dh.setEnabled(enabled)
            if not enabled:
                ui.groupbox_dh.setChecked(False)
            ui.groupbox_dh.setToolTip(None if enabled else
                                      'Energy equations must be enabled in order to specify heat of reaction')

        enabled = self.project.get_value('stiff_chemistry', default=False)
        ui.checkbox_stiff_chemistry.setChecked(bool(enabled))
        ui.groupbox_stiff_chemistry.setVisible(bool(enabled))
        self.set_widget_enabled(ui.groupbox_stiff_chemistry, enabled,
                                ignore_lock=(self.lock_level==PARTIAL_LOCK))

        ui.combobox_third_body_species.currentIndexChanged.disconnect()
        ui.combobox_third_body_species.clear()
        ui.combobox_third_body_species.addItems(self.fluid_species)
        ui.combobox_third_body_species.currentIndexChanged.connect(
            self.chemistry_handle_combobox_third_body_species)

        for w in (ui.label_des_min_pmass_frac,
                  ui.lineedit_keyword_des_min_pmass_frac):
            self.set_widget_enabled(w,
                                    self.project.solver not in [SINGLE, TFM],
                                    reason="Discrete phase models only.")


    def chemistry_update_label_n_reactions(self):
        ui = self.ui.chemistry
        n_reactions = sum(1 for (k,v) in self.project.reactions.items()
                          if v.get('chem_eq') or self.disabled_reactions.get(k))
        disabled = sum(1 for x in self.disabled_reactions if x in self.project.reactions)
        text = ("No reactions defined" if n_reactions==0
                else "1 reaction defined" if n_reactions==1
                else "%s reactions defined"%n_reactions)
        if disabled:
            text += ', %s disabled' % disabled
        ui.label_n_reactions.setText(text)

    def chemistry_extract_info(self):
        """extract additional chemistry info after loading project file"""
        for reaction in self.project.reactions.values():
            self.chemistry_check_fracdh(reaction) # sets 'phases' field


    def chemistry_event_filter(self, row, event):
        # Don't toggle reaction if we're just selecting the row
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        if event.type() not in (QEvent.MouseButtonPress,
                                QEvent.MouseButtonDblClick):
            return False
        if self.lock_level: # Don't enable/disable when paused or running
            return False
        sel = get_selected_rows(tw)
        if row not in sel: #Selection changed
            return False

        name = tw.item(row,1).text()
        chem_eq = self.project.reactions.get(name, {}).get('chem_eq')
        enabled = bool(chem_eq and chem_eq.upper()!='NONE')
        enabled = not enabled # Toggle
        cb = QCheckBox()
        cb.setChecked(enabled)
        pm = cb.grab()
        for s in sel:
            self.chemistry_toggle_reaction(s, enabled)
            tw.cellWidget(s,0).setPixmap(pm)
        return True

    def set_reaction_edited(self):
        if not self.reaction_edited:
            self.reaction_edited = True
            #import traceback
            #traceback.print_stack()

    def clear_reaction_edited(self):
        self.reaction_edited = False

    def chemistry_toggle_reaction(self, row, enabled):
        ui = self.ui.chemistry
        tw = ui.tablewidget_reactions
        name = list(self.project.reactions.keys())[row]
        reaction = self.project.reactions[name]
        chem_eq = reaction.get('chem_eq')
        self.set_widget_enabled(ui.detail_pane, enabled)
        if chem_eq == 'NONE' and not enabled:
            #Already disabled
            return
        if chem_eq and chem_eq != 'NONE' and enabled:
            #Already enabled
            return
        if enabled:
            reaction['chem_eq'] = chem_eq = self.disabled_reactions.pop(name, 'NONE')
            display_text = unicode_chem_eq(chem_eq)
            tw.item(row, COL_CHEM_EQ).setText(chem_eq)
        else:
            self.disabled_reactions[name] = reaction.get('chem_eq', 'NONE')
            display_text = "Disabled"
            reaction['chem_eq'] = 'NONE'
        tw.item(row, COL_CHEM_EQ).setText(display_text)
        self.chemistry_update_label_n_reactions()
        self.set_unsaved_flag()


    def chemistry_check_fracdh(self, reaction):
        phases = self.chemistry_reaction_phases(reaction)
        reaction['phases'] = phases
        num_phases = len(phases)
        dh = reaction.get('dh')
        key = 'fracdh'
        if dh is None:
            reaction.pop(key, None)
            return
        fracdh = reaction.get(key, {})
        for p in phases:
            if p not in fracdh:
                fracdh[p] = None

        if num_phases == 0:
            reaction.pop(key, None)
        elif num_phases == 1:
            reaction[key] = {phases[0]: 1.0}

        else:
            n_holes = list(fracdh.values()).count(None)
            if n_holes == 0:
                pass # Nothing to do
            else:
                s = sum(v or 0 for v in fracdh.values())
                for p in phases:
                    if fracdh.get(p) is None:
                        fracdh[p] = (1-s)/n_holes
            # Modify last value so they sum to 1.0
            s = sum(fracdh.values())
            if s < 1:
                fracdh[phases[-1]] += (1-s)
            idx = len(phases)-1
            s = sum(fracdh.values())
            while idx >=0  and s > 1:
                fracdh[phases[idx]] -= (s-1)
                if fracdh[phases[idx]] < 0:
                    fracdh[phases[idx]] = 0.0
                    s = sum(fracdh.values())
                    idx -= 1
                else:
                    break


            if fracdh[phases[-1]] < 0:
                fracdh[phases[-1]] = 0

        for (k,v) in fracdh.items():
            if v is not None:
                fracdh[k] = round(v, 6)


    def chemistry_to_str(self):
        # This is only for the data saved in the #!MFIX-GUI section
        # Don't save deleted disabled reactions
        disabled = {k:v for (k,v) in self.disabled_reactions.items()
                    if k in self.project.reactions}
        if disabled:
            data = {'disabled_reactions': disabled}
            return JSONEncoder().encode(data)
        return ''


    def chemistry_from_str(self, s):
        if not s:
            return
        data = JSONDecoder().decode(s)

        if data:
            RP = ReactionParser()
            val = data.get('disabled_reactions')
            if val:
                self.disabled_reactions = val
                for (name, eq) in self.disabled_reactions.items():
                    if name in self.project.reactions:
                        reaction = self.project.reactions[name]
                        reaction['reactants'], reaction['products'] = RP.parse_chem_eq(eq)
                        self.chemistry_check_fracdh(reaction)


    def reset_chemistry(self):
        ui = self.ui.chemistry
        self.current_reaction_name = None
        self.clear_reaction_edited()
        self.adding_reaction = False
        self.reaction_mass_totals = [None, None]
        self.working_reaction = None
        self.project.reactions.clear() # done in project.reset()
        self.disabled_reactions.clear()
        self.chemistry_clear_tables()
        tw = ui.tablewidget_reactions
        tw.clearSelection()
        tw.clearContents()
        tw.setRowCount(0)
        self.fixup_chemistry_table(tw, stretch_column=2)
        tw.current_row = None
        tw.setEnabled(True)
        ui.toolbutton_add_reaction.setEnabled(True)
        ui.toolbutton_delete_reaction.setEnabled(False)
        for tb in (ui.toolbutton_ok, ui.toolbutton_cancel):
            tb.setEnabled(False)
###
### Chemkin additions (Arrhenius model)
###
    def chemkin_enabled(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if reaction is None:
            return False
        key = ('arrhenius_rrates_des' if reaction_is_des(self.project, reaction)
               else 'arrhenius_rrates_fluid')
        return self.project.get_value(key, default=False)

    def chemistry_chemkin_warning(self):
        arrhenius_rrates_fluid = self.project.get_value('arrhenius_rrates_fluid',
                                                        default=False)
        arrhenius_rrates_des = self.project.get_value('arrhenius_rrates_des',
                                                        default=False)
        usr_rates_file = None
        usr_rates_des_file = None
        files = os.listdir()
        for f in files:
            if f.lower() == 'usr_rates.f':
                usr_rates_file = f
            if f.lower() == 'usr_rates_des.f':
                usr_rates_des_file = f
        if arrhenius_rrates_fluid and usr_rates_file:
            self.warning('ARRHENIUS_RRATES_FLUID is enabled and file %s is present.\n%s will not be used.' % (usr_rates_file, usr_rates_file),
                         popup=True)
        if arrhenius_rrates_des and usr_rates_des_file:
            self.warning('ARRHENIUS_RRATES_DES is enabled and file %s is present.\n%s will not be used.' % (usr_rates_des_file, usr_rates_des_file),
                         popup=True)
        fluid_reactions = []
        des_reactions = []
        for r in self.project.reactions.values():
            chem_eq = r.get('chem_eq')
            if not chem_eq:
                continue
            if chem_eq.lower() == 'none':
                continue
            if not r.get('phases'):
                phases = self.chemistry_reaction_phases(r)
                r['phases'] = phases
            if reaction_is_des(self.project, r):
                des_reactions.append(r)
            else:
                fluid_reactions.append(r)

        if fluid_reactions and not arrhenius_rrates_fluid and not usr_rates_file:
            self.warning("Fluid reactions are defined, but ARRHENIUS_RRATES_FLUID is not enabled and usr_rates.f is not present.  Reaction rates will not be calculated.", popup=True)

        if des_reactions and not arrhenius_rrates_des and not usr_rates_des_file:
            self.warning("DES reactions are defined, but ARRHENIUS_RRATES_DES is not enabled and usr_rates_des.f is not present.  Reaction rates will not be calculated.")


        def check_arrhenius(reaction):
            press_rxn_param = reaction.get('press_rxn_param')
            if press_rxn_param and 'plog' in press_rxn_param.lower():
                # PLOG reactions do not require arrhenius_coeff
                return True
            arr = reaction.get('arrhenius_coeff')
            if not arr:
                return False
            tok = arr.split()
            if len(tok) != 3 or 'None' in tok:
                return False
            return True

        if arrhenius_rrates_fluid:
            if not all(check_arrhenius(r) for r in fluid_reactions):
                self.warning("ARRHENIUS_RRATES_FLUID is enabled.  Arrhenius coefficients must be defined for all fluid reactions.", popup=True)

        if arrhenius_rrates_des:
            if not all(check_arrhenius(r) for r in des_reactions):
                self.warning("ARRHENIUS_RRATES_DES is enabled.  Arrhenius coefficients must be defined for all DES reactions.", popup=True)







    def chemistry_check_chemkin_params(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not reaction:
            return False, "No reaction"
        # Arrhenius
        if not 'PLOG' in reaction.get('press_rxn_param', ''):
            arr = reaction.get('arrhenius_coeff')
            if not arr or 'None' in arr:
                return False, "Arrhenius coefficients must be specified"

        # Reverse reaction: assume complete due to combobox

        # Reaction order:
        if ui.checkbox_reaction_order.isChecked():
            rxn_order = reaction.get("rxn_order")
            if not rxn_order or ':None' in rxn_order:
                return False, "Reaction order must be specified for all reactants"

        # That's all for DES reactions
        if reaction_is_des(self.project, reaction):
            return True, ''

        # Third-body
        if ui.checkbox_third_body.isChecked():
            gb = ui.groupbox_third_body
            if ui.combobox_third_body.currentIndex() == 0: #M_all
                if ui.combobox_third_body_add_species.currentIndex() != 0:
                    le = ui.lineedit_third_body_effcy
                    if le.text() == '':
                        return False, "Specify efficiency for " + le.key

                layout = gb.layout()
                for row in range(layout.rowCount()):
                    item = layout.itemAtPosition(row, 2)
                    if not item:
                        continue
                    le = item.widget()
                    if le.isEnabled() and le.text() == '':
                        return False, "Specify efficiency for " + le.key
            else:
                le = ui.lineedit_third_body_effcy
                if le.text() == '':
                    return False, "Specify efficiency for " + le.key

        # Pressure-dependent
        if ui.checkbox_pressure_dependent.isChecked():
            press_rxn_param = reaction.get('press_rxn_param')
            if not press_rxn_param: #?
                return False, "Specify parameters for pressure-dependent reaction"
            tok = press_rxn_param.split()
            model = tok[0].lower()
            param = self.parse_press_rxn_param(tok[1:])
            if model == 'plog':
                tw = ui.tablewidget_plog
                n_complete_rows = 0
                n_incomplete_rows = 0
                # Blank rows are OK
                for row in range(tw.rowCount()):
                    data = [None if not tw.item(row, col) else (tw.item(row,col).text() or None)
                            for col in range(1,5)]
                    if None in data:
                        if data == [None, None, None, None]:
                            continue
                        else:
                            n_incomplete_rows += 1
                    else:
                        n_complete_rows += 1
                if n_incomplete_rows or n_complete_rows==0:
                    return False, "PLOG table is incomplete"

            else: # Not PLOG
                arr = param.get('arrhenius_press')
                if not arr or None in arr:
                    return False, "Pressure Arrhenius coefficients must be specified"
                press_coeff = param.get('press_coeff')
                if model.startswith('lindemann'):
                    pass
                elif model.startswith('sri'):
                    if not press_coeff or None in press_coeff:
                        return False, "SRI coefficients must be specified"
                elif model.startswith('troe'):
                    if not press_coeff or None in press_coeff:
                        return False, "Troe coefficients must be specified"

        # Landau-Teller
        if ui.checkbox_landau_teller.isChecked():
            lt_coeff = reaction.get('lt_coeff')
            if not lt_coeff or 'None' in lt_coeff:
                return False, "Landau-Teller coefficients must be specified"

        # Rate fit
        if ui.checkbox_rate_fit.isChecked():
            fit_coeff = reaction.get('fit_coeff')
            if not fit_coeff or 'None' in fit_coeff:
                return False, "Rate-fitting coefficients must be specified"

        return True, ''

### Arrhenius
    def chemistry_handle_lineedit_arrhenius(self, widget, update_dict, args):
        reaction = self.working_reaction
        if not reaction:
            return
        ui = self.ui.chemistry
        key, val = update_dict.popitem()
        arr = reaction.get('arrhenius_coeff')
        if not arr:
            c = [None, None, None]
        else:
            c = [None if x=='None' else FloatExp(x) for x in arr.split()]
        if widget == ui.lineedit_arrhenius_A:
            c[0] = val
        elif widget == ui.lineedit_arrhenius_B:
            c[1] = val
        elif widget == ui.lineedit_arrhenius_Ea:
            c[2] = val
        reaction['arrhenius_coeff'] = '%s %s %s' % tuple(None if x in ('', None) else FloatExp(x) for x in c)
        self.set_reaction_edited()
        self.chemistry_update_buttons()


    def chemistry_update_arrhenius(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not self.chemkin_enabled():
            ui.groupbox_arrhenius.setVisible(False)
            return
        ui.groupbox_arrhenius.setVisible(True)
        key = 'arrhenius_coeff'
        val = reaction.get(key)

        if not val:
            for w in (ui.lineedit_arrhenius_A,
                      ui.lineedit_arrhenius_B,
                      ui.lineedit_arrhenius_Ea):
                w.setText('')
        else:
            A,B,Ea = (None if x is None else FloatExp(x) for x in val.split())
            ui.lineedit_arrhenius_A.updateValue(None, A)
            ui.lineedit_arrhenius_B.updateValue(None, B)
            ui.lineedit_arrhenius_Ea.updateValue(None, Ea)


### Reverse reaction
    def chemistry_handle_checkbox_reverse_reaction(self, checked):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not reaction:
            return
        combo = ui.combobox_reverse_reaction
        combo.setVisible(checked)
        ui.label_reverse_reaction.setVisible(checked)
        if not checked:
            reaction.pop('reverse_calc', None)
            return
        val = reaction.get('reverse_calc', 'fromForwardRateConstant') # This should not be set yet
        if val.lower()=='fromarrheniuscoeff':
            combo.setCurrentIndex(1)
        else:
            combo.setCurrentIndex(0)
        reaction['reverse_calc'] = val
        self.set_reaction_edited()
        self.chemistry_update_buttons()

    def chemistry_handle_combobox_reverse_reaction(self, index):
        reaction = self.working_reaction
        if not reaction:
            return
        reaction['reverse_calc'] = ('fromForwardRateConstant', 'fromArrheniusCoeff')[index]
        set_combobox_tooltip(self.ui.chemistry.combobox_reverse_reaction)
        self.set_reaction_edited()
        self.chemistry_update_buttons()

    def chemistry_update_reverse_reaction(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        check =  ui.checkbox_reverse_reaction
        combo = ui.combobox_reverse_reaction

        if not self.chemkin_enabled():
            check.setVisible(False)
            combo.setVisible(False)
            ui.label_reverse_reaction.setVisible(False)
            return

        val = reaction.get('reverse_calc')
        if val and val.lower() not in ('fromforwardrateconstant',
                               'fromarrheniuscoeff'):
            self.error("Invalid value for 'reverse_calc': %s" % val,
                       popup=True)
            reaction.pop('reverse_calc', None)
            self.set_reaction_edited()
            self.chemistry_update_buttons()
            return
        check.setVisible(True)
        check.setEnabled(True)
        if not val:
            combo.setVisible(False)
            check.setChecked(False)
            ui.label_reverse_reaction.setVisible(False)
            return
        check.setChecked(True)
        combo.setVisible(True)
        ui.label_reverse_reaction.setVisible(True)
        combo.currentIndexChanged.disconnect()
        combo.setCurrentIndex(0 if val.lower()=='fromforwardrateconstant' else 1)
        combo.currentIndexChanged.connect(
            self.chemistry_handle_combobox_reverse_reaction)

### Reaction order
    def chemistry_handle_checkbox_reaction_order(self, checked):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not reaction:
            return
        gb = ui.groupbox_reaction_order
        if not checked:
            gb.setVisible(False)
            reaction.pop('rxn_order', None)
        else:
            reactants = reaction.get('reactants', [])
            if 'rxn_order' not in reaction: # Should not be set if we get here
                reaction['rxn_order'] = ' '.join("%s:%s"%(r[0], r[1]) for r in reactants)
            self.chemistry_update_reaction_order()
        self.set_reaction_edited()
        self.chemistry_update_buttons()


    def chemistry_handle_lineedit_reaction_order(self, le, val, args):
        reaction = self.working_reaction
        if not reaction:
            return
        reaction_order = reaction.get('rxn_order')
        reactants = reaction.get('reactants', [])
        if not reaction_order:
            reaction_order = ' '.join("%s:%s"%(r[0], r[1]) for r in reactants)
        d = dict(x.split(':',1) for x in reaction_order.split())
        d.update(val)
        def _none(s):
            return 'None' if s in (None, "None", '') else s
        reaction_order = ' '.join("%s:%s"%(r[0], _none(d.get(r[0], r[1]))) for r in reactants)
        if reaction.get('rxn_order') != reaction_order:
            reaction['rxn_order'] = reaction_order
            self.set_reaction_edited()
            self.chemistry_update_buttons()


    def chemistry_update_reaction_order(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not self.chemkin_enabled():
            ui.checkbox_reaction_order.setVisible(False)
            ui.groupbox_reaction_order.setVisible(False)
            return
        ui.checkbox_reaction_order.setVisible(True)
        ui.checkbox_reaction_order.setEnabled(True)
        val = reaction.get('rxn_order')
        ui.checkbox_reaction_order.setChecked(bool(val))
        gb = ui.groupbox_reaction_order
        if not val:
            gb.setVisible(False)
            return
        # Delete all items
        layout = gb.layout()
        while layout.count() > 0:
            item = layout.itemAt(0)
            if not item:
                continue
            w = item.widget()
            if not w:
                continue
            w.hide()
            layout.removeWidget(w)
            w.setParent(None)
            w.deleteLater()

        # And make new ones
        reaction_order = reaction.get('rxn_order', '')
        d = dict(x.split(':',1) for x in reaction_order.split())
        for row,(name, coeff) in enumerate(reaction['reactants']):
            label = QLabel(name)
            le = LineEdit()
            le.key = name
            le.dtype = float
            val = d.get(name, coeff)
            le.updateValue(name, val)
            le.min = 0
            layout.addWidget(label, row, 0, 1, 1)
            layout.addWidget(le, row, 1, 1, 1)
            le.value_updated.connect(self.chemistry_handle_lineedit_reaction_order)
        gb.setVisible(True)




### Third-body reaction

    def chemistry_update_third_body(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not self.chemkin_enabled():
            ui.checkbox_third_body.setVisible(False)
            ui.groupbox_third_body.setVisible(False)
            return
        if not self.fluid_species or reaction_is_des(self.project, reaction):
            ui.checkbox_third_body.setVisible(False)
            ui.groupbox_third_body.setVisible(False)
            reaction.pop("third_body_param", None)
            return

        ui.checkbox_third_body.setVisible(True)
        ui.checkbox_third_body.setEnabled(True)
        val = reaction.get('third_body_param')
        ui.checkbox_third_body.setChecked(bool(val))
        gb = ui.groupbox_third_body
        if not val:
            gb.setVisible(False)
            return
        gb.setVisible(True)

        gbe = ui.groupbox_third_body_effcy
        # Delete all items...
        layout = gbe.layout()
        while layout.count() > 0:
            item = layout.itemAt(0)
            if not item:
                continue
            w = item.widget()
            if not w:
                continue
            w.hide()
            layout.removeWidget(w)
            if w not in (ui.combobox_third_body_add_species,
                         ui.lineedit_third_body_effcy):
                w.setParent(None)
                w.deleteLater()

        third_body_param = reaction.get('third_body_param','').strip()
        data = {}

        tok = third_body_param.split(maxsplit=1)
        if len(tok)==2:
            param = dequote(tok[1])
            if param:
                data = dict(x.split(':',1) for x in param.split())

        if third_body_param.lower().startswith('m_all'):
            ui.combobox_third_body.currentIndexChanged.disconnect()
            ui.combobox_third_body.setCurrentIndex(0)
            ui.combobox_third_body.currentIndexChanged.connect(
                self.chemistry_handle_combobox_third_body)

            ui.combobox_third_body_species.setVisible(False)
            # ... and make new ones
            names = list(self.fluid_species.keys())
            row = 0
            h = ui.combobox_third_body_species.fontMetrics().height()
            for name in self.fluid_species:
                val = data.get(name, 1.0)
                if val is None or val == 'None':
                    continue
                if float(val) == 1.0:
                    continue
                label = QLabel(name)
                le = LineEdit()
                le.key = name
                names.remove(name)
                le.updateValue(name, float(val))
                le.min = 0
                b = QToolButton()
                b.setIcon(get_icon('remove.svg'))
                # Make button obviously a button
                #b.setAutoRaise(True)
                b.clicked.connect(lambda checked, row=row:  self.chemistry_delete_third_body_row(row))
                layout.addWidget(b, row, 0, 1, 1)
                layout.addWidget(label, row, 1, 1, 1)
                layout.addWidget(le, row, 2, 1, 1)
                le.value_updated.connect(self.chemistry_handle_lineedit_third_body)
                row += 1
            if names:
                layout.addWidget(ui.combobox_third_body_add_species, row, 1, 1, 1)
                ui.combobox_third_body_add_species.show()
                ui.combobox_third_body_add_species.currentIndexChanged.disconnect()
                ui.combobox_third_body_add_species.clear()
                ui.combobox_third_body_add_species.addItem("+")
                ui.combobox_third_body_add_species.addItems(names)
                ui.combobox_third_body_add_species.currentIndexChanged.connect(
                    self.chemistry_handle_combobox_third_body_add_species)
                le = ui.lineedit_third_body_effcy
                le.value_updated.disconnect()
                le.setText('')
                le.value_updated.connect(self.chemistry_handle_lineedit_third_body)
                layout.addWidget(le, row, 2, 1, 1)
                le.show()
                le.setEnabled(False)

        else: # M_species
            if not third_body_param.lower().startswith('m_'):
                self.error("Invalid setting for third_body_param:  %s" % third_body_param, popup=True)
                reaction.pop('third_body_param', None)
                self.set_reaction_edited()
                self.chemistry_update_third_body()
                self.chemistry_update_buttons()
                return
            ui.combobox_third_body.currentIndexChanged.disconnect()
            ui.combobox_third_body.setCurrentIndex(1)
            ui.combobox_third_body.currentIndexChanged.connect(
                self.chemistry_handle_combobox_third_body)

            name = tok[0][2:]

            cb = ui.combobox_third_body_species
            cb.currentIndexChanged.disconnect()
            cb.setCurrentIndex(list(self.fluid_species.keys()).index(name))
            cb.currentIndexChanged.connect(
                self.chemistry_handle_combobox_third_body_species)
            cb.setVisible(True)
            row = 0
            label = QLabel(name)
            le = ui.lineedit_third_body_effcy
            le.setEnabled(True)
            le.key = name
            le.dtype = float
            le.min = 0
            val = data.get(name)
            le.value_updated.disconnect()
            if val:
                le.saved_value = val
                le.updateValue(name, float(val))
            else:
                le.setText('')
            le.value_updated.connect(self.chemistry_handle_lineedit_third_body)
            layout.addWidget(label, row, 0, 1, 1)
            layout.addWidget(le, row, 1, 1, 1)
            le.setVisible(True)

        gb.setVisible(True)


    def chemistry_delete_third_body_row(self, row):
        ui = self.ui.chemistry
        self.reaction_modified = True
        layout = ui.groupbox_third_body_effcy.layout()
        name = layout.itemAtPosition(row, 1).widget().text()
        for col in 2,1,0:
            item = layout.itemAtPosition(row, col)
            w = item.widget()
            layout.removeWidget(w)
            w.hide()
            w.setParent(None)
            w.deleteLater()
        self.chemistry_update_third_body_param()
        self.chemistry_update_third_body()
        self.chemistry_update_buttons()


    def chemistry_handle_checkbox_third_body(self, checked):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not reaction:
            return
        if not checked:
            reaction.pop('third_body_param', None)
        else:
            reaction['third_body_param'] = 'M_all'

        self.chemistry_update_third_body()
        self.set_reaction_edited()
        self.chemistry_update_buttons()

    def chemistry_handle_combobox_third_body(self, index):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not reaction:
            return
        if index == 0:
            # All fluid-phase reactants
            third_body_param = 'M_all '
            reaction['third_body_param'] = third_body_param
        else:
            if self.fluid_species:
                third_body_param = 'M_'+list(self.fluid_species.keys())[0]
                reaction['third_body_param'] = third_body_param
            else: # What to do?
                reaction.pop('third_body_param', None)
        self.set_reaction_edited()
        self.chemistry_update_third_body()
        self.chemistry_update_buttons()


    def chemistry_handle_combobox_third_body_species(self, index):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not reaction:
            return
        if not self.fluid_species: # Should not get here
            if reaction.pop('third_body_param', None):
                self.set_reaction_edited()
                self.chemistry_update_buttons()
            self.update_third_body()
            return
        third_body_param = reaction.get('third_body_param')
        reaction['third_body_param'] = 'M_%s' % ui.combobox_third_body_species.currentText()
        self.chemistry_update_third_body()
        ui.lineedit_third_body_effcy.setText('')
        self.set_reaction_edited()
        self.chemistry_update_buttons()


    def chemistry_handle_lineedit_third_body(self, le, val, args):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not reaction:
            return

        m_all = (ui.combobox_third_body.currentIndex() == 0)
        if m_all and le == ui.lineedit_third_body_effcy:
            if ui.combobox_third_body_add_species.currentIndex() == 0: # <Species>
                le.setText('')
                return
        self.chemistry_update_third_body_param()
        if m_all:
            self.chemistry_update_third_body()
        self.chemistry_update_buttons()


    def chemistry_update_third_body_param(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not reaction:
            return
        m_all = (ui.combobox_third_body.currentIndex() == 0)
        if m_all:
            third_body_param = 'M_all'
            gb = ui.groupbox_third_body_effcy
            layout = gb.layout()
            param = []
            for row in range(layout.rowCount()):
                item = layout.itemAtPosition(row,1)
                if not item:
                    continue
                w = item.widget()
                if isinstance(w, QComboBox):
                    text = w.currentText()
                    if text=='+':
                        continue
                else:
                    text = w.text()
                le = layout.itemAtPosition(row,2).widget()
                if le.value in ('', None):
                    value = None
                else:
                    value = le.value
                if value == 1:
                    continue
                param.append('%s:%s' % (text, value))

            if param:
                third_body_param = third_body_param + ' "' + ' '.join(param) + '"'
        else: # M_species
            val = ui.lineedit_third_body_effcy.value
            if val == '':
                return
            name = ui.combobox_third_body_species.currentText()
            third_body_param = 'M_%s "%s:%s"' % (name, name, val)
        if reaction.get('third_body_param') != third_body_param:
            reaction['third_body_param'] = third_body_param
            self.set_reaction_edited()

    def chemistry_handle_combobox_third_body_add_species(self, idx):
        ui = self.ui.chemistry
        ui.lineedit_third_body_effcy.setEnabled(idx != 0)
        ui.lineedit_third_body_effcy.key = ui.combobox_third_body_add_species.currentText()
        if idx != 0:
            self.set_reaction_edited()
        self.chemistry_update_buttons()

### Pressure-dependent reaction
    def chemistry_update_pressure(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        cb = ui.checkbox_pressure_dependent
        gb = ui.groupbox_pressure_dependent
        if not self.chemkin_enabled():
            cb.setVisible(False)
            gb.setVisible(False)
            return
        if reaction_is_des(self.project, reaction):
            reaction.pop('press_rxn_param', None)
            cb.setVisible(False)
            gb.setVisible(False)
            return

        cb.setVisible(True)
        press_rxn_param = reaction.get('press_rxn_param','')
        cb.setChecked(bool(press_rxn_param))
        gb.setVisible(bool(press_rxn_param))
        if not press_rxn_param:
            return
        tok = press_rxn_param.split()
        model = tok[0].lower()
        params = self.parse_press_rxn_param(tok[1:])
        cb = ui.combobox_pressure_model
        if model.endswith('_falloff'):
            idx = 0
        elif model.endswith('_bimo'):
            idx = 1
        elif model=='plog':
            idx = 2
        else:
            self.error('Invalid pressure model: "%s"' % tok[0],
                       popup=True)
            return
        ui.combobox_pressure_model.currentIndexChanged.disconnect()
        ui.combobox_pressure_model.setCurrentIndex(idx)
        ui.combobox_pressure_model.currentIndexChanged.connect(self.chemistry_handle_combobox_pressure_model)
        self.chemistry_handle_combobox_pressure_model(idx, in_setup=True)

        if model == 'plog':
            # Populate table
            tw = ui.tablewidget_plog
            tw.cellChanged.disconnect()
            press_coeff = params.get('press_coeff',[])
            tw.clearContents()
            tw.setRowCount(len(press_coeff) + 1) # Always leave a blank row to add more
            for (i,x) in enumerate(press_coeff):
                b = QToolButton()
                b.setIcon(get_icon('remove.svg'))
                b.clicked.connect(lambda checked, row=i: (tw.removeRow(row),
                                                     self.chemistry_update_press_rxn_param(),
                                                     self.chemistry_update_pressure(),
                                                     self.set_reaction_edited(),
                                                     self.chemistry_update_buttons()))
                # Make button obviously a button
                #b.setAutoRaise(True)
                tw.setCellWidget(i,0,b)
                item = tw.item(i,1)
                if not item:
                    item = QTableWidgetItem()
                    tw.setItem(i,1, item)
                item.setText(str(FloatExp(x)))
            i = 0
            j = 2
            arrhenius_press = params.get('arrhenius_press', [])
            while arrhenius_press:
                item = tw.item(i, j)
                if not item:
                    item = QTableWidgetItem()
                    tw.setItem(i, j, item)
                item.setText(str(FloatExp(arrhenius_press[0])))
                arrhenius_press.pop(0)
                j += 1
                if j == 5:
                    j = 2
                    i += 1
                if i == tw.rowCount():
                    break
            self.fixup_chemistry_table(tw)
            tw.cellChanged.connect(self.chemistry_handle_plog)
        else:
            if model.startswith('lindemann_'):
                idx = LINDEMANN
            elif model.startswith('troe_'):
                idx = TROE
            elif model.startswith('sri_'):
                idx = SRI
            else:
                self.error('Invalid pressure model: "%s"' % tok[0],
                           popup=True)
                return
            ui.combobox_pressure_method.currentIndexChanged.disconnect()
            ui.combobox_pressure_method.setCurrentIndex(idx)
            ui.combobox_pressure_method.currentIndexChanged.connect(
                self.chemistry_handle_combobox_pressure_method)
            self.chemistry_handle_combobox_pressure_method(idx, in_setup=True)

            arrhenius_press = params.get('arrhenius_press',[None, None, None])
            for i,s in enumerate(("A", "B", "Ea")):
                w = getattr(ui, "lineedit_pressure_arrhenius_"+s)
                w.updateValue(None, arrhenius_press[i])

            if model.startswith('troe_'):
                press_coeff = params.get('press_coeff',[None,None,None,0])
                if len(press_coeff) < 4:
                    press_coeff.extend([0.0]*(4-len(press_coeff)))
                for i,s in enumerate(('alpha', 't3star', 'tstar', 't2star')):
                    w = getattr(ui, "lineedit_troe_"+s)
                    w.updateValue(None, press_coeff[i])

            if model.startswith('sri_'):
                press_coeff = params.get('press_coeff',[None, None, None, 1, 0])
                if len(press_coeff) == 4:
                    press_coeff.append(0)
                elif len(press_coeff) == 3:
                    press_coeff.extend([1,0])
                elif len(press_coeff) < 3:
                    press_coeff.extend([None] * (3-len(press_coeff)))
                    press_coeff.extend([1,0])

                for i,s in enumerate('abcde'):
                    w = getattr(ui, "lineedit_sri_"+s)
                    w.updateValue(None, press_coeff[i])


    def chemistry_handle_checkbox_pressure_dependent(self, checked):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not reaction:
            return
        ui.groupbox_pressure_dependent.setVisible(checked)
        if checked:
            ui.combobox_pressure_model.setCurrentIndex(0)
            ui.combobox_pressure_method.setCurrentIndex(0)
            # Clear out all inputs
            for s in ('pressure_arrhenius_A', 'pressure_arrhenius_B', 'pressure_arrhenius_Ea',
                      'troe_alpha', 'troe_t3star', 'troe_tstar', 'troe_t2star',
                      'sri_a', 'sri_b', 'sri_c', 'sri_d', 'sri_e'):
                w = getattr(ui, 'lineedit_'+s)
                w.updateValue(None, 1 if s== 'sri_d' else 0 if s in ('sri_e', 'troe_t2star') else None)
            #self.chemistry_handle_combobox_pressure_model(ui.combobox_pressure_model.currentIndex(),
            #                                              in_setup=True)
            #self.chemistry_handle_combobox_pressure_method(ui.combobox_pressure_method.currentIndex(),
            #                                               in_setup=True)
            self.chemistry_update_press_rxn_param()
            self.chemistry_update_pressure()
        else:
            ui.groupbox_arrhenius.setVisible(True)
            reaction.pop('press_rxn_param', None)
        self.set_reaction_edited()
        self.chemistry_update_buttons()


    def chemistry_handle_combobox_pressure_model(self, idx, in_setup=False):
        ui = self.ui.chemistry
        set_combobox_tooltip(ui.combobox_pressure_model)
        ui.groupbox_arrhenius.setVisible(idx != PLOG)

        if idx == UNIMOL: # Fall-off
            ui.label_pressure_method.setVisible(True)
            ui.combobox_pressure_method.setVisible(True)
            ui.groupbox_pressure_arrhenius.setVisible(True)
            ui.groupbox_pressure_coeff.setVisible(True)
            ui.tablewidget_plog.setVisible(False)
        elif idx == BIMOL: # Chemically-activated
            ui.label_pressure_method.setVisible(True)
            ui.combobox_pressure_method.setVisible(True)
            ui.groupbox_pressure_arrhenius.setVisible(True)
            ui.groupbox_pressure_coeff.setVisible(True)
            ui.tablewidget_plog.setVisible(False)
        elif idx == PLOG: # PLOG
            ui.label_pressure_method.setVisible(False)
            ui.combobox_pressure_method.setVisible(False)
            ui.groupbox_pressure_arrhenius.setVisible(False)
            ui.groupbox_pressure_coeff.setVisible(False)
            ui.tablewidget_plog.setVisible(True)
        if idx != PLOG:
            self.chemistry_handle_combobox_pressure_method(
                ui.combobox_pressure_method.currentIndex(),
                in_setup=in_setup)
        if not in_setup:
            self.set_reaction_edited()
            self.chemistry_update_press_rxn_param()
            self.chemistry_update_buttons()

    def chemistry_handle_combobox_pressure_method(self, idx, in_setup=False):
        ui = self.ui.chemistry
        set_combobox_tooltip(ui.combobox_pressure_method)
        if idx == LINDEMANN:
            ui.groupbox_pressure_coeff.setVisible(False)
        elif idx == TROE:
            ui.groupbox_pressure_coeff.setVisible(True)
            ui.groupbox_pressure_troe.setVisible(True)
            ui.groupbox_pressure_sri.setVisible(False)
            if not in_setup:
                for s in ('alpha', 't3star', 'tstar', 't2star'):
                    w = getattr(ui, 'lineedit_troe_'+s)
                    w.updateValue(None, 0 if s=='t2star' else None)
        elif idx == SRI:
            ui.groupbox_pressure_coeff.setVisible(True)
            ui.groupbox_pressure_troe.setVisible(False)
            ui.groupbox_pressure_sri.setVisible(True)
            for s in ('abcde'):
                w = getattr(ui, 'lineedit_sri_'+s)
                w.updateValue(None, 1 if s=='d' else 0 if s=='e' else None)
        if not in_setup:
            self.set_reaction_edited()
            self.chemistry_update_press_rxn_param()
            self.chemistry_update_pressure()
            self.chemistry_update_buttons()


    def chemistry_handle_pressure_coeff(self, widget, update, args):
        reaction = self.working_reaction
        if not reaction:
            return
        self.chemistry_update_press_rxn_param()
        self.set_reaction_edited()
        self.chemistry_update_buttons()

    def chemistry_handle_plog(self, row, col):
        ui = self.ui.chemistry
        tw = ui.tablewidget_plog
        item = tw.item(row, col)
        if item:
            text = item.text().strip()
            if text:
                try:
                    val = float(text)
                except ValueError:
                    self.error("Invalid number: %s" % text, popup=True)
                    item.setText('')
                    return
        # TODO check all rows are complete
        empty_cell = False
        for r in range(tw.rowCount()):
            if empty_cell:
                break
            for c in range(1, 5):
                item = tw.item(r,c)
                if not item or item.text().strip()=='':
                    empty_cell = True
                    break
        if not empty_cell:
            self.chemistry_update_press_rxn_param()
            self.chemistry_update_pressure()

            #tw.insertRow(tw.rowCount())
            #header_height = tw.horizontalHeader().height()

            # Note - scrollbar status can change outside of this function.
            # Do we need to call this every time window geometry changes?
            #scrollbar_height = tw.horizontalScrollBar().isVisible() * (4+tw.horizontalScrollBar().height())
            #nrows = tw.rowCount()
            #row_height = tw.rowHeight(0)
            #tw.setMinimumHeight(header_height+scrollbar_height
            #           + nrows*row_height)
        self.set_reaction_edited()
        self.chemistry_update_buttons()

    def chemistry_update_press_rxn_param(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not reaction:
            return
        if not ui.checkbox_pressure_dependent.isChecked():
            reaction.pop('press_rxn_param', None)
            return
        method_idx = ui.combobox_pressure_method.currentIndex()
        method = ['Lindemann', 'Troe', 'SRI'][method_idx]
        model_idx = ui.combobox_pressure_model.currentIndex()
        params = []
        if model_idx == UNIMOL:
            model = method + "_falloff"
        elif model_idx == BIMOL:
            model = method + "_bimo"
        elif model_idx == PLOG:
            model = "PLOG"
        else:
            raise ValueError(model_idx)

        if model_idx in (UNIMOL, BIMOL):
            (A,B,Ea) = (None if le.text()=='' else FloatExp(le.value)
                            for le in (ui.lineedit_pressure_arrhenius_A,
                                       ui.lineedit_pressure_arrhenius_B,
                                       ui.lineedit_pressure_arrhenius_Ea))
            params.append('arrhenius_press: %s %s %s' % (A,B,Ea))

            if method_idx == TROE:
                params.append("press_coeff:")
                for s in 'alpha', 't3star', 'tstar', 't2star':
                    w = getattr(ui, 'lineedit_troe_'+s)
                    params.append(str(FloatExp(w.value)) if w.text() else ('0' if s=='t2star' else 'None'))
            elif method_idx == SRI:
                params.append("press_coeff:")
                for s in 'abcde':
                    w = getattr(ui, 'lineedit_sri_'+s)
                    params.append(str(FloatExp(w.value)) if w.text() else
                                  ('1' if s=='d' else '0' if s=='e' else 'None'))

        else: #PLOG
            tw = ui.tablewidget_plog
            nrows = tw.rowCount()
            def get_plog_val(i,j):
                item = tw.item(i,j)
                if not item:
                    return None
                text = item.text()
                if not text:
                    return None
                return FloatExp(text)

            press_coeff = [(get_plog_val(i,1),i) for i in range(nrows) if get_plog_val(i,1) is not None]
            press_coeff.sort()
            params.append("press_coeff:")
            for (v,i) in press_coeff:
                if v is not None:
                    params.append(str(v))
            params.append("arrhenius_press:")
            for (v,i) in press_coeff:
                if v is None:
                    continue
                a = get_plog_val(i,2)
                b = get_plog_val(i,3)
                c = get_plog_val(i,4)
                params.append('0' if a is None else str(a))
                params.append('0' if b is None else str(b))
                params.append('0,' if c is None else str(c)+',')
            if params[-1].endswith(','): # Trim trailing comma
                params[-1] = params[-1][:-1]
        reaction['press_rxn_param'] = '%s "%s"' % (model,
                                                   ' '.join(params))


    def parse_press_rxn_param(self, tok):
        d = {}
        key = None
        for t in tok:
            t = dequote(t)
            if t.endswith(':'):
                key = t[:-1]
                d[key] = []
            else:
                if key:
                    if t.endswith(','):
                        t = t[:-1]
                    if t == 'None':
                        val = None
                    else:
                        try:
                            val = float(t)
                        except ValueError:
                            self.error('Invalid numeric value "%s" in press_rxn_param "%s"' %
                                       (t, ' '.join(tok)),
                                       popup=True)
                            val = None
                    d[key].append(val)
                else:
                    self.error('Invalid press_rxn_param: "%s"' % ' '.join(tok),
                               popup=True)

        return d

### Landau-Teller
    def chemistry_update_landau_teller(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        cb = ui.checkbox_landau_teller
        gb = ui.groupbox_landau_teller
        if not self.chemkin_enabled():
            cb.setVisible(False)
            gb.setVisible(False)
            return
        if reaction_is_des(self.project, reaction):
            reaction.pop('lt_coeff', None)
            cb.setVisible(False)
            gb.setVisible(False)
            return
        cb.setVisible(True)
        lt_coeff = reaction.get('lt_coeff')
        cb.setChecked(bool(lt_coeff))
        gb.setVisible(bool(lt_coeff))
        if not lt_coeff:
            return
        try:
            b,c = [None if x=='None' else FloatExp(x) for x in lt_coeff.split()]
        except Exception as e:
            self.error('Invalid value for lt_coeff "%s"\n%s' %
                       (lt_coeff, str(e)), popup=True)
            reaction['lt_coeff'] = '0 0'
            b, c = 0, 0
            self.set_reaction_edited()
            self.chemistry_update_buttons()
        ui.lineedit_landau_teller_B.updateValue(None, b)
        ui.lineedit_landau_teller_C.updateValue(None, c)


    def chemistry_handle_checkbox_landau_teller(self, checked):
        ui = self.ui.chemistry
        ui.groupbox_landau_teller.setVisible(checked)
        reaction = self.working_reaction
        if not reaction:
            return
        ui.lineedit_landau_teller_B.updateValue(None, None)
        ui.lineedit_landau_teller_C.updateValue(None, None)

        if not checked:
            reaction.pop('lt_coeff', None)
        else:
            reaction['lt_coeff'] = 'None None'

        self.set_reaction_edited()
        self.chemistry_update_buttons()

    def chemistry_handle_lineedit_landau_teller(self, widget, update_dict, args):
        reaction = self.working_reaction
        if not reaction:
            return
        ui = self.ui.chemistry
        key, val = update_dict.popitem()
        arr = reaction.get('lt_coeff')
        if not arr:
            c = [None, None]
        else:
            c = [None if x=="None" else FloatExp(x) for x in arr.split()]
        if key == 'B':
            c[0] = val
        elif key == 'C':
            c[1] = val
        reaction['lt_coeff'] = '%s %s' % tuple('None' if x in ('None', None, '')
                                               else FloatExp(x) for x in c)
        self.set_reaction_edited()
        self.chemistry_update_buttons()


### Rate fitting
    def chemistry_update_rate_fit(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        cb = ui.checkbox_rate_fit
        gb = ui.groupbox_rate_fit
        if not self.chemkin_enabled():
            cb.setVisible(False)
            gb.setVisible(False)
            return
        if reaction_is_des(self.project, reaction):
            reaction.pop('fit_coeff', None)
            cb.setVisible(False)
            gb.setVisible(False)
            return

        cb.setVisible(True)
        fit_coeff = reaction.get('fit_coeff','')
        cb.setChecked(bool(fit_coeff))
        gb.setVisible(bool(fit_coeff))
        if not fit_coeff:
            return
        method = fit_coeff.split()[0].lower()
        if method == 'jan':
            expected = 9
        elif method == 'fit1':
            expected = 4
        else:
            self.error('Invalid method for fit_coeff: "%s"' % method.upper(),
                       popup=True)
            ui.combobox_rate_fit_method.currentIndexChanged.disconnect()
            ui.combobox_rate_fit_method.setCurrentIndex(JAN)
            ui.combobox_rate_fit_method.currentIndexChanged.connect(
                self.chemistry_handle_combobox_rate_fit_method)
            self.chemistry_handle_combobox_rate_fit_method(JAN, in_setup=True)
            method = 'jan'
            fit_coeff = 'JAN "coeff: ' + ' '.join(9*['0']) + '"'
            self.working_reaction['fit_coeff'] = fit_coeff
            self.set_reaction_edited()

        coeffs = dequote(fit_coeff[4:]).split()
        if coeffs[0] == 'coeff:':
            coeffs.pop(0)
        else:
            self.error('Expected "coeff:"\n"%s"' % fit_coeff,
                       popup=True)
            fit_coeff = method.upper() + ' "coeff: ' + ' '.join(expected*['0']) + '"'
            self.working_reaction['fit_coeff'] = fit_coeff
            self.set_reaction_edited()
            coeffs = expected*['0']

        if len(coeffs) != expected:
            self.error('%d parameters required for %s fitting, got %s\n"%s"',
                       (expected, method.upper(), len(coeffs), fit_coeff),
                       popup=True)
            coeffs = [0]*expected
            self.working_reaction['fit_coeff'] = method.upper() + ' "coeff: ' + ' '.join(['None'*expected]) + '"'
            self.set_reaction_edited()
        try:
            coeffs = list(map(FloatExp, coeffs))
        except Exception as e:
            self.error('Invalid value for fit_coeff: "%s"' % str(e),
                       popup=True)
            coeffs = [0]*expected
            self.working_reaction['fit_coeff'] = method.upper() + ' "coeff: ' + ' '.join(['None'*expected]) + '"'
            self.set_reaction_edited()

        idx = JAN if method=='jan' else FIT1
        ui.combobox_rate_fit_method.currentIndexChanged.disconnect()
        ui.combobox_rate_fit_method.setCurrentIndex(idx)
        ui.combobox_rate_fit_method.currentIndexChanged.connect(
            self.chemistry_handle_combobox_rate_fit_method)
        self.chemistry_handle_combobox_rate_fit_method(idx, in_setup=True)
        for i in range(expected):
            getattr(ui, 'lineedit_fit_b%s'%(i+1)).updateValue(None, coeffs[i])
        if self.reaction_edited:
            self.chemistry_update_buttons()


    def chemistry_handle_checkbox_rate_fit(self, checked):
        ui = self.ui.chemistry
        ui.groupbox_rate_fit.setVisible(checked)
        if checked:
            ui.combobox_rate_fit_method.setCurrentIndex(0)
            self.chemistry_handle_combobox_rate_fit_method(0)
        else:
            self.working_reaction.pop('fit_coeff', None)
            self.set_reaction_edited()
            self.chemistry_update_buttons()

    def chemistry_handle_combobox_rate_fit_method(self, idx, in_setup=False):
        ui = self.ui.chemistry
        set_combobox_tooltip(ui.combobox_rate_fit_method)
        for x in range(1,10):
            getattr(ui, 'lineedit_fit_b%s'%x).updateValue(None, None)
        for x in range(5, 10):
            getattr(ui, 'lineedit_fit_b%s'%x).setVisible(idx==0) # JAN
            getattr(ui, 'label_fit_b%s'%x).setVisible(idx==0)
        if not in_setup:
            self.chemistry_update_fit_coeff()

    def chemistry_handle_lineedit_rate_fit(self, *args):
        reaction = self.working_reaction
        if not reaction:
            return
        self.chemistry_update_fit_coeff()

    def chemistry_update_fit_coeff(self):
        ui = self.ui.chemistry
        reaction = self.working_reaction
        if not reaction:
            return
        fit_coeff = reaction.get('fit_coeff')
        idx = ui.combobox_rate_fit_method.currentIndex()
        if idx==JAN:
            method = 'JAN'
            expected = 9
        else:
            method = 'FIT1'
            expected = 4
        vals = []
        for i in range(expected):
            le = getattr(ui, 'lineedit_fit_b%s'%(i+1))
            vals.append('None' if le.value in (None, '') else str(le.value))

        self.working_reaction['fit_coeff'] = method.upper() + ' "coeff: ' + ' '.join(vals) + '"'
        self.set_reaction_edited()
        self.chemistry_update_buttons()
