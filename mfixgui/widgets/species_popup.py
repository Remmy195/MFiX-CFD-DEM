#!/usr/bin/env python
"""Species selector dialog for MFIX GUI, includes stand-alone test"""

# 2016-11-20  Species/alias unification
#  we will only expose 'alias' to the user.  'species' is only used
#  as a key into Burcat/THERMO_DATA, and we're going to inline all
#  of the thermodynamic data - cgw

# 2021-07-09 TODO use the Lineedit widget from mfix to do input
# validation


BURCAT, CHEMKIN, CANTERA = 1,2,3
import os
import sys
import signal
import pickle
import warnings
from copy import deepcopy

from qtpy.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog,
                            QFileDialog,  QHeaderView,
                            QLabel, QLineEdit, QPushButton,
                            QSizePolicy, QSpacerItem, QTableWidgetItem,
                            QVBoxLayout, QWidget)

from qtpy.QtGui import QDoubleValidator, QValidator
from qtpy.QtCore import Signal, QSize

from mfixgui.tools import safe_float
from mfixgui.tools.qt import (get_combobox_item, get_icon, get_selected_row,
                              get_selected_rows, set_item_noedit, set_item_enabled,
                              get_ui)
from mfixgui.widgets.plot_window import PlotWindow
from mfixgui.tools.read_thermo import (parse_section_burcat, parse_section_chemkin,
                                       expand_tabs)


def resize_column(tw, col, flags):
    tw.horizontalHeader().setSectionResizeMode(col, flags)

# Make search case-, whitespace- and punctuation-insensitive
squash_table = [c.lower() if c.isalnum() else None for c in map(chr,range(256))]
def squash(string):
    r = string.translate(squash_table)
    return r

def composition_str(c):
    return ' '.join('%s %d'%(k.title(),c[k.title()]) for k in sorted(c))

def parse_composition_str(s):
    d = {}
    tok = s.strip().split()
    if not s.strip():
        raise ValueError
    while len(tok) >= 2:
        d[tok[0].title()] = float(tok[1])
        tok = tok[2:]
    if tok:
        raise ValueError
    return d

phase_names = 'Gas', 'Liquid', 'Solid', 'Composite'
gui = None

class SpeciesPopup(QDialog):
    save = Signal()
    cancel = Signal()

    def load_burcat_pickle(self, path):
        if not os.path.exists(path):
            print("%s not found, create it by running read_thermo.py" % path)
            sys.exit(-1)
        with open(path, 'rb') as db_file:
            database = pickle.load(db_file)
        by_phase = {}

        for k, v in database.items():
            # key:  name, phase, tmin, tmax
            # value:  coeffs, mol_weight, h_f, tcom, comment
            name, phase, tmin, tmax = k
            coeffs, mol_weight, h_f, tcom, comment = v
            if phase not in by_phase:
                by_phase[phase] = {}

            key = (name, tmin, tmax)
            temps = [tmin, tcom, tmax]
            composition = None
            by_phase[phase][key] = (temps, coeffs, composition, mol_weight, h_f, comment, "BURCAT")
        self.db = by_phase
        self.db_burcat = self.db

        # build search list, lowercased
        self.haystack = []
        self.comments = {}
        for phase in 'GLSC':
            h_tmp = [((squash(k[0]), squash(v[5])), k, phase) for (k,v) in self.db[phase].items()]
            h_tmp.sort()
            self.haystack.extend(h_tmp)
            # comment fields
            self.comments[phase] = dict((k, v[5])
                                        for (k,v) in self.db[phase].items())
        self.haystack_burcat = self.haystack
        self.comments_burcat = self.comments

    def load_atomic_pickle(self, path):
        if not os.path.exists(path):
            print("%s not found, create it by running read_atomic.py" % path)
            sys.exit(-1)
        with open(path, 'rb') as db_file:
            self.atomic = pickle.load(db_file)

    def do_search(self, string=None):
        lineedit = self.ui.lineedit_search
        if string is None:
            string = lineedit.text()
        string = string.lstrip() # Don't allow leading spaces
        if string != lineedit.text():
            lineedit.setText(string)
            return

        self.ui.tablewidget_search.clearContents()
        results = []
        match_empty = True # Show all possibilities when no search string
        if match_empty or string:
            needle = squash(string)
            for (k, key, phase) in self.haystack:
                if (phase in self.phases and
                    (needle in k[0] or
                     (self.include_comments and needle in k[1]))):
                    results.append((key, phase))
        # Put exact matches & leading-substring matches first
        if string:
            results.sort(
                key=lambda x:
                (1 - x[0][0].lower().startswith(string.lower()), x))

        tw = self.ui.tablewidget_search
        nrows = len(results)
        self.ui.tablewidget_search.clearContents()
        self.ui.tablewidget_search.setRowCount(nrows)
        self.search_results = [None]*nrows

        # http://stackoverflow.com/questions/10192579/
        tw.model().blockSignals(True)
        for (i, r) in enumerate(results):
            key, phase = r
            comment = self.comments[phase][key]
            item = QTableWidgetItem(key[0])
            item.setToolTip(comment)
            set_item_noedit(item)
            tw.setItem(i, 0, item)
            item = QTableWidgetItem(phase)
            set_item_noedit(item)
            tw.setItem(i, 1, item)
            self.search_results[i] = (key, phase)
        tw.model().blockSignals(False)
        if nrows == 1: # Autoselect unique row
            tw.setCurrentCell(0,0)

    def get_species_data(self, key, phase):
        """exposes species database to external clients"""
        db = self.db.get(phase)
        if not db:
            return None
        # FIXME, this is inefficient.  remove tmin/tmax from key tuple.
        #  also, it's possible that there are multiple definitions for the
        #  same species, with different temp. ranges.  This just returns
        #  the first one
        for (keytuple, data) in db.items():
            (species, tmin, tmax) = keytuple
            if species == key:
                (coeffs, mol_weight, tcom, comment) = data
                h_f = coeffs[14]
                coeffs = coeffs[7:14] + coeffs[:7]
                return {'phase': phase,
                        'temps': [tmin, tcom, tmax],
                        'coeffs': coeffs,
                        'composition': None,
                        'mol_weight': mol_weight,
                        'h_f': h_f,
                        format: 'BURCAT'}


    def handle_search_selection(self):
        row = get_selected_row(self.ui.tablewidget_search)
        self.ui.pushbutton_import.setEnabled(row is not None)


    def handle_include_comments(self, val):
        self.include_comments = val
        self.do_search()


    def clear_species_panel(self):
        for item in self.species_panel_items:
            item.setEnabled(False)
            if hasattr(item, 'setText'):
                item.setText('')
        hdr = self.ui.hlayout_table_header
        for i in range(hdr.count()):
            item = hdr.itemAt(i)
            widget = item.widget()
            if isinstance(widget, QPushButton):
                widget.setEnabled(False)

        tw = self.ui.tablewidget_params
        for row in range(8):
            for col in range(2):
                w = tw.cellWidget(row, col)
                if w and hasattr(w,'setText'):
                    w.setText('')
                    #tw.cellWidget(i,j).setText('')

    def make_item(self, val=None, key=None):
        item = QLineEdit()
        item.setText('' if val is None else str(val))
        item.saved_value = val
        item.setValidator(QDoubleValidator(item))
        item.setFrame(False)
        if key:
            item.editingFinished.connect(self.make_handler(item=item, key=key))
            #inputRejected is only available in Qt 5.12 and above
            #item.inputRejected.connect(make_reject_handler(item=item, key=key))
            self.make_reject_handler(item=item, key=key)
            self.make_foe(item, key)
        return item

    def make_foe(self, item, key):
        item._focusOutEvent = item.focusOutEvent
        def foe(ev, item=item):
            key = item.key
            if key == 'max_temp':
                key = ('temps', -1)
            txt = item.text()
            try:
                if key == 'h_f' and not txt: #1662 allow blank Hf/R
                    data = self.defined_species[self.current_species]
                    data[key] = None
                    self.check_data()
                    item.saved_value = None
                else:
                    if key=='composition':
                        val = parse_composition_str(txt)
                    else:
                        val = float(txt)
                    item.setText('' if val is None else composition_str(val) if key=='composition'
                                 else str(val))
                    item.saved_value = val
            except ValueError:
                val = getattr(item, 'saved_value', None)
                item.setText('' if val is None else composition_str(val) if key=='composition'
                             else str(val))
                item.reject_handler(item)
                self.check_data()
            return item._focusOutEvent(ev)
        item.key = key
        item.focusOutEvent = foe


    def make_handler(self, item, key):
        def handler(item=item):
            ui = self.ui
            if not self.current_species:
                return
            key = item.key
            if key == 'max_temp':
                key = ('temps', -1)
            val = item.text()
            try:
                data = self.defined_species[self.current_species]
                if key == 'h_f' and not val: # 1662
                    val = None
                elif key== 'composition':
                    val = parse_composition_str(val)
                else:
                    val = float(val)
                if isinstance(key, tuple):
                    data[key[0]][key[1]] = val
                else:
                    data[key] = val

            except ValueError:
                # should not get here, field has been validated
                #print("VE")
                pass
            item.saved_value = val
            format = data['format']
            if format=='CHEMKIN':
                if isinstance(key, tuple) and key[0]=='coeffs' and key[1]<7:
                    self.recompute_h_f()
                if key=='composition':
                    try:
                        composition = data['composition']
                        mol_weight = sum(v*self.atomic[k][2] for (k,v) in composition.items())
                        mol_weight = round(mol_weight, 8)
                        data['mol_weight'] = mol_weight
                        ui.lineedit_mol_weight.setText(str(mol_weight))
                        ui.lineedit_composition.setText(composition_str(composition))
                    except:
                        data['mol_weight'] = None
                        ui.lineedit_mol_weight.setText('')
                        ui.lineedit_mol_weight.saved_value = None
            self.check_data()
        item.key = key
        return handler

    def make_reject_handler(self, item, key):
        def handler(item=item):
            if not self.current_species:
                return
            key = item.key
            if key == 'max_temp':
                key = ('temps', -1)
            data = self.defined_species[self.current_species]
            val = getattr(item, 'saved_value', 0.0)
            if isinstance(key, tuple):
                data[key[0]][key[1]] = val
            elif key == 'composition':
                if val:
                    data[key] = val
            else:
                data[key] = val

        item.key = key
        item.reject_handler = handler
        return handler


    def enable_species_panel(self):
        ui = self.ui
        self.in_setup = True
        for item in self.species_panel_items:
            item.setEnabled(True)
        hdr = self.ui.hlayout_table_header
        buttons = []
        for i in range(hdr.count()):
            item = hdr.itemAt(i)
            widget = item.widget()
            if isinstance(widget, QPushButton):
                buttons.append(widget)
        for b in buttons:
            # Don't allow delete if only 2 columns
            b.setEnabled(len(buttons) > 3)
        buttons[-1].setEnabled(True)

        species = self.current_species
        data = self.defined_species.get(species)
        ui.lineedit_alias.setText(data['alias'])
        ui.combobox_phase.setEnabled(False)

        i = 'GLSC'.index(data['phase'])
        ui.combobox_phase.setCurrentIndex(i)
        ui.combobox_phase.setToolTip(phase_names[i])
        format = data['format']
        temps = data['temps']
        ui.lineedit_mol_weight.setReadOnly(format=='CHEMKIN')
        ui.lineedit_h_f.setReadOnly(format=='CHEMKIN')
        if format=='CHEMKIN':
            ui.lineedit_mol_weight.setToolTip("Molecular weight is computed from composition.")
            ui.lineedit_h_f.setToolTip("Heat of formation is computed from coeficcinents.")
        else:
            ui.lineedit_mol_weight.setToolTip(None)
            ui.lineedit_h_f.setToolTip(None)

        composition = data.get('composition')
        if composition:
            ui.lineedit_composition.setText(composition_str(composition))
            ui.lineedit_composition.saved_value = composition
        else:
            ui.lineedit_composition.setText(None)
        ui.lineedit_mol_weight.setText(str(data['mol_weight']))
        ui.lineedit_mol_weight.saved_value = data['mol_weight']
        ui.lineedit_ljsig.saved_value = data.get('ljsig')
        ui.lineedit_ljeps.saved_value = data.get('ljeps')

        ui.lineedit_h_f.setText('' if data['h_f'] is None else str(data['h_f']))
        if self.density_enabled:
            density = data.get('density')
            ui.lineedit_density.setText('' if density is None else str(density))
            handler = self.make_handler(ui.lineedit_density, 'density')
            ui.lineedit_density.editingFinished.connect(handler)
        tw = ui.tablewidget_params
        temps = data['temps']
        cb = ui.combobox_data_format
        if len(temps) > 3:
            cb.setCurrentIndex(1)
            get_combobox_item(cb,0).setEnabled(False)
            get_combobox_item(cb,0).setToolTip("CHEMKIN format required for more than 2 temperature ranges")
        else:
            get_combobox_item(cb,0).setEnabled(True)
            get_combobox_item(cb,0).setToolTip(None)
        n_ranges = len(temps) - 1
        layout = ui.hlayout_table_header
        while tw.columnCount() < n_ranges:
            self.add_column()
        while tw.columnCount() > n_ranges:
            pb = layout.itemAt(layout.count()-2).widget()
            self.remove_column(pb)

        for i in range(n_ranges):
            tw.setCellWidget(0, i, self.make_item(temps[i], key=('temps',i)))

        ui.lineedit_max_temp.setText(str(temps[-1]))

        for i in range(n_ranges):
            for (j,x) in enumerate(data['coeffs'][i*7:i*7+7]):
                tw.setCellWidget(j+1, i, self.make_item(x, key=('coeffs', i*7+j)))
        item = get_combobox_item(ui.combobox_data_format, 0) # Burcat
        if format=='BURCAT' or len(temps) < 4:
            set_item_enabled(item, True)
            item.setToolTip(None)
        else:
            set_item_enabled(item, False)
            item.setToolTip("Burcat format only supports 2 temperature ranges")
        pb = ui.pushbutton_add_column
        if format=='BURCAT':
            pb.setEnabled(False)
            pb.setToolTip("Burcat format only supports 2 temperature ranges")
        else:
            pb.setEnabled(True)
            pb.setToolTip(None)
        ui.lineedit_composition.setEnabled(data['format']=='CHEMKIN')
        ui.combobox_data_format.setCurrentIndex(0 if data['format']=='BURCAT'
                                                else 1)
        for k in ('config_transport', 'alpha_transport', 'mu_transport', 'zrot_transport',
                  'ljsig', 'ljeps'):
            val = data.get(k)
            if k=='config_transport':
                ui.combobox_config_transport.setCurrentIndex(0 if val is None else val+1)
            else:
                if val is None and k.endswith('_transport'):
                    val = 0
                getattr(ui, 'lineedit_'+k).setText(str(val) if val is not None else '')
        self.in_setup = False


    def add_column(self):
        ui = self.ui
        tw = ui.tablewidget_params
        layout = ui.hlayout_table_header
        ncols = tw.columnCount()
        pb = QPushButton("")
        pb.setIcon(get_icon("remove"))
        pb.setToolTip("Remove column %s"%ncols)
        pb.setAutoDefault(False)
        pb.column = ncols
        layout.insertWidget(ncols, pb)
        pb.clicked.connect(lambda arg, pb=pb: self.remove_column(pb))
        tw = ui.tablewidget_params
        tw.setColumnCount(ncols + 1)
        resize_column(tw, ncols, QHeaderView.Stretch)
        if self.current_species:
            data = self.defined_species[self.current_species]
            coeffs = data['coeffs']
            temps = data['temps']
            if not self.in_setup:
                if len(coeffs) != 7*(ncols+1):
                    coeffs.extend([None]*7)
                if len(temps) != ncols + 2:
                    temps.insert(-1,None) # Keep max temp at end

            tw.setCellWidget(0, ncols, self.make_item(temps[ncols], key=('temps',ncols)))
            for i in range(7):
                tw.setCellWidget(i+1, ncols, self.make_item(coeffs[ncols*7+i], key=('coeffs',ncols*7+i)))
        self.setup_column_buttons()
        self.resize_table_header()
        if not self.in_setup:
            self.check_data()

    def remove_column(self, pb):
        ui = self.ui
        tw = ui.tablewidget_params
        col = pb.column
        tw.removeColumn(col-1)
        ui.hlayout_table_header.removeWidget(pb)
        pb.hide()
        pb.deleteLater()
        self.setup_column_buttons()
        self.resize_table_header()
        if self.current_species and not self.in_setup:
            data = self.defined_species[self.current_species]
            del data['coeffs'][(col-1)*7: col*7]
            del data['temps'][(col-1)]
            self.recompute_h_f()

        for j in range(col-1, tw.columnCount()):
            for i in range(tw.rowCount()):
                item = tw.cellWidget(i,j)
                if item:
                    key = item.key
                    if isinstance(key, tuple):
                        key, idx = key
                        #print(j,i,item.key, end=" ")
                        if key=='coeffs':
                            idx -= 7
                        else:
                            idx -= 1
                        item.key = (key, idx)
                        #print("-> ", item.key)
        if not self.in_setup:
            self.check_data()

    def setup_column_buttons(self):
        hdr = self.ui.hlayout_table_header
        buttons = []
        for i in range(hdr.count()):
            item = hdr.itemAt(i)
            widget = item.widget()
            if isinstance(widget, QPushButton):
                buttons.append(widget)
        for i,b in enumerate(buttons[:-1],1):
            b.setToolTip("Remove column %s"%(i))
            b.column = i
            b.setEnabled(len(buttons)>3)
        buttons[-1].setEnabled(True)
        hdr.update()


    def check_composition_str(self, s):
        if not self.atomic:
            return False, "Error: atomic data not found"
        if not s.strip():
            return False, "Composition may not be empty"
        tok = s.strip().split()
        while len(tok)>=2:
            if tok[0].title() not in self.atomic:
                return False, "Invalid element %s in composition" % tok[0]
            try:
                val = float(tok[1])
            except ValueError:
                return False, "Invalid number %s in composition" % tok[0]
            if val <= 0:
                return False, "Atom counts must be > 0 in composition"
            tok = tok[2:]
        if tok:
            return False, "Odd number of terms in composition"
        return True, "Composition OK"


    def check_data(self):
        ui = self.ui
        if not self.current_species:
            self.set_done_button(False)
            return
        data = self.defined_species[self.current_species]
        self.data_ok = True
        mw = data.get('mol_weight', 0)
        temps = data.get('temps')
        h_f = data.get('h_f')
        # Avoid zero-divide in plotting
        ui.pushbutton_plot.setEnabled(bool(mw and mw >0))

        # Check user inputs
        coeffs = data['coeffs']
        if data['format'] == 'CHEMKIN':
            # data checks for Chemkin format
            ok, msg = self.check_composition_str(ui.lineedit_composition.text())
            if not ok:
                self.data_ok = False
                self.set_done_button(False, msg)
            elif None in temps:
                self.data_ok = False
                self.set_done_button(False, "All temperatures must be defined")
            elif temps != sorted(temps):
                self.data_ok = False
                self.set_done_button(False, "Temperatures must be ascending")
            elif None in coeffs:
                self.data_ok = False
                self.set_done_button(False, "All coefficients must be defined")
            else:
                for x in range(len(coeffs)//7):
                    if all(coeffs[i]==0 for i in range(x*7, x*7+5)):
                        self.data_ok = False
                        self.set_done_button(False, "All coefficients may not be 0 within a range")

        if data['format'] == 'BURCAT':
            tmin, tcom, tmax = temps
            if mw <= 0:
                self.data_ok = False
                self.set_done_button(False, "Molecular weight must be > 0")
            elif tmin is None:
                self.data_ok = False
                self.set_done_button(False, "Low-temp range limit must be defined")
            elif tmin < 0:
                self.data_ok = False
                self.set_done_button(False, "Low-temp range limit must be > 0K")
            elif tcom is None:
                self.set_done_button(False, "High temp range start must be defined")
            elif tmin >= tcom:
                self.data_ok = False
                self.set_done_button(False, "Low-temp range limit must be < %sK" % tcom)
                return
            elif tmax is None:
                self.data_ok = False
                self.set_done_button(False, "Maximum temperature must be defined")
                return
            elif tmax <= tmin:
                self.data_ok = False
                self.set_done_button(False, "Maximum temperature must be > %s"%tmin)
            #elif h_f is None:
            #    self.data_ok = False
            #    self.set_done_button(False, "Heat of formation must be defined")
            elif (tmin < tcom and all(x==0 for x in coeffs[:5])):
                self.data_ok = False
                self.set_done_button(False,
                    "Low-temp coefficients can not all be 0 when t_low < tcom")
            elif (tmax > tcom and all(x==0 for x in coeffs[7:12])):
                self.data_ok = False
                self.set_done_button(False,
                    "High-temp coefficients can not all be 0 when t_high > tcom")

        if self.cantera_enabled:
            for (name, data) in self.defined_species.items():
                for (name, data) in self.defined_species.items():
                    if data.get('config_transport') not in (0,1,2):
                        self.data_ok = False
                        self.set_done_button(False,
                                             "Molecule shape must be defined for "+name)
                        break
                    elif 'ljeps' not in data or 'ljsig' not in data:
                        self.data_ok = False
                        self.set_done_button(False,
                                             "Lennard-Jones parameters must be defined for " + name)
                        break


        elif self.lennard_jones_enabled:
            for (name, data) in self.defined_species.items():
                if 'ljsig' not in data or 'ljeps' not in data:
                    self.data_ok = False
                    self.set_done_button(False,
                                         "Lennard-Jones parameters must be defined for " + name)
                    break

        if self.data_ok and self.alias_ok:
            self.set_done_button(True)


    def enable_density(self, enable):
        self.density_enabled = enable
        ui = self.ui
        if not enable:
            if ui.lineedit_density in self.species_panel_items:
                self.species_panel_items.remove(ui.lineedit_density)
                ui.lineedit_density.setEnabled(False)
            ui.lineedit_density.clear()
        else:
            ui.lineedit_density.setEnabled(True)
            self.species_panel_items.append(ui.lineedit_density)

    def enable_lennard_jones(self, enable):
        self.lennard_jones_enabled = enable
        self.ui.groupbox_lennard_jones.setVisible(enable)
        self.setup_combobox_plot_type()


    def enable_cantera(self, enable):
        self.cantera_enabled = enable
        self.ui.groupbox_cantera.setVisible(enable)
        self.setup_combobox_plot_type()

    def setup_combobox_plot_type(self):
        ui = self.ui
        cb = ui.combobox_plot_type
        for (idx, key) in ((1, 'kg_model'),
                           (2, 'mu_g_model')):
            item = get_combobox_item(cb, idx)
            enable = gui.project.get_value(key) in ("CANTERA_POLY", "LENNARD_JONES")
            set_item_enabled(item, enable)
            item.setToolTip("Requires Lennard-Jones or Cantera polynomial model"
                            if not enable else None)

        if not get_combobox_item(cb, cb.currentIndex()).isEnabled():
            cb.setCurrentIndex(0)

    def handle_defined_species_selection(self):
        #self.ui.tablewidget_search.clearSelection()
        tw = self.ui.tablewidget_defined_species
        row = get_selected_row(tw)

        if row is None:
            self.current_species = None
            self.clear_species_panel()
            self.ui.pushbutton_copy.setEnabled(False)
            self.ui.combobox_phase.setEnabled(False)
            self.ui.lineedit_density.setEnabled(False)
            self.ui.combobox_config_transport.setEnabled(False)
        else:
            self.ui.pushbutton_copy.setEnabled(True)
            self.current_species = tw.item(row, 0).text()
            self.enable_species_panel()
            if self.density_enabled:
                self.ui.lineedit_density.setEnabled(True)
            self.ui.combobox_config_transport.setEnabled(True)

    def make_alias(self, name):
        #Aliases must be unique.
        #Aliases are limited to 32 characters and must follow FORTRAN variable
        #naming conventions (i.e., alphanumeric combinations with a letter as the first
        #character).
        #Aliases are not case sensitive.
        #Aliases cannot conflict with existing MFiX variable names (e.g., a species
        # alias of MU_g will cause an error when compiling MFiX).

        for (pat, repl) in [(' ', '_'),
                            ('(', '_'),
                            (')', ''),
                            ('-', '_'),
                            ('__', '_')]:
            name = name.replace(pat, repl)

        alias = ''.join([c for c in name if c.isalnum() or c=='_'])

        while alias and not alias[0].isalpha(): # strip leading _ and digits
            alias = alias[1:]

        if len(alias) > 32:
            alias = alias[:32]

        if alias.lower() not in self.reserved_aliases:
            return alias

        count = 1
        base = alias[:28] # leave room for _nn
        # Strip _nn suffix
        if '_' in alias:
            i = alias.rindex('_')
            if i > 0 and alias[i+1:].isdigit():
                count = 1 + int(alias[i+1:])
                base = alias[:i][:28]

        while alias.lower() in self.reserved_aliases:
            alias = '%s_%s' % (base, count)
            count += 1

        return alias


    def make_user_species_name(self):
        n = 1
        while "Species_%d" % n in self.user_species_names:
            n += 1
        name = "Species_%d" % n
        self.user_species_names.add(name)
        return name


    def do_import(self):
        tw = self.ui.tablewidget_search
        rows = get_selected_rows(tw)
        for row in rows:
            self.do_import_row(row)
            tw.item(row,0).setSelected(False)
            tw.item(row,1).setSelected(False)
            if not(self.data_ok and self.alias_ok):
                tw.scrollToItem(tw.item(row,0), 1)# Scroll to top
                if self.err_msg:
                    err_msg = self.err_msg.replace('<', '&lt;').replace('>', '&gt;')
                    self.parent.error(err_msg, popup=True)
                break

    def do_import_row(self, row):
        self.ui.combobox_phase.setEnabled(False)
        rowdata = self.search_results[row]
        key, phase = rowdata
        data = self.db[phase][key]
        (species, tmin, tmax) = key
        (temps, coeffs, composition, mol_weight, h_f, comment, format) = data

        alias = self.make_alias(species)

        species_data = {
            'phase': phase,
            'alias': alias,
            'species': species,
            'temps': temps,
            'coeffs': coeffs,
            'composition': composition,
            'mol_weight': mol_weight,
            'h_f': h_f,
            'format': format,
        }

        if self.density_enabled:
            species_data['density'] = None # ? where do we get this?
        self.defined_species[alias] = species_data
        self.add_defined_species_row(alias, select=True)
        self.check_data()


    def update_defined_species(self):
        self.ui.tablewidget_defined_species.clearSelection()
        self.ui.tablewidget_defined_species.setRowCount(0)
        for species_key in self.defined_species.keys():
            self.add_defined_species_row(species_key, select=False)


    def add_defined_species_row(self, alias, select=False):
        species_data = self.defined_species[alias]
        ui = self.ui
        tw = ui.tablewidget_defined_species
        nrows = tw.rowCount()
        tw.setRowCount(nrows+1)
        phase = species_data['phase']
        item = QTableWidgetItem(alias)
        set_item_noedit(item)
        tw.setItem(nrows, 0, item)
        item = QTableWidgetItem(phase)
        set_item_noedit(item)
        tw.setItem(nrows, 1, item)
        if select:
            tw.setCurrentCell(nrows, 0) # select new row


    def handle_copy(self):
        tw = self.ui.tablewidget_defined_species
        row = get_selected_row(tw)
        if row is None:
            return
        species = tw.item(row, 0).text()
        alias = self.make_alias(species)
        if species not in self.defined_species:
            return
        species_data = deepcopy(self.defined_species[species])
        species_data['alias'] = alias
        species_data['species'] = species
        species = alias
        self.defined_species[species] = species_data
        self.current_species = species
        self.enable_species_panel()
        self.add_defined_species_row(alias, select=True)
        lineedit = self.ui.lineedit_alias
        lineedit.selectAll()
        lineedit.setFocus(0)
        self.ui.combobox_phase.setEnabled(True)
        self.check_data()

    def handle_new(self):
        phase = self.default_phase
        alias = species = self.make_user_species_name()
        mol_weight = 0
        density = None
        h_f = None
        tmin = 200.0
        tmax = 6000.0
        tcom = 1000.0
        coeffs = [0.0]*14
        species_data = {'phase': phase,
                        'alias': alias,
                        'species': species,
                        'temps': [tmin, tcom, tmax],
                        'coeffs': coeffs,
                        'compostion': None,
                        'mol_weight': mol_weight,
                        'h_f': h_f,
                        'density': density,
                        'format': 'BURCAT'
                        }

        self.defined_species[species] = species_data
        self.current_species = species
        self.enable_species_panel()
        self.add_defined_species_row(alias, select=True)
        lineedit = self.ui.lineedit_alias
        lineedit.selectAll()
        lineedit.setFocus(0)
        self.ui.combobox_phase.setEnabled(True)
        self.check_data()

    def handle_alias(self):
        val = self.ui.lineedit_alias.text() # Already validated (?)
        tw = self.ui.tablewidget_defined_species
        row = get_selected_row(tw)
        if row is None: # No selection
            return
        #note, making a new item here, instead of changing item inplace
        item = QTableWidgetItem(val)
        set_item_noedit(item)
        tw.setItem(row, 0, item)
        defined_species = {}
        for (key, data) in self.defined_species.items():
            if key == self.current_species:
                key = val
                data['alias'] = val
            defined_species[key] = data
        self.current_species = val
        self.defined_species = defined_species

    def set_done_button(self, state, msg=''):
        self.err_msg = msg
        self.ui.pushbutton_done.setEnabled(state)
        self.ui.label_status.setText(msg)
        self.ui.label_status.setStyleSheet("color: red;" if state is False
                                           else "color: blue;")

    def handle_combobox_phase(self, index):
        phase = 'GLSC'[index]
        if not self.current_species:
            return
        species = self.defined_species[self.current_species]
        species['phase'] = phase
        self.ui.combobox_phase.setToolTip(phase_names[index])

    def handle_combobox_config_transport(self, index):
        if not self.current_species:
            return
        key = 'config_transport'
        species = self.defined_species[self.current_species]
        if index == 0:
            species.pop(key, None)
        else:
            if species.get(key) == index-1:
                return
            species[key] = index-1
        if self.in_setup:
            return
        self.check_data()


    def reset_signals(self):
        for sig in (self.cancel, self.save):
            try:
                sig.disconnect()
            except:
                pass

    def handle_phase(self):
        phases = ''
        for phase in 'GLSC':
            button = getattr(self.ui, 'pushbutton_%s' % phase)
            if button.isChecked():
                phases += phase
        if phases == self.phases:
            return
        self.phases = phases
        self.default_phase = phases[0] if phases else ''
        self.do_search()


    def __init__(self, parent=None, phases='GLCS'):
        global gui
        super(SpeciesPopup, self).__init__(parent)
        gui = self.parent = parent
        self.include_comments = True
        self.default_phase = phases[0] if phases else ''
        self.density_enabled = True
        self.err_msg = ''
        self.current_species = None
        self.in_setup = False
        self.lennard_jones_enabled = False
        self.cantera_enabled = False

        datadir = os.path.abspath(os.path.dirname(__file__))
        self.load_burcat_pickle(os.path.join(datadir, 'burcat.pickle'))
        self.load_atomic_pickle(os.path.join(datadir, 'atomic.pickle'))

        self.format = BURCAT

        ui = self.ui = get_ui('species_popup.ui', self)
        # key=alias, val=data tuple.  can add phase to key if needed
        self.defined_species = {}
        self.reserved_aliases = set() # To support enforcing uniqueness
        if parent:
            self.mfix_keywords = set(k.lower()
                                     for k in parent.keyword_doc.keys())
        else:
            self.mfix_keywords = set()

        self.search_results = []
        self.user_species_names = set()
        self.plot_window_cp = None
        self.plot_window_k = None
        self.plot_window_mu = None

        # Set up UI
        ui.combobox_data_source.currentIndexChanged.connect(
            self.set_data_source)
        ui.combobox_data_format.currentIndexChanged.connect(self.handle_data_format)

        ui.lineedit_search.textChanged.connect(self.do_search)
        ui.pushbutton_import.clicked.connect(self.do_import)
        ui.pushbutton_import.setEnabled(False)
        ui.tablewidget_search.itemSelectionChanged.connect(
            self.handle_search_selection)
        ui.tablewidget_defined_species.itemSelectionChanged.connect(
            self.handle_defined_species_selection)

        ui.pushbutton_new.clicked.connect(self.handle_new)
        ui.pushbutton_copy.clicked.connect(self.handle_copy)
        ui.checkbox_include_comments.setChecked(True)
        ui.checkbox_include_comments.clicked.connect(self.handle_include_comments)
        ui.pushbutton_choose_file.setVisible(False)
        ui.pushbutton_choose_file.clicked.connect(self.choose_file)

        self.phases = ''
        if phases is not None:
            self.set_phases(phases)

        for phase in 'GLSC':
            button = getattr(self.ui, 'pushbutton_%s' % phase)
            button.clicked.connect(self.handle_phase)

        cb = ui.combobox_phase
        cb.currentIndexChanged.connect(self.handle_combobox_phase)
        for i,t in enumerate(phase_names):
            get_combobox_item(cb, i).setToolTip(t)

        cb = ui.combobox_config_transport
        cb.currentIndexChanged.connect(self.handle_combobox_config_transport)
        self.parent.add_tooltip(cb, 'config_transport')
        for i in range(3):
            self.parent.add_tooltip(get_combobox_item(cb, 1+i),
                                    'config_transport', value=i)

        #http://stackoverflow.com/questions/15845487/how-do-i-prevent-the-enter-key-from-closing-my-qdialog-qt-4-8-1
        # Do not use buttonbox.  https://mfix.netl.doe.gov/gitlab/develop/mfix/issues/101
        buttons = (ui.pushbutton_done, ui.pushbutton_cancel)
        buttons[0].clicked.connect(lambda: (self.save.emit(), self.close()))
        buttons[1].clicked.connect(lambda: (self.cancel.emit(), self.close()))

        # Fake header for tablewidget_params
        pb = ui.pushbutton_add_column
        pb.setIcon(get_icon("add"))
        #pb.setFlat(True)
        pb.setToolTip("Add new column")
        pb.clicked.connect(lambda arg: self.add_column())

        #ss = """QPushButton { border-radius: 0px; }
        #        QPushButton:pressed { background-color: #d0d0d0; }"""
        #pb.setStyleSheet(ss)

        layout = ui.hlayout_table_header
        for i in range(1,3):
            pb = QPushButton("")
            pb.setAutoDefault(False)
            pb.setIcon(get_icon("remove"))
            pb.setToolTip("Remove column %s"%(i))
            pb.column = i
            pb.clicked.connect(lambda arg, pb=pb: self.remove_column(pb))
            pb.setEnabled(False)
            #pb.setStyleSheet(ss)
            #pb.setFlat(True)
            layout.insertWidget(i,pb)

        tw = ui.tablewidget_params
        tw._resizeEvent = tw.resizeEvent
        tw.resizeEvent = self.resize_table_header

        self.alias_ok = False
        self.data_ok = False
        class AliasValidator(QValidator):
            # Make sure aliases are unique
            def __init__(self, parent=None):
                super(AliasValidator, self).__init__()
                self.parent = parent

            def validate(self, text, pos):
                # 'parent' here is the popup, not the mfix gui
                self.parent.alias_ok = False
                if text == "":
                    self.parent.set_done_button(False)
                    return (QValidator.Intermediate, text, pos)
                if (text[0].isdigit()
                    or text[0] == '_'
                    or not all (c=='_' or c.isalnum() for c in text)):
                    return (QValidator.Invalid, text, pos)
                tlower = text.lower()
                if tlower in self.parent.mfix_keywords:
                    self.parent.set_done_button(False, '%s is an MFiX keyword'%text)
                    return (QValidator.Intermediate, text, pos)
                # This generates extra warnings when we are reviewing
                # species definitions
                #if tlower in self.parent.reserved_aliases:
                #    self.parent.set_done_button(False, 'Alias must be unique')
                #    return (QValidator.Intermediate, text, pos)
                self.parent.alias_ok = True
                #self.parent.check_data()
                return (QValidator.Acceptable, text, pos)

        lineedit = ui.lineedit_alias
        lineedit.setValidator(AliasValidator(parent=self))
        lineedit.editingFinished.connect(self.handle_alias)
        lineedit._focusOutEvent = lineedit.focusOutEvent
        def handle_alias_foe(ev):
            if lineedit.text() == '':
                lineedit.setText(self.current_species or '')
            return lineedit._focusOutEvent(ev)

        lineedit.focusOutEvent = handle_alias_foe
        lineedit.setMaxLength(32)

        for line_edit in (ui.lineedit_mol_weight,
                          ui.lineedit_h_f,
                          ui.lineedit_density,
                          ui.lineedit_max_temp):
            line_edit.setValidator(QDoubleValidator())

        self.species_panel_items = [ui.lineedit_alias,
                                    ui.combobox_data_format,
                                    ui.lineedit_composition,
                                    ui.lineedit_mol_weight,
                                    ui.lineedit_h_f,
                                    ui.lineedit_max_temp,
                                    ui.lineedit_density,
                                    ui.tablewidget_params,
                                    ui.lineedit_alpha_transport,
                                    ui.lineedit_mu_transport,
                                    ui.lineedit_zrot_transport,
                                    ui.lineedit_ljsig,
                                    ui.lineedit_ljeps,
                                    ui.widget_plot_options]

        hv = QHeaderView
        for tw in (self.ui.tablewidget_search,
                   self.ui.tablewidget_defined_species):
            resize_column(tw, 0, hv.Stretch)
            resize_column(tw, 1, hv.ResizeToContents)
        tw = self.ui.tablewidget_params
        for i in (0, 1):
            resize_column(tw, i, hv.Stretch)

        # plot button
        b = ui.pushbutton_plot
        b.clicked.connect(self.do_plot)

        for k in ('mol_weight', 'h_f', 'max_temp', 'composition',
                  'alpha_transport', 'mu_transport', 'zrot_transport',
                  'ljeps', 'ljsig'):
            le = getattr(ui, 'lineedit_'+k)
            handler = self.make_handler(le, k)
            le.editingFinished.connect(handler)
            self.make_reject_handler(le, k)
            self.make_foe(le, k)
            if k.endswith('transport') or k.startswith('lj'):
                self.parent.add_tooltip(le, k)
                self.parent.add_tooltip(getattr(ui, 'label_'+k), k)

        self.set_done_button(False) # nothing to accept
        self.clear_species_panel()

    def set_phases(self, phases):
        if phases == self.phases:
            return
        self.phases = phases
        for phase in 'GLSC':
            button = getattr(self.ui, 'pushbutton_%s' % phase)
            button.setChecked(phase in phases)
        self.default_phase = phases[0] if phases else ''
        self.do_search()

    def popup(self):
        self.raise_()
        self.show()
        self.resize_table_header()
        self.combobox_data_source.setCurrentIndex(0)
        self.combobox_data_format.setCurrentIndex(0)
        self.activateWindow()
        self.check_data()

    def resize_table_header(self, ev=None):
        # Fancy fake header for tablewidget_params
        ui = self.ui
        hl = ui.hlayout_table_header
        tw = ui.tablewidget_params
        vh = tw.verticalHeader()
        spacer = hl.itemAt(0)
        width = vh.width()
        line_width = 1 # lines between columns in table
        spacer.changeSize(width+line_width, 0, QSizePolicy.Fixed, QSizePolicy.Fixed)
        button_add = ui.pushbutton_add_column
        button_add.setFixedWidth(button_add.height())
        # Header has 2 more items than table columns:  spacer at left and
        # add button at right
        ncols = tw.columnCount()
        for i in range(1, ncols+1):
            b = hl.itemAt(i).widget()
            w = tw.columnWidth(i-1) + line_width
            if i == ncols:
                w -= button_add.width()
                if tw.verticalScrollBar().isVisible():
                    w += tw.verticalScrollBar().width() + 1
            b.setFixedWidth(w)
        if ev:
            return tw._resizeEvent(ev)

    def recompute_h_f(self):
        ui = self.ui
        if not self.current_species:
            return
        data = self.defined_species[self.current_species]
        try:
            T = 298.15
            a = data['coeffs']
            h_f = a[5] + sum((a[i-1]*T**i)/i for i in range(1,6))
            data['h_f'] = h_f
            ui.lineedit_h_f.setText(str(h_f))
        except:
            data['h_f'] = None
            ui.lineedit_h_f.setText('')


    def do_plot(self):
        ui = self.ui
        species = self.current_species
        if species is None:
            return
        t_range = sorted([safe_float(ui.lineedit_plt_temp_from.text(), 250),
                   safe_float(ui.lineedit_plt_temp_to.text(), 4000)])
        self.t_range = t_range
        plot_type = self.ui.combobox_plot_type.currentIndex()
        if plot_type == 0:
            if (plot_window := self.plot_window_cp) is None:
                plot_window = self.plot_window_cp = PlotWindow(self, 0)
        elif plot_type == 1:
            if (plot_window := self.plot_window_k) is None:
                plot_window = self.plot_window_k = PlotWindow(self, 1)
        elif plot_type == 2:
            if (plot_window := self.plot_window_mu) is None:
                plot_window = self.plot_window_mu = PlotWindow(self, 2)
        plot_window.do_plot(species)
        plot_window.raise_()


    def set_data_source(self):
        ui = self.ui
        cb = ui.combobox_data_source
        idx = cb.currentIndex()
        ui.pushbutton_choose_file.setVisible(idx > 0)
        if idx==0:  # Burcat
            self.db = self.db_burcat
            self.haystack = self.haystack_burcat
            self.comments = self.comments_burcat
            self.do_search()
        else: # idx==1:
            self.format = idx # 1=BURCAT, 2=CHEMKIN
            self.haystack = []
            self.comments = []
            self.db = {}
            self.do_search()

    def choose_file(self):
        dialog = QFileDialog()
        fname, filter = dialog.getOpenFileName()
        if not fname:
            return
        try:
            self.load_from_file(fname, self.format)
        except Exception as e:
            self.parent.error(str(e), popup=True)
            raise
            return

    def load_from_file(self, fname, format=1):
        if format==BURCAT:
            return self.load_burcat_format(fname)
        elif format==CHEMKIN:
            return self.load_chemkin_format(fname)
        elif format==CANTERA:
            return self.load_cantera_transport(fname)
        else:
            raise ValueError(format)

    def handle_data_format(self, idx):
        if not self.current_species:
            return
        self.defined_species[self.current_species]['format'] = ['BURCAT','CHEMKIN'][idx]
        if not self.in_setup:
            self.enable_species_panel()


    def load_burcat_format(self, fname):
        # see read_thermo.py
        by_phase = {}
        lines = []
        with open(fname, encoding='utf8') as f:
            for line in f.readlines():
                line = line.split('!')[0].rstrip()
                line = expand_tabs(line)
                lines.append(line)
        with warnings.catch_warnings(record=True) as ws:
            for (name, phase, coeffs, temps, composition, mol_weight, h_f, comment) in parse_section_burcat(lines):
                if phase not in by_phase:
                    by_phase[phase] = {}
                tmin, tcom, tmax = temps
                key = (name, tmin, tmax)
                composition = None
                if key in by_phase[phase]:
                    # uniquify names? not now, just skip
                    print("duplicate key %s, skipping" % str(key))
                else:
                    by_phase[phase][key] = (temps, coeffs, composition, mol_weight, h_f, comment, "BURCAT")
        for w in ws:
            self.parent.warn(str(w.message), popup=True)

        self.db = by_phase
        self.haystack = []
        self.comments = {}
        for phase in 'GLSC':
            h_tmp = [((squash(k[0]), squash(v[5])), k, phase) for (k,v) in self.db.get(phase,{}).items()]
            h_tmp.sort()
            self.haystack.extend(h_tmp)
            # comment fields
            self.comments[phase] = dict((k, v[5])
                                        for (k,v) in self.db.get(phase,{}).items())
        self.do_search()


    def load_chemkin_format(self, fname):
        by_phase = {}
        lines = []
        with open(fname, encoding='utf8') as f:
            for line in f.readlines():
                line = line.split('!')[0].rstrip()
                line = expand_tabs(line)
                lines.append(line)
        with warnings.catch_warnings(record=True) as ws:
            for (name, phase, temps, coeffs, composition, mol_weight, h_f, comment) in parse_section_chemkin(
                    lines, atomic=self.atomic):
                if phase not in by_phase:
                    by_phase[phase] = {}
                tmin, tmax = temps[0], temps[-1]
                key = (name, tmin, tmax)
                if key in by_phase[phase]:
                    # uniquify names? not now, just skip
                    print("duplicate key %s, skipping" % str(key))
                else:
                    by_phase[phase][key] = (temps, coeffs, composition, mol_weight, h_f, comment, "CHEMKIN")

        for w in ws:
            self.parent.warn(str(w.message), popup=True)

        self.db = by_phase
        self.haystack = []
        self.comments = {}
        for phase in 'GLSC':
            h_tmp = [((squash(k[0]), squash(v[5])), k, phase) for (k,v) in self.db.get(phase,{}).items()]
            h_tmp.sort()
            self.haystack.extend(h_tmp)
            # comment fields
            self.comments[phase] = dict((k, v[5])
                                        for (k,v) in self.db.get(phase,{}).items())
        self.do_search()

    def load_cantera_transport(self, fname):
        try:
            with open(fname, 'r') as f:
                for line in f:
                    data = {}
                    line = line.strip()
                    if not line:
                        continue
                    line = line.split('!')[0].strip()
                    if (line == ''
                        or line.startswith('!')
                        or not line[0].isalpha()):
                        continue
                    tokens = line.split()
                    if len(tokens) < 7:
                        continue
                    name = tokens[0]
                    data['config_transport'] = int(tokens[1])
                    for k,t in zip(('ljeps', 'ljsig', 'mu_transport', 'alpha_transport', 'zrot_transport'),
                                   tokens[2:8]):
                        data[k] = float(t)

                    # Make names case-insentitive
                    for key in self.defined_species.keys():
                        if key.upper() == name.upper():
                            self.defined_species[key].update(data)
                            break
        except Exception as e:
            self.parent.error("Error loading %s: %s" % (fname, str(e)),
                              popup=True)

        self.handle_defined_species_selection()
        self.check_data()
