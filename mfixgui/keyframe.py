# -*- coding: utf-8 -*-
""" Keyframe pane"""

import csv
import os
import shutil

from json import JSONDecoder, JSONEncoder

from matplotlib import pyplot as plt

from qtpy.QtCore import Qt
UserRole = Qt.UserRole

from qtpy.QtWidgets import (QAbstractItemView, QComboBox,  QFileDialog,
                            QGridLayout, QGroupBox, QHeaderView, QLabel,
                            QMessageBox, QPushButton, QRadioButton,
                            QTableWidgetItem, QVBoxLayout, QWidget)

from qtpy.QtCore import QFileSystemWatcher

from mfixgui.constants import *

from mfixgui.widgets.base import LineEdit

from mfixgui.tools.qt import (get_combobox_item, get_selected_row,
                              set_combobox_tooltip, get_pixmap,
                              set_item_noedit, sub_icon_height)

kf_keys = {
    "Unassigned": {None: None},

    "Fluid velocity": {"x":"bc_u_g",
                       "y":"bc_v_g",
                       "z":"bc_w_g"},

    "Solids velocity": {"x":"bc_u_s",
                        "y":"bc_v_s",
                        "z":"bc_w_s"},

    "Boundary wall velocity" : {"x":"bc_wall_vel",
                                "y":"bc_wall_vel",
                                "z":"bc_wall_vel"},

    "Boundary wall rotation center" : {"x":"bc_wall_rot_center",
                                       "y":"bc_wall_rot_center",
                                       "z":"bc_wall_rot_center"},

    "Boundary wall rotation speed" : {"x": "bc_wall_omega",
                                      "y": "bc_wall_omega",
                                      "z": "bc_wall_omega"},

    "Internal surface velocity" : {"x":"is_wall_vel",
                                   "y":"is_wall_vel",
                                   "z":"is_wall_vel"},

    "Internal surface rotation center" : {"x":"is_wall_rot_center",
                                          "y":"is_wall_rot_center",
                                          "z":"is_wall_rot_center"},

    "Internal surface rotation speed" : {"x": "is_wall_omega",
                                         "y": "is_wall_omega",
                                         "z": "is_wall_omega"},

    "DES rigid motion velocity" : {"x":"des_rigid_motion_vel",
                                   "y":"des_rigid_motion_vel",
                                   "z":"des_rigid_motion_vel"},

    "DES rigid rotation center" : {"x":"des_rigid_motion_rot_center",
                                   "y":"des_rigid_motion_rot_center",
                                   "z":"des_rigid_motion_rot_center"},

    "DES rigid rotation speed" : {"x": "des_rigid_motion_omega",
                                  "y": "des_rigid_motion_omega",
                                  "z": "des_rigid_motion_omega"}
}

def dequote(txt):
    txt = txt.strip()
    while txt and txt[0]==txt[-1] and txt[0] in ('"',"'"):
        txt = txt[1:-1].strip()
    return txt

class Keyframe:
    #Keyframe Task Pane Window

    def init_keyframe(self):
        ui = self.ui.keyframe

        self.kf = {} # key: index.  value: data dictionary for internal surface
        self.kf_current_index = None
        self.kf_file_watcher = QFileSystemWatcher()
        self.kf_file_watcher.fileChanged.connect(self.kf_file_changed)
        self.kf_last_plotted_file = None
        self.kf_last_plotted_column = None

        #  The top of the task pane is where users define/select keyframes.
        ui.toolbutton_add.clicked.connect(self.kf_add_file)
        ui.toolbutton_add.key = "keyframe" # locatability
        ui.toolbutton_delete.clicked.connect(self.kf_delete)
        ui.toolbutton_delete.setEnabled(False) # Need a selection

        tw = ui.tablewidget_keyframes
        tw.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        tw.setSelectionBehavior(QAbstractItemView.SelectRows)
        tw.itemSelectionChanged.connect(self.handle_kf_selection)
        self.in_setup_keyframe = False
        cb = ui.combobox_keyword_omega_unit
        for (i, tip) in enumerate(("Radians/second", "Degrees/second", "Rotations/second",
                                   "Radians/minute", "Degrees/minute", "Rotations/minute")):
            get_combobox_item(cb,i).setToolTip(tip)

    def setup_keyframe(self, allow_disabled_tab=False):
        ui = self.ui.keyframe
        self.fixup_kf_table()
        if self.lock_level:
            for w in (ui.toolbutton_add, ui.toolbutton_delete):
                self.disable_widget(w)

        if self.in_setup_keyframe:
            return

        self.kf_clear_bottom_pane()
        self.in_setup_keyframe = True # Prevent keyword updates

        tw = ui.tablewidget_keyframes
        row = get_selected_row(tw)
        # Autoselect if only 1 row
        if row is None and tw.rowCount() == 1:
            row = 0
            self.in_setup_keyframe = False
            tw.selectRow(row)
            return

        if self.kf_current_index is None:
            self.in_setup_keyframe = False
            return

        idx = self.kf_current_index
        fname = self.kf.get(idx,{}).get("filename")
        if not fname:
            self.in_setup_keyframe = False
            return

        with open(fname) as f:
            try:
                csvreader = csv.reader(f)
                headers = next(csvreader)
            except Exception as e:
                self.error("Bad CSV file %s %s" %(fname, str(e)),
                           popup=True)
                self.in_setup_keyframe = False
                return
            for (col, header) in enumerate(headers, 1):
                label = "Column %s: %s" % (col, dequote(header))
                gb = QGroupBox(label)
                gb.setStyleSheet("QGroupBox {font-weight: bold;}")
                gb.setFlat(True)
                gb.key = "kf_keyword"
                ui.bottom_pane.addWidget(gb, col-1, 0)
                gb.setLayout(QGridLayout())
                gb.layout().setContentsMargins(5,5,5,5)
                layout = gb.layout()
                if col == 1:
                    layout.addWidget(QLabel("Independent variable (time)"), 0, 0)
                else:
                    layout.addWidget(QLabel("Dependent variable"), 0, 0)
                    cb = QComboBox()
                    cb.addItems(kf_keys.keys())
                    for i in range(cb.count()):
                        item = get_combobox_item(cb, i)
                        if item.text().startswith("Boundary"):
                            item.setEnabled(bool(self.bcs))
                        if item.text().startswith("Internal"):
                            item.setEnabled(bool(self.iss))
                        if item.text().startswith("DES"):
                            item.setEnabled(bool(self.solids))

                    layout.addWidget(cb, 0, 1)
                    cb.currentIndexChanged.connect(lambda val, col=col: self.kf_handle_combobox(col, val))
                    val = self.project.get_value("kf_keyword", default=None, args=[idx, col-1])
                    if val is None:
                        cb.setCurrentIndex(0) # Unassigned
                        self.kf_handle_combobox(col,0) # Display plot button
                    else:
                        vlower = val.split('(')[0].lower()
                        for (j, data) in enumerate(kf_keys.items()):
                            if vlower in data[1].values():
                                cb.setCurrentIndex(j)
                                break
                        else:
                            self.error("Invalid value '%s' for KF_KEYWORD(%s,%s)" % (val,idx,col), popup=True)
                            cb.setCurrentIndex(0)

                        if '(' in val:  #keyword has args
                            args = val.split('(')[1].split(')')[0].split(',')
                            cb2 = layout.itemAt(3)
                            if cb2:
                                cb2 = cb2.widget()
                                try:
                                    target = int(args[0])
                                except ValueError as e:
                                    self.error("Invalid arguments '%s': %s" % (val, str(e)), popup=True)
                                    target = None
                                for i in range(cb2.count()):
                                    indices = get_combobox_item(cb2, i).data(UserRole)
                                    if target in indices:
                                        cb2.setCurrentIndex(i)
                                        break
                                else:
                                    if target is not None:
                                        self.warning("Invalid index %s in %s" % (args[0], val), popup=True)
                            if len(args) > 1:
                                try:
                                    target = int(args[1]) - 1
                                except ValueError as e:
                                    self.error("Invalid index %s in %s" % (args[1], val), popup=True)
                                    target = None
                                if target is not None:
                                    cb3 = layout.itemAt(5)
                                    if cb3:
                                        cb3 = cb3.widget()
                                        cb3.setCurrentIndex(target)
                        for (i,x) in enumerate(('_u_', '_v_', '_w_')):
                            if x in val.lower():
                                cb3 = layout.itemAt(5)
                                if cb3:
                                    cb3 = cb3.widget()
                                    cb3.setCurrentIndex(i)
                                break

        # Interpolation
        key = "kf_interp"
        gb = QGroupBox("Interpolation")
        gb.setStyleSheet("QGroupBox {font-weight: bold;}")
        gb.setFlat(True)
        gb.key = key
        ui.bottom_pane.addWidget(gb, col, 0)
        gb.setLayout(QGridLayout())
        gb.layout().setContentsMargins(5,5,5,5)
        rb_linear = QRadioButton("Linear")
        rb_linear.clicked.connect(lambda _: self.set_kf_interp("LINEAR")) # TODO update plot
        gb.layout().addWidget(rb_linear,0,0)
        rb_step = QRadioButton("Step")
        rb_step.clicked.connect(lambda _: self.set_kf_interp("STEP"))  # TODO update plot
        gb.layout().addWidget(rb_step,0,1)

        vals = ("LINEAR", "STEP")
        args = [self.kf_current_index]
        val = self.project.get_value(key, default="LINEAR", args=args)
        if val.upper() not in vals:
            self.warning("Invalid KF_INTERP(%s) '%s', setting to LINEAR"%
                         (self.kf_current_index, val),
                         popup=True)
            val = "LINEAR"
            self.update_keyword(key, val, args=args)
        val = val.upper()
        rb_linear.setChecked(val=="LINEAR")
        rb_step.setChecked(val=="STEP")

        self.in_setup_keyframe = False
        self.fixup_kf_table()


    def copy_to_project(self, fname):
        basename = os.path.basename(fname)
        if QMessageBox.question(self,
                                   "Copy file?",
                                   "Copy CSV file to project directory?",
                                   QMessageBox.Yes | QMessageBox.No,
                                   QMessageBox.Yes) != QMessageBox.Yes:
            return False
        if os.path.exists(os.path.basename(fname)):
            if QMessageBox.question(self,
                                    "Replace?",
                                    "'%s' exists in project directory, replace it?" % basename,
                                    QMessageBox.Yes | QMessageBox.No,
                                    QMessageBox.Yes) != QMessageBox.Yes:
                return False
        try:
            shutil.copy2(fname,
                         os.path.join(os.getcwd(), basename))
        except Exception as e:
            self.error("File not copied: %s" % str(e),
                       popup=True)
            return False
        return True

    def kf_cancel_add(self):
        ui = self.ui.keyframe

        for item in (ui.toolbutton_add,
                     ui.tablewidget_keyframes):
            item.setEnabled(True)

        if get_selected_row(ui.tablewidget_keyframes) is not None:
            for item in (ui.bottom_frame,
                         ui.toolbutton_delete):
                item.setEnabled(True)

    def kf_extract_keyframes(self):
        if self.kf:
            # We assume that keyframe definitions have been initialized correctly
            # from mfix_gui_comments
            return
        # TODO handle non-GUI cases
        kwlist = list(self.project.keywordItems())
        if any (k.key.startswith(('kf_', 'keyframe')) for k in kwlist):
            self.warning("Keyframe keys are set but keyframe files not found.", popup=True)


    def save_keyframe_files(self, index=None):
        # Save all files if index==None
        if not self.kf:
            if os.path.exists('.KF'):
                shutil.rmtree('.KF')
            return True
        try:
            os.makedirs(".KF", exist_ok=True)
        except Exception as e:
            self.error("Cannot create directory .KF: %s" % str(e),
                       popup=True)
            return False
        for (idx, data) in self.kf.items():
            if index is not None and idx!=index:
                continue
            try:
                target = ".KF/keyframe_%04d.csv"%idx
                shutil.copyfile(data["filename"], target)
            except Exception as e:
                print(e)
                self.error("Cannot create file %s: %s" % (target, str(e)),
                           popup=True)
                return False
        return True

    def kf_add_file(self):
        # Interactively add files regions to define keyframes
        ui = self.ui.keyframe
        tw = ui.tablewidget_keyframes
        self.kf_cancel_add() # Re-enable input widgets
        fname, filter = QFileDialog.getOpenFileName(self, "Choose CSV file",
                                                    os.getcwd(),
                                                    "CSV files (*.csv)")
        if not fname:
            return
        fname = os.path.relpath(fname, os.getcwd())

        if '/.KF/' in fname or '\\.KF\\' in fname or fname.startswith(('.KF/','.KF\\')):
            self.error("Files in the '.KF' subdirectory are for internal use by the mfix solver and may be deleted or renamed without warning.",
                       popup=True)
            return

        if (os.path.realpath(os.path.dirname(fname)) !=
            os.path.realpath(os.getcwd())):
            if self.copy_to_project(fname):
                fname = os.path.basename(fname)

        index = 1 + len(self.kf)
        self.kf_add_file_1(fname, index, select=True)


    def kf_add_file_1(self, fname, index, select=False):
        if not os.path.exists(fname):
            self.error("File %s not found" % fname)
            self.unset_keyword("keyframe", args=[index])
            self.unset_keyword("kf_keyword", args=[index])
            self.unset_keyword("kf_interp", args=[index])
            return
        self.kf_file_watcher.addPath(fname)
        ui = self.ui.keyframe
        tw = ui.tablewidget_keyframes
        n = tw.rowCount() + 1
        self.kf[index] = {"filename": fname}
        tw.setRowCount(n)
        def make_item(val):
            item = QTableWidgetItem('' if val is None else str(val))
            set_item_noedit(item)
            return item
        item = make_item(os.path.basename(fname))
        item.setData(UserRole, fname)
        tw.setItem(n-1, 0, item)
        item = make_item(n)
        tw.setItem(n-1, 1, item)
        self.update_keyword("keyframe", True, args=[index])
        if select:
            tw.selectRow(n-1)

    def handle_kf_selection(self):
        ui = self.ui.keyframe
        tw = ui.tablewidget_keyframes
        row = get_selected_row(tw)
        self.enable_widget(ui.toolbutton_add)
        if row is None:
            self.kf_current_index = None
            self.kf_clear_bottom_pane()
            self.disable_widget(ui.toolbutton_delete)
        else:
            self.kf_current_index = row + 1
            self.setup_keyframe()
            self.enable_widget(ui.toolbutton_delete)


    def kf_handle_combobox(self, col, val):
        # 'col' is CSV file column
        ui = self.ui.keyframe
        kf_idx = self.kf_current_index
        if kf_idx is None: #?
            return
        gb = ui.bottom_pane.itemAt(col-1).widget()
        while True: # Clear bottom of combobox
            item = gb.layout().takeAt(2)
            if not item:
                break
            if item.widget():
                item.widget().deleteLater()
        gb.layout().update()
        row = 1
        if val == 0: # Unassigned
            self.unset_keyword('kf_keyword', args=[kf_idx, col]) # FIXME multi-region BCs and ICs
            key_data = None
        else:
            key_data = kf_keys[list(kf_keys.keys())[val]]
        if key_data:
            key = key_data[list(key_data.keys())[0]].lower()
            tw = None

            cb = QComboBox()
            if key.startswith('bc_'):
                tw = self.ui.boundary_conditions.tablewidget_regions # yuk
                gb.layout().addWidget(QLabel("Boundary"), row, 0)
            elif key.startswith('is_'):
                tw = self.ui.internal_surfaces.tablewidget_regions # yuk
                gb.layout().addWidget(QLabel("Surface"), row, 0)
            elif key.startswith('des_'):
                gb.layout().addWidget(QLabel("Solids phase"), row, 0)
                for (idx, name) in enumerate(self.solids, 1):
                    cb.addItem(name)
                    get_combobox_item(cb, idx-1).setData((idx,), UserRole)
            cb.currentIndexChanged.connect(lambda val, col=col: self.kf_handle_combobox_2(col, val))
            if tw:
                for trow in range(tw.rowCount()):
                    name = tw.item(trow,0).text()
                    indices, regions = tw.item(trow,0).data(UserRole)
                    cb.addItem(name)
                    # Why are args reversed for setData on combobox?
                    get_combobox_item(cb, trow).setData(indices,UserRole)
            gb.layout().addWidget(cb, row, 1)
            row += 1
            if len(key_data) > 1:
                gb.layout().addWidget(QLabel("Component"), row, 0)
                cb = QComboBox()
                cb.currentIndexChanged.connect(lambda val, col=col: self.kf_handle_combobox_3(col, val))
                for k in key_data.keys():
                    cb.addItem(k.upper())
                gb.layout().addWidget(cb, row, 1)
                row += 1
        pb = QPushButton("Plot")
        pb.setAutoDefault(False)
        pb.clicked.connect(lambda checked, fname=self.kf[kf_idx].get('filename'), col=col: self.kf_plot(fname, col))
        gb.layout().addWidget(pb, row, 0, 1, 2)
        gb.updateGeometry()
        self.set_kf_keyword(col)

    def kf_handle_combobox_2(self, col, val):
        kf_idx = self.kf_current_index
        if kf_idx is None: #?
            return
        self.set_kf_keyword(col)

    def kf_handle_combobox_3(self, col, val):
        kf_idx = self.kf_current_index
        if kf_idx is None: #?
            return
        self.set_kf_keyword(col)

    def set_kf_keyword(self, col):
        if self.in_setup_keyframe:
            return
        kf_idx = self.kf_current_index
        if kf_idx is None:
            return

        ui = self.ui.keyframe
        gb = ui.bottom_pane.itemAt(col-1).widget()
        cb = gb.layout().itemAt(1).widget()
        key_type = cb.currentText()
        key_data = kf_keys.get(key_type, {})

        if key_type == "Unassigned":
            self.unset_keyword("kf_keyword", args=[kf_idx])
            return
        cb2 = gb.layout().itemAt(3)
        if cb2 is None: #UI not set up
            return
        cb2 = cb2.widget()
        indices = get_combobox_item(cb2, cb2.currentIndex()).data(UserRole)
        region_index = indices[0] ## FIXME multi-region BC/IC
        if len(key_data) == 1:
            key = base_key = list(key_data.values())[0]
            if region_index:
                key += "(%s)"%region_index
            self.update_keyword("kf_keyword", key, args=[kf_idx, col-1])
        else:
            cb3 = gb.layout().itemAt(5)
            if cb3 is None: #UI not set up
                return
            cb3 = cb3.widget()
            comp = cb3.currentText()
            key = base_key = key_data[comp.lower()]
            args = [kf_idx, col-1]
            key += "(%s"%region_index
            if base_key.endswith(("_vel", "_center", "_omega")): # inconsistency, sigh
                key += ",%s"%(1+cb3.currentIndex())
            key += ")"
            self.update_keyword("kf_keyword", key, args=args)


    def set_kf_interp(self, val):
        kf_idx = self.kf_current_index
        if not kf_idx:
            return
        self.update_keyword("kf_interp", val, args=[kf_idx])
        fname = self.kf[kf_idx].get('filename')
        if fname == self.kf_last_plotted_file and plt.get_fignums():
            # update plot window if open and showing this file
            self.kf_plot(self.kf_last_plotted_file, self.kf_last_plotted_column)


    def kf_plot(self, fname, col):
        if not fname:
            self.error("No filename specified", popup=True)
            return
        try:
            with open(fname) as f:
                csvreader = csv.reader(f)
                header = next(csvreader)
                xs = []
                ys = []
                for row in csvreader:
                    if not (row and row[0] and len(row)>col-1 and row[col-1]):
                        # Blank or short line
                        continue
                    # TODO popup error on bad data.
                    xs.append(float(row[0]))
                    ys.append(float(row[col-1]))
            interp = self.project.get_value('kf_interp', args=[self.kf_current_index], default='LINEAR')
            linewidth = 2.5
            plt.clf()
            if interp.lower() == 'step':
                plt.step(xs, ys, 'o-', linewidth=linewidth, where='post')
            else:
                plt.plot(xs,ys, 'o-', linewidth=linewidth)
            title = dequote(header[col-1])
            plt.xlabel("Simulation time (s)")
            plt.ion()
            plt.grid(True)
            plt.title(title)
            plt.get_current_fig_manager().set_window_title(os.path.split(fname)[-1])
            plt.show()
            self.kf_last_plotted_file = fname
            self.kf_last_plotted_column = col

        except Exception as e:
            self.error("Cannot plot %s: %s" % (fname, str(e)),
                       popup=True)
            return


    def kf_clear_bottom_pane(self):
        ui = self.ui.keyframe
        grid = ui.bottom_pane
        while True:
            item = grid.takeAt(0)
            if not item:
                break
            if item.widget():
                item.widget().deleteLater()
        grid.update()

    def reset_keyframe(self):
        self.kf.clear()
        self.kf_current_index = None
        for f in self.kf_file_watcher.files():
            self.kf_file_watcher.removePath(f)
        ui = self.ui.keyframe
        ui.tablewidget_keyframes.clearContents()
        ui.tablewidget_keyframes.setRowCount(0)
        self.in_setup_keyframe = False

    def kf_delete(self):
        ui = self.ui.keyframe
        tw = ui.tablewidget_keyframes
        row = get_selected_row(tw)
        if row is None: # No selection
            return
        fname = tw.item(row,0).text()
        self.kf_file_watcher.removePath(fname)
        for idx in range(self.kf_current_index, max(self.kf.keys())):
            self.update_keyword("keyframe", self.project.get_value("keyframe", args=[idx+1]),
                                args=[idx])
            kwlist = list(self.project.keywordItems())
            for kw in kwlist:
                key, args = kw.key, kw.args
                if key.startswith('kf_') and args and args[0] == idx+1:
                    val = self.project.get_value(key, args=args)
                    self.update_keyword(key, val, args=[idx]+args[1:])
                kwlist = list(self.project.keywordItems())
                for kw in kwlist:
                    key, args = kw.key, kw.args
                    if key.startswith('kf_') and args and args[0] == idx:
                        a = self.project.get_value(key, args=args)
                        b = self.project.get_value(key, args=[idx+1]+args[1:])
                        if a != b:
                            self.unset_keyword(key, args=args)

            self.kf[idx] = self.kf[idx+1]

        idx = max(self.kf.keys())
        del self.kf[idx]
        self.unset_keyword("keyframe", args=[idx])
        kwlist = list(self.project.keywordItems())
        for kw in kwlist:
            key, args = kw.key, kw.args
            if key.startswith('kf_') and args and args[0] == idx:
                self.unset_keyword(key, args=args)

        self.kf_current_index = None
        tw.removeRow(row)
        for i in range(tw.rowCount()):
            tw.item(i,1).setText(str(i+1))
        self.fixup_kf_table()
        self.setup_keyframe()
        self.handle_kf_selection() # refresh bottom pane


    def fixup_kf_table(self):
        ui = self.ui.keyframe
        hv = QHeaderView
        tw = ui.tablewidget_keyframes
        resize = tw.horizontalHeader().setSectionResizeMode
        ncols = tw.columnCount()
        for n in range(0, ncols):
            resize(n, hv.Stretch if n == 0 else hv.ResizeToContents)

        # trim excess vertical space - can't figure out how to do this in designer
        header_height = tw.horizontalHeader().height()

        # Note - scrollbar status can change outside of this function.
        # Do we need to call this every time window geometry changes?
        scrollbar_height = (tw.horizontalScrollBar().isVisible()
                            * (4+tw.horizontalScrollBar().height()))
        nrows = tw.rowCount()
        if nrows == 0:
            row_height = 0
            height = header_height+scrollbar_height
        else:
            row_height = tw.rowHeight(0)
            height = (header_height+scrollbar_height
                      + nrows*row_height + 4) # extra to avoid unneeded scrollbar

        icon_height = sub_icon_height() + 8
        ui.top_frame.setMaximumHeight(height+icon_height)
        ui.top_frame.setMinimumHeight(header_height+icon_height+row_height*min(nrows, 5))
        ui.top_frame.updateGeometry()
        tw.setMaximumHeight(height)
        tw.setMinimumHeight(header_height)
        tw.updateGeometry() #? needed?

    def kf_file_changed(self, fname):
        if os.path.exists(fname):
            self.kf_file_watcher.addPath(fname)
        if fname == self.kf_last_plotted_file and plt.get_fignums():
            # update plot window if open and showing this file
            self.kf_plot(self.kf_last_plotted_file, self.kf_last_plotted_column)

        for (i, data) in self.kf.items():
            if data.get('filename') == fname:
                self.save_keyframe_files(i)


    def kf_to_str(self):
        return JSONEncoder().encode(self.kf)

    def kf_from_str(self, s):
        ui = self.ui.keyframe
        if not s:
            return
        d = JSONDecoder().decode(s)
        for (index, data) in d.items():
            self.kf_add_file_1(data.get("filename", "???"), int(index))
            self.fixup_kf_table()
