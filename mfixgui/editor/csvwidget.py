#!/usr/bin/env python

# Written by cgw (Charles G Waldman) 2024-09-29
### TODO:
# Input type verification (are values always numeric?)
# Drag/drop rows and columns
# Clean up separation of signals between CSVWidget and CSVEditor

import sys, os
from pandas import read_csv
from qtpy.QtWidgets import (QAction, QWidget, QLineEdit, QInputDialog,
                            QMessageBox, QTableWidget, QTableWidgetItem,
                            QGridLayout, QFileDialog, QMenu, QInputDialog,
                            QAbstractItemView, QHeaderView, QPushButton)

import warnings

from qtpy.QtCore import Qt, Signal
class CSVWidget(QTableWidget):

    def __init__(self, parent=None):
        QTableWidget.__init__(self, parent)
        self.parent = parent
        header = self.horizontalHeader()
        le = self.lineedit = QLineEdit(parent=header.viewport())
        le.editingFinished.connect(self.set_header_data)
        le.hide()
        header.sectionDoubleClicked.connect(self.edit_header)
        #header.setCascadingSectionResizes(True)
        #header.setSectionsMovable(True) - reorder cols?
        col_menu = self.col_menu = QMenu()
        self.col_actions = {name: col_menu.addAction(name)
                            for name in ("Add column right",
                                         "Add column left",
                                         "Delete column",
                                         "Save",
                                         "Set column to constant",
                                         "Rename column (double-click)")}

        col_menu.triggered.connect(self.handle_col_menu)
        header.setContextMenuPolicy(Qt.CustomContextMenu)
        header.customContextMenuRequested.connect(self.show_col_menu)

        row_menu = self.row_menu = QMenu()
        self.row_actions = {name: row_menu.addAction(name)
                            for name in ("Add row above",
                                         "Add row below",
                                         "Delete row")}

        row_menu.triggered.connect(self.handle_row_menu)
        vheader = self.verticalHeader()
        vheader.setContextMenuPolicy(Qt.CustomContextMenu)
        vheader.customContextMenuRequested.connect(self.show_row_menu)
        self.cellChanged.connect(self.handle_cell_changed)
        self.default_triggers = self.editTriggers()
        self.unsaved_flag = False
        self.read_only = False
        self.set_no_data()


    def set_no_data(self):
        self.data = None
        self.num_rows = self.num_cols = 1
        self.setRowCount(self.num_rows)
        self.setColumnCount(self.num_cols)
        self.column_names = ['1']
        self.setHorizontalHeaderLabels(self.column_names)


    def handle_cell_changed(self, row, col):
        if self.loading:
            return
        self.unsaved_flag = True
        self.needsSaved.emit()

    def warning(self, str):
        print("Warning: ", str)

    def error(self, str):
        print("Error: ", str)

    def load_file(self, fname):
        def cleanup(s): # Strip and dequote
            s = str(s).strip()
            if not s:
                return ''
            while len(s) > 1 and s[0]==s[-1] and s[0] in ('"', "'"):
                s = s[1:-1].strip()
            return s

        self.loading = True
        # Don't fail on empty file
        if os.stat(fname).st_size == 0:
            self.num_rows = self.num_cols = 0
            self.column_names = ['1']
            self.num_cols = 1
            self.num_rows = 0
        else:
            try:
                with warnings.catch_warnings(record=True) as ws:
                    df = read_csv(fname, skip_blank_lines=True,
                                  header='infer',
                                  comment='#',
                                  on_bad_lines='warn')
                for w in ws:
                    self.warning(str(w.message))
            except Exception as e:
                self.error("Error loading %s: %s" % (fname, e))
                self.set_no_data()
                self.loading = False
                return False
            self.num_rows, self.num_cols = df.shape
            self.setRowCount(max(self.num_rows, 1))
            self.setColumnCount(max(self.num_cols, 1))
            self.column_names = [cleanup(x) for x in df.columns]

        self.setHorizontalHeaderLabels(self.column_names)
        for i in range(self.num_rows):
            for j in range(self.num_cols):
                self.setItem(i, j, QTableWidgetItem(cleanup(str(df.iat[i, j]))))
        self.loading = False
        self.set_column_width()
        #self.resize_table()
        return True

    def set_column_width(self):
        header = self.horizontalHeader()
        # This doesn't allow the user to resize columns
        #for i in range(self.num_cols-1):
        #    header.setSectionResizeMode(i, QHeaderView.ResizeToContents)
        #if self.num_cols: # Stretch last column
        #    header.setSectionResizeMode(self.num_cols-1, QHeaderView.Stretch)
        fm = self.fontMetrics()
        for j in range(self.num_cols-1):
            header.setSectionResizeMode(j, QHeaderView.Interactive)
        if self.num_cols: # Stretch last column
            header.setSectionResizeMode(self.num_cols-1, QHeaderView.Stretch)
        for j in range(self.num_cols):
            w = max(fm.boundingRect(t).width()
                    for t in [self.column_names[j]] +
                    [self.item(i,j).text() if self.item(i,j) else ''
                     for i in range(self.num_rows)])
            header.resizeSection(j, w+20)

    def save_file(self, fname):
        if not fname:
            return False
        tmpfile = fname + '.tmp'
        with open(tmpfile, 'w') as f:
            cv = self.parent.comment_view
            comments = cv.toPlainText()
            if comments:
                f.write(comments)
                if not comments.endswith(os.linesep):
                    f.write(os.linesep)
            def quote(s):
                return '"'+s+'"' if ',' in s else s
            header = ','.join(quote(x.strip()) for x in self.column_names)
            f.write(header.rstrip() + os.linesep)
            for i in range(self.num_rows):
                row = ','.join(str(self.item(i,j).text()
                                    if self.item(i,j) else '').strip()
                                for j in range(self.num_cols))
                f.write(row.rstrip() + os.linesep)

        os.replace(tmpfile, fname)

    def resize_table(self):
        return #Not ready for prime time
        h = 0
        vheader = self.verticalHeader()
        for i in range(vheader.count()):
            print(i)
            if not vheader.isSectionHidden(i):
                h += vheader.sectionSize(i)
                print(h)
            #if not self.horizontalScrollBar().isHidden():
            #    h += self.horizontalScrollBar().height()
            #    print("SB", h)
            if not self.horizontalHeader().isHidden():
                print("HH", h)
                h += self.horizontalHeader().height()
            print(h)
        self.setMaximumHeight(h)
        self.setMinimumHeight(h)

    def show_col_menu(self, pos):
        header = self.horizontalHeader()
        self.menu_column = header.logicalIndexAt(pos)
        self.col_menu.exec(header.mapToGlobal(pos))
        self.clearSelection()


    def handle_col_menu(self, action):
        for name in self.col_actions:
            if self.col_actions[name] == action:
                if name == 'Add column right':
                    self.insertColumn(self.menu_column+1)
                    self.column_names.insert(self.menu_column+1, str(self.menu_column+2))
                    self.setHorizontalHeaderLabels(self.column_names)
                    self.set_column_width()
                    self.num_cols += 1
                elif name == 'Add column left':
                    self.insertColumn(self.menu_column)
                    self.column_names.insert(self.menu_column, str(self.menu_column+1))
                    self.setHorizontalHeaderLabels(self.column_names)
                    self.set_column_width()
                    self.num_cols += 1
                elif name == 'Delete column':
                    self.removeColumn(self.menu_column)
                    del self.column_names[self.menu_column]
                    self.num_cols -= 1
                elif name.startswith("Rename"):
                    self.edit_header(self.menu_column)
                elif name.startswith("Set"):
                    reply, ok = QInputDialog.getText(self, "Set constant value",
                                               "Enter constant value:")
                    if ok: # False on cancel, (val, True) on OK
                        val = str(reply.strip())
                        for i in range(self.num_rows):
                            if item := self.item(i, self.menu_column):
                                item.setText(val)
                            else:
                                self.setItem(i, self.menu_column, QTableWidgetItem(val))
                break

    def show_row_menu(self, pos):
        vheader = self.verticalHeader()
        self.menu_row = vheader.logicalIndexAt(pos)
        self.row_menu.exec(vheader.mapToGlobal(pos))
        self.clearSelection()

    def handle_row_menu(self, action):
        for name in self.row_actions:
            if self.row_actions[name] == action:
                if name == 'Add row below':
                    self.insertRow(self.menu_row+1)
                    self.num_rows += 1
                elif name == 'Add row above':
                    self.insertRow(self.menu_row)
                    self.num_rows += 1
                elif name == 'Delete row':
                    self.removeRow(self.menu_row)
                    self.num_rows -= 1
                else:
                    break
        #self.resize_table()

    def set_header_data(self):
        le = self.lineedit
        text = le.text().strip()
        self.lineedit.hide()
        if not text:
            return
        self.column_names[self.edit_section] = text
        self.setHorizontalHeaderLabels(self.column_names)
        self.set_column_width()
        self.clearSelection()
        self.unsaved_flag = True
        self.needsSaved.emit()

    def edit_header(self,section):
        self.cursor_changed.emit(0, section+1)
        self.edit_section = section
        le = self.lineedit
        header = self.horizontalHeader()
        geo = le.geometry()
        geo.setWidth(header.sectionSize(section))
        geo.setHeight(header.height())
        geo.moveLeft(header.sectionViewportPosition(section))
        le.setGeometry(geo)
        le.setText('')
        le.setFocus()
        le.show()

    def setReadOnly(self, readonly):
        self.read_only = readonly
        self.setEditTriggers(QAbstractItemView.NoEditTriggers
                             if readonly else self.default_triggers)

    def isReadOnly(self):
        return self.read_only

    def goto(self, line):
        # Subtract 1 for 1-based line count, another for header line
        self.selectRow(max(line-2, 0))
