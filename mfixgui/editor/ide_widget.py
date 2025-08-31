"""
Integrated Development Environment for MFiX.
"""

import os
import shutil

from qtpy.QtWidgets import (QWidget, QSplitter, QTabWidget, QHBoxLayout,
                            QTabBar, QFileDialog, QShortcut, QMessageBox,
                            QInputDialog)
from qtpy.QtCore import Qt, Signal
from qtpy.QtGui import QKeySequence

from mfixgui.editor.tree import FileWidget
from mfixgui.editor.code_editor import CodeEditorWidget
from mfixgui.editor.csv_editor import CSVEditorWidget
from mfixgui.tools import get_mfix_src, get_full_path
from mfixgui.tools.qt import get_icon, sub_icon_size


class IDEWidget(QWidget):
    unsaved_signal = Signal()
    file_changed = Signal(str)
    help_request = Signal() # Should we include the string?
    def __init__(self, parent=None, mfixgui=None):
        QWidget.__init__(self, parent)

        self.mfixgui = mfixgui

        self.current_directory = '.'
        self.project_dir = None
        self.project_file = None

        layout = QHBoxLayout(self)
        layout.setContentsMargins(5, 5, 5, 5)

        splitter = QSplitter(self)
        layout.addWidget(splitter)

        # Tree
        self.tree = FileWidget(self)
        self.tree.mfix_source_file_clicked.connect(self.open_from_source)
        self.tree.project_file_clicked.connect(self.open)
        self.tree.new_file.connect(self.new_file)
        splitter.addWidget(self.tree)

        # tab widget
        self.tabs = QTabWidget()
        # configure tab widget
        self.tabs.setTabsClosable(True)
        self.tabs.tabCloseRequested.connect(self.remove_tab)
        tab_bar = self.tabs.tabBar()
        tab_bar.setElideMode(Qt.ElideRight)

        # create a tab at the end to act like "new" button
        tab_bar.setSelectionBehaviorOnRemove(QTabBar.SelectLeftTab)
        new_tab = QWidget()
        self.tabs.addTab(new_tab, '')
        idx = self.tabs.count()-1
        self.tabs.setTabIcon(idx, get_icon('add.svg'))
        self.tabs.setIconSize(sub_icon_size())
        tab_bar.tabBarClicked.connect(self.handle_tab_clicked)

        # remove close btn
        right_btn = self.tabs.tabBar().tabButton(idx, QTabBar.RightSide)
        if right_btn:
            right_btn.resize(0, 0)
        left_btn = self.tabs.tabBar().tabButton(idx, QTabBar.LeftSide)
        if left_btn:
            left_btn.resize(0, 0)

        splitter.addWidget(self.tabs)

        splitter.setSizes([120, 600])

        # keyboard shortcuts
        shortcuts = [
            ('Ctrl+O', self.open),
            #('Ctrl+N', self.new_tab),
            ]

        for key, callback in shortcuts:
            q = QShortcut(QKeySequence(key), self)
            q.activated.connect(callback)

    def handle_tab_clicked(self, idx):
        # This is a little crazy, it would be nice if we could attach a callback directly
        #  to the tab with the "+" icon
        if (idx >=0
            and idx == self.tabs.count()-1
            and self.tabs.tabText(idx)==''
            and not self.tabs.tabIcon(idx).isNull()):
            self.new_file()

    def new_file(self):
        fname, ok = QInputDialog.getText(self, "New file", "Enter name for new file")
        if not ok:
            return
        if not fname:
            return
        if os.path.exists(fname):
            self.mfixgui.error("File %s exists"%fname, popup=True)
            return
        # Create the file
        try:
            d = os.path.dirname(fname)
            if d:
                os.makedirs(os.path.dirname(fname), exist_ok=True)
            f = open(fname, 'w')
            f.close()
        except Exception as e:
            self.mfixgui.error("Error creating %s: %s" % (fname, e), popup=True)
            return
        self.open(fname)


    def show_location(self, path, lineno, column=0, error=False):
        fullpath = get_full_path(self.project_dir, path)
        read_only = not fullpath.startswith(os.path.realpath(self.project_dir))
        open = self.open_from_source if read_only else self.open
        open(fullpath, read_only=read_only).editor.goto(lineno, column=column, error=error)

    def open_from_source(self, fname, read_only=False):
        editor = self.open(fname, read_only=True)
        editor.copy_project_widget.setVisible(True)
        return editor

    def copy_to_project(self, fname):
        # copy file to project directory
        if self.project_dir is not None:
            basename = os.path.basename(fname)
            cfilename = os.path.join(self.project_dir, basename)

            # check if file exists in project
            if os.path.exists(cfilename):
                rsp = QMessageBox.warning(
                    self, 'Overwrite file?',
                    'Warning: File exists in the project directory. Overwrite file?',
                    QMessageBox.Yes | QMessageBox.Cancel)

                if rsp == QMessageBox.Yes:
                    shutil.copy(fname, cfilename)
            else:
                shutil.copy(fname, cfilename)

            # remove read only tab
            open_index = self.is_open(fname)
            self.remove_tab(open_index)

            # open file in project
            editor = self.open(cfilename)
            editor.editor.setReadOnly(False)

            # clean up tree selection
            self.tree.project_tree.select(cfilename)
            self.tree.mfix_tree.clearSelection()
            return editor

    def open(self, fname=None, read_only=False, goto_tab=True, closeable=True,
             project_file=False):
        if fname is None:
            fname = QFileDialog.getOpenFileName(self)[0]
            if not fname:
                return

        if project_file:
            self.project_file = fname

        open_index = self.is_open(fname)
        if open_index > -1:
            if goto_tab:
                self.tabs.setCurrentIndex(open_index)
            editor = self.tabs.widget(open_index)
        else:
            tab_name = os.path.basename(fname)
            if read_only:
                tab_name += ' [read only]'
            editor = self.new_tab(fname=tab_name, closeable=closeable)
            editor.open(fname, read_only)
        return editor

    def open_project(self, prj):
        self.project_dir = prj
        self.tree.set_project_tree(prj)

    def new_tab(self, fname=None, closeable=True):

        if fname and fname.lower().endswith('.csv'):
            w = CSVEditorWidget(self, self.mfixgui)
        else:
            w = CodeEditorWidget(self, self.mfixgui)
        index = self.tabs.count()-1
        w.editor.needsSaved.connect(lambda w=w: self.unsaved(w))
        w.editor.help_request.connect(self.help_request.emit)
        w.editor.saved.connect(lambda w=w: self.handle_saved(w))
        self.tabs.insertTab(index, w, fname if fname is not None else 'untitled')
        self.tabs.setCurrentIndex(index)
        if not closeable:
            self.hide_close_button(index)
        return w

    def hide_close_button(self, idx=None):
        right_btn = self.tabs.tabBar().tabButton(idx, QTabBar.RightSide)
        if right_btn:
            right_btn.resize(0, 0)
        left_btn = self.tabs.tabBar().tabButton(idx, QTabBar.LeftSide)
        if left_btn:
            left_btn.resize(0, 0)

    def is_open(self, fname):
        """ see if fname is already open, return index else -1 """
        fname = os.path.abspath(fname)

        for index in range(self.tabs.count()):
            tab = self.tabs.widget(index)
            if hasattr(tab, 'editor') and tab.editor.file_name is not None and fname == os.path.abspath(tab.editor.file_name):
                return index
        return -1

    def remove_tab(self, index, check_unsaved=True):
        """ handle signal to close tab """
        tab = self.tabs.widget(index)

        if check_unsaved and tab.editor.unsaved_flag:
            message_box = QMessageBox(self)
            message_box.setWindowTitle("Save?")
            message_box.setIcon(QMessageBox.Warning)
            message_box.setText("Save changes before closing?")
            for b in [QMessageBox.Yes, QMessageBox.No, QMessageBox.Cancel]:
                message_box.addButton(b)
            message_box.setDefaultButton(QMessageBox.Cancel)
            ret = message_box.exec_()

            if ret == QMessageBox.Yes:
                success = tab.save()
                if not success:
                    return False
            elif ret == QMessageBox.Cancel:
                return False

        self.tabs.removeTab(index)
        tab.deleteLater()
        return True

    def remove_untitled_tabs(self):
        for index in range(self.tabs.count()):
            tab = self.tabs.widget(index)
            if hasattr(tab, 'editor') and tab.editor.file_name is None:
                self.remove_tab(index, check_unsaved=False)

    def reset(self):
        for i in range(self.tabs.count()-1):
            self.remove_tab(0)

    def unsaved(self, tab):
        index = self.tabs.indexOf(tab)
        t = self.tabs.tabText(index)
        if not t.startswith('*'):
            self.tabs.setTabText(index, '*' + t)
        self.unsaved_signal.emit()

    def check_unsaved(self):
        for index in range(self.tabs.count()):
            tab = self.tabs.widget(index)
            if hasattr(tab, 'editor'):
                if tab.editor.unsaved_flag:
                    return True
        return False

    def handle_saved(self, tab):
        index = self.tabs.indexOf(tab)
        self.tabs.setTabText(index, os.path.basename(tab.editor.file_name))
        self.file_changed.emit(tab.editor.file_name)
