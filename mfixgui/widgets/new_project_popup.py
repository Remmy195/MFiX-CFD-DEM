# -*- coding: utf-8 -*-
import os

from qtpy.QtWidgets import QDialog, QFileDialog
from qtpy.QtGui import QRegExpValidator
from qtpy.QtCore import QRegExp

from mfixgui.tools.qt import get_icon, get_ui, SETTINGS
from mfixgui.regexes import re_valid_run_name_qt

class NewProjectDialog(QDialog):
    """ Dialog for selecting the RUN_NAME & path when creating a new project
    (from a template) """

    def __init__(self, parent, run_name):
        QDialog.__init__(self, parent)
        ui = self.ui = get_ui('new_project_popup.ui', self)
        self.setWindowTitle('Create a new project')
        ui.lineedit_project_name.setText(run_name)
        self.ok_button = ui.buttonBox.button(ui.buttonBox.Ok)

        ui.toolbutton_browse.clicked.connect(self.browse)
        ui.toolbutton_browse.setIcon(get_icon('folder.svg'))

        ui.lineedit_project_name.setValidator(QRegExpValidator(
            QRegExp(re_valid_run_name_qt)))
        ui.lineedit_project_name.textChanged.connect(lambda text:
                                                     self.update_ok_button(text))
        ui.combobox_location.editTextChanged.connect(lambda loc:
            ui.buttonBox.button(ui.buttonBox.Ok).setEnabled(os.path.isdir(loc)))
        locs = [loc.rstrip(os.path.sep) for loc in SETTINGS.value('project_locations', '').split(',')
                if os.path.isdir(loc)] or [os.path.expanduser('~').rstrip(os.path.sep)]
        # filter dups resulting from bug 1036
        seen = set()
        tmp = []
        for l in locs:
            if l not in seen:
                seen.add(l)
                tmp.append(l)
        self.ui.combobox_location.addItems(tmp)
        self.update_ok_button(run_name)

    def update_ok_button(self, text):
        self.ok_button.setEnabled(bool(text))

    def get_name_and_location(self):
        """ Returns tuple (RUN_NAME, project_dir), or None if user cancels """

        if self.exec_() != QDialog.Accepted:
            return None

        cb = self.ui.combobox_location
        le = self.ui.lineedit_project_name
        locs = [cb.itemText(i) for i in range(cb.count())]
        # User may have typed in a directory name instead of used 'browse', so
        #  make current project dir default for next time.
        loc = cb.currentText().rstrip(os.path.sep)
        #New one to the head of the list
        if loc in locs:
            locs.remove(loc)
        locs.insert(0, loc)
        # only save 5 most recent ones
        locs = locs[:5]
        # What if dirname has a comma in it?
        SETTINGS.setValue('project_locations', ','.join(locs))
        SETTINGS.sync()
        cb = self.ui.combobox_location
        cb.clear()
        cb.addItems(locs)
        cb.setCurrentIndex(0)
        return (le.text(), loc)


    def browse(self):
        cb = self.ui.combobox_location
        loc = QFileDialog.getExistingDirectory(self, 'Location', cb.currentText())
        if not loc:
            return
        loc = loc.rstrip(os.path.sep)
        locs = [cb.itemText(i) for i in range(cb.count())]
        #New one to the head of the list
        if loc in locs:
            locs.remove(loc)
        locs.insert(0, loc)
        #Keep most recent 5
        locs = locs[:5]
        cb.clear()
        cb.addItems(locs)
        cb.setCurrentIndex(0)
        SETTINGS.setValue('project_locations', ','.join(locs))
        SETTINGS.sync()
