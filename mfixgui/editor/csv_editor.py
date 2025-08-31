'''
Simple CSV editor
'''
import os
import copy

### TODO
# add/remove column/row needs to set unsaved_flag
#  unsaved_flag needs to be split into project/other files

from qtpy.QtWidgets import (QAction, QFileDialog, QHBoxLayout,
    QInputDialog, QLabel, QListWidget, QMenu, QPlainTextEdit,
    QPushButton, QShortcut, QShortcut, QTableWidget, QToolButton,
    QVBoxLayout, QWidget)

from qtpy.QtGui import QKeySequence, QIcon, QAction
from qtpy.QtCore import (Signal, QTimer, QFileSystemWatcher)

from mfixgui.tools.qt import get_ui, get_icon, get_pixmap, sub_icon_size, SETTINGS
from mfixgui.editor.constants import DEFAULT_COLOR_SCHEME
from mfixgui.editor.syntax_highlighter import get_syntax_highlighter
from mfixgui.editor.completer import get_completer
from mfixgui.editor.widgets import (FindReplaceWidget, ReloadWidget,
    InfoBarWidget, CopyToProjectWidget)

from mfixgui.editor.csvwidget import CSVWidget


class CSVEditorWidget(QWidget):

    def __init__(self, parent=None, mfixgui=None):
        QWidget.__init__(self, parent)

        self.ide_widget = parent
        self.mfixgui = mfixgui
        self.file_name = None
        self.editor = CSVEditor(self, mfixgui)
        self.editor.changed_on_disk.connect(self.show_reload)
        self.find_widget = FindReplaceWidget(self, self.editor)
        self.find_widget.setVisible(False)
        self.reload_widget = ReloadWidget(self, self.editor)
        self.reload_widget.setVisible(False)
        self.info_bar = InfoBarWidget(self, self.editor)
        self.copy_project_widget = CopyToProjectWidget(self)
        self.copy_project_widget.setVisible(False)

        self.comments = []
        self.comment_view = QPlainTextEdit()
        self.comment_view.setVisible(False)
        self.comment_view.document().contentsChanged.connect(self.comments_edited)
        self.inserting_comments = False

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self.reload_widget)
        layout.addWidget(self.copy_project_widget)
        layout.addWidget(self.comment_view)
        layout.addWidget(self.editor)
        # "Add row" button
        #b = self.add_row_button = QToolButton()
        #b.setIcon(get_icon("add.svg"))
        #b.setIconSize(sub_icon_size())
        #b.setAutoRaise(True)
        #w = QWidget()
        #l2 = QHBoxLayout(w)
        #l2.setContentsMargins(0, 0, 0, 0)
        #w.setLayout(l2)
        #l2.addWidget(b)
        #l2.addStretch(1)
        #layout.addWidget(w)
        #layout.addStretch(1)
        layout.addWidget(self.find_widget)
        layout.addWidget(self.info_bar)
        self.editor.cursor_changed.connect(self.info_bar.set_position)

        # keyboard shortcuts
        shortcuts = [
            ('Ctrl+S', self.save), # FIXME does not work if enabled in main app
            ('Ctrl+G', self.goto),
            ('Esc', self.close_dialogs)]

        for key, callback in shortcuts:
            q = QShortcut(QKeySequence(key), self)
            q.activated.connect(callback)


    def contextMenuEvent(self, event):
        menu = QMenu()
        # &Save does not work, but it does in code_editor (?)
        save_action = QAction(QIcon(""), "Save", menu)
        save_action.setEnabled(self.editor.unsaved_flag)
        save_action.triggered.connect(lambda checked: self.save(None))
        menu.addAction(save_action)
        goto_action = QAction(QIcon(""), "Goto line", menu)
        goto_action.triggered.connect(self.goto)
        menu.addAction(goto_action)
        menu.exec_(event.globalPos())

    def open(self, fname=None, read_only=False):
        # Look for comment block at top of file
        self.comments = []
        try:
            with open(fname) as f:
                for line in f:
                    if line.startswith('#'):
                        self.comments.append(line)
                    else:
                        break
        except Exception as e:
            self.mfixgui.error("Error %s loading file" % e)
            return
        if self.comments:
            self.comment_view.setVisible(True)
            self.comment_view.setReadOnly(read_only)
            self.inserting_comments = True
            self.comment_view.clear()
            for line in self.comments:
                self.comment_view.appendPlainText(line.rstrip())
            self.inserting_comments = False
            self.adjust_comment_view()
        else:
            self.inserting_comments = True
            self.comment_view.clear()
            self.inserting_comments = False
            self.comment_view.setVisible(False)

        self.editor.open(fname, read_only)
        self.info_bar.set_path(os.path.abspath(fname))
        self.info_bar.set_language("CSV")

    def save(self, fname=None):
        self.reload_widget.hide()
        # Emitting the signal here results in a double reload.  Our own
        # file watcher catches the file changed event, even when it's us
        # who is writing the file - cgw
        #self.ide_widget.file_changed.emit(self.editor.file_name)
        ok =  self.editor.save(fname)
        if ok:
            self.mfixgui.print_internal("Saved %s" % self.editor.file_name,
                                        color="blue")
        return ok

    def reload(self):
        self.editor.open(self.editor.file_name, self.editor.isReadOnly())
        self.ide_widget.file_changed.emit(self.editor.file_name)
        self.error_line = None
        self.reload_widget.hide()

    def copy_to_project(self):
        # Preserve current cursor position, error highlight, selection and scrollbar
        editor = self.ide_widget.copy_to_project(self.editor.file_name).editor
        editor.error_line = self.editor.error_line
        # editor.setTextCursor(self.editor.textCursor())  # causing issue #1788
        editor.verticalScrollBar().setValue(
            self.editor.verticalScrollBar().value())


    def goto(self):
        l, ok = QInputDialog.getText(self, "Goto line", "Enter line number")
        if not ok:
            return
        try:
            l = int(l)
        except ValueError:
            return
        self.editor.goto(l)

    def show_reload(self):
        if not os.path.exists(self.editor.file_name):
            return
        if self.editor.unsaved_flag:
            self.reload_widget.show()
        else:
            self.reload()

    def comments_edited(self):
        if self.inserting_comments:
            return
        self.editor.unsaved_flag = True
        self.editor.needsSaved.emit()
        self.adjust_comment_view()

    def close_dialogs(self):

        if self.editor.completer.isVisible():
            self.editor.completer.hide()

        else:
            for w in [self.find_widget]:
                if w.isVisible():
                    w.hide()


    def adjust_comment_view(self):
        #https://stackoverflow.com/questions/45028105/get-the-exact-height-of-qtextdocument-in-pixels
        doc = self.comment_view.document()
        layout = doc.documentLayout()
        h = 0
        b = doc.begin()
        while b != doc.end():
            h += layout.blockBoundingRect(b).height()
            b = b.next()
        self.comment_view.setFixedHeight(int(h
                                             + doc.documentMargin()
                                             + 2 * self.comment_view.frameWidth()
                                             + 1))


class CSVEditor(CSVWidget):
    needsSaved = Signal()
    saved = Signal()
    changed_on_disk = Signal()
    cursor_changed = Signal(int, int)
    help_request = Signal()

    def __init__(self, parent=None, mfixgui=None):
        CSVWidget.__init__(self, parent)
        self.parent = parent
        self.mfixgui = mfixgui
        self.unsaved_flag = False
        self.loading = False
        self.file_name = None
        self.watch_file_on_disk = True

        self.error_line = None
        # colors
        self.background_color = DEFAULT_COLOR_SCHEME.get('background')
        self.sidearea_color = DEFAULT_COLOR_SCHEME.get('sideareas')
        self.currentline_color = DEFAULT_COLOR_SCHEME.get('currentline')

        # signals
        self.currentCellChanged.connect(self.cursor_moved)

        # file system watcher
        self.file_watcher = QFileSystemWatcher()
        self.file_watcher.fileChanged.connect(self.file_changed_on_disk)
        self.file_watcher.directoryChanged.connect(self.dir_changed)

    def cursor_moved(self, row, col, prev_row, prev_col):
        self.cursor_changed.emit(row+1, col+1)

    def modified(self):
        self.unsaved_flag = True
        if not self.loading:
            self.needsSaved.emit()


    def open(self, fname, read_only=False):
        self.loading = True
        self.file_name = fname

        if not self.load_file(fname):
            return
        self.setReadOnly(read_only)

        # add path to file watcher
        if os.path.exists(fname):
            self.file_watcher.addPath(fname)
        self.file_watcher.addPath(os.path.dirname(fname))

        self.loading = False
        self.unsaved_flag = False

    def file_changed_on_disk(self, fname):
        if self.watch_file_on_disk:
            self.changed_on_disk.emit()

    def dir_changed(self, dirname):
        if self.file_name and os.path.exists(self.file_name):
            self.file_watcher.addPath(self.file_name)

    def error(self, msg):
        if self.mfixgui:
            self.mfixgui.error(msg, popup=True)
        else:
            print(msg)

    def warning(self, msg):
        if self.mfixgui:
            self.mfixgui.warning(msg, popup=True)
        else:
            print(msg)

    def save(self, fname=None):
        if fname is None:
            fname = self.file_name
        if fname is None:
            fname = QFileDialog.getSaveFileName(self)[0]
            if not fname:
                return False
            self.file_name = fname

        # update file watcher
        watch_paths = self.file_watcher.files()
        if watch_paths:
            self.file_watcher.removePaths(watch_paths)

        try:
            self.save_file(fname)
        except Exception as e:
            self.error("Error saving file: %s" % e)
            return False

        self.saved.emit()
        self.unsaved_flag = False

        # update file watcher
        if os.path.exists(fname):
            self.file_watcher.addPath(fname)
        self.file_watcher.addPath(os.path.dirname(fname))
        return True
