# -*- coding: utf-8 -*-
"""
Dialog to build the mfixsolver from the GUI

"""

import os
from os.path import dirname, join
import shutil
import sys

from qtpy.QtCore import QProcess, QProcessEnvironment, QTimer, QSize, Qt, QUrl
from qtpy.QtGui import QColor, QFont, QFontMetrics, QMouseEvent
from qtpy.QtWidgets import QApplication, QDialog

from mfixgui.tools import plural, safe_int
from mfixgui.tools.qt import get_ui, SETTINGS
import mfixgui.build_mfixsolver


WINDOWS = sys.platform.startswith('win')

# Build types
RELWITHDEBINFO, DEBUG, CUSTOM = range(3)
BuildTypes = ['RELWITHDEBINFO', 'DEBUG', 'CUSTOM']
BuildTypeNames = ["Release with debug information",
                  "Debug",
                  "Custom (specify compiler flags)"]

PYMFIX, BATCH = range(2)
ServerTypes = ['PYMFIX', 'BATCH']
ServerTypeNames = ["Enabled",
                   "None (batch)"]


class BuildPopup(QDialog):
    """ Dialog box for running build_mfixsolver to build the solver """

    def __init__(self, parent):
        QDialog.__init__(self, parent)
        gui = self.gui = parent
        self.cwd = os.getcwd()
        self.build_log = None
        ui = self.ui = get_ui("build_popup.ui", self)
        ui.sizeHint = self.sizeHint
        tb = ui.textbrowser
        le = ui.lineedit_build_command
        m = QFontMetrics(tb.font())
        ui.min_width = m.width(" " * 120)
        font_height = m.lineSpacing() # not .height!
        tb.setMinimumSize(QSize(ui.min_width, 20*font_height + 10))
        le.setMinimumSize(QSize(m.width(" "*40), font_height+10))
        self.buffer = []
        tb.contextMenuEvent = self.handle_context_menu
        tb.anchorClicked.connect(self.handle_anchor)
        ui.button_first_error.clicked.connect(self.goto_first_error)
        ui.button_next_error.clicked.connect(self.goto_next_error)
        ui.button_prev_error.clicked.connect(self.goto_prev_error)
        self.finished.connect(gui.handle_compile_finished) # NB compile_finished != finished_building  FIXME
        self.build_proc = None #QProcess()
        #ui.setModal(True)
        ui.setWindowTitle("Build solver")
        ui.combobox_server.addItems(ServerTypeNames)
        ui.combobox_buildtype.addItems(BuildTypeNames)
        self.timer = QTimer()
        self.timer.timeout.connect(self.flush)
        self.timer.setInterval(250) # 4 Hz
        #ui.adjustSize()
        self.init_handlers()

        #desktop = QApplication.instance().desktop()
        #screen = desktop.screenNumber(desktop.cursor().pos())
        #geo = self.ui.frameGeometry()
        #geo.moveCenter(desktop.screenGeometry(screen).center())
        #self.ui.move(geo.center())


    def popup(self):
        # only show server option in dev mode
        self.gui.ui.toolbutton_compile.setEnabled(False)
        self.ui.progressbar.setRange(0,100)
        self.ui.progressbar.setValue(0)
        dev_mode = int(SETTINGS.value("developer_mode", 0))
        self.init_visibility(dev_mode)
        if self.gui.project:
            self.load_settings(self.gui.project.mfix_gui_comments)

        self.set_output_visible(False)
        self.update_build_cmd()
        self.set_build_enabled(True)
        self.update_compiler_box()
        self.ui.checkbox_dmp.stateChanged.connect(self.update_compiler_box)

        self.cwd = os.getcwd()
        self.show()
        self.raise_()
        self.activateWindow()
        #desktop = QApplication.instance().desktop()
        #screen = desktop.screenNumber(desktop.cursor().pos())
        #geo = self.ui.frameGeometry()
        #geo.moveCenter(desktop.screenGeometry(screen).center())
        #self.ui.move(geo.center())

    def sizeHint(self):
        # calling 'geometry' causes subsequent height and width to be correct (?)
        geo = self.geometry()
        ui = self.ui
        height = ui.height()
        width = ui.width()
        show_output = ui.textbrowser.isVisible()

        if not show_output:
            height -= (ui.textbrowser.height() +
                       ui.layout_error_buttons.geometry().height() +
                       10) # spacing + margin

            # would be nicer to put these things into a common frame that
            # we could show/hide
        return QSize(width, height)



    def handle_context_menu(self, event):
        # Avoid exposing internal urls
        menu = self.ui.textbrowser.createStandardContextMenu(event.globalPos())
        for ac in menu.actions():
            if ac.text() == 'Copy &Link Location':
                menu.removeAction(ac)
        menu.exec_(event.globalPos())


    def goto_first_error(self):
        self.goto_error(1)

    def goto_prev_error(self):
        self.goto_error(self.current_error-1)

    def goto_next_error(self):
        self.goto_error(self.current_error+1)

    def update_compiler_box(self):
        self.ui.combobox_compiler.clear()
        compilers = (
            [] if self.ui.checkbox_dmp.isChecked() else ["gfortran", "ifort"]
        ) + ["mpifort", "mpiifort", "mpif90"]
        self.ui.combobox_compiler.addItems([fc for fc in compilers if shutil.which(fc)])


    def intercept_close(self):
        self.timer.stop()
        if self.build_proc and self.build_proc.state() > 0:
            self.kill_build()
        if self.gui.project:
            self.save_build_settings(self.gui.project.mfix_gui_comments)
        self.close()

    def handle_anchor(self, url):
        url = url.url() # it's a QUrl
        n = int(url.split('#error_',1)[1])
        self.current_error = n
        self.goto_error(n)


    def goto_error(self, n):
        if n < 1:
            n = 1
        elif n > self.error_count:
            n = self.error_count
        self.current_error = n
        self.ui.textbrowser.setSource(QUrl("#error_%s"%n))
        fname, line, column = self.errors[n-1]
        tab = self.gui.editor_widget.open(fname)
        tab.editor.goto(line, column, error=True)
        self.gui.change_mode('editor')
        self.update_error_buttons()
        #self.gui.raise_()


    def cancel(self):
        """ if build is running, stop build
            if build is not running, close dialog box """
        self.timer.stop()
        self.flush()
        if self.build_proc and self.build_proc.state() > 0:
            self.kill_build()
            self.set_build_enabled(True)
            self.ui.label_progress.setText("Build canceled.")
        else:
            if self.gui.project:
                self.save_build_settings(self.gui.project.mfix_gui_comments)
            self.ui.close()

    def kill_build(self):
        self.build_proc.canceling = True
        if self.build_proc and self.build_proc.state() > 0:
            self.build_proc.kill()
            self.build_proc.waitForFinished(300)
        self.timer.stop()
        self.flush()

    def toggle_output(self, show=False):
        """ hide or show the build output TextBrowser """
        if not self.ui.textbrowser.isVisible() or show:
            self.set_output_visible(True)
            self.scroll_to_end()
        else:
            self.set_output_visible(False)

    def scroll_to_end(self):
        sb = self.ui.textbrowser.verticalScrollBar()
        sb.setValue(sb.maximum())

    def set_output_visible(self, show):
        ui = self.ui
        for b in (ui.button_first_error, ui.button_prev_error, ui.button_next_error):
            b.setVisible(show)
        ui.textbrowser.setVisible(show)
        if show:
            ui.pushbutton_show_out.setText("Hide build output")
        else:
            ui.pushbutton_show_out.setText("Show build output")

        ui.adjustSize()



    def update_build_cmd(self):
        ui = self.ui
        def quote_if_needed(x):
            return ('"' + x + '"' if ' ' in x else x)
        cmd = ' '.join(map(quote_if_needed, self.get_build_cmd()))
        ui.lineedit_build_command.setText(cmd)
        can_build = bool(ui.combobox_compiler.currentText()
                         or not ui.checkbox_dmp.isChecked())
        status = ('Press "Build solver" to compile.' if can_build
                  else "Specify compiler [wrapper] for DMP.")

        self.pushbutton_build.setEnabled(can_build)
        self.label_progress.setText(status)
        self.label_progress.setStyleSheet("color: None;")

        for w in (ui.label_flags, ui.lineedit_compiler_flags):
            w.setEnabled(self.ui.combobox_buildtype.currentIndex()==CUSTOM)

    def clean(self):
        ui = self.ui
        ui.textbrowser.clear()
        ui.label_progress.setText("Cleaning build directory...")
        ui.label_progress.setStyleSheet("color: None;")

        ui.progressbar.setValue(0)
        ui.progressbar.setRange(0, 0)

        self.print_to_output("Cleaning build directory: %s\n" % self.cwd, color='blue')
        self.set_build_enabled(False)

        removed_stuff = mfixgui.build_mfixsolver.do_clean()

        ui.progressbar.setValue(0)
        self.print_to_output(removed_stuff)
        self.flush()
        ui.label_progress.setText("Cleaning build directory... done")
        self.set_build_enabled(True)

        ui.progressbar.setValue(0)
        ui.progressbar.setRange(0, 100)

    def build(self):
        """ Start a build_mfixsolver process """
        ui = self.ui
        tb = ui.textbrowser
        tb.clear()

        self.timer.start()
        self.error_count = 0
        self.errors = []
        self.current_error = None
        self.update_error_buttons()
        self.scroll_to_end()
        ui.label_progress.setText("Checking requirements...")
        ui.label_progress.setStyleSheet("color: None;")
        ui.progressbar.setRange(0, 0)
        self.build_log = open("build.log", 'w')
        if self.gui.unsaved_flag:
            confirm = self.gui.message(text="Save project before building?",
                                   title="Save?",
                                   icon='question',
                                   buttons=['yes', 'no'],
                                   default='no')
            if confirm == 'yes':
                self.gui.save_project()

        if self.gui.editor_widget.check_unsaved():
            confirm = self.gui.message(text="Save editor files before building?",
                                   title="Save?",
                                   icon='question',
                                   buttons=['yes', 'no'],
                                   default='no')
            if confirm == 'yes':
                self.gui.save_editor_files()

        self.build_proc = QProcess()
        self.build_proc.canceling = False
        self.build_proc.error.connect(self.error)
        self.build_proc.readyReadStandardOutput.connect(self.read_out)
        self.build_proc.readyReadStandardError.connect(self.read_err)
        self.build_proc.finished.connect(self.finished_building)
        self.build_proc.setWorkingDirectory(str(self.cwd))
        cmd = self.get_build_cmd()
        self.print_to_output("Running %s\n" % ' '.join(cmd), color='blue')
        self.flush()
        # Set up environment
        env = QProcessEnvironment.systemEnvironment()
        env.insert("PYTHONUNBUFFERED", "y")
        env.remove("TERM")
        self.build_proc.setProcessEnvironment(env)

        self.set_build_enabled(False)
        self.build_proc.start(cmd[0], cmd[1:])


    def update_error_buttons(self):
        ui = self.ui
        if self.error_count == 0:
            for b in (ui.button_first_error, ui.button_next_error, ui.button_prev_error):
                b.setEnabled(False)
        else:
            ui.button_first_error.setEnabled(True)
            ui.button_next_error.setEnabled(self.current_error is not None
                                            and self.current_error < self.error_count)
            ui.button_prev_error.setEnabled(self.current_error is not None
                                            and self.current_error > 1)



    def update_gui_comment(self, key, val):
        if not self.gui.project:
            return
        comments = self.gui.project.mfix_gui_comments
        if comments.get(key) != val:
            comments[key] = val
            self.gui.set_unsaved_flag()

    def finished_building(self, exit_code, exit_status):
        self.timer.stop()
        self.flush()
        if self.build_log:
            self.build_log.close()
            self.build_log = None
        if exit_code == 0 and exit_status == 0:
            self.ui.progressbar.setRange(0, 100)
            self.ui.progressbar.setValue(100)
            self.ui.label_progress.setText("Build succeeded.")
            self.ui.label_progress.setStyleSheet("color: blue;")
            exe = os.path.join(os.getcwd(), 'mfixsolver')
            if self.dmp:
                exe += '_dmp'
            if self.smp:
                exe += '_smp'
            if WINDOWS:
                exe += '.exe' if self.batch else '.bat'
            self.update_gui_comment("mfix_exe", exe)
            if self.gui.unsaved_flag:
                self.gui.save_project()
        else:
            if self.error_count > 0:
                msg = "Build failed with %s." % plural(self.error_count, 'error')
            else:
                msg = "Build failed."
            #self.ui.progressbar.reset()
            self.ui.progressbar.setRange(0,100)
            self.ui.label_progress.setText(msg)
            self.ui.label_progress.setStyleSheet("color: red;")

            if not self.build_proc.canceling:
                self.toggle_output(show=True)
        self.set_build_enabled(True)
        self.ui.pushbutton_show_out.setEnabled(True)

    def set_build_enabled(self, enabled):
        self.ui.pushbutton_cancel.setText("Close" if enabled else "Cancel")
        self.ui.pushbutton_build.setEnabled(enabled)
        self.ui.pushbutton_clean.setEnabled(enabled)
        self.ui.groupbox_compiler.setEnabled(enabled)
        self.ui.groupbox_options.setEnabled(enabled)

    def error(self, error):
        exit_status = self.build_proc.ExitStatus()
        exit_code = 0 if exit_status==0 else 1
        cmd = " ".join(self.get_build_cmd())
        if error == QProcess.FailedToStart:
            msg = "Process failed to start: " + cmd
        elif error == QProcess.Crashed:
            msg = "Process exit: " + cmd
        elif error == QProcess.Timedout:
            msg = "Process timeout: " + cmd
        elif error in (QProcess.WriteError, QProcess.ReadError):
            msg = "Process communication error " + cmd
        else:
            msg = "Unknown error: " + cmd

        self.ui.label_progress.setText("Process error")
        self.finished_building(exit_code, exit_status)
        self.print_to_output(msg, color='red')
        self.flush()
        if not self.build_proc.canceling:
            self.toggle_output(show=True)
        self.finished_building(1,1)


    def read_out(self):
        if self.build_proc is None:
            self.timer.stop()
            return
        output = str(self.build_proc.readAllStandardOutput(),
                     encoding='ascii', errors='ignore')
        if output:
            if self.build_log:
                self.build_log.write(output)
                self.build_log.flush()
            self.check_progress(output, is_err=False)

    def read_err(self):
        if self.build_proc is None:
            self.timer.stop()
            return
        output = str(self.build_proc.readAllStandardError(),
                     encoding='ascii', errors='ignore')
        if output:
            if self.build_log:
                self.build_log.write(output)
                self.build_log.flush()
            self.check_progress(output, is_err=True)

    def check_progress(self, output, is_err):
        for line in output.splitlines():
            self.check_progress1(line, is_err)

    def check_progress1(self, line, is_err):
        ui = self.ui
        color = 'black'
        err = None
        line = line.rstrip()
        ltmp = line
        drive = ''
        if is_err:
            color = 'red'
            # Look for filename:line:col:  or filename:line:
            #  This is slightly tricky because filename may contain : on Posix
            if line.endswith(':'):
                if WINDOWS: # Most Windows paths will have a drive letter, unless they are net mounts
                    if line[1:2] == ':':
                        drive = line[:2]
                        ltmp = line[2:]
                parts = ltmp.rsplit(':', 3)[:-1] # delete final empty item
                if parts[-1].isdigit(): # Looks good, ends with a number
                    if len(parts) == 3:
                        if parts[1].isdigit(): # line:column
                            fname = drive + parts[0]
                            lno = int(parts[1])
                            col = int(parts[2])
                            err = (fname, lno, col)
                        else: # colon was part of filename?
                            fname = drive + parts[0] + ':' + parts[1]
                            lno = int(parts[2])
                            err = (fname, lno, 0)
                    elif len(parts) == 2:
                        fname = drive + parts[0]
                        lno = int(parts[1])
                        err = (fname, lno, 0)
            if err:
                self.errors.append(err)
                text = '%s:%s:%s:' % err
                self.error_count += 1
                anchor = 'error_%s' % self.error_count
                self.print_to_output('<a href=#%s name=%s>%s</a>'
                                     %(anchor, anchor, text))
                line = '' # Nothing to print, it's been turned into the above link
                self.update_error_buttons()


        else:
            percent = None
            if line[:1] == '[' and line[4:6] == '%]':
                try:
                    percent = int(line[1:4])
                except ValueError:
                    pass
            if percent:
                ui.label_progress.setText("Compiling...")
                ui.progressbar.setRange(0, 100)
                ui.progressbar.setValue(percent)
            elif self.ui.progressbar.value() == 100:
                ui.label_progress.setText("Linking...")
                ui.progressbar.setRange(0, 0)
            if 'build succ' in line.lower():
                color = 'blue'
            elif 'build fail' in line.lower():
                color = 'red'
        #if not line.endswith('\n'):
        #    line += '\n'
        line = line.replace('&', '&amp;')
        line = line.replace('>', '&gt;')
        line = line.replace('<', '&lt;')
        self.print_to_output(line+'\n', color)


    def print_to_output(self, msg, color='black'):
        """ display message in the build output TextBrowser """

        if "<a href=" in msg:
            html = msg # Already html, don't change a thing
        else:
            msg = msg.replace('\r', '') # Windows output from build_mfixsolver
            msg = msg.replace('\n', '<br>')
            msg = msg.replace(' ', '&nbsp;')
            html = '<span style="color:%s">%s</span>' % (color, msg) # add color
        self.buffer.append(html)

    def flush(self):
        if not self.buffer:
            return
        ui = self.ui
        tb = ui.textbrowser

        cursor = tb.textCursor()
        cursor.movePosition(cursor.End)

        sb = tb.verticalScrollBar()
        at_end = sb.value() == sb.maximum()
        cursor.insertHtml(''.join(self.buffer))
        self.buffer.clear()
        if at_end:
            self.scroll_to_end()


    def init_visibility(self, dev_mode):
        self.ui.textbrowser.setVisible(False)
        self.ui.checkbox_dmp.setVisible(not WINDOWS)
        if not WINDOWS:
            self.ui.checkbox_dmp.setChecked(False)
        devmode_only = [
            self.ui.combobox_server,
            self.ui.label_server,
        ]
        for widget in devmode_only:
            widget.setVisible(dev_mode)

    def init_handlers(self):
        ui = self.ui
        for sig in (ui.checkbox_dmp.toggled,
                    ui.checkbox_parallel.toggled,
                    ui.checkbox_smp.toggled,
                    ui.checkbox_verbose.toggled,
                    ui.combobox_buildtype.currentIndexChanged,
                    ui.combobox_compiler.currentIndexChanged,
                    ui.combobox_compiler.editTextChanged,
                    ui.combobox_server.currentIndexChanged,
                    ui.combobox_server.editTextChanged,
                    ui.lineedit_compiler_flags.textChanged):
            sig.connect(self.update_build_cmd)

        ui.finished.connect(self.intercept_close)
        ui.pushbutton_build.clicked.connect(self.build)
        ui.pushbutton_cancel.clicked.connect(self.cancel)
        ui.pushbutton_clean.clicked.connect(self.clean)
        ui.pushbutton_show_out.clicked.connect(self.toggle_output)


    def get_build_cmd(self):
        """ return command argument list for building the solver """
        return ["build_mfixsolver"] + self.get_args()

    def get_args(self):
        args = []
        idx = self.ui.combobox_buildtype.currentIndex()
        if idx == DEBUG:
            args += ["-DCMAKE_BUILD_TYPE=" + BuildTypes[idx]]
        elif idx == CUSTOM:
            args += ["-DCMAKE_Fortran_FLAGS=" + self.ui.lineedit_compiler_flags.text()]

        if self.ui.checkbox_parallel.isChecked():
            args += ["-j"]

        if self.ui.checkbox_smp.isChecked():
            args += ["--smp"]
            self.smp = True
        else:
            self.smp = False

        compiler = self.ui.combobox_compiler.currentText()
        if self.ui.checkbox_dmp.isChecked():
            self.dmp = True
            args += ["--dmp"]
        else:
            self.dmp = False

        if compiler:
            args += ["-DCMAKE_Fortran_COMPILER=" + compiler]

        self.batch = False
        if int(SETTINGS.value("developer_mode", 0)):
            if self.ui.combobox_server.currentIndex() == BATCH:
                args += ["--batch"]
                self.batch = True

        if self.ui.checkbox_verbose.isChecked():
            args += ["--verbose"]

        return args

    def load_settings(self, comments):
        ui = self.ui
        ui.checkbox_dmp.setChecked(safe_int(comments.get("BUILD_DMP", "0")) and not WINDOWS)
        ui.checkbox_parallel.setChecked(safe_int(comments.get("BUILD_PARALLEL", "1")))
        ui.checkbox_smp.setChecked(safe_int(comments.get("BUILD_SMP", "0")))
        ui.checkbox_verbose.setChecked(safe_int(comments.get("BUILD_VERBOSE", "0")))
        ui.combobox_compiler.setEditText(comments.get("BUILD_FC", ""))
        ui.lineedit_compiler_flags.setText(comments.get("BUILD_FC_FLAGS", ""))

        val = comments.get("BUILD_TYPE", "RELWITHDEBINFO")
        if val in BuildTypes:
            self.ui.combobox_buildtype.setCurrentIndex(BuildTypes.index(val))

        val = comments.get("BUILD_INTERACTIVE", "PYMFIX")
        if val in ServerTypes:
            self.ui.combobox_server.setCurrentIndex(ServerTypes.index(val))

    def save_build_settings(self, comments):
        """ Save Build dialog settings in gui_comments """
        ui = self.ui
        def update(key, val):
            old = comments.get(key)
            if isinstance(val, int):
                if old is not None:
                    old = int(old)
            if val != old:
                comments[key] = val
                self.gui.set_unsaved_flag()
        for (key, val) in (("BUILD_DMP", int((ui.checkbox_dmp.isChecked()) and not WINDOWS)),
                           ("BUILD_INTERACTIVE", ServerTypes[ui.combobox_server.currentIndex()]),
                           ("BUILD_PARALLEL", int(ui.checkbox_parallel.isChecked())),
                           ("BUILD_SMP", int(ui.checkbox_smp.isChecked())),
                           ("BUILD_TYPE", BuildTypes[ui.combobox_buildtype.currentIndex()])):
            update(key, val)
        flags = ui.lineedit_compiler_flags.text().strip()
        if flags:
            update("BUILD_FC_FLAGS", flags)
        else:
            if comments.pop("BUILD_FC_FLAGS", None):
                self.gui.set_unsaved_flag()
        if self.gui.unsaved_flag:
            self.gui.save_project()
