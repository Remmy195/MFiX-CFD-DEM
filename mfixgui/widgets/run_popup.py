# TODO clean up the messy and confusing hiding of widgets

import os
import shutil
import json
import multiprocessing
import subprocess

from os.path import join
from glob import glob

from qtpy.QtCore import Qt, QObject, Signal

from qtpy.QtWidgets import (QComboBox, QFileDialog, QLabel,
                            QSpinBox, QListWidget, QLineEdit)

from mfixgui.tools import plural, safe_int
from mfixgui.tools.qt import (get_icon, sub_icon_size,
                              get_ui, widget_iter)

from mfixgui.template_manager import (QueueTemplateManager, QueueTemplate,
                                      init_template_manager)

from mfixgui.tools.qt import SETTINGS
from mfixgui.solver.manager import SolverManager, GOOD, PENDING, ERROR


SPX_GLOB = ["*.sp?"]
RESTART_2_REMOVE_GLOB = ["*.sp?", "*.pvd", "*.vtp", "*.vtu"]

def open_run_popup(mfixgui):
    mfixgui.chemistry_chemkin_warning()
    if not check_custom_solver_if_udfs(mfixgui): #rename
        return
    # Comments
    # 1) Why are we creating a new RunPopup each time, and why are we sticking it in
    # the parent's namespace?
    # 2) We need better way of getting hold of global singleton mfixgui.
    # No need to be passing it as a parameter (twice here even)
    # 3) is it a Dialog or Popup?
    mfixgui.run_popup = RunPopup(mfixgui, mfixgui)
    mfixgui.run_popup.popup()


def check_custom_solver_if_udfs(mfixgui): # rename
    project_dir = mfixgui.get_project_dir()
    mfixsolver = glob(join(project_dir, "mfixsolver*"))
    if mfixsolver:
        return True

    project_dir = mfixgui.get_project_dir()
    if not list(glob(join(project_dir, "*.f"))): # TODO add .f90
        # TODO check solver newer than sources
        return True
    return confirm_ignore_udfs(mfixgui,
                               (
                                   "Fortran source files exist for this project, but "
                                   "there is no custom solver in the project directory. "
                                   "This case will not run correctly unless this "
                                   "project's custom mfixsolver is selected. Proceed anyway?"
                               ))


def confirm_ignore_udfs(mfixgui, udf_msg):
    project_dir = mfixgui.get_project_dir()
    if not list(glob(join(project_dir, "*.f"))):
        return True

    response = mfixgui.message(
        title="Warning",
        icon="question",
        text=udf_msg,
        buttons=["yes", "no"],
        default="no",
    )
    return response == "yes"


n_cpus = multiprocessing.cpu_count()

class RunPopup(QObject):
    templates_updated = Signal()
    sig_solvers_updated = Signal()

    def __init__(self, parent, mfixgui):
        super(RunPopup, self).__init__()
        self.mfixgui = mfixgui
        self.gui_comments = self.mfixgui.project.mfix_gui_comments
        self.solver_manager = None
        self.template_manager = None

        # load ui
        ui = self.ui = get_ui("run_popup.ui")
        ui.setParent(parent, Qt.Dialog)
        ui.setModal(True)
        ui.closeEvent = self.closeEvent
        ui.layout.setSizeConstraint(ui.layout.SetFixedSize)

        ui.toolbutton_browse.clicked.connect(self.handle_browse_exe)
        ui.toolbutton_browse.setIcon(get_icon("add.svg"))
        ui.toolbutton_browse.setIconSize(sub_icon_size())

        ui.toolbutton_remove.clicked.connect(self.handle_remove_exe)
        ui.toolbutton_remove.setIcon(get_icon("remove.svg"))
        ui.toolbutton_remove.setIconSize(sub_icon_size())

        ui.toolbutton_view_error.clicked.connect(self.show_solver_error)

        ui.listwidget_solver_list.itemSelectionChanged.connect(self.slot_select_solver)

        ui.groupbox_queue.setChecked(False)
        ui.widget_queue.setVisible(False)
        ui.groupbox_queue.toggled.connect(self.show_hide_queue)

        ui.button_run.clicked.connect(self.handle_run)
        ui.button_cancel.clicked.connect(self.ui.close)

        self.smp = self.dmp = False
        txt = plural(n_cpus, "core")
        ui.groupbox_smp_options.setTitle("SMP options (%s available)" % txt)
        ui.groupbox_dmp_options.setTitle("DMP options (%s available)" % txt)
        ui.groupbox_queue.toggled.connect(self.toggle_run_btn_text)
        ui.combobox_template.setEnabled(True)

        ui.checkbox_keyword_bdist_io.clicked.connect(
            lambda x, mfixgui=mfixgui: mfixgui.update_keyword('bdist_io', x))
        mfixgui.add_tooltip(ui.checkbox_keyword_bdist_io,
                            key="bdist_io")
        mfixgui.add_tooltip(ui.spinbox_nodesi, key='nodesi')
        mfixgui.add_tooltip(ui.spinbox_nodesj, key='nodesj')
        mfixgui.add_tooltip(ui.spinbox_nodesk, key='nodesk')
        mfixgui.add_tooltip(ui.lineedit_mpirun_flags, key=None,
                            description="Additional flags for mpirun command.")
        ui.checkbox_use_hwthread.setToolTip("Add <b>‑‑use‑hwthread‑cpus</b> flag to mpirun command (if supported).")
        ui.checkbox_use_hwthread.tooltip0 = ui.checkbox_use_hwthread.toolTip()
        self.init_ui()
        self.init_templates()

        self.update_dialog()
        self.ui.listwidget_solver_list.setCurrentRow(0)
        # TODO: simplify
        self.sig_solvers_updated.connect(self.update_solver_listbox)
        self.solver_manager = SolverManager(self.mfixgui, self.sig_solvers_updated)

    def slot_select_solver(self):
        self.update_dialog()
        self.save_selected_exe()

    @property
    def solver(self):
        """The currently selected solver"""
        item = self.ui.listwidget_solver_list.currentItem()
        if item is None:
            return None

        return self.solver_manager.solver_dict.get(item.text())

    def project_dir(self):
        projdir = self.mfixgui.get_project_dir()
        if projdir is None:
            return None
        return projdir

    def toggle_run_btn_text(self):
        txt = "Submit" if self.ui.groupbox_queue.isChecked() else "Run"
        self.ui.button_run.setText(txt)

    def update_gui_comment(self, key, val):
        if str(self.gui_comments.get(key)) != str(val):
            self.gui_comments[key] = val
            self.mfixgui.set_unsaved_flag()

    # UI update functions
    def init_ui(self):
        self.ui.setWindowTitle("Run MFiX solver")
        self.init_restart()
        self.init_smp()
        self.init_dmp()
        self.init_solver_list()
        self.ui.groupbox_queue.setChecked(
            bool(int(self.gui_comments.get("submit_to_queue", False))))
        self.ui.checkbox_keyword_bdist_io.setChecked(
            bool(int(self.mfixgui.project.get_value('bdist_io', False))))

    def init_restart(self):
        ui = self.ui
        spx_files = self.mfixgui.get_output_files(SPX_GLOB)
        res_files = self.mfixgui.get_res_files()

        # Should we select restart_1 based on run_type in file, or
        #  tstop being extended?
        b = ui.button_restart_0
        b.setChecked(True)

        b = ui.button_restart_1
        b.setEnabled(bool(spx_files and res_files))
        if not b.isEnabled():
            b.setToolTip("No restart files found")
            b.setChecked(False)
        else:
            b.setToolTip("")

        b = ui.button_restart_2
        b.setEnabled(bool(res_files))
        if not b.isEnabled():
            b.setToolTip("No restart files found")
            b.setChecked(False)
        else:
            b.setToolTip("")


    def init_smp(self):
        project_threads = self.gui_comments.get("OMP_NUM_THREADS", "1")
        n_threads = os.environ.get("OMP_NUM_THREADS", project_threads)
        n_threads = safe_int(n_threads, default=1)
        sb = self.ui.spinbox_threads
        sb.setValue(n_threads)
        sb.valueChanged.connect(self.update_total_smp)
        self.update_total_smp(n_threads)

    def init_dmp(self):
        nodes = [(self.ui.spinbox_nodesi, "nodesi"),
                 (self.ui.spinbox_nodesj, "nodesj"),
                 (self.ui.spinbox_nodesk, "nodesk")]
        for (sb, kw) in nodes:
            val = self.mfixgui.project.get_value(kw,
                      default=self.mfixgui.get_retained_keyword(kw, default=1))
            sb.setValue(val)
            sb.valueChanged.connect(self.update_total_dmp)
        self.update_total_dmp()
        mpirun_flags = self.gui_comments.get('mpirun_flags','').strip()
        use_hwthread = self.gui_comments.get('use_hwthread', 'True')
        use_hwthread = (use_hwthread=="True") # Stored value is str not bool
        if '--use-hwthread-cpus' in mpirun_flags:
            use_hwthread = True
            mpirun_flags = mpirun_flags.replace('--use-hwthread-cpus', '')
            mpirun_flags = mpirun_flags.replace('  ', ' ')
            mpirun_flags = mpirun_flags.strip()
        # Test whether '--use-hwthread-cpus' is valid
        arg = '--use-hwthread-cpus'
        status, output = subprocess.getstatusoutput('mpirun ' + arg)
        if 'unrecognized' in output.lower():
            supported = False
        else:
            supported = True
        self.ui.checkbox_use_hwthread.setEnabled(supported)

        if not supported:
            self.ui.checkbox_use_hwthread.setToolTip("This version of mpirun does not support the <b>‑‑use‑hwthread‑cpus</b> flag.")
            self.ui.checkbox_use_hwthread.tooltip0 = self.ui.checkbox_use_hwthread.toolTip()
            use_hwthread = False
        self.ui.checkbox_use_hwthread.setChecked(use_hwthread)
            # These were added to work around error messages we're not seeing any more
            #"-mca", "mpi_warn_on_fork", "0",  #https://www.open-mpi.org/faq/?category=tuning#fork-warning
            #"-mca", "mca_base_component_show_load_errors", "0", #https://github.com/open-mpi/ompi/issues/7752
        if mpirun_flags == 'None':
            mpirun_flags = ''
        if mpirun_flags:
            self.ui.lineedit_mpirun_flags.setText(mpirun_flags)

    def init_solver_list(self):
        self.update_solver_listbox()
        self.ui.listwidget_solver_list.setCurrentRow(0)

    def update_total_dmp(self):
        val =(self.ui.spinbox_nodesi.value()
            * self.ui.spinbox_nodesj.value()
            * self.ui.spinbox_nodesk.value())
        label = self.ui.label_total_dmp
        txt = str(val)
        if self.dmp and val > n_cpus and not self.ui.groupbox_queue.isChecked():
            txt += ' ⚠'
            label.setStyleSheet("color: red;")
            label.setToolTip("Number of cores requested exceeds cores available.")
        else:
            label.setStyleSheet("")
            label.setToolTip("")
        label.setText(txt)
        self.set_num_cores()

    def update_total_smp(self, val):
        sb = self.ui.spinbox_threads
        for w in widget_iter(sb):
            if isinstance(w, QLineEdit):
                w.setStyleSheet("color: red;" if self.smp and val > n_cpus else "")
        l = self.ui.label_smp_warning
        if self.smp and val > n_cpus:
            l.setText("⚠")
            l.setStyleSheet("color: red;")
            l.setToolTip("Number of threads requested exceeds available cores.")
        else:
            l.setText('')
            l.setStyleSheet('')
            l.setToolTip("")

    def init_templates(self):
        self.template_manager = QueueTemplateManager(self.templates_updated)
        init_template_manager(self.template_manager)

        self.templates_updated.connect(self.update_template_combobox)
        self.ui.browse_template_toolButton.clicked.connect(self._handle_browse_template)
        self.ui.browse_template_toolButton.setIcon(get_icon("add.svg"))
        self.ui.browse_template_toolButton.setIconSize(sub_icon_size())

        self.ui.combobox_template.currentIndexChanged.connect(self.select_template)

        self.ui.delete_template_toolButton.clicked.connect(
            lambda: self.template_manager.remove( self.ui.combobox_template.currentText()))

        self.ui.delete_template_toolButton.setIcon(get_icon("remove.svg"))
        self.ui.delete_template_toolButton.setIconSize(sub_icon_size())

        curr_template = "Joule"
        template_values = None
        temp = self.gui_comments.get("queue_template")
        if temp:
            template_values = json.loads(temp)
            t_name = template_values.get("template")
            if t_name:
                curr_template = t_name
        self.update_template_combobox(curr_template)
        if template_values and self.current_template():
            self.current_template().init_values(template_values)

    def current_template(self):
        ct = self.get_template(self.ui.combobox_template.currentText())
        return ct


    def get_template(self, template_name):
        if template_name in self.template_manager.template_keys():
            return self.template_manager[template_name]
        return None

    def show_hide_queue(self, val):
        self.ui.widget_queue.setVisible(val)
        self.update_total_dmp() # update warning

    def update_template_combobox(self, name="Joule"):
        """" Update when template is added or removed """
        self.ui.combobox_template.clear()
        self.ui.combobox_template.addItems(self.template_manager.template_keys())
        cb = self.ui.combobox_template
        for itm in range(cb.count()):
            if str(name).lower() == str(cb.itemText(itm)).lower():
                cb.setCurrentIndex(itm)
                break
        self.update_queue_widgets()

    def select_template(self):
        self.ui.combobox_template.setToolTip(self.ui.combobox_template.currentText())
        self.update_queue_widgets()

    def update_queue_widgets(self):
        """" Update widgets when selected template changes """
        queue_template = self.current_template()

        if queue_template is None:
            return
        self.ui.delete_template_toolButton.setEnabled(not queue_template.is_builtin)

        if not self._check_submit_command_exists(queue_template):
            return

        layout = self.ui.groupbox_queue_options_gridlayout

        # remove previous queue widgets
        while layout.count():
            item = layout.takeAt(0)
            if item is None:
                continue
            widget = item.widget()
            if widget is None:
                pass #?
                #layout.removeItem(item)
            else:
                # layout.removeWidget(widget)?
                widget.setVisible(False)

        # add the selected queue widgets
        layout.addWidget(queue_template.widget, 0, 0)
        queue_template.widget.setVisible(True)
        self.set_job_name()
        self.set_num_cores()

    def set_job_name(self):
        if not self.template_manager:
            return
        queue_template = self.current_template()
        if not queue_template:
            return
        for w in widget_iter(queue_template.widget):
            if w.objectName() == 'JOB_NAME':
                w.setText(self.mfixgui.project.get_value(
                    "run_name", "${PROJECT_NAME}"))

    def set_num_cores(self):
        if not self.template_manager:
            return
        queue_template = self.current_template()
        if not queue_template:
            return
        if self.dmp:
            nodesi = self.ui.spinbox_nodesi.value()
            nodesj = self.ui.spinbox_nodesj.value()
            nodesk = self.ui.spinbox_nodesk.value()
            cores = nodesi*nodesj*nodesk
        elif self.smp:
            cores = self.ui.spinbox_threads.value()
        else:
            cores = 1
        for w in widget_iter(queue_template.widget):
            if w.objectName() == 'CORES':
                w.setValue(cores)
            if w.objectName() == 'JOB_NAME':
                w.setText(self.mfixgui.project.get_value(
                    "run_name", "${PROJECT_NAME}"))



    def _check_submit_command_exists(self, queue_template):
        """ Display a warning and return False if the submit command doesn't exist, otherwise return True """

        cmd = queue_template.submit().split()[0]
        if shutil.which(cmd):
            self.ui.button_run.setEnabled(True)
            return True

        label = QLabel(
            'The submission command "{}" does not exist in '
            "the current environment. Please select another "
            "template, edit the template, and/or check your "
            "environment.".format(cmd)
        )
        label.setStyleSheet("color:red")
        label.setWordWrap(True)
        layout =  self.ui.groupbox_queue_options.layout()

        while layout.count():
            item = layout.takeAt(0)

            if item is None:
                continue
            widget = item.widget()
            if widget is None:
                pass #?
                #layout.removeItem(item)
            else:
                # layout.removeWidget(widget)?
                widget.setVisible(False)

        layout.addWidget(label, 0, 0)
        label.setVisible(True)
        self.ui.button_run.setEnabled(False)

        return False


    def update_dialog(self):
        """ Enable or disable options based on self.solver features,
        local or remote settings """

        ui = self.ui

        if self.solver is None:
            ui.button_run.setEnabled(False)
            ui.groupbox_queue.setEnabled(False)
            ui.groupbox_smp_options.setEnabled(False)
            ui.groupbox_dmp_options.setEnabled(False)
            ui.label_solver_error.setText("")
            ui.toolbutton_view_error.setVisible(False)
            ui.toolbutton_remove.setEnabled(False)
            return

        ui.toolbutton_remove.setEnabled(True)
        error = self.solver.status == ERROR
        good = self.solver.status == GOOD
        self.smp = smp = self.solver.smp_enabled()
        self.dmp = dmp = self.solver.dmp_enabled()
        no_k = self.mfixgui.project.get_value("no_k")
        ui.toolbutton_view_error.setVisible(error)
        ui.button_run.setEnabled(good)
        ui.groupbox_queue.setEnabled(good)

        no_smp_msg = "Selected solver does not support SMP"
        gb = ui.groupbox_smp_options
        gb.setEnabled(smp)
        gb.setToolTip(None if smp else no_smp_msg)
        self.update_total_smp(ui.spinbox_threads.value()) # update warning

        no_dmp_msg = "Selected solver does not support DMP"
        gb = ui.groupbox_dmp_options
        gb.setEnabled(dmp)
        gb.setToolTip(None if dmp else no_dmp_msg)

        for w in (ui.spinbox_nodesi,
                  ui.spinbox_nodesj,
                  ui.checkbox_keyword_bdist_io,
                  ui.checkbox_use_hwthread,
                  ui.lineedit_mpirun_flags):
            w.setToolTip(w.tooltip0 if dmp else no_dmp_msg)
        sb = ui.spinbox_nodesk
        sb.setEnabled(dmp and not no_k)
        sb.setToolTip(sb.tooltip0 if dmp and not no_k
                      else no_dmp_msg if not dmp
                      else "Simulation is 2-dimensional")
        self.update_total_dmp() # update warning


        status_text = self.solver.get_status_text()
        python = self.solver.python_enabled()
        if good and not python:
            status_text = "Cannot control selected solver (not built with interactive support)"
        ui.label_solver_error.setText(status_text)

    def show_solver_error(self):
        error_msg = self.solver.get_error()
        self.mfixgui.message(text="The solver test failed with the following error:",
                             info_text=error_msg)

    def popup(self):
        self.ui.show()
        self.ui.raise_()
        self.ui.activateWindow()

    def closeEvent(self, _event):
        """save information on close"""
        # save solver list
        self.save_selected_exe()
        self.update_gui_comment("mpirun_flags", self.ui.lineedit_mpirun_flags.text())

        # queue
        self.save_template()
        self.update_gui_comment("submit_to_queue", int(self.ui.groupbox_queue.isChecked()))
        self.update_gui_comment("OMP_NUM_THREADS", str(self.ui.spinbox_threads.value()))

    def save_template(self):
        """Save the current template data"""
        template_txt = self.ui.combobox_template.currentText()
        current_template = self.current_template()
        if current_template is None:
            return
        widget_values = self.current_template().widget_values()
        widget_values["template"] = template_txt
        self.update_gui_comment("queue_template", json.dumps(widget_values))

    def handle_run(self):
        if self.solver is None:
            return
        if self.finish_with_dialog():
            use_queue = self.ui.groupbox_queue.isChecked()
            template = self.current_template() if use_queue else None
            omp_num_threads = self.ui.spinbox_threads.value()

            # collect nodes[ijk] from project to guarantee that mpirun matches
            nodesi = self.mfixgui.project.get_value("nodesi", 1)
            nodesj = self.mfixgui.project.get_value("nodesj", 1)
            nodesk = self.mfixgui.project.get_value("nodesk", 1)
            mpirun_flags = self.ui.lineedit_mpirun_flags.text()
            if self.ui.checkbox_use_hwthread.isChecked() and '--use-hwthread-cpus' not in mpirun_flags:
                mpirun_flags = '--use-hwthread-cpus ' + mpirun_flags
            self.mfixgui.process_manager.start_solver(self.solver,
                                                      template,
                                                      omp_num_threads,
                                                      (nodesi, nodesj, nodesk),
                                                      mpirun_flags=mpirun_flags)

        self.mfixgui.slot_update_runbuttons()

    def finish_with_dialog(self):
        """ save run options in project file, then emit run signal """
        ui = self.ui
        if not self.check_custom_solver_selected():
            return False

        self.template_manager.save_settings()
        self.save_template()

        if ui.button_restart_0.isChecked():
            if not self.confirm_new_run():
                return False
        elif ui.button_restart_1.isChecked():
            if not self.confirm_restart_1():
                return False
        elif ui.button_restart_2.isChecked():
            if not self.confirm_restart_2():
                return False
        else:
            return False # One of the above must be selected
        if self.smp:
            thread_count = str(self.ui.spinbox_threads.value())
            self.update_gui_comment("OMP_NUM_THREADS", thread_count)
        if self.dmp:
            mpirun_flags =  self.ui.lineedit_mpirun_flags.text().strip()
            self.update_gui_comment("mpirun_flags", mpirun_flags)
            use_hwthread = self.ui.checkbox_use_hwthread.isChecked()
            self.update_gui_comment("use_hwthread", use_hwthread)
            self.save_dmp_keywords()
        else:
            nodesi = self.mfixgui.project.get_value("nodesi", 1)
            nodesj = self.mfixgui.project.get_value("nodesj", 1)
            nodesk = self.mfixgui.project.get_value("nodesk", 1)
            if nodesi != 1:
                self.mfixgui.retain_keyword("nodesi")
            if nodesj != 1:
                self.mfixgui.retain_keyword("nodesj")
            if nodesk != 1:
                self.mfixgui.retain_keyword("nodesk")
            self.mfixgui.unset_keyword("nodesi")
            self.mfixgui.unset_keyword("nodesj")
            self.mfixgui.unset_keyword("nodesk")


        if self.mfixgui.unsaved_flag:
            # run_type keyword updated and/or nodesi/nodesj/nodesk
            self.mfixgui.save_project()
        else:
            stl = join(self.project_dir(), "geometry.stl")  # is this needed?
            self.mfixgui.vtkwidget.export_stl(str(stl))

        self.ui.close()
        self.mfixgui.signal_update_runbuttons.emit("") # Why?
        return True


    def confirm_new_run(self):
        self.mfixgui.update_keyword("run_type", "new")
        output_files = self.mfixgui.get_output_files()
        if not output_files:
            return True

        confirm_delete = self.mfixgui.remove_output_files(
            output_files,
            message_text="Starting a new run requires deleting the following files from the run directory:",
            force_remove=True,
        )
        if confirm_delete:
            return True

        return False

    def confirm_restart_1(self):
        self.mfixgui.update_keyword("run_type", "restart_1")
        return True

    def confirm_restart_2(self):
        self.mfixgui.update_keyword("run_type", "restart_2")
        spx_files = self.mfixgui.get_output_files(RESTART_2_REMOVE_GLOB)
        if self.mfixgui.remove_output_files(spx_files, force_remove=True):
            return True
        return False


    def save_dmp_keywords(self):
        # collect nodes[ijk]
        dmp = self.solver is not None and self.solver.dmp_enabled()
        nodesi = self.ui.spinbox_nodesi.value() if dmp else 1
        nodesj = self.ui.spinbox_nodesj.value() if dmp else 1
        nodesk = self.ui.spinbox_nodesk.value() if dmp else 1

        no_k = self.mfixgui.project.get_value("no_k")

        # write the correct nodes[ijk] to project file
        self.mfixgui.update_keyword("nodesi", nodesi)
        self.mfixgui.update_keyword("nodesj", nodesj)
        self.mfixgui.update_keyword("nodesk", 1 if no_k else nodesk)

    def handle_remove_exe(self):
        item = self.ui.listwidget_solver_list.currentItem()
        if item:
            self.solver_manager.remove(item.text())

    def handle_browse_exe(self):
        """ Handle file open dialog for user specified exe """
        new_exe, _ = QFileDialog.getOpenFileName(
            self.ui,
            "Select executable",
            directory=str(self.project_dir()),
            options=QFileDialog.DontResolveSymlinks,
        )

        if not new_exe:
            return

        self.add_solver(new_exe)

    def add_solver(self, new_exe):
        """ Add the string representing a new solver to the SolverManager"""
        if not os.path.exists(new_exe):
            return
        display_name = self.solver_manager.display_name(new_exe)
        if display_name in self.solver_manager.solver_dict:
            self.mfixgui.message(text="The selected solver is already in the list of available solvers.")

        try:
            solver_path = new_exe
            self.solver_manager.add(solver_path)
            self.update_solver_listbox()
            self.ui.listwidget_solver_list.setCurrentRow(0)
        except ValueError as val_err:
            self.mfixgui.message(text=val_err)

    def update_solver_listbox(self):
        lw = self.ui.listwidget_solver_list
        selected_item = lw.currentItem()
        selected_item_text = selected_item.text() if selected_item else ""
        lw.clear()
        if self.solver_manager:
            lw.addItems(list(self.solver_manager.solver_dict.keys()))
            for i, (name, solver) in enumerate(
                    self.solver_manager.solver_dict.items()):
                item = lw.item(i)
                item.setIcon(get_icon(solver.get_icon()))
                if name == selected_item_text:
                    lw.setCurrentItem(item)
        if lw.count() > 0 and lw.currentItem() is None:
            lw.setCurrentRow(0)
        self.update_dialog()

    def _handle_browse_template(self):
        new_temp, _ = QFileDialog.getOpenFileName(self.ui,
                                                  "Select a queue template file",
                                                  directory=str(self.project_dir()))
        if new_temp:
            try:
                template = QueueTemplate(new_temp)
                self.template_manager.add(template)
                self.update_template_combobox(str(template.display_name()))
            except ValueError as err:
                self.mfixgui.message(text=f"Invalid template: {err}")

    def save_selected_exe(self):
        """ add new executable to recent list, save in project file and config,
        send signal(s) """
        if self.solver is None:
            self.mfixgui.warn("No solver selected")
            return

        solver_path = str(self.solver.path)
        self.update_gui_comment("mfix_exe", solver_path)


    def check_custom_solver_selected(self):
        project_dir = self.mfixgui.get_project_dir()
        if self.solver.path.startswith(project_dir):
            os.path.relpath(self.solver.path, project_dir)
            return True

        return confirm_ignore_udfs(
            self.mfixgui,
            (
                "Fortran source files exist for this project, but "
                "the selected mfixsolver is not in the project directory. "
                "This case probably won't run correctly unless this "
                "project's custom mfixsolver is selected. Proceed anyway?"
            ),
        )
