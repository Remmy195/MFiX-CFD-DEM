import os
import sys
import subprocess
import shutil

from enum import Enum
from threading import Thread

from mfixgui.tools import SCRIPT_DIR
from mfixgui.tools.qt import SETTINGS

WINDOWS = sys.platform.startswith('win')

if WINDOWS:
    SOLVER_NAMES = ["mfixsolver" + smp + ext
                    for smp in ('', '_smp')
                    for ext in ('.bat', '.exe')]
else:
    # mfixsolver.sh is symlinked to mfixsolver on Linux
    SOLVER_NAMES = ["mfixsolver" + dmp + smp
                    for dmp in ('', '_dmp')
                    for smp in ('', '_smp')]

PENDING, GOOD, ERROR = range(3)

def is_executable(fname):
    r = fname and os.path.isfile(fname) and os.access(fname, os.X_OK)
    return r


class Solver:
    """ Represents a builtin or custom solver """

    def __init__(self, path, sig_solvers_updated):
        self.name = str(path) # Should not be needed unless we still have pathlib re
        self.path = path # TODO remove pathlib remnants, name == path
        self.flags = ""
        self.stderror = None
        self.status = PENDING
        self.sig_solvers_updated = sig_solvers_updated
        if sig_solvers_updated:
            self._thread = Thread(target=self.run_print_flags)
            self._thread.start()

    def run_print_flags(self):
        try:
            self.flags = subprocess.check_output(
                [os.path.abspath(self.name), "--print-flags"],
                stderr=subprocess.STDOUT).decode("utf-8", "ignore")
            self.status = GOOD
        except subprocess.CalledProcessError as proc_err:
            self.stderror = proc_err.output
            self.status = ERROR
        finally:
            if self.sig_solvers_updated:
                self.sig_solvers_updated.emit()

    def dmp_enabled(self):
        return "dmp" in self.flags

    def smp_enabled(self):
        return "smp" in self.flags

    def python_enabled(self):
        return "python" in self.flags

    def get_error(self):
        return self.stderror.decode("utf-8", "replace")

    def get_status(self): # This is so annoying
        return self.status

    def get_status_text(self):
        return ["Checking solver", "", "Solver error:"][self.status]

    def get_icon(self):
        return ["timelapse.svg", "check_outline.svg", "error_outline.svg"][self.status]


class SolverManager:
    """ Represent the set of solvers available.  After constructor runs, it has no solvers """

    def __init__(self, mfixgui, sig_solvers_updated):
        self.mfixgui = mfixgui
        self.project_dir = self.mfixgui.get_project_dir()
        self.solver_dict = {}
        self.sig_solvers_updated = sig_solvers_updated
        conda_prefix = os.getenv('CONDA_PREFIX')
        last_solver = mfixgui.project.mfix_gui_comments.get("mfix_exe")

        # Don't advertise default solvers from previous MFiX versions
        ## TODO embed mfix version in solver script as a comment
        ## or in binary as variable (--print-mfix-version flag)
        ## and suggest rebuilding solvers if not current
        if last_solver and is_executable(last_solver):
            if conda_prefix:
                if (last_solver.startswith(os.path.dirname(conda_prefix))
                    and not last_solver.startswith(conda_prefix)):
                    pass
                else:
                    self.add(last_solver)
            else:
                self.add(last_solver)
        if conda_prefix:
            if WINDOWS:
                d = conda_prefix + '/Scripts/'
            else:
                d = conda_prefix + '/bin/'
            for s in SOLVER_NAMES:
                if is_executable(d+s):
                    self.add(d+s)
        # If you don't install with Conda you don't have default solvers.
        # But let's check to see if somehow it got installed into $PATH (e.g. pip/wheel)
        for s in SOLVER_NAMES: # Current directory which is project dir
            x = self.project_dir + '/' + s
            if is_executable(x):
                self.add(x)

        for s in SOLVER_NAMES:
            if x:=shutil.which(s):
                self.add(x)


    def add(self, solver_path):
        """ Add a new solver if path is not already present"""
        if not os.path.exists(solver_path):
            return
        key = self.display_name(solver_path)
        if key in self.solver_dict:
            return
        else:
            solver = Solver(solver_path, self.sig_solvers_updated)
            self.solver_dict[key] = solver
            self.solvers_updated()


    def remove(self, key):
        """ Remove solver """
        if key in self.solver_dict:
            del self.solver_dict[key]
            self.solvers_updated()

    def solvers_updated(self):
        if self.sig_solvers_updated:
            self.sig_solvers_updated.emit()

    def display_name(self, solver_path):
        """ Given a solver_path, show the string to be displayed in Run dialog """
        real_path = os.path.abspath(solver_path)
        project_dir = os.path.abspath(self.project_dir)
        if real_path.startswith(project_dir):
            return  "[project]/" + os.path.basename(solver_path)
        default_dir = os.path.dirname(os.path.dirname(sys.executable))
        conda_prefix = os.getenv('CONDA_PREFIX')

        if (real_path.startswith(default_dir) or
            (conda_prefix and real_path.startswith(conda_prefix))):
            return "[default]/" +  os.path.basename(solver_path)

        return real_path

def get_default_solver():
    if conda_prefix := os.getenv('CONDA_PREFIX'):
        s = conda_prefix + ('/Scripts/mfixsolver.bat' if WINDOWS
                            else '/bin/mfixsolver')
        if is_executable(s):
            return s
    if s := shutil.which('mfixsolver.bat' if WINDOWS else 'mfixsolver'):
        return s
