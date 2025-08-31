"""" Module for showing dialog box for generating a bug report as a zipfile,
whenever an uncaught exception occurs in MFiX"""

import shutil
import subprocess
import traceback

import sys
import re
import os
from datetime import datetime
from threading import Thread
from zipfile import ZipFile

from mfixgui.tools import SCRIPT_DIR
from mfixgui.tools.qt import SETTINGS
from qtpy import QtWidgets
from qtpy.QtCore import Signal, QObject

from mfixgui.version_info import get_version_info

try:
    import git
    have_git = True
except:
    have_git = False

if sys.platform == 'linux':
    have_glxinfo = bool(shutil.which('glxinfo'))
else:
    have_glxinfo = False

if sys.platform.startswith("win32"):
    conda_cmd = "conda.bat"
else:
    conda_cmd = "conda"

have_conda = bool(shutil.which(conda_cmd))

# Workaround for menuinst problem where condabin is not properly in $PATH
# "C:\Users\cgw\AppData\Local\Temp\_MEI220042\condabin" (where is that _MEI coming from?)
# This only happens when starting by clicking desktop icon on Windows

conda_prefix = os.environ.get("CONDA_PREFIX")
if conda_prefix and not have_conda:
    d = os.path.dirname(conda_prefix)
    while len(d) > 3:
        x = os.path.join(d,'condabin/', conda_cmd)
        if os.path.exists(x):
            conda_cmd = x
            have_conda = True
            break
        d = os.path.dirname(d)



FORUM_URL = "https://mfix.netl.doe.gov/forum/c/mfix"

popup = None
# Signals have to be class variables and they have to be
#  connected to a QObject
class SigWrapper(QObject):
    sig = Signal(str)

sig_report_ready = SigWrapper()

GUI = None

def excepthook(gui, devmode, etype, exc, tback):
    # set as sys.excepthook during MFiX startup
    if devmode or gui is None:
        traceback.print_exception(etype, exc, tback)
        if gui is None:
            return

    if gui.message_box and gui.message_box.isVisible():
        # Avoid cascading dialog boxes
        return

    info = (
        f'An error occurred while running MFiX. Please report this error at\
        <a href="{FORUM_URL}">{FORUM_URL}</a>'
        "<p>If you\
    continue running, the application may become unstable. Consider\
    saving your work now.")

    traceback_text = cleanup_traceback(tback, "Error: %s\n" % exc)

    ret = gui.message(
        text=info,
        buttons=["yes", "no"],
        default="no" if int(SETTINGS.value("developer_mode", 0)) else "ok",
        traceback_text=traceback_text,
        post_traceback_text="Save diagnostic report?")

    if ret == "yes":
        save_bug_report(gui, traceback_text)


def save_bug_report(gui, traceback_text=None):
    global GUI, popup
    GUI = gui
    project_file = gui.get_project_file()
    if project_file is None:
        zdirname = os.getcwd()
        stem = "mfix"
    else:
        zdirname = os.path.dirname(project_file)
        stem = gui.project.get_value("run_name", default="mfix").lower()

    timestamp = datetime.now().isoformat().replace(":", "")
    zip_filename = os.path.join(zdirname, "%s_%s.zip" % (stem,timestamp))

    popup = QtWidgets.QDialog(gui)
    popup.setMinimumWidth(200)
    popup.setMinimumHeight(50)
    layout = QtWidgets.QVBoxLayout()
    layout.addWidget(QtWidgets.QLabel("Creating MFiX bug report..."))
    popup.setLayout(layout)
    popup.setWindowTitle("Please wait")
    popup.setModal(False)
    popup.show()
    popup.raise_()
    popup.activateWindow()

    t = Thread(target=create_bug_report, args=(project_file, traceback_text, zip_filename))
    t.start()

# These two are separate functions just so we can do a non-blocking popup notification -
#  create_bug_report is slow b/c of `conda env list`

def create_bug_report(project_filename, traceback_text, zip_filename):
    b = BugReport(project_filename,traceback_text, zip_filename)
    sig_report_ready.sig.emit(zip_filename)

def review_bug_report(zip_filename):
    if popup:
        popup.close()
    zipf = ZipFile(zip_filename)
    filenames = "\n".join([str(f) for f in zipf.namelist()])

    GUI.message(
        text=f'Bug report saved to:\n\
    <a href="file://{zip_filename}">{zip_filename}</a>\n\
    You can open the zip file and see what information is being disclosed. If you consent to sharing this data,\
    post the zipfile at <a href="{FORUM_URL}">{FORUM_URL}</a>.\n\
    Press "Show Details" to see the list of file names.',
        title="Bug report saved",
        icon="info",
        detailed_text=filenames,
        print_console=False)

sig_report_ready.sig.connect(review_bug_report)


class BugReport:
    def __init__(self, project_filename, traceback_text, zip_filename):
        self.project_filename = project_filename
        self.project_dir = os.path.dirname(self.project_filename) if self.project_filename else os.getcwd()
        self.traceback_text = traceback_text
        self.zip_filename = zip_filename
        self.zipdir =  os.path.splitext(os.path.basename(zip_filename))[0]
        self.zipf = ZipFile(self.zip_filename, mode="w")
        try:
            self.write_version()
        except:
            pass
        if have_conda:
            try:
                self.write_conda_config()
            except:
                pass
            try:
                self.write_conda_list()
            except:
                pass
        try:
            self.write_traceback()
        except:
            pass
        try:
            self.write_project_files()
        except:
            pass
        if have_git:
            try:
                self.write_git_history()
            except:
                pass
        if sys.platform == 'linux' and have_glxinfo:
            try:
                self.write_glxinfo()
            except:
                pass
        self.zipf.close()

    def write_git_history(self):
        subprocess.run(['git', 'bundle', 'create', 'git-history', '--all'],
                       stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT,
                       check=False)
        if os.path.exists('git-history'):
            self.zipf.write("git-history", arcname=self.zipdir+"/git-history")

    def write_project_files(self):
        to_write = set()
        skip = {'conda_env.yml',
                'conda_list.txt',
                'git-history',
                'mfixversioninfo.txt',
                'core',
                'makefile',
                'traceback.txt'} # already handled
        skip_ext = {'.vtk', '.vtp', '.vtu', '.pvd', '.res', '.zip', '.mp4', '.dll', '.a'}
        for root, dirs, files in os.walk(self.project_dir, topdown=True):
            for d in ('build', 'build_dmp', 'build_smp', 'build_dmp',
                      'build_dmp_smp', '.git', 'mod'):
                if d in dirs:
                    dirs.remove(d)
            for d in dirs[:]:
                if d.lower() == 'cmakefiles':
                    dirs.remove(d)
            for f in files:
                if f.lower() in skip:
                    continue
                if f.lower().startswith(('core', 'mfixsolver',
                                         'libmfix', 'cmake')):
                    continue
                if f.endswith('~'): # backup file
                    continue
                if os.path.splitext(f)[-1].lower() in skip_ext:
                    continue
                to_write.add(os.path.join(root, f))
        for filepath in sorted(to_write):
            path = self.zipdir + "/" + os.path.relpath(filepath, self.project_dir)
            self.zipf.write(filepath, arcname=path)

    def write_version(self):
        path = self.zipdir + "/mfixversioninfo.txt"
        self.zipf.writestr(path, '\n'.join(get_version_info()) + '\n')

    def write_traceback(self):
        if self.traceback_text is not None:
            path = self.zipdir + "/traceback.txt"
            self.zipf.writestr(path, self.traceback_text)

    def write_glxinfo(self):
        # Linux
        glxinfo = subprocess.run(
            ["glxinfo"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False).stdout.decode('utf-8', errors='ignore')
        path = self.zipdir + "/glxinfo.txt"
        self.zipf.writestr(path, glxinfo)

    def write_dxdiag(self):
        # Windows, but is this useful?
        dxdiag = subprocess.run(
            ["dxdiag", "/dontskip", "/whql:off", "/64bit", "-t", "dxdiag.txt"],
            check=False)
        if os.path.exists("dxdiag.txt"):
            self.zipf.write("dxdiag.txt", arcname=self.zipdir+"/dxdiag.txt")


    def write_conda_list(self):
        conda_env = subprocess.run(
            [conda_cmd, "list"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False).stdout.decode('utf-8', errors='ignore')
        conda_re = re.compile("https://mfix.*")
        conda_env = conda_re.sub("mfix", conda_env)
        path = self.zipdir + "/conda_list.txt"
        self.zipf.writestr(path, conda_env)

    def write_conda_config(self):
        conda_env = subprocess.run(
            [conda_cmd, "config", "--show"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False).stdout.decode('utf-8', errors='ignore')
        conda_re = re.compile("https://mfix.*")
        conda_env = conda_re.sub("mfix", conda_env)
        path = self.zipdir + "/conda_config.txt"
        self.zipf.writestr(path, conda_env)

def cleanup_traceback(tback, err_str):
    # Remove some leading whitespace
    tb_list = [
        line[2:] if line.startswith("  ") else line
        for line in traceback.format_tb(tback)
    ]

    # Don't let the traceback be too long (e.g. "recursion too deep")
    tb_list = tb_list[-20:]

    # Shorten long pathnames
    tb_list = [line.replace(SCRIPT_DIR, "...") for line in tb_list]

    def fix_html(txt):
        return txt.replace('<', '&lt;').replace('>','&gt;')
    tb_list = [fix_html(line) for line in tb_list]

    err_str = fix_html(err_str)
    details = [err_str] + tb_list

    return "".join(details)
