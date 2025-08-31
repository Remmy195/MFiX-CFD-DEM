import os
import errno
import pprint
import time

from qtpy.QtCore import QTimer
from mfixgui.constants import LOCKED, UNLOCKED, PARTIAL_LOCK

class StatusManager:
    """ This class is responsible for updating the UI based on changes in job status """

    # TODO merge with job manager

    def __init__(self, mfixgui):
        self.job_state = None
        self.mfixgui = mfixgui
        self.last_elapsed = 0.0
        self.first_run_msg_printed = False
        self.last_status_msg = ""

    def update_status_if_changed(self):
        new_job_state = "stopped" # Stopped status is never reported
        if self.mfixgui.job_manager and self.mfixgui.job_manager.job:
            new_job_state = self.update_status()
        if new_job_state != self.job_state:
            prev_state = self.job_state
            self.job_state = new_job_state
            self.update_runbuttons(print_in_console=(prev_state is not None))

    def update_status(self):
        status = self.mfixgui.job_manager.job.status # oh god
        jobman = self.mfixgui.job_manager

        self.update_residuals_pane()

        response_text = self.mfixgui.ui.responses.toPlainText().lower()
        if not response_text:
            self.mfixgui.ui.responses.setText("ok")

        if not self.mfixgui.job_manager.job.is_paused():
            self.mfixgui.update_plots(status)

        self.mfixgui.update_dashboard(status)

        # Handle logging of dt, nit, residuals
        self.mfixgui.update_logger_status(status)

        self.update_progress_bar(status)
        self.mfixgui.job_manager.simtime = status.get("time", 0)
        if self.mfixgui.job_manager.job.is_paused():
            self.mfixgui.job_manager.pausing = False
            self.mfixgui.job_manager.stopping = False

        new_job_state = ("stopped" if not self.mfixgui.job_manager.job # never reported
                         else "stopping" if self.mfixgui.job_manager.stopping
                         else "pausing" if self.mfixgui.job_manager.pausing
                         else "paused" if self.mfixgui.job_manager.job.is_paused()
                         else "running")
        self.update_status_message(status, new_job_state)

        # If the job has finished, collect final output from the webserver
        #  and then tell the server it can go away

        if status.get("crashed") or status.get("finished"):
            #self.mfixgui.job_manager.exit_mfix()
            #self.mfixgui.job_manager.stopping = False
            #self.mfixgui.job_manager.job_finished()
            if not self.mfixgui.job_manager.last_rites:
                QTimer.singleShot(1000, self.mfixgui.job_manager.exit_mfix)
                self.mfixgui.job_manager.last_rites = True
            #self.mfixgui.process_manager.slot_finish(status) #XXX
            self.mfixgui.job_manager.job_finished()
            self.mfixgui.slot_rundir_timer()
            if status.get('crashed'):
                msg = 'MFiX process has crashed'
            else:
                msg = 'MFiX process has stopped'
            self.mfixgui.signal_update_runbuttons.emit(msg)

            self.mfixgui.remove_mesh_tempfiles()
            self.mfixgui.close_log_files()  #  XXX issues/1698
            if self.mfixgui.job_manager.job:
                self.mfixgui.job_manager.job.cleanup_and_exit()
                self.mfixgui.job_manager.job = None
            if self.mfixgui.job_manager.pidfile:
                try:
                    os.unlink(self.mfixgui.job_manager.pidfile)
                    self.mfixgui.job_manager.pidfile = None
                except OSError as err:
                    if err.errno != errno.ENOENT:
                        raise

        return new_job_state


    def update_residuals_pane(self): #This is the 'mfix_status' pane
        status_text = pprint.pformat(self.mfixgui.job_manager.job.status)
        header = "<html><body><pre>"
        footer = "</pre></body></html>"
        status_text = "%s%s%s" % (header, status_text, footer)
        self.mfixgui.ui.residuals.setText(status_text)

    def update_progress_bar(self, status):
        if not status.get("running"):
            self.mfixgui.progress_bar.hide()
            return

        percent = None
        if status.get("finished"):
            percent = 100
        else:
            t = status.get("time", None)
            ts = status.get("tstop", None)
            if t is not None and ts is not None and ts>0:
                percent = 100 * t/ts
        if percent is not None:
            percent = int(percent)
            percent = 100 if percent > 100 else 0 if percent < 0 else percent
            self.mfixgui.progress_bar.setValue(percent)
        self.mfixgui.progress_bar.show()

    def update_status_message(self, status, new_job_state):
        w = status.get("elapsed_time", 0)
        m, s = divmod(w, 60)
        h, m = divmod(m, 60)
        elapsed_msg = '%d:%02d:%02d' % (h, m, s)

        t = status.get("time", 0)
        if 0 <= t < 0.1:
            simtime_msg = '%0.3fms' % (1000*t)
        elif 0.1 <= t < 60:
            simtime_msg = '%0.3fs' % t
        elif 60 <= t < 3600:
            simtime_msg = '%d:%0.3f' % divmod(t, 60)
        else:
            m, s = divmod(t, 60)
            h, m = divmod(m, 60)
            simtime_msg = '%d:%02d:%0.3f' % (h, m, s)
        if self.mfixgui.job_manager.job:
            # Clean up this hack!  It's ugly but at worst some messages
            # get skipped or duplicated
            print_in_console = (new_job_state != 'running' or
                                w-self.last_elapsed >= 60
                                or not self.first_run_msg_printed)
            if print_in_console and new_job_state=='running':
                self.first_run_msg_printed = True
            if (w-self.last_elapsed > 60):
                self.last_elapsed = 60 * int(w/60)
            self.status_message(
                "MFiX %s, simulation time: %s elapsed time: %s" % (new_job_state, simtime_msg, elapsed_msg),
                print_in_console=print_in_console)
            if new_job_state != 'running':
                self.last_elapsed = 0
                self.first_run_msg_printed = False

    def update_runbuttons(self, message=None, print_in_console=True):
        # This is the main state-transition handler
        if message is not None:
            if message != self.last_status_msg:
                if not "running" in message.lower():
                    # highlight for visibility, this is an important state change
                    self.mfixgui.print_internal(message, color="blue")
                self.last_status_msg = message

        # TODO: set this in __init__ or another early setup method
        # assemble list of available executables

        project_file = os.path.basename(self.mfixgui.get_project_file() or "")
        project_open = bool(project_file and self.mfixgui.open_succeeded)
        pending = self.mfixgui.job_manager.is_job_pending()
        # why both paused and unpaused states?
        paused = bool(self.mfixgui.job_manager.job and self.mfixgui.job_manager.job.is_paused())
        pausing = self.mfixgui.job_manager.pausing
        if paused:
            pausing = False

        unpaused = bool(self.mfixgui.job_manager.job and not paused)
        resumable = (bool(self.mfixgui.get_res_files())
#                     and self.mfixgui.job_manager.simtime > 0
                     and not self.mfixgui.job_manager.job)

        stopping = self.mfixgui.job_manager.stopping
        self.mfixgui.update_window_title()  # put run state in window titlebar

        if not project_open:
            lock_level = LOCKED
        else:
            if paused or resumable:
                lock_level = PARTIAL_LOCK
            elif pending or unpaused:
                lock_level = LOCKED
            else:
                lock_level = UNLOCKED
        self.mfixgui.enable_input(lock_level)

        if pending:
            self.set_run_button(enabled=False)

        elif unpaused:
            self.set_run_button(enabled=False)

        elif paused:
            self.set_run_button(text="Unpause", enabled=True)
            if not "paused" in self.mfixgui.ui.label_status.text().lower():
                self.status_message("MFiX paused", print_in_console=False) # Want solver to print message

        elif resumable:
            self.set_run_button(text="Resume", enabled=True)
            self.status_message("Previous MFiX run is resumable.")

        else:  # Not running
            self.set_run_button(text="Run", enabled=project_open)
            if self.mfixgui.sms_mode and not self.mfixgui.mesh_accepted:
                self.status_message("Generate, inspect and accept mesh to proceed.")
            else:
                self.status_message("Ready", print_in_console=print_in_console)

        self.set_reset_button(enabled=resumable)
        self.set_pause_button(text="Pause",
                              enabled=unpaused and not (stopping or pausing))

        self.set_stop_button(enabled=pending or unpaused or paused)
        self.set_build_button(enabled=not (pending or unpaused or paused))

    def set_run_button(self, text=None, enabled=None):
        b = self.mfixgui.ui.toolbutton_run_mfix
        if self.mfixgui.sms_mode:
            if not self.mfixgui.mesh_accepted:
                b.setEnabled(False)
                b.setToolTip('Mesh must be generated and accepted before starting simulation')
            else:
                if enabled is not None:
                    b.setEnabled(enabled)
                    self.set_reset_button(enabled)
                b.setToolTip('Run MFiX simulation')
            return

        if text is not None:
            b.setToolTip("Resume previous MFiX run, or start new run" if text == "Resume"
                         else text + " MFiX simulation")
        if enabled is not None:
            b.setEnabled(enabled)
            # disable reset (delete) button while a simulation is running (issues/403)
            self.set_reset_button(enabled)

    def set_pause_button(self, text=None, enabled=None):
        b = self.mfixgui.ui.toolbutton_pause_mfix
        if enabled is not None:
            b.setEnabled(enabled)
        if text is not None:
            b.setToolTip(text + " MFIX")

    def set_stop_button(self, enabled):
        b = self.mfixgui.ui.toolbutton_stop_mfix
        b.setEnabled(enabled)
        # tooltip?

    def set_reset_button(self, enabled):
        files = False
        if enabled:
            files = self.mfixgui.get_output_files()
        b = self.mfixgui.ui.toolbutton_reset_mfix
        b.setEnabled(enabled and bool(files))

    def set_build_button(self, enabled):
        b = self.mfixgui.ui.toolbutton_compile
        b.setEnabled(enabled)

    def status_message(self, message="", print_in_console=True):
        """set the status text"""
        if message.strip() == self.mfixgui.ui.label_status.text().strip():
            return
        # pad text (why?)
        message += "  "
        self.mfixgui.ui.label_status.setText(message)

        message_prefix = message.split(":", 1)[0]
        if print_in_console:
            self.mfixgui.print_internal(message, color="blue")
