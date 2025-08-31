# -*- coding: utf-8 -*-

import errno
import simplejson as json
import subprocess
import os
import re
import signal
import socket
import time

import psutil

from qtpy.QtCore import QObject, QTimer, QUrl, Signal
from qtpy.QtNetwork import QNetworkAccessManager, QNetworkRequest, QNetworkReply

from mfixgui.tools import replace_with_dict


class JobManager(QObject):

    sig_update_job_status = Signal()
    sig_change_run_state = Signal()

    def __init__(self, parent):
        super(JobManager, self).__init__()
        self.job = None
        self.parent = parent 
        self.pidfile = None
        self.simtime = 0
        self.warning = parent.warning
        self.error = parent.error
        self.stopping = False
        self.pausing = False
        self.last_rites = False
        self.traceback_reported = False
        #self.fswatcher = QFileSystemWatcher() #Does not work, use run_dir_timer

    def cleanup_stale_pidfile(self, pidfile):
        """Delete pidfile if it refers to a job that is no longer running.
        Return True if file was stale/deleted"""
        try:
            process_info = get_process_info(pidfile)
        except Exception as e:
            self.error("Removing file %s due to exception %s" % (pidfile, e))
            try:
                os.unlink(pidfile)
            except OSError as e:
                if e.errno != errno.ENOENT:
                    print(e)
            return True

        pid = process_info.get('pid')
        job_id = process_info.get('job_id')
        url = process_info.get('url')
        host = process_info.get('host')

        if host not in ('127.0.0.1', 'localhost', socket.gethostname()):
            # TODO check status of remote job using qstat
            return False # Assume remote job is alive

        if pid is not None and psutil.pid_exists(pid):
            return False

        try:
            #self.warning("Removing %s" % pidfile)
            os.unlink(pidfile)
        except OSError as e:
            if e.errno != errno.ENOENT:
                self.error('%s: %s' % (pidfile, e))
        return True


    def try_connect(self, pidfile):
        if self.job: # Already connected
            return
        if pidfile is None:
            return
        if not os.path.isfile(pidfile):
            return

        if self.cleanup_stale_pidfile(pidfile):
            return

        self.traceback_reported = False
        self.last_rites = False
        self.parent.setup_watcher(truncate=False)

        self.job = Job(self, pidfile)

        self.pidfile = pidfile
        # connect Job signals
        self.job.sig_job_exit.connect(self.job_finished)
        self.job.sig_update_job_status.connect(self.sig_update_job_status)
        self.job.sig_change_run_state.connect(self.sig_change_run_state)

    def is_job_pending(self):
        return self.job and self.job.pending

    def is_job_paused(self):
        return self.job and self.job.is_paused()

    def pause_job(self):
        if self.job:
            self.pausing = True
            self.job.pause()

    def submit_command(self,
                       qscript,
                       submit_cmd,
                       delete_cmd, # XXX UNUSED
                       status_cmd,
                       job_id_regex,
                       replace_dict):

        project_dir = self.parent.get_project_dir()
        script_name = 'qsubmit_script'

        # write the script
        with open(os.path.join(project_dir, script_name), 'w', encoding='utf-8', errors='replace') as f:
            f.write(qscript)

        replace_dict['SCRIPT'] = script_name

        submit_cmd = replace_with_dict(submit_cmd, replace_dict)

        # submit the job
        self.parent.print_internal("Job submit command: {}".format(submit_cmd),
                                   color='blue')

        proc = subprocess.Popen(submit_cmd,
                     shell=True,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE,
                     cwd=self.parent.get_project_dir(),
                     env=dict(os.environ, LD_PRELOAD=""))
        out, err = proc.communicate()
        if job_id_regex is not None and out:
            job_id = re.findall(job_id_regex, str(out))
        else:
            job_id = []
        if job_id:
            job_id = job_id[0]
            self.parent.print_internal("Job successfully submitted with job id: {}".format(job_id),
                                       color='blue')
        else:
            self.parent.error('Could not determine job id')
            job_id = None
        if err:
            self.parent.error('Error with submission:\n{}'.format(err.decode('utf-8',errors='ignore')),
                              popup=True)

        #TODO:
        # use status_cmd to see if it queued, running, etc.

    def reset_mtime(self):
        if self.job is not None:
            self.job.mtime = 0

    def exit_mfix(self):
        if self.job is not None:
            self.job.mtime = 0
            self.job.exit_mfix()
        # XXX needed to stop reattached jobs, where the process
        # signals are not connected.
        if not self.parent.process_manager.mfix_process:
            self.parent.process_manager.slot_finish(None)


    def stop_mfix(self):
        if self.job is not None:
            self.stopping = True
            self.job.mtime = 0
            self.job.stop_mfix()
            if not os.path.isfile(self.pidfile):
                self.force_stop_mfix()


    def force_stop_mfix(self):
        if self.job is not None:
            self.stopping = True
            self.job.mtime = 0
            self.job.force_stop_mfix()
            self.job_finished()


    def disconnect(self):
        """disconnect job"""
        if self.job:
            self.job.stop_timers() # cleanup_and_exit recurses
        self.job = None
        self.sig_update_job_status.emit()
        self.sig_change_run_state.emit()

    def job_finished(self):
        """Job ended or exited. Destroy Job object and remove pidfile"""
        # TODO: this handles local only, not queued jobs
        if self.job:
            self.job.stop_timers() # cleanup_and_exit recurses
            self.job = None
            self.sig_update_job_status.emit()
            self.sig_change_run_state.emit()

        pidfile = self.pidfile
        if pidfile:
            try:
                process_info = get_process_info(pidfile)
                pid = process_info.get('pid')
                job_id = process_info.get('qjobid')
            except OSError as e:
                if e.errno != errno.ENOENT:
                    self.error("Cannot get PID: %s" % e)
                pid = job_id = None
            try:
                os.unlink(pidfile)
            except OSError as e:
                if e.errno != errno.ENOENT:
                    self.error("Error removing file: %s" % e)
        self.pidfile = None

class Job(QObject):
    """Class for managing and monitoring an MFiX job. This class contains
    methods for issuing requests to and handling responses from the MFiX API.

    pidfile: Name of file containing MFiX API connection information.
    """

    # Signals
    sig_job_exit = Signal()
    sig_response = Signal(str, dict)
    sig_update_job_status = Signal()
    sig_change_run_state = Signal()

    def __init__(self, parent, pidfile):
        super(Job, self).__init__()
        self.parent = parent # job manager
        self.warning = parent.warning
        self.error = parent.error
        self.mfix_pid = None
        self.job_id = None
        self.mtime = time.time() # Job start time
        self.status = {}
        self.pidfile = pidfile
        self.api = None
        self.pending = True
        self.api = PymfixAPI(self)

        self.mfix_pid = self.api.mfix_pid
        self.job_id = self.api.job_id
        self.run_cmd = self.api.run_cmd

        self.status_timer = QTimer()
        self.status_timer.timeout.connect(self.update_job_status)
        self.status_timer.start(100)

    def update_job_status(self):
        self.api.get_status()
        # TODO track state and emit this signal when it changes
        #self.parent.status_manager.update_runbuttons()
        #self.sig_update_job_status.emit()

    def stop_timers(self):
        # FIXME: this causes "RuntimeError: wrapped C/C++ object of type QTimer has been deleted" at exit
        if self.status_timer and self.status_timer.isActive():
            self.status_timer.stop()

    def cleanup_and_exit(self):
        try:
            self.stop_timers()
        except Exception as e:
            print(e)
        self.sig_job_exit.emit()


    def reinit(self, project_file_contents, autostart=False):
        """reinitialize job. Sanity checks (paused, gui.unsaved_flag, etc)
        must be handled by caller"""

        data = json.dumps({'mfix_dat': project_file_contents,
                           'autostart': autostart})

        # TODO move unicode conversion into 'post' method
        data = data.encode('utf-8')
        headers = {'content-type': 'application/json'}
        self.api.post('reinit', data=data, headers=headers)
        # TODO:  handlers for reinit success/failure (don't unpause job if reinit failed)
        # TODO restore GUI values to values matching running job if reinit failed (how?)


    def is_paused(self):
        return self.status.get('paused')

    def pause(self):
        self.api.get('pause')

    def unpause(self):
        self.api.get('unpause')

    def stop_mfix(self):
        self.api.get('stop')

    def exit_mfix(self):
        self.api.get('exit')

    def force_stop_mfix(self):
        """Send stop request"""
        if self.job_id:

            # Queued job; TODO qkill
            return

        if not self.mfix_pid:
            # Nothing to kill
            return

        try:
            p = psutil.Process(self.mfix_pid)
        except Exception as e:
            if isinstance(e, psutil.NoSuchProcess):
                return # Nothing to do
            else:
                raise

        is_mpi = False
        kill_pid = self.mfix_pid
        if self.run_cmd and 'mpirun' in self.run_cmd:
            while p:
                is_mpi = any('mpirun' in x for x in p.cmdline()) # cmdline is list of strings
                if is_mpi:
                    kill_pid = p.pid
                    break
                p = p.parent()
            else:
                return

        if is_mpi:
            sig = signal.SIGTERM
        else:
            sig = signal.SIGTERM
        try:
            os.kill(kill_pid, sig)
        except Exception as e:
            self.error("force_kill: %s" % e)


    def handle_stop_mfix(self, req_id, response):
        """Handler for responses to stop requests"""
        self.status_timer.stop()
        self.handle_response(req_id, response)
        self.sig_change_run_state.emit()
        self.cleanup_and_exit()

def get_process_info(filename):
    """Read contents of provided MFiX job pid file and set supported
    key-value pairs as dictionary members
    """
    d = {}
    try:
        with open(filename, encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                tok = [x.strip() for x in line.split('=', 1)]
                if len(tok) != 2:
                    raise ValueError(line)
                k, v = tok
                if v.isdigit():
                    v = int(v)
                else:
                    if v == 'None':
                        v = None
                d[k] = v
    except OSError as e:
        if e.errno != errno.ENOENT: # Return empty dict if no file
            raise
    return d


class PymfixAPI(QNetworkAccessManager):

    def __init__(self, parent):

        super(PymfixAPI, self).__init__()
        self.parent = parent #Job
        self.warning = parent.warning
        self.error = parent.error
        self.pidfile = parent.pidfile
        self.process_info = get_process_info(self.pidfile)
        self.mfix_pid = self.process_info.get('pid')
        self.job_id = self.process_info.get('job_id')
        self.run_cmd = self.process_info.get('run_cmd')

        self.methods = {'put': super(PymfixAPI, self).put,
                        'get': super(PymfixAPI, self).get,
                        'post': super(PymfixAPI, self).post,
                        'delete': super(PymfixAPI, self).deleteResource}

    def get_status(self):
        url = self.process_info.get('url')
        if not url:
            return
        host = self.process_info.get('host') # could get from url
        if not url:
            errmsg = "cannot connect, no URL found for mfix process"
            self.error(errmsg)
            self.parent.status = {'error': errmsg}
            #self.get_status_error()
            self.parent.parent.parent.slot_update_status()
            self.parent.status_manager.update_runbuttons()

        # issues/412 circumvent network check if job is local
        if host in ('127.0.0.1', 'localhost'):
            self.setNetworkAccessible(True)
        request = QNetworkRequest(QUrl(url + '/status'))
        reply = self.methods['get'](request) # ugh
        reply.errorOccurred.connect(lambda error, reply=reply, self=self: self.get_status_error(error, reply)) # pylint: disable=undefined-variable
        reply.finished.connect(lambda reply=reply, self=self: self.get_status_finished(reply)) # pylint: disable=undefined-variable


    def get_status_finished(self, reply):
        try:
            content_type = reply.header(QNetworkRequest.ContentTypeHeader)
            if not content_type:
                return
            response = reply.readAll()
            if not response:
                return

            self.parent.pending = False

            data = response.data().decode('utf8', errors='ignore')
            if 'application/json' in content_type:
                data = json.loads(data)

                mfix_status = data.get('mfix_status',
                                       {'error': 'no status returned'})

                self.parent.status = mfix_status # FIXME: leads to stale status
                self.parent.parent.parent.slot_update_status() # WTF
                self.parent.parent.parent.slot_update_runbuttons()
                #     ^     ^      ^
                #     \job  \mgr   \GUI
                # FIXME FIXME FIXME
                # let's change to passing always-fresh status as a parameter
                #self.parent.parent.parent.slot_update_status(mfix_status)
            elif 'text/html' in content_type:
                #self.parent.parent.parent.error(data, popup=False)
                #no, too many popups!
                self.show_response(data)
            else:
                print("Unknown content type: %s", content_type)
        except Exception as e:
            print("JSON decode error:", e)

        finally:
            reply.deleteLater()

    def show_response(self, text): #
        try:
            rsp = self.parent.parent.parent.ui.responses
            if hasattr(rsp, 'setText'):
                rsp.setText(text) # TODO clean this up
            if not ("200" in text
                    or "2 Connection closed" in text
                    or "1 Connection refused" in text):
                self.parent.parent.parent.error(text)

        except Exception as e:
            print(e, text)


    def get_status_error(self, error_code, reply):
        try:
            if error_code == QNetworkReply.ConnectionRefusedError:
                ## todo,  increment error count
                return
            self.show_response("%s %s" % (error_code, reply.errorString()))
        finally:
            reply.deleteLater()



    def request(self, method, endpoint, data=None, headers=None):
        if not headers:
            headers = {}
        url = self.process_info.get('url')
        token = self.process_info.get('token')
        host = self.process_info.get('host')
        if host in ('127.0.0.1', 'localhost'):
            self.setNetworkAccessible(True)
        req = QNetworkRequest(QUrl('%s/%s' % (url, endpoint)))
        method = self.methods[str(method).lower()]

        # assemble headers
        if 'content-type' not in [x.lower() for x in headers.keys() if isinstance(x, str)]:
            headers['content-type'] = 'text/plain'
        headers['content-length'] = str(len(data if data else ''))
        if token:
            (name, value) = token.split(':')
            headers[name] = value
        for k, v in headers.items():
            kb = k.encode('utf-8') # why?
            vb = v.encode('utf-8')
            req.setRawHeader(kb, vb)

        if data is not None:
            reply = method(req, data)
        else:
            reply = method(req)
        reply.errorOccurred.connect(lambda error, reply=reply, self=self: # pylint: disable=undefined-variable
                            self.handle_error(error, reply))
        reply.finished.connect(lambda reply=reply, self=self: # pylint: disable=undefined-variable
                               self.print_reply(reply))
        return reply

    def print_reply(self, reply):
        try:
            content_type = reply.header(QNetworkRequest.ContentTypeHeader)
            if not content_type:
                return
            response = reply.readAll()
            if not response:
                return
            data =  response.data().decode('utf-8', errors='ignore')

            if '<title>' not in data: # don't mess up preformmated page
                http_code = reply.attribute(QNetworkRequest.HttpStatusCodeAttribute)
                data = "%s %s" % ( http_code, data)
            if 'text/plain' in content_type:
                self.show_response(data)
            elif 'text/html' in content_type:
                self.show_response(data)
            else:
                print("Unknown content_type: %s" % content_type)

        finally:
            reply.deleteLater() # Will python gc these?

    def handle_error(self, error_code, reply):
        try:
            if error_code == QNetworkReply.ConnectionRefusedError:
                if self.parent.parent.stopping:
                    return
            self.show_response("%s %s" % (error_code, reply.errorString()))
        finally:
            reply.deleteLater() # Will python gc these?

    def get(self, endpoint, handlers=None):
        return self.request('get', endpoint)

    def put(self, endpoint, data=b'', handlers=None):
        return self.request('put', endpoint, data=data)

    def post(self, endpoint, data=b'', headers=None):
        return self.request('post', endpoint, data=data, headers=headers)
