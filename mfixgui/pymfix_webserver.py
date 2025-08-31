#!/usr/bin/env python

# This module is only intended to be used with 'pymfix' since it depends
#  on pymfix to start and manage the mfix thread, etc.
# This code in a separate module so that pymfix can optionally run without
#  the web interface, in which case we skip loading this module

import simplejson # make sure Flask uses 'simplejson' backend
from flask import Flask, jsonify, request
from werkzeug.serving import ThreadedWSGIServer
from threading import Thread
import os
import sys
import random
import string
import ctypes
from time import sleep

import flask.cli
flask.cli.show_server_banner = lambda *_: None

# Prevent logging every request to stdout
import logging
log = logging.getLogger('werkzeug')
log.setLevel(logging.ERROR)

#Globals
F = Flask(__name__)
W = None # WSGIServer

run_name = 'unset' # set by pymfix
mfix_thread = None # set by pymfix
M = None # set by pymfix
job_is_remote = False #set by pymfix
crashed = False # set by pymfix
error = None # set by pymfix
do_time_step_summary = False

# Views
@F.route('/')
def default():
    return "hello"

@F.route('/exit')
def exit():
    def timebomb():
        sleep(1)
        W.shutdown()
    Thread(target=timebomb).start()
    return "goodbye"

@F.route('/stop')
def stop():
    exit_flag.value = True
    pause_flag.value = False
    # Write an 'MFIX.STOP' file?  kill the process?
    return "stop OK"

@F.route('/pause')
def pause():
    pause_flag.value = True
    return "pause OK"

@F.route('/unpause')
def unpause():
    pause_flag.value = False
    return "unpause OK"

@F.route('/status')
def status():
    running = mfix_thread and mfix_thread.is_alive()
    elapsed_time = M.__time_mod_MOD_get_elapsed_time()
    paused_time = M.__time_mod_MOD_get_paused_time()
    remaining_time = M.__time_mod_MOD_get_remaining_time()
    io_time = M.__time_mod_MOD_get_io_time()
    finished = bool(exit_flag) # or time.value + 0.1*dt.value >= tstop.value #keep this check for backward compat

    if version.value is None:
        ver = 'Unknown'
    else:
        ver = version.value.decode('utf8', errors='ignore').strip()

    d = {
        'crashed': crashed,
        'dt': dt.value,
        'elapsed_time': elapsed_time,
        'error': error,
        'finished': finished,
        'io_time': io_time,
        'nit': nit.value,
        'paused': pause_flag.value,
        'paused_time': paused_time,
        'pid': os.getpid(),
        'remaining_time': remaining_time,
        'run_name': run_name,
        'running': running,
        'time': time.value,
#       'tstart': tstart.value,
        'tstop': tstop.value,
        'version': ver,
    }

    if group_resid.value:
        get_count, get_name, get_val = (M.__residual_pub_MOD_get_resid_grp_string_count,
                                        M.__residual_pub_MOD_get_resid_grp_string,
                                        M.__residual_pub_MOD_get_resid_grp)
    else:
        get_count, get_name, get_val = (M.__residual_pub_MOD_get_resid_string_count,
                                        M.__residual_pub_MOD_get_resid_string,
                                        M.__residual_pub_MOD_get_resid)

    seen = set()
    residuals = []
    buf = ctypes.create_string_buffer(8)
    for i in range(1,1+get_count()):
        arg = ctypes.byref(ctypes.c_int(i))
        get_name(arg, buf, ctypes.c_int(8))
        name = buf.value.strip()
        if isinstance(name, bytes):
            name = name.decode('utf8', errors='ignore')
        if not name:
            continue
        if name in seen:
            continue
        seen.add(name)
        val = get_val(arg)
        residuals.append((name, val))

    d['residuals'] = residuals

    if do_time_step_summary:
        # Time step summary
        tssh_str = tssh.value.decode('utf8', errors='ignore')
        tssi_str = tssi.value.decode('utf8', errors='ignore')

        tss = []
        # Time step summary (header and info) can have up to 100 entries, 12 characters each
        for i in range(0,100):
            name = tssh_str[i*12:(i+1)*12].strip()
            val  = tssi_str[i*12:(i+1)*12].strip()
            if name != '' and val != '' and val !='-':
            # if name is not '' and val is not '':
                tss.append((name, val))

        d['time step summary'] = tss

    output_bytes = {'stdout_bytes':  os.fstat(sys.stdout.fileno()).st_size,
                    'stderr_bytes':  os.fstat(sys.stderr.fileno()).st_size}

    return jsonify(mfix_status=d,
                   output_bytes=output_bytes)


@F.route('/reinit', methods=['POST'])
def reinit():
    if not pause_flag:
        return("model not paused", 403)

    try:
        data = request.get_json(force=True)
        mfix_dat = data.get('mfix_dat')
        autostart = data.get('autostart')
        # It's really too bad that we allowed utf-8 into mfix.dat,
        # can we undo that decision?
        mfix_dat = mfix_dat.encode('utf-8', 'ignore')
        #M.pause.set_reinit_data(mfix_dat)
        M.__pause_MOD_set_reinit_data(ctypes.create_string_buffer(mfix_dat),
                                      ctypes.byref(ctypes.c_int(len(mfix_dat))))
        autostart_flag.value = bool(autostart)
        reinit_flag.value = True
        return "reinit OK"
    except Exception as e:
        return(str(e), 500)


@F.route('/help')
def help():
    views = [v for v in F.view_functions.keys()
             if v not in ('default', 'static')]
    return '<br>'.join(views)


# placate pylint
c_bool = c_double = c_int = lambda x,y: None

def init(M_, mfix_thread_, run_name_, job_is_remote_):
    global M, autostart_flag, dt, exit_flag, group_resid, job_is_remote
    global mfix_thread, nit, pause_flag, reinit_flag, run_name, time, time_start
    global tstop, version, wall0, wall_paused, wall_start
    global tssh, tssi

    M = M_
    mfix_thread = mfix_thread_
    run_name = run_name_
    job_is_remote = job_is_remote_

    def mangle(mod, var):
        return f'__{mod}_MOD_{var}'

    # Define useful functions.  pylint doesn't understand this
    for t in 'int', 'double', 'bool', 'char_p':
        globals()['c_'+t] = lambda mod, var,t=t: getattr(ctypes, 'c_'+t).in_dll(M, f'__{mod}_MOD_{var}')

    # Avoid doing the attribute lookups on every status request
    autostart_flag = c_bool('pause', 'autostart_flag')
    dt =  c_double('run', 'dt')
    exit_flag = c_bool('exit', 'exit_flag')
    group_resid = c_bool('residual_pub', 'group_resid')
    nit = c_int('iterate', 'nit')
    pause_flag = c_bool('pause', 'pause_flag')
    reinit_flag = c_bool('pause', 'reinit_flag')
    time =  c_double('run', 'time')
    tstart = c_double('run', 'tstart')
    tstop = c_double('run', 'tstop')
    version = (ctypes.c_char*80).in_dll(M, mangle('run', 'project_version'))

# Time step summary
    tssh = (ctypes.c_char*1200).in_dll(M, mangle('output', 'tssh'))
    tssi = (ctypes.c_char*1200).in_dll(M, mangle('output', 'tssi'))

    M.__time_mod_MOD_get_elapsed_time.restype = ctypes.c_double
    M.__time_mod_MOD_get_remaining_time.restype = ctypes.c_double
    M.__time_mod_MOD_get_paused_time.restype = ctypes.c_double
    M.__time_mod_MOD_get_io_time.restype = ctypes.c_double
    M.__residual_pub_MOD_get_resid_grp.restype = ctypes.c_double
    M.__residual_pub_MOD_get_resid.restype = ctypes.c_double
    M.__residual_pub_MOD_get_resid_string.restype = None
    M.__residual_pub_MOD_get_resid_grp_string.restype = None

def start(port=None):
    global W
    W = ThreadedWSGIServer(
        host='0.0.0.0',
        port=port,
        app=F)

    t = Thread(target=W.serve_forever)
    t.start()


#    F.run(host='0.0.0.0',
#          port=port,
#          use_reloader=False)
