#!/usr/bin/env python

import argparse
import glob
import os
import errno
import signal
import socket
import sys
import ctypes
try:
    import setproctitle
except ImportError:
    setproctitle = None


from threading import Thread

port = 0 # OS finds free port
job_id = None
job_is_remote = False

run_cmd = os.environ.get('MFIX_RUN_CMD', ' '.join(sys.argv))

run_name = None # gets set when *.mfix file is loaded

WINDOWS = sys.platform.startswith('win')

HERE = os.path.dirname(os.path.abspath(__file__))

class MfixThread(Thread):
    def __init__(self, filename, loglevel):
        Thread.__init__(self)
        self.filename = filename
        self.loglevel = loglevel

    def run(self):
        r = run_mfix(self.filename, self.loglevel)
        return r

def parse_args():
    parser = argparse.ArgumentParser(description='Welcome to PYMFiX')
    ARG = parser.add_argument

    ARG("-f", "--file",
        default="mfix.dat",
        help="specify an input file (*.mfx or *.dat)")

    ARG("-m", "--mfixsolver",
        default=None,
        action="store",
        help="full path to mfixsolver module (mfixsolver.so, mfixsolver_dmp.so, etc.)")

    ARG("-p", "--print-flags",
        action="store_true",
        help="display the solver configuration flags")

    ARG("-l", "--log",
        metavar="LOG",
        action="store",
        default="INFO",
        choices=["ERROR", "WARNING", "STATUS", "INFO"],
        help="set logging level (ERROR, WARNING, STATUS, INFO)")

    ARG("-c", "--clean",
        action="store_true",
        help="clean output files before starting run")

    ARG("-s", "--server",
        action="store_true",
        help="start HTTP server")

    return parser.parse_args()


def clean():
    """Remove output files from previous solver run. Must be called AFTER
    parallel_init() so that mype is initialized.
    """
    if mype.value != 0: # DMP
        return
    # For SMP we might still have multiple threads, so treat ENOENT as non-fatal
    for ext in 'vt?', 'sp?', 'out', 'res', 'pvd', 'res', 'log', 'pid', 'env', 'msh':
        for ext in ext, ext.upper():
            for f in glob.glob('*.' + ext):
                try:
                    os.remove(os.path.normcase(os.path.realpath(f)))
                    #print("removing %s" % f, flush=True)
                except OSError as e:
                    if e.errno == errno.ENOENT: #File deleted already
                        pass
                    else:
                        raise


def get_run_name(filename):
    global run_name
    for line in open(filename, encoding='utf-8', errors='replace'):
        line = line.strip()
        if line.lower().startswith('run_name') and '=' in line:
            tok = line.split('=', 1)[1].strip()
            # Remove quotes if present.
            # NB: Due to a bug, in some files, run_name may lack quotes
            if tok.startswith('"') or tok.startswith("'"):
                tok = tok[1:]
            if tok.endswith('"') or tok.endswith("'"):
                tok = tok[:-1]
            run_name = tok.strip()
            break



def write_pidfile():
    # Avoid importing Flask in global namespace, to keep -p etc fast
    pid = os.getpid()
    protocol = 'http'
    if job_is_remote:
        hostname = socket.gethostname()
    else:
        # issues/407 just use loopback for local jobs
        hostname = '127.0.0.1'
    with open(run_name + '.pid', 'w', encoding='utf-8', errors='replace') as f:
        f.write('host=%s\n' % hostname)
        f.write('pid=%s\n' % pid)
        f.write('url=%s://%s:%s\n' %
                (protocol, hostname, port))
        f.write('job_id=%s\n' % job_id)
        f.write('run_cmd=%s\n' % run_cmd)
    # Save environment for debugging
    with open(run_name+'.env', 'w', encoding='utf-8', errors='replace') as f:
        for (k,v) in sorted(os.environ.items()):
            f.write('%s=%s\n' % (k,v))


def find_free_port():
    # Note that this is racy - the port could get used by
    # another program before the webserver binds to it
    global port
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.bind(("", 0))
    port = sock.getsockname()[1]
    sock.close()


def redirect_stdio():
    # I do not know why these are reversed on Windows
    stdout_file = os.open(run_name+'.stdout', os.O_RDWR|os.O_CREAT|os.O_TRUNC)
    stderr_file = os.open(run_name+'.stderr', os.O_RDWR|os.O_CREAT|os.O_TRUNC)
    os.dup2(stdout_file, sys.stdout.fileno())
    os.dup2(stderr_file, sys.stderr.fileno())

server = None
def start_webserver():
    # Avoid importing pymfix_webserver in global namespace,
    # because we don't always need it and loading Flask is slow
    global server
    try:
        import mfixgui.pymfix_webserver as W
    except ImportError:
        import pymfix_webserver as W
    server = W

    # Set globals in W
    W.init(M, mfix_thread, run_name, job_is_remote)
    W.start(port=port)

def stop_webserver():
    if server:
        server.exit()

def main():
    global M, DBG, mype, mfix_thread, options, job_id, job_is_remote

    options = parse_args()

    if not options.print_flags and not os.path.exists(options.file):
        print("Error: file not found: %s"%options.file, file=sys.stderr)
        sys.exit(-1)

    job_id = (
        os.environ.get("JOB_ID")  # GridEngine
        or os.environ.get("PBS_JOBID")  # Torque
        or os.environ.get("SLURM_JOB_ID")  # SLURM
        or None  # local
    )
    job_is_remote = (job_id is not None) or ("qsub" in run_cmd) or ("sbatch" in run_cmd)

    loglevels = {"ERROR": 1, "WARNING": 2, "STATUS": 3, "INFO": 4}
    loglevel = loglevels.get(options.log, 3)

    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)



    if not options.mfixsolver:
        print("-m/--mfixsolver must be specified", flush=True)
        sys.exit(-1)

    so_path = options.mfixsolver

    if not os.path.isabs(so_path):
        so_path =  os.getcwd() + '/' + so_path

    if not os.path.exists(so_path):
        print("Cannot find mfixsolver module", flush=True)
        if so_path:
            print("Tried", so_path, flush=True)
        sys.exit(-1)

    dirname, filename = os.path.split(options.file)
    if dirname:
        os.chdir(dirname)

    name = 'mfixsolver'
    for x in '_dmp', '_smp':
        if x in so_path.lower():
            name += x
    if setproctitle and not options.print_flags:
        setproctitle.setproctitle(name + ": " + os.path.basename(options.file))

    # CDLL does not accept relative path
    if not os.path.isabs(so_path):
        so_path = os.path.join(HERE, so_path)

    # set the global run_name variable
    # We do this here instead of getting it from the solver,
    #  because the solver may fail to load the input file

    if options.server:
        get_run_name(filename)
        if not run_name:
            if mype.value == 0:
                print("ERROR, run_name not set in %s" % filename, flush=True)
            sys.exit(1)
        redirect_stdio()

    M = DBG = None
    if WINDOWS:
        parent = os.path.dirname(HERE)
        grandparent = os.path.dirname(parent)
        greatgrandparent = os.path.dirname(grandparent)
        d3_path = None
        # Ugh, install location keeps changing
        for d in (os.path.join(greatgrandparent,'Scripts'),
                  HERE,
                  parent,
                  os.path.join(parent, 'mfixgui'),
                  os.path.join(parent, 'mfix_solver'),
                  grandparent,
                  greatgrandparent):
            tmp = d + '\\d3.dll'
            if os.path.exists(tmp):
                d3_path = tmp
                break
        if d3_path:
            try:
                DBG = ctypes.CDLL(d3_path)
                DBG.get_error_code.restype = ctypes.c_uint64
                DBG.get_error_msg.restype = ctypes.c_char_p
                DBG.load_module.restype = ctypes.c_uint64
                DBG.load_module.argtypes = (ctypes.c_char_p,)
                dh = DBG.load_module(bytes(so_path, 'utf8'))
                if dh == 0:
                    print("Error loading %s from d3" % so_path, flush=True)
                else:
                    M = ctypes.CDLL('mfixsolver.so', handle=dh)
                    print("Using %s"%d3_path, flush=True)
                sys.stdout.flush()
            except Exception as e:
                print("Error loading d3.dll", e, flush=True)
                M = DBG = None

    if M is None:
        sys.stdout.flush()
        M = ctypes.CDLL(so_path)

    if options.print_flags:
        M.__main_MOD_print_flags()
        sys.exit(0)

    M.__parallel_mpi_MOD_parallel_init()

    mype = ctypes.c_int.in_dll(M, '__compar_MOD_mype')
    sys.stdout.flush()
    if mype.value == 0:
        print("Using", so_path, flush=True)

    # TODO check if pidfile exists & process is running.  If so, refuse
    #  to start
    sys.stdout.flush()

    if options.clean:
        clean()

    if options.server:
        get_run_name(filename)
        if not run_name:
            if mype.value == 0:
                print("ERROR, run_name not set in %s" % filename, flush=True)
            sys.exit(1)
        if mype.value == 0:
            find_free_port()
            write_pidfile()
        mfix_thread = MfixThread(filename, loglevel)
        mfix_thread.start()

        # start the Flask server on rank 0
        if mype.value == 0:
            start_webserver()
        mfix_thread.join()
        stop_webserver()
        try:
            os.unlink(run_name + '.pid')
        except:
            pass

    else: # --server not specified, don't start webserver
        # don't need separate thread for mfix
        r = run_mfix(filename, loglevel)
    M.__parallel_mpi_MOD_parallel_fin()

def run_mfix(filename, loglevel):
    filename_buf = ctypes.create_string_buffer(bytes(filename, 'utf8', 'ignore'))
    assert M.gfortran_init_() == 57 # must be done per-thread
    try:
        if DBG:
            sys.stdout.flush()
            DBG.run_mfix_dbg(filename_buf, ctypes.c_int(loglevel))
            if DBG.get_error_code():
                errmsg = DBG.get_error_msg().decode('utf8', 'ignore')
                loc = errmsg.splitlines()[-1].split('at', 1)[1].strip()
                if server:
                    server.crashed = True
                    server.error = errmsg
                sys.stderr.write(' >>>>>-----\n')
                sys.stderr.write("Error from %s\n" % loc)
                sys.stderr.write("Error: Solver crash!\n")
                sys.stderr.write("The MFiX solver has terminated unexpectedly\n")
                sys.stderr.write("Error information:\n")
                sys.stderr.write(errmsg)
                sys.stderr.write('<<<<<-----\n')
                sys.stderr.flush()
                return False
            else:
                return True
        else: # no d3.dll
            run_mfix0 = M.__main_MOD_run_mfix0
            run_mfix0(filename_buf,
                      ctypes.byref(ctypes.c_int(len(filename))),
                      ctypes.byref(ctypes.c_int(loglevel)))
            return True

    # On Windows, if we are not using the d3 wrapper,
    # the built-in signal handers for SIGSEGV and SIGFPE turn into Python exceptions
    # which we trap here.  These do not contain information about the solver stack.
    except OSError as e:
        if server:
            server.crashed = True
            server.error = str(e)
        sys.stderr.write(' >>>>>-----\n')
        sys.stderr.write("ERROR: Solver crash!\n")
        sys.stderr.write("The MFiX solver has terminated unexpectedly\n")
        sys.stderr.write("Error information:\n")
        sys.stderr.write(str(e) + '\n')
        sys.stderr.write('<<<<<-----\n')
        sys.stderr.flush()
        return False

if __name__ == '__main__':
    main()
