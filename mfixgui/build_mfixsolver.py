#!/usr/bin/env python

import argparse
import platform
import shutil
import subprocess
import sys
import os
import errno

from mfixgui.tools import get_mfix_src

WINDOWS = sys.platform.startswith('win')

def main():
    try:
        err = main_args(*sys.argv[1:])
        sys.exit(err)
    except subprocess.CalledProcessError as err:
        print(74*"=",
              "                     BUILD FAILED",
              "",
              sep=os.linesep)
        if err.stdout:
            print("Output:")
            print(err.stdout)
        if err.stderr:
            print("Error messages:")
            print(err.stderr)
        sys.exit(err.returncode)

def main_args(*args):
    global config
    config, cmake_args = parse_args(args)
    name = 'mfixsolver'
    build_dir = 'build'
    if config.dmp:
        name += '_dmp'
        build_dir += '_dmp'
    if config.smp:
        name += '_smp'
        build_dir += '_smp'
    if config.postmfix:
        build_dir += '_postmfix'

    if config.clean:
        msg = do_clean()
        print(msg)
        return

    run_dir = os.getcwd()
    # CMake does not like '\'
    run_dir = run_dir.replace('\\', '/')
    check_run_dir(run_dir)

    check_for_mfx_file()

    if not shutil.which("cmake"):
        print("CMake not found. To build MFiX, install CMake 3.16 or newer.")
        sys.exit(-1)

    env = os.environ

    env.pop("LD_PRELOAD", None) # issues/1594

    # Set FC and CC if they are unset
    if not env.get('FC'):
        fc = shutil.which('gfortran')
        if fc:
            #print("Setting FC to %s" % fc)
            env['FC'] = fc
    if not env.get('CC'):
        cc = shutil.which('gcc')
        if cc:
            #print("Setting CC to %s" % cc)
            env['CC'] = cc

    cwd = os.getcwd()
    if not os.path.exists(build_dir):
        os.makedirs(build_dir, exist_ok=True)
    os.chdir(build_dir)

    # Should we set the build dir in CMakeLists.txt instead?  For users who
    #  don't use "build_mfixsolver"?

    # Clean up CMakeCache to prevent confusion
    try:
        os.unlink("CMakeCache.txt")
    except Exception as e:
        if e.errno != errno.ENOENT:
            raise

    err = build_target(run_dir, cmake_args)
    if err:
        if config.quiet:
            sys.exit(err)
        else:
            return err

    os.chdir(cwd)

    if not (config.batch or config.postmfix):
        if not WINDOWS:  # symlink requires extra privileges on Win
            try:
                try:
                    os.unlink(name)
                except OSError as e:
                    if e.errno != errno.ENOENT:
                        raise
                os.symlink(name+'.sh', name)
            except OSError as e:
                print("Error creating symbolic link:", e)
                if config.quiet:
                    sys.exit(1)
                return 1
    sys.stdout.flush()
    return 0

def do_clean():
    removed = False
    msg = ''
    # This is heavy-handed, should we only remove dir corresponding
    # to current settings?
    for d in  ("build", "build_dmp", "build_smp", "build_dmp_smp",
               "build_postmfix", ".build", "lib"):
        if not os.path.isdir(d):
            continue
        removed = True
        msg += ("Removing directory: %s\n" % d)
        shutil.rmtree(d, ignore_errors=True)
#    for f in ('mfixsolver' + dmp + smp + ext
#             for dmp in ('', '_dmp')
#             for smp in ('', '_smp')
#             for ext in ('', '.sh', '.bat', '.exe', '.so')):
#        if os.path.exists(f):
#            msg += ("Removing: %s\n" % f)
#            removed = True
#            os.remove(f)
    if not removed:
        msg = 'Nothing to clean.'
    return msg

def solver_name():
    desc = ''
    if config.dmp and config.smp: desc='DMP+SMP'
    elif config.dmp: desc='DMP'
    elif config.smp: desc='SMP'
    else: desc='serial'
    if config.batch: desc += ' batch'
    else: desc += ' interactive'
    return ' '.join((desc, 'solver'))

def check_for_mfx_file():
    #Guarantee that there is not more than one project file in the run dir
    # See issues/1216

    mfx_files = [f for f in os.listdir('.')
                 if f.lower().endswith('.mfx')
                 or f.lower() == 'mfix.dat']
    if len(mfx_files) > 1:
        print("ERROR, multiple .mfx files found:")
        for mfx_file in mfx_files:
            print(mfx_file)
        sys.exit(-1)

    if not config.quiet:
        if config.postmfix:
            print("Building postmfix")
        else:
            if mfx_files:
                print("Building custom %s for %s" % (solver_name(),(mfx_files[0])))
            else:
                print("Building generic %s" % solver_name())
        sys.stdout.flush()


def check_run_dir(run_dir):
    #Check for problem characters in the build_dir path
    ok = True
    for c in '&=':
        if c in run_dir:
            print("Unable to build solver because the build path:")
            print(run_dir)
            print("contains the character '%s'. Rename your working directory."%c)
            ok = False
    if not ok:
        sys.exit(-1)


def parse_args(cmdline_args):

    parser = argparse.ArgumentParser(
        description="Build MFiX solver binary"
    )

    parser.add_argument("--batch",
                       action="store_true",
                       help="Build solver without support for GUI interaction")

    parser.add_argument("--clean",
                        action="store_true",
                        help="Remove results of previous build")

    parser.add_argument("--dmp", "--enable-dmp",
                        action="store_true",
                        help="Build with DMP (MPI) support")

    parser.add_argument("-j",
                        "--jobs",
                        nargs="?",
                        type=int,
                        const=min(os.cpu_count(), 10),  # not much improvement past -j10
                        default=1,
                        action="store",
                        help="Run make with multiple parallel jobs")

    parser.add_argument("--lib-dir",
                        action="store",
                        help="Build accelerator for running test suite",
                        type=str)


    parser.add_argument("--postmfix",
                        action="store_true",
                        help="Build postmfix instead of mfixsolver")

    parser.add_argument("--quiet", "-q",
                        action="store_true",
                        help='Quiet output')

    parser.add_argument("--smp", "--enable-smp",
                        action="store_true",
                        help="Build with SMP (OpenMP) support")

    parser.add_argument("--verbose", "-v",
                        action="store_true",
                        help='Build with "make VERBOSE=1"')




    config, args = parser.parse_known_args(cmdline_args)
    if "clean" in args:
        config.clean = True
        return config, []

    check_deprecated_args(args)

    fortran_compiler = "-DMPI_Fortran_COMPILER=" if config.dmp else "-DCMAKE_Fortran_COMPILER="

    cmake_args = []
    for arg in args:
        if arg.startswith("FC="):
            arg = arg.replace("FC=", fortran_compiler)
        if arg == 'postmfix': # User did not specify -- in front of postmfix
            config.postmfix = True
            continue
        if arg.startswith("FCFLAGS="):
            arg = arg.replace("FCFLAGS=", "-DCMAKE_Fortran_FLAGS=")
        if arg.startswith("CC="):
            arg = arg.replace("CC=", "-DCMAKE_C_COMPILER=")
        cmake_args.append(arg)

    if config.postmfix:
        config.dmp = False
        config.smp = False
        config.batch = True

    if not config.batch:
        cmake_args.append("-DENABLE_PYMFIX=ON")

    if config.dmp:
        cmake_args.append("-DENABLE_MPI=ON")

    if config.smp:
        cmake_args.append("-DENABLE_OpenMP=ON")

    if config.postmfix:
        cmake_args.append("-DENABLE_POSTMFIX=ON")

    return config, cmake_args


def generate_makefile(args, run_dir):
    #Run CMake to generate Makefile(s)
    if sys.platform.startswith('win'):
        args.extend(["-G", "MinGW Makefiles"])
    args.append("-DCMAKE_INSTALL_PREFIX="+run_dir)
    if not config.postmfix:
        args.append("-DUDF_DIR="+run_dir)
    mfix_src = get_mfix_src()
    mfix_src = mfix_src.replace('\\', '/')
    args.append(mfix_src)

    cmd = ["cmake"] + args
    if not config.quiet:
        print("CMake command:", end=' ')
        for arg in cmd:
            if ' ' in arg:
                print('"'+arg+'"', end=' ')
            else:
                print(arg, end=' ')
        print()
        sys.stdout.flush()

    sp=subprocess.run(cmd,
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE,
                      text=True)
    if not config.quiet:
        print(sp.stdout)
    if sp.stderr:
        print(sp.stderr, file=sys.stderr)
        return(-1)
    if sp.returncode:
        print("Build failed with error code ", sp.returncode)
        return(sp.returncode)


def build_target(run_dir, cmake_args):
    # Run make
    err = generate_makefile(cmake_args, run_dir)
    if err:
        if config.quiet:
            sys.exit(err)
        else:
            return err
    jobs = config.jobs
    cmd = ["cmake", "--build", ".", "--target", "install"]
    if jobs > 1:
        cmd.extend(["-j", str(jobs)])
    if config.verbose:
        cmd.append("--verbose")
    if not config.quiet:
        print("Build command:", *cmd)
        sys.stdout.flush()
    if config.quiet:
        sp = subprocess.run(cmd,
                            stdout=subprocess.DEVNULL,  # Suppress standard output
                            stderr=subprocess.PIPE,
                            text=True)
        if sp.returncode:
            print("Build failed with error code", sp.returncode)
            print(sp.stderr)
            return sp.returncode

    else:
        subprocess.check_call(cmd)

    if not config.quiet:
        print(74*"=",
              "                     BUILD SUCCESSFUL",
              "",
              sep=os.linesep)


def check_deprecated_args(cmake_args):
    #To avoid user confusion, stop if more than one of FC, CMAKE_Fortran_COMPILER, or MPI_Fortran_COMPILER are specified on the command line.

    fc_args = ["-DMPI_Fortran_COMPILER", "-DCMAKE_Fortran_COMPILER", "FC"]
    fc_args_used = [
        fc_arg
        for arg in cmake_args
        for fc_arg in fc_args
        if arg.startswith(f"{fc_arg}")
    ]
    if len(fc_args_used) > 1:
        print(f"Only one of {fc_args_used} can be used.")
        print()
        print("To specify a compiler and build with MPI support:")
        print("> build_mfixsolver FC=/path/to/fortran/compiler --dmp")
        print()
        print("To specify a compiler and build without MPI support:")
        print("> build_mfixsolver FC=/path/to/fortran/compiler")
        sys.exit(-1)


if __name__ == "__main__":
    main()
