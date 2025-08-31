import sys,os
import signal
import argparse
try:
    import setproctitle
except ImportError:
    setproctitle = None

from qtpy import QtGui, QtCore, QtWidgets
Qt = QtCore.Qt

from .version import get_version

from .tools import find_project_file
from .tools.qt import (get_pixmap, get_preferred_style, SETTINGS)


# NB not all imports are at the top of this file.  Ignore style
# guides, this is for a reason.  We want to make the splash screen
# appear as early as possible, and the imports can take significant
# time, especially the first time the app runs.  First impressions matter.

def main():
    main_args(sys.argv[1:])

def main_args(sys_argv):
    global gui
    gui = None
    # 'gui' is initialized here instead of at module-level
    #  so that interpreter.py can do 'from main import gui' without
    #  clobbering reference

    # handle command-line arguments
    styles_avail = [s.lower() for s in QtWidgets.QStyleFactory.keys()]
    parser = argparse.ArgumentParser(description='MFiX GUI')
    ARG = parser.add_argument
    ARG('project', action='store', nargs='?', default=None,
        help="open mfix.dat or <RUN_NAME>.mfx project file or search a specified directory for project files")
    ARG('-s', '--style', metavar='STYLE', action='store', default=None,
        choices=styles_avail,
        help='specify app style %s'%styles_avail)
    ARG('-n', '--noload', action='store_true',
        help='do not autoload previous project')
    ARG('-w', '--nonodeworks', action='store_false',
        help='do not load the nodeworks environment')
    ARG('-k', '--novtk', action='store_false',
        help='do not load vtk')
    ARG('-g', '--default_geo', action='store_true',
        help="use default geometry, don't restore previous state.")
    ARG('-d', '--developer', action='store_true',
        help="enable developer mode.")
    ARG('-c', '--clear', action='store_true',
        help="clear all saved settings.")
    ARG('-t', '--test', action='store_true',
        help="enable test mode.")
    ARG('--save', action='store_true',
        help="save the project in test mode.")
    ARG('-v', '--version', action='version', version=get_version())

    args = parser.parse_args(sys_argv)

    if args.clear:
        print("Clearing all MFIX settings from ", SETTINGS.fileName())
        SETTINGS.clear()
        SETTINGS.sync()
        return

    project_file = args.project
    if project_file and os.path.isdir(project_file):
        project_file = find_project_file(project_file)
        if project_file is None:
            print("Can't find *.mfx or mfix.dat in directory: %s" % project_file)
            parser.print_help()
            return

    elif project_file and not os.path.isfile(project_file):
        print("%s: is not a file " % project_file)
        parser.print_help()
        return

    if setproctitle:
        if project_file:
            txt = os.path.basename(project_file)
        else:
            txt = "ready"
        setproctitle.setproctitle('mfix: ' + txt)

    # Disable Qt scaling, it looks bad on most displays
    #QtCore.QCoreApplication.setAttribute(Qt.AA_Use96Dpi, False)
    #QtCore.QCoreApplication.setAttribute(Qt.AA_DisableHighDpiScaling, False)
    #QtCore.QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)


    #https://bugreports.qt.io/browse/QTBUG-70431
    QtCore.QCoreApplication.setAttribute(Qt.AA_DisableWindowContextHelpButton)

    # create the QApplication
    qapp = QtWidgets.QApplication([])
    QFont = QtGui.QFont
    font_size = SETTINGS.value('font_size', 12)
    if not font_size:
        font_size = 12
    font_size = int(font_size)
    use_system_font = int(SETTINGS.value('use_system_font', 0))
    # splash screen
    splash = None
    if not args.test:
        splash_font_size = 14
        fname = 'splash.png'
        # Rounded corners are ugly on Joule
        if os.environ.get('VGL_ISACTIVE','').startswith(('1','Y')):
            fname = 'splash-vgl.png'
        pm = get_pixmap(fname, 600, 350)
        splash = QtWidgets.QSplashScreen(pm)
        #splash.setMask(pm.mask())  # jaggies or worse
        f = QFont("Arial", splash_font_size, QFont.Bold)
        f.setPointSize(splash_font_size)
        splash.setFont(f)
        splash.setWindowFlags(Qt.SplashScreen |
                              Qt.BypassWindowManagerHint |
                              Qt.FramelessWindowHint |
                              Qt.WindowStaysOnTopHint)
        #splash.setWindowIcon(get_icon('mfix.png'))
        splash.show()
        splash.activateWindow() # issues/1651 (OSX)
        splash.raise_()         # issues/1651

    def set_splash_text(text):
        if splash is None:
            return
        splash.showMessage(text,
                           int(Qt.AlignHCenter|Qt.AlignBottom),
                           QtGui.QColor('#101010'))
        qapp.processEvents()

    set_splash_text('Starting...')

    # set style
    current_style = qapp.style().objectName().lower()
    use_style = get_preferred_style(args.style, current_style)
    if use_style:
        qapp.setStyle(use_style.lower())

    # Set exception handler early, so we catch any errors in initialization
    test_mode = args.test
    if not args.test:
        from mfixgui.bug_report import excepthook
        def hook(etype, exc, tb):
            """ avoid exceptions when catching exceptions """
            try:
                excepthook(gui, args.developer, etype, exc, tb)
            except Exception as e:
                print("Error displaying bug report")
                print(e)
        sys.excepthook = hook

    # create the gui
    set_splash_text('Loading MFiX')
    # Please do not move this import to the top of the file!
    #  it's here for a reason, to make the splash screen appear as early as possible
    from mfixgui.gui import MfixGui

    gui = MfixGui(qapp,
                  project_file=project_file,
                  loadnodeworks=args.nonodeworks,
                  loadvtk=args.novtk,
                  set_splash_text=set_splash_text,
                  style=use_style,
                  test_mode=test_mode)

    gui.font_size = font_size

    #fd = QtWidgets.QFontDialog()
    #font = fd.getFont()[0]
    #print(font.family(), font.pointSize())

    if int(SETTINGS.value('use_system_font', 0)):
        pass
    else:
        font = QFont("Arial", font_size, QFont.Normal)
        #print("USING ARIAL", font_size)
        font.setPointSize(font_size)
        qapp.setFont(font)
        gui.change_font('', font_size) # fixes up buttons and nav tree widths

    geo = SETTINGS.value('geometry')
    if geo is not None and not args.default_geo:
        # set previous geometry
        gui.restoreGeometry(geo)
        left_right = SETTINGS.value('splitter_left_right')
        if left_right is not None:
            gui.ui.splitter_left_right.setSizes([int(num) for num in left_right])

        cmd_output = SETTINGS.value('splitter_graphics_cmd_output')
        if cmd_output is not None:
            gui.ui.splitter_graphics_cmd_output.setSizes([int(num) for num in cmd_output])
    else:
        # default geometry
        geo = gui.frameGeometry()
        screen = qapp.desktop().screenNumber(qapp.desktop().cursor().pos())
        centerPoint = qapp.desktop().screenGeometry(screen).center()
        geo.moveCenter(centerPoint)
        gui.move(geo.topLeft())

    # set developer mode
    gui.enable_developer_mode(int(SETTINGS.value('developer_mode', 0)) or args.developer)

    # show the gui, unless disabled
    if not args.test:
        gui.show()

    # close the splash
    if splash is not None:
        splash.close()

    # Start update check
    # This seems to block GUI updates. Move this to a QThread.
    import threading
    threading.Thread(target=gui.check_update_available,
                     args=(get_version(),)).start()

    # Handle SMS mode (issues/1032)
    sms_enabled = SETTINGS.value('SMS', 0)
    if sms_enabled in ('false', 'False'):
        sms_enabled = 0
    elif sms_enabled in ('true', 'True'):
        sms_enabled = 1
    gui.set_sms_enabled(int(sms_enabled))

    last_project = None
    if project_file is None and not args.noload:
        # autoload last project
        last = SETTINGS.value('project_file')
        last_project = project_file = os.path.normpath(last) if last else None

    if project_file and not args.noload and os.path.exists(project_file):
        gui.open_project(project_file, interactive=(not args.test))
    else:
        gui.set_no_project()
        gui.open_file_menu()

    # have to initialize vtk after the widget is visible!
    if not (args.novtk or args.test):
        gui.vtkwidget.vtkiren.Initialize()

    # exit with Ctrl-C at the terminal
    signal.signal(signal.SIGINT, signal.SIG_DFL)

    if not args.test:
        qapp.exec_()

    else:  # Run internal test suite
        gui.navigate_all()

        # save project
        if args.save:
            gui.handle_save()


if __name__ == '__main__':
    main()
