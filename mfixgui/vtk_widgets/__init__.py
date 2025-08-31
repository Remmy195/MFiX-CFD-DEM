import logging
import traceback
LOG = logging.getLogger(__name__)

VTK_AVAILABLE = True
VTK_IMPORT_INFO = []
try:
    import vtk

    # try to get the version
    try:
        from vtk import vtkVersion
    except ImportError:
        # VTK 9+
        from vtkmodules.vtkCommonCore import vtkVersion

    VTK_VERSION_STRING = vtkVersion.GetVTKVersion()
    LOG.info('VTK version: %s', VTK_VERSION_STRING)
    VTK_IMPORT_INFO.append('Successfuly imported vtk version: %s' % VTK_VERSION_STRING)
    VTK_MAJOR_VERSION = vtkVersion.GetVTKMajorVersion()

except ImportError:
    VTK_AVAILABLE = False
    vtk = None
    VTK_MAJOR_VERSION = None
    ex = traceback.format_exc()
    LOG.info("can't import vtk:\n%s", ex)
    VTK_IMPORT_INFO.append('Could not import VTK, please check your installation.\n'
                           'traceback:\n{}'.format(ex))
