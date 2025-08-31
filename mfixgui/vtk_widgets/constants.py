"""
Constants for vtk
"""
import vtk
from qtpy.QtGui import QColor

DEFAULT_MESH_NAME = 'mesh.vtu'

CELL_TYPE_ENUM = {
    0:  'empty_cell',
    1:  'vertex',
    2:  'poly_vertex',
    3:  'line',
    4:  'poly_line',
    5:  'triangle',
    6:  'triangle_strip',
    7:  'polygon',
    8:  'pixel',
    9:  'quad',
    10: 'tetra',
    11: 'voxel',
    12: 'hexahedron',
    13: 'wedge',
    14: 'pyramid',
    15: 'pentagonal_prism',
    16: 'hexagonal_prism',
    21: 'quadratic_edge',
    22: 'quadratic_triangle',
    23: 'quadratic_quad',
    24: 'quadratic_tetra',
    25: 'quadratic_hexahedron',
    26: 'quadratic_wedge',
    27: 'quadratic_pyramid',
    28: 'biquadratic_quad',
    29: 'triquadratic_hexahedron',
    30: 'quadratic_linear_quad',
    31: 'quadratic_linear_wedge',
    32: 'biquadratic_quadratic_wedge',
    33: 'biquadratic_quadratic_hexahedron',
    34: 'biquadratic_triangle',
    35: 'cubic_line',
    36: 'quadratic_polygon',
    }

# Procedural
PROCEDURAL_DICT = {
    'cylinder':            None,
    'bend':                None,
    'circle_to_rectangle': None,
    'body_of_revolution':  None,

}
DEFAULT_PROCEDURAL_PARAMS = {
    'radius':                  0.5,
    'radiusresolution':        20,
    'frontradius':             0.5,
    'backradius':              0.5,
    'bendmajorradius':         0.5,
    'bendminorradius':         0.5,
    'width':                   1.0,
    'widthresolution':         20,
    'height':                  1.0,
    'heightresolution':        20,
    'length':                  1.0,
    'lengthresolution':        20,
    'frontlength':             1.0,
    'frontresolution':         20,
    'backlength':              1.0,
    'backresolution':          20,
    'circumferencestartangle': 0.0,
    'bendstartangle':          0.0,
    'circumferenceendangle':   360.0,
    'bendendangle':            90.0,
    'circumferenceresolution': 32,
    'bendresolution':          32,
    'bottomcap':               False,
    'bottomresolution':        20,
    'topcap':                  False,
    'topresolution':           20,
    'offset':                  0,
    'offsetx':                 0,
    'offsety':                 0,
    'offsetz':                 0,
    'geo_type':                'procedural',
    'type':                    'cylinder',
    'visible':                 True,
    'centerx':                 0.0,
    'centery':                 0.0,
    'centerz':                 0.0,
    'lengthx':                 1.0,
    'lengthy':                 1.0,
    'lengthz':                 1.0,
    'rotationx':               0.0,
    'rotationy':               0.0,
    'rotationz':               0.0,
    'startangle':              0.0,
    'endangle':                360.0,
    'segmentdict':            {},
}

# Implicits
IMPLICIT_DEFAULT_RES = 40
IMPLICIT_DICT = {
    'sphere':       vtk.vtkSphere,
    'box':          vtk.vtkBox,
    'cylinder':     vtk.vtkCylinder,
    'cone':         vtk.vtkCone,
    'quadric':      vtk.vtkQuadric,
    'superquadric': vtk.vtkSuperquadric,
    }

DEFAULT_IMPLICIT_PARAMS = {
    'bounds':    [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0],
    'minx':      -1.0,
    'maxx':      1.0,
    'miny':      -1.0,
    'maxy':      1.0,
    'minz':      -1.0,
    'maxz':      1.0,
    'centerx':   0.0,
    'centery':   0.0,
    'centerz':   0.0,
    'rotationx': 0.0,
    'rotationy': 0.0,
    'rotationz': 0.0,
    'radius':    1.0,
    'lengthx':   1.0,
    'lengthy':   1.0,
    'lengthz':   1.0,
    'height':    1.0,
    'a0':        0.5,
    'a1':        1.0,
    'a2':        0.2,
    'a3':        0.0,
    'a4':        0.1,
    'a5':        0.0,
    'a6':        0.0,
    'a7':        0.2,
    'a8':        0.0,
    'a9':        0.0,
    'phi':       1.0,
    'theta':     1.0,
    'thickness': 0.3333,
    'toroidal':  True,
    'geo_type':  'implicit',
    'type':      '',
    'visible':   True,
    'sample':    None,
    }

# Primitives
PRIMITIVE_DICT = dict([
    ('sphere',   vtk.vtkSphereSource),
    ('box',      vtk.vtkCubeSource),
    ('cylinder', vtk.vtkCylinderSource),
    ('cone',     vtk.vtkConeSource),
    ('superquadric', vtk.vtkSuperquadricSource)
    ])

DEFAULT_PRIMITIVE_PARAMS = {
    'centerx':         0.0,
    'centery':         0.0,
    'centerz':         0.0,
    'rotationx':       0.0,
    'rotationy':       0.0,
    'rotationz':       0.0,
    'radius':          1.0,
    'directionx':      1.0,
    'directiony':      0.0,
    'directionz':      0.0,
    'lengthx':         1.0,
    'lengthy':         1.0,
    'lengthz':         1.0,
    'height':          1.0,
    'resolution':      10,
    'thetaresolution': 10,
    'phiresolution':   10,
    'theta':           1.0,
    'phi':             1.0,
    'size':            0.5,
    'scalex':          1.0,
    'scaley':          1.0,
    'scalez':          1.0,
    'toroidal':        0,
    'thickness':       0.3333,
    'visible':         True,
    'geo_type':        'primitive',
    'type':            '',
    }

# Parametrics
PARAMETRIC_DICT = dict([
    ('torus',           vtk.vtkParametricTorus),
    ('boy',             vtk.vtkParametricBoy),
    ('conic_spiral',    vtk.vtkParametricConicSpiral),
    ('cross_cap',       vtk.vtkParametricCrossCap),
    ('dini',            vtk.vtkParametricDini),
    ('ellipsoid',       vtk.vtkParametricEllipsoid),
    ('enneper',         vtk.vtkParametricEnneper),
    ('figure_8_klein',  vtk.vtkParametricFigure8Klein),
    ('klein',           vtk.vtkParametricKlein),
    ('mobius',          vtk.vtkParametricMobius),
    ('random_hills',    vtk.vtkParametricRandomHills),
    ('roman',           vtk.vtkParametricRoman),
    ('super_ellipsoid', vtk.vtkParametricSuperEllipsoid),
    ('super_toroid',    vtk.vtkParametricSuperToroid),
    ])

DEFAULT_PARAMETRIC_PARAMS = {
    'translationx':       0.0,
    'translationy':       0.0,
    'translationz':       0.0,
    'centerx':            0.0,
    'centery':            0.0,
    'centerz':            0.0,
    'rotationx':          0.0,
    'rotationy':          0.0,
    'rotationz':          0.0,
    'radius':             1.0,
    'radiusx':            1.0,
    'radiusy':            1.0,
    'radiusz':            1.0,
    'ringradius':         1.0,
    'crosssectionradius': 0.5,
    'zscale':             0.125,
    'ascale':             0.8,
    'bscale':             0.2,
    'bfunc':              1.0,
    'cfunc':              0.1,
    'nfunc':              2.0,
    'nhills':             30,
    'variancex':          2.5,
    'scalex':             0.3,
    'variancey':          2.5,
    'scaley':             0.3,
    'amplitude':          2.0,
    'scaleamplitude':     0.3,
    'allowrandom':        True,
    'n1':                 1.0,
    'n2':                 1.0,
    'visible':            True,
    'geo_type':           'parametric',
    'type':               ''}

# Filters
FILTER_DICT = dict([
    ('sample_implicit',       vtk.vtkContourFilter),
    ('transform',             vtk.vtkTransformPolyDataFilter),
    ('flip_normals',          vtk.vtkReverseSense),
    ('clean',                 vtk.vtkCleanPolyData),
    ('fill_holes',            vtk.vtkFillHolesFilter),
    ('triangle',              vtk.vtkTriangleFilter),
    ('decimate',              vtk.vtkDecimatePro),
    ('quadric_decimation',    vtk.vtkQuadricDecimation),
    ('quadric_clustering',    vtk.vtkQuadricClustering),
    ('linear_subdivision',    vtk.vtkLinearSubdivisionFilter),
    ('loop_subdivision',      vtk.vtkLoopSubdivisionFilter),
    ('butterfly_subdivision', vtk.vtkButterflySubdivisionFilter),
    ('smooth',                vtk.vtkSmoothPolyDataFilter),
    ('windowed_sinc',         vtk.vtkWindowedSincPolyDataFilter),
    ('reverse_sense',         vtk.vtkReverseSense),
    ])

DEFAULT_FILTER_PARAMS = {
    'linestopoints':        True,
    'polystolines':         True,
    'stripstopolys':        True,
    'maximumholesize':      1.0,
    'processvertices':      True,
    'processlines':         True,
    'targetreduction':      0.2,
    'preservetopology':     True,
    'splitmesh':            True,
    'deletevertices':       False,
    'divisionsx':           10,
    'divisionsy':           10,
    'divisionsz':           10,
    'autoadjustdivisions':  True,
    'visible':              True,
    'relaxation':           0.01,
    'iterations':           20,
    'boundarysmoothing':    True,
    'featureangle':         45.0,
    'featureedgesmoothing': False,
    'edgeangle':            15.0,
    'passband':             0.1,
    'manifoldsmoothing':    False,
    'normalize':            False,
    'reversecells':         False,
    'reversenormals':       True,
    'flipnormals':          True,
    'rotationx':            0.0,
    'rotationy':            0.0,
    'rotationz':            0.0,
    'scalex':               1.0,
    'scaley':               1.0,
    'scalez':               1.0,
    'translatex':           0.0,
    'translatey':           0.0,
    'translatez':           0.0,
    'samplesx':             40,
    'samplesy':             40,
    'samplesz':             40,
    'sample':               None,
    'geo_type':             'filter',
    'type':                 ''}

DEFAULT_BOOLEAN_PARAMS = {
    'children': [],
    'visible':  True,
    'type':     '',
    'geo_type': 'boolean',
    'bounds':   [-1, 1, -1, 1, -1, 1],
    'sample':   None,
    }

DEFAULT_STL_PARAMS = {
    'type':            'stl',
    'filename':        None,
    'centerx':         0.0,
    'centery':         0.0,
    'centerz':         0.0,
    'rotationx':       0.0,
    'rotationy':       0.0,
    'rotationz':       0.0,
    'translationx':    0.0,
    'translationy':    0.0,
    'translationz':    0.0,
    'visible':         True,
    'extentxmin':      0.0,
    'extentxmax':      0.0,
    'extentymin':      0.0,
    'extentymax':      0.0,
    'extentzmin':      0.0,
    'extentzmax':      0.0,
    'scale':           1.0,
    'aboutorigin':     False,
    'units':           'm',
    'geo_type':        'stl',
    'flipnormals':     False,
    }

DEFAULT_PARAMS = {
    'primitive':        DEFAULT_PRIMITIVE_PARAMS,
    'parametric':       DEFAULT_PARAMETRIC_PARAMS,
    'procedural':       DEFAULT_PROCEDURAL_PARAMS,
    'filter':           DEFAULT_FILTER_PARAMS,
    'boolean':          DEFAULT_BOOLEAN_PARAMS,
    'boolean_implicit': DEFAULT_BOOLEAN_PARAMS,
    'stl':              DEFAULT_STL_PARAMS,
    'implicit':         DEFAULT_IMPLICIT_PARAMS,
}

DEFAULT_TEXT_COLOR = QColor(0, 0, 0)
DEFAULT_GEO_COLOR = QColor(224, 224, 224)

DEFAULT_VISUAL_PROPS = {
    'mesh':            {'color': QColor(244,  67,  40), 'visible': False,
                        'opacity': 1, 'rep': 'solid', 'color_by': 'Volume'},
    'background_mesh': {'color': QColor(100, 182, 247), 'visible': True,
                        'opacity': 1, 'rep': 'wire'},
    'geometry':        {'color': DEFAULT_GEO_COLOR, 'visible': True,
                        'opacity': 0.3, 'rep': 'solid'},
    'regions':         {'color': QColor(224, 224, 224), 'visible': False,
                        'opacity': 0.5, 'rep': 'solid', 'linewidth': 1.0},
    'normals':         {'color': QColor(0, 0, 224), 'visible': False,
                        'scale': 0.1, 'count': 1000},
    'axes':            {'color': QColor(0, 0, 224), 'visible': True},
    'boundary':        {'color': QColor(0, 0, 224), 'opacity': 1,
                        'rep': 'edges', 'visible': False, 'color_by': 'bc_id'},
    'particles':       {'color': QColor(214, 215, 176), 'visible': True, 'opacity':1},
    }

# add edge color
for geo in list(DEFAULT_VISUAL_PROPS.keys()):
    DEFAULT_VISUAL_PROPS[geo]['edge'] = DEFAULT_VISUAL_PROPS[geo]['color'].darker()

# Default color_bar options
DEFAULT_COLOR_BAR = {
    'position': 'right',
    'n_labels': 11,
    'label_fmt': '%.2f',
    'color': DEFAULT_TEXT_COLOR,
    'shadow': False,
    'opacity': 1.0,
    'title_size': 12,
    'label_size': 12,
}


# Splat shader code for the vtkPointGaussianMapper
# from Paraview:
# https://gitlab.kitware.com/paraview/paraview/blob/master/ParaViewCore/ClientServerCore/Rendering/vtkPointGaussianRepresentation.cxx
SPLAT_SHADERS = {
    'sphere':
        "//VTK::Color::Impl\n"
        "float dist = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);\n"
        "if (dist > 1.0) {\n"
        "  discard;\n"
        "} else {\n"
        "  float scale = (1.0 - dist);\n"
        "  ambientColor *= scale;\n"
        "  diffuseColor *= scale;\n"
        "}\n",
    'black_edge_circle':
        "//VTK::Color::Impl\n"
        "float dist = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);\n"
        "if (dist > 1.0) {\n"
        "  discard;\n"
        "} else if (dist > 0.8) {\n"
        "  ambientColor = vec3(0.0, 0.0, 0.0);\n"
        "  diffuseColor = vec3(0.0, 0.0, 0.0);\n"
        "}\n",
    'circle':
        "//VTK::Color::Impl\n"
        "float dist = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);\n"
        "if (dist > 1.0) {\n"
        "  discard;\n"
        "};\n",
    'triangle':
        "//VTK::Color::Impl\n",
    'square_outline':
        "//VTK::Color::Impl\n"
        "if (abs(offsetVCVSOutput.x) > 2.2 || abs(offsetVCVSOutput.y) > 2.2) {\n"
        "  discard;\n"
        "}\n"
        "if (abs(offsetVCVSOutput.x) < 1.5 && abs(offsetVCVSOutput.y) < 1.5) {\n"
        "  discard;\n"
        "}\n",
    'square':
        "//VTK::Color::Impl\n"
        "  if (abs(offsetVCVSOutput.x) > 1.0 || abs(offsetVCVSOutput.y) > 1.0) { discard; }\n",
}
