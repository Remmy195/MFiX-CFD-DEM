from pathlib import Path

from mfixgui.vtk_widgets.pvd_reader import PvdReader, read_griddata, read_partdata
from mfixgui.vtk_widgets.histogram_viewer import HistogramViewer

TEST_VTU = Path(__file__).parent / 'test.vtu'
TEST_VTP = Path(__file__).parent / 'test.vtp'
TEST_PVD = Path(__file__).parent / 'test.pvd'


def test_grid_reader():
    gd = read_griddata(TEST_VTU.as_posix())

    assert gd.array_info.keys() == set(
        ['Volume', 'Relative_Volume', 'Aspect_Ratio', 'Delta_h']
    )
    assert gd.array_info['Volume'].components == 1
    assert gd.array_info['Volume'].magnitude == (
        4.770579948853992e-07,
        1.7313333955826238e-05,
    )
    assert gd.array_info['Volume'].range_ == [
        (4.770579948853992e-07, 1.7313333955826238e-05)
    ]

    assert gd.array_info['Relative_Volume'].components == 1
    assert gd.array_info['Relative_Volume'].range_ == [(0.02755437046289444, 1.0)]
    assert gd.array_info['Relative_Volume'].magnitude == (0.02755437046289444, 1.0)
    assert gd.array_info['Aspect_Ratio'].components == 1
    assert gd.array_info['Aspect_Ratio'].range_ == [
        (1.2231849431991577, 3.852201223373413)
    ]
    assert gd.array_info['Aspect_Ratio'].magnitude == (
        1.2231849431991577,
        3.852201223373413,
    )
    assert gd.array_info['Delta_h'].components == 1
    assert gd.array_info['Delta_h'].range_ == [(0.0, 0.022193606942892075)]
    assert gd.array_info['Delta_h'].magnitude == (0.0, 0.022193606942892075)


def test_get_partdata():
    pd = read_partdata(TEST_VTP.as_posix())

    assert pd.array_info.keys() == set(['Diameter', 'Velocity'])
    assert pd.array_info['Diameter'].magnitude == (
        0.009999999776482582,
        0.009999999776482582,
    )
    assert pd.array_info['Diameter'].number_of_tuples == 1977
    assert pd.array_info['Diameter'].range_ == [
        (0.009999999776482582, 0.009999999776482582)
    ]

    assert pd.array_info['Velocity'].magnitude == (
        0.4543472760838399,
        0.7534805024072301,
    )
    assert pd.array_info['Velocity'].number_of_tuples == 1977
    assert pd.array_info['Velocity'].range_ == [
        (-0.19630415737628937, 0.1985410749912262),
        (-0.7436310648918152, -0.4476633667945862),
        (-8.981057542621629e-13, 6.18555987680347e-13),
    ]


def test_find_pvds(qtbot):
    pvd_dir = PvdReader()
    assert len(pvd_dir.pvd_filenames()) == 0

    pvd_dir.refresh(TEST_PVD.parent)
    assert len(pvd_dir.pvd_filenames()) == 1
    assert pvd_dir.pvd('test').type_ == "vtp"
    assert pvd_dir.pvd('test').series[0] == "HOPPER_P_0000.vtp"
    assert pvd_dir.pvd('test').series[0.01044444] == "HOPPER_P_0001.vtp"
    assert pvd_dir.pvd('test').series[0.02007407] == "HOPPER_P_0002.vtp"
