from dataclasses import dataclass
from pathlib import Path

from mfixgui.tools import safe_int
from mfixgui.vtk_widgets.tools import parse_pvd_file
from qtpy.QtCore import QTimer, Qt, QObject, Signal

from vtk.util.numpy_support import vtk_to_numpy
import vtk
import numpy as np

def find_nearest(array, value):
    array = np.asarray(array)
    return (np.abs(array - value)).argmin()


class PvdReader(QObject):
    """ Mapping of pvd basenames -> pvd files in project_dir """
    signal_new_files = Signal(list)
    signal_removed_files = Signal(list)
    def __init__(self):
        QObject.__init__(self)
        self._project_dir = None
        self._pvd_files = {}
        self._pvd_index = -1

        # look for new vtk files
        self.file_timer = QTimer()
        self.file_timer.timeout.connect(self.refresh)

    def __iter__(self):
        return self

    def __next__(self):
        if self._pvd_index < -1 or self._pvd_index >= len(self._pvd_files)-1:
            self._pvd_index = -1
            raise StopIteration
        self._pvd_index = self._pvd_index + 1
        k = list(self._pvd_files.keys())[self._pvd_index]
        return self._pvd_files[k]

    def set_project_dir(self, project_dir):
        self._project_dir = Path(project_dir)
        self._pvd_files = {}
        self.refresh()
        if not self.file_timer.isActive():
            self.file_timer.start(1000)

    def pvd_filenames(self):
        return self._pvd_files.keys()

    def pvd(self, filename):
        return self._pvd_files.get(filename, None)

    def refresh(self):
        new_pvd_files = []
        if self._project_dir is not None:
            for pvd_filename in self._project_dir.glob('*.pvd'):
                if pvd_filename.stem not in self._pvd_files:
                    # we could read this file before it is completely written,
                    # catch all errors
                    try:
                        self._pvd_files[pvd_filename.stem] = self._make_pvdfile(pvd_filename)
                        new_pvd_files.append(pvd_filename.stem)
                    except:
                        continue
                else:
                    u = self._pvd_files[pvd_filename.stem].refresh()
                    if u:
                        new_pvd_files.append(pvd_filename.stem)

        # check for deleted files
        rm_list = []
        for key, value in self._pvd_files.items():
            if not value.pvdfilename.exists():
                rm_list.append(key)
        for key in rm_list:
            self._pvd_files.pop(key)

        if new_pvd_files:
            self.signal_new_files.emit(new_pvd_files)
        if rm_list:
            self.signal_removed_files.emit(rm_list)

    def _make_pvdfile(self, pvdfilename):
        times_to_fnames = parse_pvd_file(pvdfilename)
        type_ = None
        fnames = times_to_fnames.values()
        if all(fname.endswith('vtu') for fname in fnames):
            type_ = 'vtu'
        elif all(fname.endswith('vtp') for fname in fnames):
            type_ = 'vtp'
        else:
            print("Invalid PVD file!  Not all VTU or all VTP: %s" % pvdfilename)
        return PvdFile(self._project_dir, pvdfilename, times_to_fnames, type_)


class PvdFile:
    """ Mapping of time -> (VTU/VTP) filename """

    def __init__(self, project_dir, pvdfilename, series, type_):
        self.project_dir = project_dir
        self.pvdfilename = pvdfilename
        self.series = series
        self.type_ = type_
        self.array_info_cache = {}
        self.time_series = list(self.series.keys())
        ts = np.asarray(self.time_series)
        self.dt = np.mean(ts[1:] - ts[:-1])

        # init the reader
        time = self.time_series[0]
        path = self.project_dir / self.series[time]
        if self.type_ == 'vtp':
            self.data = read_partdata(path)

        if self.type_ == 'vtu':
            self.data = read_griddata(path)

    @property
    def reader(self):
        return self.data.reader

    def refresh(self):
        update = False
        self.series = parse_pvd_file(self.pvdfilename)
        if len(self.series) != len(self.time_series):
            update = True
        self.time_series = list(self.series.keys())
        ts = np.asarray(self.time_series)
        self.dt = np.mean(ts[1:] - ts[:-1])
        return update

    def change_time(self, time):
        if not self.time_series:
            return None
        time = find_nearest(self.time_series, time)
        return self.change_index(time)

    def change_index(self, index):
        time = list(self.series.keys())[index]

        path = self.project_dir / self.series[time]
        if not path.exists():
            #print("File does not exist: %s" % path)
            return None
        self.data.reader.SetFileName(str(path))
        self.data.reader.Update()
        return self.data

    def get_arrayinfo(self):
        return self.data.get_arrayinfo()

    def max_index(self):
        return len(self.series) - 1

    def histogram(self, var, component, bins_value):
        return self.data.histogram(var, component, bins_value)


def read_griddata(path):
    ugrid_reader = (
        vtk.vtkXMLPUnstructuredGridReader()
        if path.suffix == "pvtu"
        else vtk.vtkXMLUnstructuredGridReader()
    )
    ugrid_reader.SetFileName(path.as_posix())
    ugrid_reader.Update()
    return GridData(ugrid_reader)


def read_partdata(path):
    particle_reader = (
        vtk.vtkXMLPPolyDataReader()
        if path.suffix == ".pvtp"
        else vtk.vtkXMLPolyDataReader()
    )
    particle_reader.SetFileName(path.as_posix())
    particle_reader.Update()
    return PartData(particle_reader)


@dataclass(frozen=True)
class Histogram:
    counts: dict
    edges: dict


@dataclass(frozen=True)
class PartData:
    reader: vtk.vtkObject

    def histogram(self, var, component, bins_value):
        data = self.reader.GetOutput()
        point_data = data.GetPointData()
        data = point_data.GetAbstractArray(var)
        return make_histogram(data, component, bins_value)

    def get_arrayinfo(self):
        data = self.reader.GetOutput()
        point_data = data.GetPointData()
        array_info = {}
        n_tuples = None
        for i in range(point_data.GetNumberOfArrays()):
            array = point_data.GetArray(i)
            n_comp = array.GetNumberOfComponents()
            n_tuples = array.GetNumberOfTuples()

            name = point_data.GetArrayName(i)
            array_info[name] = ParticleArray(
                number_of_tuples=n_tuples,
                components=n_comp,
                range_=[array.GetRange(i) for i in range(n_comp)],
                magnitude=array.GetRange(-1),
            )

        # read field data
        field_data_info = {}
        field_data = data.GetFieldData()
        for i in range(field_data.GetNumberOfArrays()):
            name = field_data.GetArrayName(i)
            if 'SQP' in name:
                name = safe_int(name.split('_')[-1], None)
                if name is None:
                    continue
            array = field_data.GetArray(i)
            field_data_info[name] = vtk_to_numpy(array)
        array_info['field_data'] = field_data_info
        return array_info


@dataclass(frozen=True)
class GridData:
    reader: vtk.vtkObject

    def histogram(self, var, component, bins_value):
        data = self.reader.GetOutput()
        cell_data = data.GetCellData()
        data = cell_data.GetAbstractArray(var)
        return make_histogram(data, component, bins_value)

    def get_arrayinfo(self):
        data = self.reader.GetOutput()

        cell_data = data.GetCellData()
        array_info = {}
        for i in range(cell_data.GetNumberOfArrays()):
            array = cell_data.GetArray(i)
            n_comp = array.GetNumberOfComponents()
            name = cell_data.GetArrayName(i)
            array_info[name] = GridArray(
                magnitude=array.GetRange(-1),
                components=n_comp,
                range_=[array.GetRange(i) for i in range(n_comp)],
            )
        return array_info


@dataclass(frozen=True)
class GridArray:
    magnitude: list
    components: int
    range_: list


def make_histogram(data, component, bins_value):
    if data is None:
        return Histogram(np.asarray([0]), np.asarray([0, 1]))
    array = vtk_to_numpy(data)
    if len(array.shape) > 1:
        if component == 'mag':
            array = np.sqrt((array * array).sum(axis=1))
        else:
            array = array[:, ['x', 'y', 'z'].index(component)]
    counts, edges = np.histogram(array, bins_value)
    return Histogram(counts, edges)


@dataclass(frozen=True)
class ParticleArray:
    number_of_tuples: int
    components: int
    range_: list
    magnitude: list


PVD_READER = None
def get_pvd_reader():
    # construct a singleton pvd reader to act as a server
    global PVD_READER
    if not PVD_READER:
        PVD_READER = PvdReader()
    return PVD_READER
