import math
import trimesh
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import vtk
from vtk.util import numpy_support


def main():

# Examples of particle seeding

# Example 1: STL file (adjust file file name and path)
    lattice_df1 = seed_particles(geo_type='stl', geo_input='~/suzanne.stl',
                                lattice=simple_cubic,
                                reverse = True,
                                diameter=0.2,
                                max_wall_overlap = 0.0,
                                shift=[0.0, 0.0, 0.0],
                                spacing=[0.2, 0.2, 0.2])

# Example 2 Superquadric shape
    lattice_df2 = seed_particles(geo_type='sqp', geo_input=[1.0, 2.0, 3.0, 2, 2],
                                lattice=body_centered_cubic,
                                reverse = False,
                                diameter=0.5,
                                max_wall_overlap=0.0,
                                shift=[0.0, 0.0, 0.0],
                                spacing=[0.5, 0.5, 0.5])

# Assign one of the lattices from the examples above
    lattice_df = lattice_df1


    plot_lattice(lattice_df) # scatter plot (optional)

    PolyData = df_to_polydata(lattice_df)

    write_polydata_to_vtp(PolyData, 'out.vtp', binary=True)




def simple_cubic(origin=[0.0, 0.0, 0.0], end=[1.0, 1.0, 1.0], spacing=[1.0, 1.0, 1.0]):
    """ Generates points arranged as a simple cubic lattice.
        Returns a numpy array of points, each point is a [x,y,z] array

    Arguments:
        origin       : Coordinates [x,y,z] of lattice lower bound corner
        end          : Coordinates [x,y,z] of lattice upper bound corner
        spacing      : Spacing in [x,y,z] direction (lattice unit cell size)

    Returns:
        Numpy array [[x,y,z]] coordinates
    """

    lattice = []

#Number of seeds in x,y,z direction
    seeds = np.ceil((end - origin)/spacing).astype(int) + 1

# Outer loop is j, which usually correspond to the vertical direction in MFiX
# TODO: Add option to permute i,j,k

    for j in range(0, seeds[1]):
        for i in range(0, seeds[0]):
            for k in range(0, seeds[2]):
                x = origin[0] + i * spacing[0]
                y = origin[1] + j * spacing[1]
                z = origin[2] + k * spacing[2]

                lattice.append((x, y, z))

    return np.array(lattice)




def body_centered_cubic(origin=[0.0, 0.0, 0.0], end=[1.0, 1.0, 1.0], spacing=[1.0, 1.0, 1.0]):
    """ Generates points arranged as a body-centered cubic lattice.
        Returns a numpy array of points, each point is a [x,y,z] array

    Arguments:
        origin       : Coordinates [x,y,z] of lattice lower bound corner
        end          : Coordinates [x,y,z] of lattice upper bound corner
        spacing      : Spacing in [x,y,z] direction (lattice unit cell size)

    Returns:
        Numpy array [[x,y,z]] coordinates
    """

    lattice = []

#Number of seeds in x,y,z direction
    seeds = np.ceil((end - origin)/spacing).astype(int) + 1

# Outer loop is j, which usually correspond to the vertical direction in MFiX
# TODO: Add option to permute i,j,k
    for j in range(0, seeds[1]):
        for i in range(0, seeds[0]):
            for k in range(0, seeds[2]):
                x = origin[0] + i * spacing[0]
                y = origin[1] + j * spacing[1]
                z = origin[2] + k * spacing[2]

                lattice.append((x, y, z))

                x = origin[0] + (i + 0.5) * spacing[0]
                y = origin[1] + (j + 0.5) * spacing[1]
                z = origin[2] + (k + 0.5) * spacing[2]

                lattice.append((x, y, z))


    return np.array(lattice)

def hexagonal(origin=[0.0, 0.0, 0.0], end=[1.0, 1.0, 1.0], spacing=[1.0, 1.0, 1.0]):
    """ Generates points arranged as an hexagonal lattice.
        Returns a numpy array of points, each point is a [x,y,z] array

    Arguments:
        origin       : Coordinates [x,y,z] of lattice lower bound corner
        end          : Coordinates [x,y,z] of lattice upper bound corner
        spacing      : Spacing in [x,y,z] direction (lattice unit cell size)

    Returns:
        Numpy array [[x,y,z]] coordinates
    """

    lattice = []
    seeds = [1,1,1]

    sq3 = math.sqrt(3.0)
    sq6 = math.sqrt(6.0)
    sq6o3 = sq6/3.0
    sq3o2 = sq3/2.0

#Number of seeds in x,y,z direction
    span = end - origin

    seeds[0] = math.floor(span[0]/(spacing[0]*sq3))*2+1
    seeds[1] = math.floor(span[1]/(spacing[1]*sq3o2))+1
    seeds[2] = math.floor(span[2]/(spacing[2]*1.5))*2+1

# Outer loop is j, which usually correspond to the vertical direction in MFiX
# TODO: Add option to permute i,j,k
    for j in range(0, seeds[1]):
        for k in range(0, seeds[2]):
            for i in range(0, seeds[0]):

                x = origin[0] + (2*i +(j+k)%2)      * spacing[0]/2
                y = origin[1] + (2*j*sq6o3)         * spacing[1]/2
                z = origin[2] + (k + (j%2)/3.0)*sq3 * spacing[2]/2

                lattice.append((x, y, z))

    return np.array(lattice)


def seed_particles(geo_type='sqp', geo_input=[1.0, 2.0, 3.0, 2, 2], lattice=simple_cubic, reverse = False, diameter=0.1, max_wall_overlap=0.0, shift=[0.0, 0.0, 0.0], spacing=[1.0, 1.0, 1.0]):
    """ Seeds points in a Superquadric shape.
        Returns a Panda dataframe

    Arguments:
        geo_type           : Geometry input type:
                                 - 'stl' for an STL file
                                 - 'sqp' for a superquadric shape
        geo_input          : Geometry input:
                                 - stl filename if geo_type='stl'
                                 - Superquadric parameters [a, b, c, m, n] if geo_type='sqp' (same definition as in MFiX)
        lattice            : Type of lattice arrangement
                                 - simple_cubic
                                 - body_centered_cubic
                                 - hexagonal
        reverse            : Option to seed particles outside the geometry (reverse= True) instead of
                             seeding inside the geometry (reverse= False, default value)
        diameter           : Particle diameter
        max_wall_overlap   : Maximum overlap between particle's surface and the geometry surface (wall)
                             max_wall_overlap=0.0: Particles do not overlap the geometry wall, but can potentially touch the wall
                             max_wall_overlap>0.0: Particles can overlap the wall as much as max_wall_overlap
                             max_wall_overlap<0.0: Particles do not overlap the wall and are at least max_wall_overlap away from it
        shift              : Shift in [x, y, z] direction applied to origin of the lattice. This can be used to nudge the particles location for better results
        spacing            : Spacing in [x,y,z] direction (lattice unit cell size)

    Returns:
        Panda dataframe
    """

# Setup corner points of the lattice
    if geo_type=='stl':
        mesh   = trimesh.load(geo_input)
        origin = mesh.bounding_box.bounds[0] - max_wall_overlap + shift
        end    = mesh.bounding_box.bounds[1] + max_wall_overlap + shift

    elif geo_type=='sqp':
        origin = -np.array(geo_input[0:3]) - max_wall_overlap + shift
        end    =  np.array(geo_input[0:3]) + max_wall_overlap + shift

    else:
        raise ValueError("Unknown geometry input type:",geo_type," . Must be either 'stl' or 'sqp' ")


# Generate lattice
    xyz = lattice(origin=origin, end=end, spacing=spacing)

# Determine if point is inside or outside geometry
    if geo_type=='stl':
        sdf = trimesh.proximity.signed_distance(mesh,xyz)
    elif geo_type=='sqp':
        sdf = sqp_eval(xyz, geo_input)

# build a dataframe
    df = pd.DataFrame(xyz, columns = ['X','Y','Z'])
    df['SDF'] = sdf
    df['Diameter'] = diameter

# Keep points either inside or outside superquadric shape
# TODO: Need to modify the test for sqp
    if reverse:
        df_new = df[df['SDF'] > diameter/2 - max_wall_overlap]
    else:
        df_new = df[df['SDF'] < max_wall_overlap - diameter/2]

    return df_new


def sqp_eval(xyz, sqp_param):
    """ Evaluate superquadric at Sweeps a set of points along a curve.
        Returns a list of values (sqpf).
        A negative value of sqpf means the point is inside the superquadric shape.
        A positive value of sqpf means the point is outside the superquadric shape.
        A zero value of sqpf means the point is exactly on the surface of the superquadric shape.
        Note: This is not a dimensional signed distance

    Arguments:
        xyz          : Numpy array of points, each point is a [x,y,z] array
        sqp_param    : Superquadric parameters [a, b, c, m, n]

    Returns:
        Numpy array of the sqp function at each xyz location
    """

    a, b, c, m, n = sqp_param

    sqpf = []
    for x,y,z in xyz:
        s = pow(pow(x/a,m) + pow(y/b,m),n/m) + pow(z/c,n) - 1
        sqpf.append(s)

    return sqpf


def df_to_polydata(df):
    """ Converts a Panda dataframe into vtk PolyData

    Arguments:
        df  : A Panda dataframe, assumed to have at least 'X', 'Y', and 'Z' columns

    Returns:
        A PolyData object
    """

    # Create a vtk object so we can save the file
    Points = vtk.vtkPoints()
    PolyData  = vtk.vtkPolyData()

# point coordinates
    for x,y,z in zip(df['X'],df['Y'],df['Z']):
        Points.InsertNextPoint([x,y,z])

    PolyData.SetPoints(Points)

# Add point data from each column of the dataframe, except the coordinates
    for column in df:
        if column not in ['X', 'Y', 'Z']:
            vtk_array = numpy_support.numpy_to_vtk(df[column])
            vtk_array.SetName(column)
            PolyData.GetPointData().AddArray(vtk_array)

    return PolyData


#Save the file (vtp format)
def write_polydata_to_vtp(PolyData, filename, binary=True):
    """ Writes a PolyData object into a vtp file

    Arguments:
        PolyData  : A PolyData object
        filename  : Name of the vtp file
        binary:   : Option to save file in binary format (binary=True) or ASCII format (binary=False)

    Returns:
        Nothing. Saved file to disc
    """

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    if binary:
        writer.SetDataModeToBinary()
    else:
        writer.SetDataModeToAscii()
    writer.SetInputData(PolyData)
    writer.Update()


# Simple scatter plot
def plot_lattice(df):

    plt.rcParams["figure.figsize"] = [7.00, 3.50]
    plt.rcParams["figure.autolayout"] = True
    fig, ax = plt.subplots(1, 1, subplot_kw={'projection': '3d'})

    ax.scatter(df['X'], df['Y'], df['Z'], c=df['SDF'], alpha=1)
    plt.show()


if __name__ == '__main__':
    main()
