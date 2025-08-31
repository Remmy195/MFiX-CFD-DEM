#!/usr/bin/env python

"""Glued sphere particle generator"""

import errno
import math
import numpy as np
import os
import pandas as pd
import random
import shutil
import sys
import vtk
from vtk.util import numpy_support

verbose = False
def dprint(*args):
    if verbose:
        print(*args)

progress_callback = None
def set_progress_callback(func):
    global progress_callback
    progress_callback = func

stop_flag = False
def set_stop_flag(val=True):
    global stop_flag
    stop_flag = val

def is_cube_inside_stl(cubePolyData, stlPolyData):
    # Create a vtkSelectEnclosedPoints filter
    enclosedPoints = vtk.vtkSelectEnclosedPoints()
    enclosedPoints.SetInputData(cubePolyData)
    enclosedPoints.SetSurfaceData(stlPolyData)
    enclosedPoints.Update()

    # Check if all points of the cube are inside the STL
    points = cubePolyData.GetPoints()
    for i in range(points.GetNumberOfPoints()):
        if not enclosedPoints.IsInside(i):
            return False
    return True

def is_cube_outside_stl(cubePolyData, stlPolyData):
    # Create a vtkSelectEnclosedPoints filter
    enclosedPoints = vtk.vtkSelectEnclosedPoints()
    enclosedPoints.SetInputData(cubePolyData)
    enclosedPoints.SetSurfaceData(stlPolyData)
    enclosedPoints.Update()

    # Check if all points of the cube are outside the STL
    points = cubePolyData.GetPoints()
    for i in range(points.GetNumberOfPoints()):
        if enclosedPoints.IsInside(i):
            return False
    return True

def decimation_vtk(input_path, output_path, reduction_rate):
    readerSTL = vtk.vtkSTLReader()
    readerSTL.SetFileName(input_path)
    readerSTL.Update()

    inputPoly = readerSTL.GetOutput()
    dprint("Before decimation\n"
          "-----------------\n"
          "There are " + str(inputPoly.GetNumberOfPoints()) + "nodes.\n"
          "There are " + str(inputPoly.GetNumberOfPolys()) + "elements.\n")

    decimate = vtk.vtkQuadricDecimation()
    decimate.SetInputData(inputPoly)
    decimate.SetTargetReduction(1-reduction_rate)
    decimate.Update()

    decimatedPoly = vtk.vtkPolyData()
    decimatedPoly.ShallowCopy(decimate.GetOutput())

    dprint("After decimation \n"
          "-----------------\n"
          "There are " + str(decimatedPoly.GetNumberOfPoints()) + "nodes.\n"
          "There are " + str(decimatedPoly.GetNumberOfPolys()) + "elements.\n")

    stlWriter = vtk.vtkSTLWriter()
    stlWriter.SetFileName(output_path)
    stlWriter.SetFileTypeToBinary()
    stlWriter.SetInputConnection(decimate.GetOutputPort())
    stlWriter.Write()

def calculate_signed_distance(clonedPolyData, Points):
# calculate total signed distance of a surface from each point in surface
#   if <0 means it is inside the geo
#      >0 means outside of the geo
#      =0 means the point is on the geo surface

    total_sdf = 0
    sdt = 0
    vertexsdf_list = []
    implicitPolyDataDistance = vtk.vtkImplicitPolyDataDistance()
    implicitPolyDataDistance.SetInput(clonedPolyData)
    for point in Points:
        sdt = implicitPolyDataDistance.EvaluateFunction(point)
        vertexsdf_list.append(sdt)
        total_sdf += sdt
    return total_sdf, np.array(vertexsdf_list)

# create dir if nonexist, otherwise empty dir and rewrite
def dir_create(dir_path):
    try:
        shutil.rmtree(dir_path)
    except Exception as e:
        if e.errno != errno.ENOENT:
            raise
    os.makedirs(dir_path)


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
            # Ensure the column data type is supported by VTK
            if np.issubdtype(df[column].dtype, np.number):
                # Convert the column to numpy array if it's a supported numeric type
                vtk_array = numpy_support.numpy_to_vtk(df[column].values)
            else:
                # Otherwise, convert the column to a supported float type
                vtk_array = numpy_support.numpy_to_vtk(df[column].astype(float).values)
            vtk_array.SetName(column)
            PolyData.GetPointData().AddArray(vtk_array)

    return PolyData

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
def create_gsp_config_header(filename, dimension, phase_id, spheres_per_particle):
    """
    Creates a .dat file with a specified format.

    :param filename: Name of the file to be created.
    :param dimension: Dimension of the model (2D or 3D).
    :param phase_id: Phase ID of the solid phases.
    :param spheres_per_particle: Number of spheres in each glued-sphere particle for each phase.
    """
    header = [
        "Version: 1.0\n",
        "======================================================================================================\n",
        "Instructions:\n",
        "Dimension: Enter \"2\" for 2D, \"3\" for 3D (Integer)\n",
        "Phase_ID: ID of solid phases using glued sphere model. All others are assumed sphere.\n",
        "Spheres_per_particle: Number of spheres in each glued-sphere particle.\n",
        "Coordinates are (X,Y) in 2D, (X,Y,Z) in 3D. These are the relative positions of spheres to its host \n",
        "non-spherical particle. Also, they are not the absolute distance and are the ratio between absolute \n",
        "distance and diameter of non-spherical particle.\n",
        "Phase_ID and Relative Diameter are scalars\n",
        "Relative Diameter is the ratio of the sphere diameter to the diameter of its host non-spherical particle\n",
        "Data starts from all the spheres in the first phase, then the second phase. \n",
        "Data starts at line 21. Each column corresponds to a variable.\n",
        "======================================================================================================\n",
        f"Dimension: {dimension}\n",
        f"Phase_ID: {' '.join(map(str, phase_id))}\n",
        f"Spheres_per_particle: {' '.join(map(str, spheres_per_particle))}\n",
        "======================================================================================================\n",
        "| X | Y | Z | Relative Diameter | surfaceArea | W | E | S | N | Top | Bot | Wa | Ea | Sa | Na | Ta | Ba \n",
        "======================================================================================================\n"
    ]

    with open(filename, 'w') as file:
        file.writelines(header)

def create_dat_file_with_header(filename, dimension, phase_id, spheres_per_particle, TotalComponentSpheres, isTp = False):
    """
    Creates a .dat file with a specified format.

    :param filename: Name of the file to be created.
    :param dimension: Dimension of the model (2D or 3D).
    :param phase_id: Phase ID of the solid phases.
    :param spheres_per_particle: Number of spheres in each non-spherical particle for each phase.
    """
    header = [
        "Version 2.0\n",
        "======================================================================================================\n",
        "Instructions: \n",
        "Dimension: Enter \"2\" for 2D, \"3\" for 3D (Integer) \n",
        "Particles: Number of particles to read from this file (Integer) \n",
        "For each variable, enter whether it will be read from the file \n",
        "(\"T\" to be read (True), \"F\" to not be read (False)). When not read, enter the \n",
        "default values assign to all particles. \n",
        "Coordinates are always read (X,Y in 2D, X,Y,Z in 3D) \n",
        "Phase_ID, Diameter, Density, Temperature are scalars \n",
        "Velocity requires U,V in 2D, U,V,W in 3D \n",
        "Temperature is only read/set if the energy equation is solved \n",
        "Species are only read/set if the species equations are solved, and requires \n",
        "all species to be set (if phase 1 has 2 species and phase 2 has 4 species, \n",
        "all particles must have 4 entries, even phase 1 particles (pad list with zeros) \n",
        "User scalars need DES_USR_VAR_SIZE values \n",
        "Data starts at line 35. Each column corresponds to a variable. \n",
        "======================================================================================================\n",
        f"Dimension:                                              {dimension}\n",
        f"TotalComponentSpheres:                                  {TotalComponentSpheres}\n",
        "======================================================================================================\n",
        "Variable      Read(T/F)     Default value (when Read=F) \n",
        "======================================================================================================\n",
        "Coordinates:     T          Must always be T \n", # vector
        "Phase_ID:        T          1.0 \n", # scalar
        "Diameter:        T          0.0 \n", # scalar
        "Density:         T          0.0 \n", # scalar
        "Velocity:        T          0.0 0.0 0.0 \n", # vector
        "Quaternion:      T          1.0 0.0 0.0 0.0 \n", #vector
        "Temperature:     T          (Ignored if energy eq. is not solved) \n", # scalar
        "Species:         T          (Ignored if species eq. are not solved) \n", # vecotor
        "User_Scalar:     T          (Ignored if no user scalars are defined) \n", # vector
        "======================================================================================================\n",
        "| X | Y | Z | Phase_ID | Bounding_Diameter | Density | U | V | W | q1 | q2 | q3 | q4\n",
        "======================================================================================================\n"
    ]
    if isTp:
        header[-2] = header[-2].strip() + " | Temperature \n"
    with open(filename, 'w') as file:
        file.writelines(header)

def append_to_dat_file(filename, data):
    """
    Appends data to an existing .dat file.

    :param filename: Name of the file to append data to.
    :param data: Data to append (DataFrame, list of strings, or a single string).
    """
    with open(filename, 'a') as file:
        if isinstance(data, pd.DataFrame):
            data_string = data.to_string(header=False, index=False)
            file.write(data_string)
            if not data_string.endswith(os.linesep):
                file.write(os.linesep)
        elif isinstance(data, list):
            file.writelines(data)
        else:
            file.write(data)
            if not data.endswith(os.linesep):
                file.write(os.linesep)


def remove_tiny_component_spheres(df_usr0, df_particle, critical_vol=0.01, critical_totalvol=0.99):
    dprint(f"Start to remove tiny component spheres with volume fraction lower than {critical_vol}")
    dprint(f"The total volume fraction remained should not be lower than {critical_totalvol}")
    Nbefore = df_usr0.shape[0]
    df_particle["v_frac"] = 0.0
    lastcol = df_particle.shape[1]
    for i in range(0, df_particle.shape[0]):
        cur_r = df_particle.iloc[i,3]
        df_particle.iloc[i, lastcol-1] = (4/3)*math.pi*(cur_r**3)
    total_volume = df_particle["v_frac"].sum()
    max_vlolume = df_particle["v_frac"].max()
    for i in range(0, df_particle.shape[0]):
        df_particle.iloc[i, lastcol-1] = (df_particle.iloc[i, lastcol-1]/max_vlolume)
    indices = df_particle.index[df_particle['v_frac'] < critical_vol].tolist()

    # sort the indices based on the value of volume fraction, always wish to remove the most tiny sphere first
    sorted_indices = df_particle.loc[indices, 'v_frac'].sort_values().index.tolist()
    if sorted_indices:
        for k in range(len(sorted_indices)):
            indice = sorted_indices[k]
            # first, transfer the neigh_sa into surfaceArea
            # by loop over 2~7 column the 'indice' row
            for j in range(2,8):
                if df_usr0.iloc[indice,j] > 0:
                    cur_neighbor = df_usr0.iloc[indice,j]
                    cur_neighbor_sa = df_usr0.iloc[indice,j+6]
                    cur_neighbor_indice = cur_neighbor - 1
                    df_usr0.iloc[cur_neighbor_indice,0] = df_usr0.iloc[cur_neighbor_indice,0] + cur_neighbor_sa
                    cur_sub_row = df_usr0.iloc[cur_neighbor_indice,2:8].tolist()
                    try:
                        column_index = cur_sub_row.index(indice+1)+2
                    except ValueError:
                        sys.exit()

                    df_usr0.iloc[cur_neighbor_indice,column_index] = -1
                    df_usr0.iloc[cur_neighbor_indice,column_index+6] = 0.0

            # rebuilt the gid list. if a neigh_sid > indice+1, neigh_sid = neigh_sid + 1
            df_usr0.iloc[:, 2:8] = df_usr0.iloc[:, 2:8].map(lambda x: x - 1 if x > (indice + 1) else x)
            df_usr0 = df_usr0.drop(index=indice, axis=0)
            df_particle = df_particle.drop(index=indice, axis=0)
            df_usr0 = df_usr0.reset_index(drop=True)
            df_particle = df_particle.reset_index(drop=True)
            after_reduce_volume_frac = df_particle["v_frac"].sum()
            after_reduce_volume_frac *= max_vlolume
            after_reduce_volume_frac /= total_volume
            # if the total volume remained less than 0.99, stop removing
            if after_reduce_volume_frac <= critical_totalvol:
                break
            # after remove a specific sphere is done, update the indices in sorted_indics list then start the next removal
            sorted_indices = [x-1 if x>indice else x for x in sorted_indices]
    else:
        dprint("No components spheres are removed")
        return df_usr0, df_particle, df_usr0.shape[0]
    dprint(f"After removing tiny spheres, the total volume fraction become {after_reduce_volume_frac}")
    nsphere_after_removal = df_usr0.shape[0]
    dprint(f"The number of component spheres changes from {Nbefore} to {nsphere_after_removal}")
    return df_usr0, df_particle, nsphere_after_removal

def combined_seed(df_usr0, df_particle, df_glueposcopy, bounding_diameter,
                  Nphase, rho, isUsrVel, usrvel, isTp, Tp, isRandomQ,
                  isRemoveTinyVolume, critical_vol_fraction, critical_totalvol_fraction,
                  save_dir='example1'):
    os.makedirs(save_dir, exist_ok=True)
    df_usr0 = df_usr0.reset_index(drop=True)
    df_particle = df_particle.reset_index(drop=True)

    # if not removal flag
    nsphere_after = df_usr0.shape[0]
    # if removal flag
    if isRemoveTinyVolume:
        df_usr0, df_particle, nsphere_after = remove_tiny_component_spheres(
            df_usr0, df_particle, critical_vol_fraction, critical_totalvol_fraction)

    # format usr0_input for auto seed
    first_column = df_usr0.iloc[:, 0]
    first_column = first_column / bounding_diameter / bounding_diameter
    int_columns = df_usr0.iloc[:, 1:8] # only read 1~7 columns
    float_columns = df_usr0.iloc[:, 8:] # read the remaining columns, start from column 8
    # not include effective heat transfer distance here
    float_columns = float_columns / bounding_diameter / bounding_diameter
    first_column = first_column.map(lambda x: f"{x:.8e}")
    float_columns = float_columns.map(lambda x: f"{x:.8e}")
    df_usr0 = pd.concat([first_column, int_columns, float_columns], axis=1)

    # format relative values to bounding diameters
    df_particle["Radius"] = df_particle["Radius"] * 2
    df_particle.columns = ['diameter' if x == 'Radius' else x for x in df_particle.columns]
    df_xyzd = df_particle.iloc[:, :4]

    df_xyzd = df_xyzd / bounding_diameter
    df_rest = df_particle.iloc[:, 4:]

    df_glueposcopy['phase_id'] = Nphase
    df_glueposcopy['diameter'] = bounding_diameter
    df_glueposcopy['density'] = rho
    if isUsrVel:
        df_glueposcopy['U'] = usrvel[0]
        df_glueposcopy['V'] = usrvel[1]
        df_glueposcopy['W'] = usrvel[2]
    else:
        df_glueposcopy['U'] = 0.0
        df_glueposcopy['V'] = 0.0
        df_glueposcopy['W'] = 0.0
        for i in range(0,df_glueposcopy.shape[0]):
            xrand = random.uniform(-1,1)
            yrand = random.uniform(-1,1)
            zrand = random.uniform(-1,1)
            df_glueposcopy.iloc[i,6] = df_glueposcopy.iloc[i,6] + xrand
            df_glueposcopy.iloc[i,7] = df_glueposcopy.iloc[i,7] + yrand
            df_glueposcopy.iloc[i,8] = df_glueposcopy.iloc[i,8] + zrand

    df_glueposcopy['gsp_q1'] = 1.0
    df_glueposcopy['gsp_q2'] = 0.0
    df_glueposcopy['gsp_q3'] = 0.0
    df_glueposcopy['gsp_q4'] = 0.0

    if isRandomQ:
        par = df_glueposcopy.shape[0]
        random.uniform(-1,1)
        for i in range(0,par):
            q1 = random.uniform(-1,1)
            q2 = random.uniform(-1,1)
            q3 = random.uniform(-1,1)
            q4 = random.uniform(-1,1)

            qnorm = q1**2 + q2**2 + q3**2 + q4**2
            qnorm = math.sqrt(qnorm)
            q1 = q1/qnorm
            q2 = q2/qnorm
            q3 = q3/qnorm
            q4 = q4/qnorm

            df_glueposcopy.iloc[i,9] = q1
            df_glueposcopy.iloc[i,10] = q2
            df_glueposcopy.iloc[i,11] = q3
            df_glueposcopy.iloc[i,12] = q4
    if isTp:
        df_glueposcopy['t_s'] = Tp

    df1 = pd.concat([df_xyzd, df_usr0], axis=1)
    df2 = df_glueposcopy
    df1 = df1.drop('gid', axis=1)
    df1['pid'] = Nphase

    # write df to polydata to check the config in each phase
    dfwrite = df1.copy()
    # scale back to show the realistic config
    dfwrite['X'] *= bounding_diameter
    dfwrite['Y'] *= bounding_diameter
    dfwrite['Z'] *= bounding_diameter
    dfwrite['diameter'] *= bounding_diameter
    polydata_dfwrite= df_to_polydata(dfwrite)
    vtp_filename = os.path.join(save_dir, f'phase_{Nphase}.vtp')
    write_polydata_to_vtp(polydata_dfwrite, vtp_filename, binary=True)

    return df1, df2, nsphere_after

def write_df(global_df_config,
             global_df_particle,
             phase_id,
             spheres_per_particle,
             isTp,
             save_dir='example1'):
    dimension = 3
    particles = global_df_particle.shape[0] # this is actually the number of glued particles
    TotalComponentSpheres = 0
    for i in range(0,particles):
        tmpphase = global_df_particle.iloc[i,3]
        curN = spheres_per_particle[tmpphase-1]
        TotalComponentSpheres = TotalComponentSpheres + curN

    os.makedirs(save_dir, exist_ok=True)
    # no matter auto seeding or not, always write global_df_config
    # create_gsp_config_header(filename, dimension, phase_id, spheres_per_particle)
    global_df_config.columns = ['x','y','z','rel_d','rel_sa','w','e','s','n','t','b',\
                                'w_na','e_na','s_na','n_na','t_na','b_na','pid']
    with open(os.path.join(save_dir, 'gsp_config.csv'), 'w') as f:
        f.write("#version = 1\n")
        global_df_config.to_csv(f, index=False, header=True, sep=",")


def superquadratic_surface(a, b, c, m, n, u_samples=30, v_samples=30):
    u = np.linspace(-np.pi, np.pi, u_samples)
    v = np.linspace(-np.pi / 2, np.pi / 2, v_samples)
    u, v = np.meshgrid(u, v)

    cos_v = np.cos(v)
    sin_v = np.sin(v)
    cos_u = np.cos(u)
    sin_u = np.sin(u)

    x = a * np.sign(cos_v) * np.abs(cos_v)**(2/n) * np.sign(cos_u) * np.abs(cos_u)**(2/m)
    y = b * np.sign(cos_v) * np.abs(cos_v)**(2/n) * np.sign(sin_u) * np.abs(sin_u)**(2/m)
    z = c * np.sign(sin_v) * np.abs(sin_v)**(2/n)

    return np.vstack((x.flatten(), y.flatten(), z.flatten())).T

def is_nested_list_empty(nested_list):
    if not nested_list:
        return True
    for item in nested_list:
        if isinstance(item, list):
            if not is_nested_list_empty(item):
                return False
        else:
            return False
    return True
# ==========================================================================================================
# Main seeding function starts here
# ==========================================================================================================
# generate the dataframe
def seed_particles_gsp(file_type='stl',
                       sq_input = (0.01, 0.01, 0.01, 2, 2),
                       q = (1,0,0,0),
                       center = (0,0,0),
                       fname='gsp_sample_data/gsp.stl',
                       discretization = (2, 2, 2),
                       wall_extend_ratio = 0.0,
                       mesh_decimation = False,
                       r_ratio = 0.1,
                       isVTK = False,
                       isUsrVel = True,
                       usrvel = (0.0, 0.0, 0.0),
                       rho_gp = 1000,
                       unit = 1,
                       sq_resolution = 30,
                       save_dir='example1'):

    if file_type == 'stl':
        # Read the STL file
        stlReader = vtk.vtkSTLReader()
        stlReader.SetFileName(fname)
        stlReader.Update()
        stlPolyData = stlReader.GetOutput()
    elif file_type == 'step':
        occReader = vtk.vtkOCCReader()
        occReader.SetFileName("your_step_file.step")
        occReader.Update()
        stlPolyData = occReader.GetOutput()
    elif file_type == 'sq':
        asq, bsq, csq, msq, nsq = sq_input
        # vtk sq input is different with what we have, do a transform
        epsilon_1 = 2.0/nsq
        epsilon_2 = 2.0/msq
        asq *= 2
        bsq *= 2
        csq *= 2

        sq_source = vtk.vtkSuperquadricSource()
        sq_source.SetAxisOfSymmetry(2)
        # set center and semi axis
        sq_source.SetCenter(0.0, 0.0, 0.0)
        sq_source.SetScale(asq, bsq, csq)
        # shape resolution
        sq_source.SetThetaResolution(sq_resolution)
        sq_source.SetPhiResolution(sq_resolution)
        # sq roundness
        sq_source.SetPhiRoundness(epsilon_1)
        sq_source.SetThetaRoundness(epsilon_2)
        sq_source.Update()

        # Save the stl for future reference
        stlWriter = vtk.vtkSTLWriter()
        stlWriter.SetFileName("sq2stl.stl")
        stlWriter.SetInputConnection(sq_source.GetOutputPort())
        stlWriter.Write()

        dprint("finish creating sqpoints")

        fname = 'sq2stl.stl'

        stlReader = vtk.vtkSTLReader()
        stlReader.SetFileName(fname)
        stlReader.Update()
        stlPolyData = stlReader.GetOutput()
    else:
        raise ValueError(file_type)

    if mesh_decimation:
        # Decimate the mesh
        dprint("Before decimation\n"
            "-----------------\n"
            "There are " + str(stlPolyData.GetNumberOfPoints()) + "nodes.\n"
            "There are " + str(stlPolyData.GetNumberOfPolys()) + "elements.\n")
        decimate = vtk.vtkQuadricDecimation()
        decimate.SetInputData(stlPolyData)
        decimate.SetTargetReduction(1-r_ratio)
        decimate.Update()
        stlPolyData = decimate.GetOutput()
        # stlPolyData.ShallowCopy(decimate.GetOutput())
        dprint("After decimation \n"
            "-----------------\n"
            "There are " + str(stlPolyData.GetNumberOfPoints()) + "nodes.\n"
            "There are " + str(stlPolyData.GetNumberOfPolys()) + "elements.\n")

    M = None
    if tuple(q) != (1,0,0,0):
        M = q_to_m(q)
        # Due to the surface resolution, the center of mass filter might not be accurate enough.
        # Maybe a weighted average method should be used here.
        centerOfMassFilter = vtk.vtkCenterOfMass()
        centerOfMassFilter.SetInputData(stlPolyData)
        centerOfMassFilter.SetUseScalarsAsWeights(False)
        centerOfMassFilter.Update()

        stl_centroid = centerOfMassFilter.GetCenter()

        stl_points = stlPolyData.GetPoints()
        num_points = stl_points.GetNumberOfPoints()

        for i in range(num_points):
            if M is not None:
                stl_point = np.array(stl_points.GetPoint(i))
                rotated_point = M @ (stl_point - stl_centroid) + stl_centroid
                stl_points.SetPoint(i, rotated_point)

        stlPolyData.Modified()

        # Save the stl for future reference
        stlWriter = vtk.vtkSTLWriter()
        stlWriter.SetFileName("after_rotation.stl")
        stlWriter.SetInputData(stlPolyData)
        stlWriter.Write()
    fname = 'after_rotation.stl'

    # Calculate the geo info
    massProperties = vtk.vtkMassProperties()
    massProperties.SetInputData(stlPolyData)
    geo_volume = massProperties.GetVolume() * unit ** 3
    geo_area = massProperties.GetSurfaceArea() * unit ** 2
    equiv_radius = (geo_volume * 3 / (4 * np.pi)) ** (1.0/3.0)
    equiv_diameter = 2 * equiv_radius

    # Obtain the bounding box, based on the wall extend ratio to make sure the bounding box is include the whole geometry.
    xmin, xmax, ymin, ymax, zmin, zmax = stlPolyData.GetBounds()
    stl_origin = np.array([xmin, ymin, zmin])
    stl_end = np.array([xmax, ymax, zmax])

    dlength = stl_end - stl_origin

    if wall_extend_ratio != 0.0:
        eps = dlength*wall_extend_ratio/2
        stl_origin -= eps
        stl_end += eps

    spacing = (stl_end-stl_origin)/discretization

    # Container initialization, used in the iteration
    nba = []
    sArea = []
    scx = []
    scy = []
    scz = []
    sradius = []
    vertices_sds = []

    # variables initialization, used in the iteration
    gid = 1
    sphere_cnt = 0
    cur_iteration = 0
    total_iteration = discretization[0] * discretization[1] * discretization[2]
    totalV = 0.0
    partition_id = 0
    sphere_id_dict = dict()

    # Some vtk filters will change the polydata implicitly (ex. vtk.vtkImplicitPolyDataDistance()).
    # For safe purpose, deep copy the original stlPolyData and use it whenever needed.
    clonedPolyData = vtk.vtkPolyData()
    clonedPolyData.DeepCopy(stlPolyData)

    for j in range(0, discretization[1]):
        for i in range(0, discretization[0]):
            for k in range(0, discretization[2]):

                cur_iteration = cur_iteration + 1
                cube_start = stl_origin + np.array([i, j, k]) * spacing
                cube_end = cube_start + spacing
                ijk = f"{i}{j}{k}"

                ret_intersect = stl_cube_clipclosesruface(ijk, cube_start,
                                                      cube_end, spacing,
                                                      stlPolyData, clonedPolyData,
                                                      isVTK,
                                                      discretization,
                                                      save_dir=save_dir)

                is_valid_intersect, component_vol, component_centroid, neigh_area, particle_surface_area= ret_intersect
                # Create local sphere id only if it is a valid intersection
                if not is_valid_intersect:
                    percentage = (cur_iteration / total_iteration) * 100
                    # dprint(f"Intersecting---->{percentage:.1f}%")
                    if stop_flag:
                        return
                    if progress_callback: progress_callback(percentage)
                    continue
                else:
                    sphere_cnt += 1
                    partition_id = partition_id + 1
                    current_list = [i,j,k]
                    id_key = tuple(current_list)
                    if id_key in sphere_id_dict:
                        # duplicated key is not allowed
                        raise KeyError(id_key)
                    else:
                        sphere_id_dict[tuple(current_list)] = partition_id

                    component_radius = ((3 * component_vol) / (4 * math.pi)) ** (1.0/3.0)
                    scx.append(component_centroid[0] * unit)
                    scy.append(component_centroid[1] * unit)
                    scz.append(component_centroid[2] * unit)
                    sradius.append(component_radius * unit)
                    nba.append(neigh_area)
                    sArea.append(particle_surface_area)

                    totalV += component_vol

                percentage = (cur_iteration / total_iteration) * 100
                # dprint(f"Intersecting---->{percentage:.1f}%")
                if stop_flag:
                    return
                if progress_callback: progress_callback(percentage)

    dprint("")
    dprint("Geo operation is done, saving output files")
    # Sphere info related to usr0.input
    # Create a empty dataframe to store neighbors id using dict container
    columns = ['left', 'right', 'front', 'back', 'top', 'bot']
    df = pd.DataFrame(columns=columns)

    for locijk, sid in sphere_id_dict.items():
        i,j,k = locijk
        #find neighbors based on i,j,k only
        neighbors = {
            'left': sphere_id_dict.get(tuple([i-1,j,k]), -1),
            'right': sphere_id_dict.get(tuple([i+1, j, k]), -1),
            'front': sphere_id_dict.get(tuple([i, j-1, k]), -1),
            'back': sphere_id_dict.get(tuple([i, j+1, k]), -1),
            'top': sphere_id_dict.get(tuple([i, j, k+1]), -1),
            'bot': sphere_id_dict.get(tuple([i, j, k-1]), -1),
        }
        df.loc[sid] = neighbors

    # Dump the neighbor area to a dataframe
    leftda = []
    rightda = []
    frontda = []
    backda = []
    tda = []
    bda = []
    for sp in nba:
        leftd,rightd,frontd,backd,td,bd = sp
        leftda.append(leftd)
        rightda.append(rightd)
        frontda.append(frontd)
        backda.append(backd)
        tda.append(td)
        bda.append(bd)

    # Construct df for gsp_config.dat
    df.insert(0,"gid",gid,True)
    df.insert(0,"surfaceArea",sArea,True)

    df["leftA"] = leftda
    df["rightA"] = rightda
    df["frontA"] = frontda
    df["backA"] = backda
    df["topA"] = tda
    df["botA"] = bda

    df["surfaceArea"] = df["surfaceArea"] * unit * unit
    total_sa = df["surfaceArea"].sum()
    df["leftA"] = df["leftA"] * unit * unit
    df["rightA"] = df["rightA"] * unit * unit
    df["frontA"] = df["frontA"] * unit * unit
    df["backA"] = df["backA"] * unit * unit
    df["topA"] = df["topA"] * unit * unit
    df["botA"] = df["botA"] * unit * unit

    # Construct df for particle_input.csv
    df_particleinput = pd.DataFrame()
    df_particleinput["X"] = scx
    df_particleinput["Y"] = scy
    df_particleinput["Z"] = scz
    df_particleinput["Radius"] = sradius
    df_particleinput["density"] = rho_gp

    # Assign random velocity or usr-specific velocity
    if isUsrVel:
        df_particleinput["U"] = usrvel[0]
        df_particleinput["V"] = usrvel[1]
        df_particleinput["W"] = usrvel[2]
    else:
        xrand = random.random()
        yrand = random.random()
        zrand = random.random()
        df_particleinput["U"] = xrand
        df_particleinput["V"] = yrand
        df_particleinput["W"] = zrand

    # Calculate the glued particle center of mass
    r3 = np.array(sradius) ** 3
    mass = np.mean(r3)
    gp_cx = np.mean(scx * r3)/mass
    gp_cy = np.mean(scy * r3)/mass
    gp_cz = np.mean(scz * r3)/mass

    # Make glued particle center of mass as the origin
    df_particleinput["X"] -= gp_cx
    df_particleinput["Y"] -= gp_cy
    df_particleinput["Z"] -= gp_cz

    # Calculate the bounding sphere diameter, pick the max{distance + radius}
    maxR = 0
    for _, row in df_particleinput.iterrows():
        tmpR = math.sqrt(row['X']**2 + row['Y']**2 + row['Z']**2) + row['Radius']
        maxR = max(maxR, tmpR)
    bounding_diameter = 2 * maxR

    nspheres = sphere_cnt
    totalV = totalV * unit ** 3
    surface_devia = abs(geo_area - total_sa)/geo_area * 100
    vol_devia = abs(geo_volume - totalV)/geo_volume * 100
    dprint("")
    dprint(f"Surface area:{geo_area:.4e} m2")
    dprint(f"Component total surface area:{total_sa:.4e}")
    dprint(f"The deviation on surface is around {surface_devia:.4f}%")
    dprint(f"Volume from original shape:{geo_volume:.4e} m3")
    dprint(f"Volume from the sum of intersection volumes:{totalV:.4e} m3")
    dprint(f"The deviation on volume is around {vol_devia:.4f}%")
    dprint(f"Volume equivalent diameter:{equiv_diameter:.4e} m")
    dprint(f'Bounding sphere diameter: {bounding_diameter:.4e} m')
    dprint(f"Component spheres in current glued particle: {nspheres}")
    dprint("")

    return nspheres, df, df_particleinput, equiv_diameter, bounding_diameter, geo_area, geo_volume, (gp_cx,gp_cy,gp_cz), surface_devia, vol_devia

def stl_cube_clipclosesruface(ijk,
                          cube_start, cube_end, spacing,
                          stlPolyData, clonedPolyData, isVTK, discretization, save_dir='example1'):


    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputData(stlPolyData)
    cleaner.Update()
    stlPolyData = cleaner.GetOutput()

    if isVTK:
        dir_name_stl = 'intersections_files'
        dir_name_stl = os.path.join(save_dir, dir_name_stl)
        if not os.path.exists(dir_name_stl):
            os.makedirs(dir_name_stl)

    particle_surface_area = 0.0
    area_left, area_right, area_front, area_back, area_top, area_bot = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    neigh_area = []

    is_valid_intersection = False
    component_vol = 0.0
    component_centroid = np.array([0.0, 0.0, 0.0])

    # Create a cube source but only to determine if a cube is completely inside or outside
    cube_center = (cube_start + cube_end) / 2.0
    cubeSource = vtk.vtkCubeSource()
    cubeSource.SetXLength(spacing[0])
    cubeSource.SetYLength(spacing[1])
    cubeSource.SetZLength(spacing[2])
    cubeSource.SetCenter(cube_center[0], cube_center[1], cube_center[2])
    cubeSource.Update()

    # Trianglize the cube
    triangleFilter = vtk.vtkTriangleFilter()
    triangleFilter.SetInputConnection(cubeSource.GetOutputPort())
    triangleFilter.Update()

    cubePolyData = triangleFilter.GetOutput()

    # Save cube for future reference
    if isVTK:
        cubeStlWriter1 = vtk.vtkSTLWriter()
        tmp_file_path = os.path.join(dir_name_stl, f"cube_{ijk}.stl")
        cubeStlWriter1.SetFileName(tmp_file_path)
        cubeStlWriter1.SetInputData(cubePolyData)
        cubeStlWriter1.Write()

    # Check if the current intersection is total inside of the stl
    if is_cube_inside_stl(cubePolyData, stlPolyData):
        massProperties = vtk.vtkMassProperties()
        massProperties.SetInputData(cubePolyData)
        component_vol = massProperties.GetVolume()

        centerOfMassFilter = vtk.vtkCenterOfMass()
        centerOfMassFilter.SetInputData(cubePolyData)
        centerOfMassFilter.SetUseScalarsAsWeights(False)
        centerOfMassFilter.Update()
        component_centroid = centerOfMassFilter.GetCenter()
        is_valid_intersection = True
        area_left = spacing[0] * spacing[2]
        area_right = area_left
        area_front = spacing[1] * spacing[2]
        area_back = area_front
        area_top = spacing[0] * spacing[1]
        area_bot = area_top

        neigh_area = [area_left, area_right, area_front, area_back, area_top, area_bot]
        particle_surface_area = 0.0
        # dprint(f"The cube {ijk} is completely inside the STL, skip the current intersection.")
        if isVTK:
            cubeStlWriter1 = vtk.vtkSTLWriter()
            tmp_file_path = os.path.join(dir_name_stl, f"cube_{ijk}.stl")
            cubeStlWriter1.SetFileName(tmp_file_path)
            cubeStlWriter1.SetInputData(cubePolyData)
            cubeStlWriter1.Write()
        return is_valid_intersection, component_vol, component_centroid, neigh_area, particle_surface_area

    # Check if the current clip is total outside of the stl.
    # When discretization in any direction equal to one, the clip will also produce zero intersection.
    if is_cube_outside_stl(cubePolyData, stlPolyData) and all(d > 1 for d in discretization):
        is_valid_intersection = False
        # dprint(f"The cube {ijk} is completely outside of the STL, skip the current intersection.")
        return is_valid_intersection, component_vol, component_centroid, neigh_area, particle_surface_area

    planes = vtk.vtkPlaneCollection()

    # if discretization[0] > 1:
    planeyz_start = vtk.vtkPlane()
    planeyz_start.SetOrigin(cube_start[0], 0, 0)
    planeyz_start.SetNormal(1, 0, 0)
    planes.AddItem(planeyz_start)

    planeyz_end = vtk.vtkPlane()
    planeyz_end.SetOrigin(cube_end[0], 0, 0)
    planeyz_end.SetNormal(-1, 0, 0)
    planes.AddItem(planeyz_end)

    # if discretization[1] > 1:
    planexz_start = vtk.vtkPlane()
    planexz_start.SetOrigin(0, cube_start[1], 0)
    planexz_start.SetNormal(0, 1, 0)
    planes.AddItem(planexz_start)

    planexz_end = vtk.vtkPlane()
    planexz_end.SetOrigin(0, cube_end[1], 0)
    planexz_end.SetNormal(0, -1, 0)
    planes.AddItem(planexz_end)

    # if discretization[2] > 1:
    planexy_start = vtk.vtkPlane()
    planexy_start.SetOrigin(0, 0, cube_start[2])
    planexy_start.SetNormal(0, 0, 1)
    planes.AddItem(planexy_start)

    planexy_end = vtk.vtkPlane()
    planexy_end.SetOrigin(0, 0, cube_end[2])
    planexy_end.SetNormal(0, 0, -1)
    planes.AddItem(planexy_end)

    # Clipping the stl and forming the closed surface
    clipClosedSurface = vtk.vtkClipClosedSurface()
    clipClosedSurface.SetInputData(stlPolyData)
    clipClosedSurface.SetClippingPlanes(planes)
    clipClosedSurface.GenerateFacesOn()
    clipClosedSurface.Update()

    clippedPolyData = clipClosedSurface.GetOutput()

    if clippedPolyData.GetNumberOfPoints() == 0:
        is_valid_intersection = False
        # dprint(f"The cube {ijk} is completely outside of the STL, skip the current intersection.")
        return is_valid_intersection, component_vol, component_centroid, neigh_area, particle_surface_area

    is_valid_intersection = True

    # Use vtkConnectivityFilter to identify distinct regions (separate volumes)
    connectivityFilter = vtk.vtkConnectivityFilter()
    connectivityFilter.SetInputData(clippedPolyData)
    # connectivityFilter.SetExtractionModeToLargestRegion()
    connectivityFilter.Update()

    num_regions = connectivityFilter.GetNumberOfExtractedRegions()
    # if num_regions > 1:
    #     dprint(f" {ijk} encounter {num_regions} separate volumes, special treatment to merge them")

    # Save the intersection and the cube to an STL file
    if isVTK:
        stlWriter = vtk.vtkSTLWriter()
        tmp_file_path = os.path.join(dir_name_stl, f"intersection_{ijk}.stl")
        stlWriter.SetFileName(tmp_file_path)
        stlWriter.SetInputData(clippedPolyData)
        stlWriter.Write()

    triangles_cnt = 0

    # Iterate over each region and process them.
    # If separate volumes exist, their surface area and neighbor area should be added up together,
    # but need to process the tol_volume and mass_center separately.
    for region_id in range(num_regions):
        total_area = 0.0
        center = [0.0, 0.0, 0.0]
        connectivityFilter.SetExtractionModeToSpecifiedRegions()
        connectivityFilter.InitializeSpecifiedRegionList()
        connectivityFilter.AddSpecifiedRegion(region_id)
        connectivityFilter.Update()
        regionPolyData = connectivityFilter.GetOutput()

        # Calculate the surface area of the intersection
        massProperties = vtk.vtkMassProperties()
        massProperties.SetInputData(regionPolyData)
        cur_volume = massProperties.GetVolume()
        component_vol += cur_volume

        polys = regionPolyData.GetPolys()
        polys.InitTraversal()

        idList = vtk.vtkIdList()

        while polys.GetNextCell(idList):
            if idList.GetNumberOfIds() == 3:  # Ensure it's a triangle
                is_particle_boundary = True
                triangle_nodes = []
                triangles_cnt += 1
                points = vtk.vtkPoints()
                for i in range(3):
                    pointId = idList.GetId(i)
                    point = regionPolyData.GetPoint(pointId)
                    triangle_nodes.append(point)
                    points.InsertNextPoint(point)
                triangle = vtk.vtkTriangle()
                triangle.GetPoints().SetPoint(0, points.GetPoint(0))
                triangle.GetPoints().SetPoint(1, points.GetPoint(1))
                triangle.GetPoints().SetPoint(2, points.GetPoint(2))
                triangle_area = triangle.TriangleArea(points.GetPoint(0), points.GetPoint(1), points.GetPoint(2))

                centroid = [(points.GetPoint(0)[i] + points.GetPoint(1)[i] + points.GetPoint(2)[i]) / 3.0 for i in range(3)]
                center[0] += triangle_area * centroid[0]
                center[1] += triangle_area * centroid[1]
                center[2] += triangle_area * centroid[2]
                total_area += triangle_area
                triangle_nodes.append(centroid)

                triangle_nodes_np = np.array(triangle_nodes)
                total_sdf, vertices_sds= calculate_signed_distance(stlPolyData, triangle_nodes_np)
                vertices_sdf_tolerance = 1e-6
                for vertices_sd in vertices_sds:
                    if abs(vertices_sd) > vertices_sdf_tolerance:
                        is_particle_boundary = False

                x_coords, y_coords, z_coords = np.transpose(triangle_nodes_np)

                if np.all(np.isclose(x_coords, cube_start[0])) and not is_particle_boundary:
                    area_left = area_left + triangle_area
                elif np.all(np.isclose(x_coords, cube_end[0])) and not is_particle_boundary:
                    area_right = area_right + triangle_area
                elif np.all(np.isclose(y_coords, cube_start[1])) and not is_particle_boundary:
                    area_front = area_front + triangle_area
                elif np.all(np.isclose(y_coords, cube_end[1])) and not is_particle_boundary:
                    area_back = area_back + triangle_area
                elif np.all(np.isclose(z_coords, cube_start[2])) and not is_particle_boundary:
                    area_bot = area_bot + triangle_area
                elif np.all(np.isclose(z_coords, cube_end[2])) and not is_particle_boundary:
                    area_top = area_top + triangle_area
                else:
                    particle_surface_area = particle_surface_area + triangle_area

        center[0] /= total_area
        center[1] /= total_area
        center[2] /= total_area
        # collect the sum of current intersection center multiply by current volume
        component_centroid[0] += center[0] * cur_volume
        component_centroid[1] += center[1] * cur_volume
        component_centroid[2] += center[2] * cur_volume

    # Weighted average of the centroid if there are multiple intersects at one time
    component_centroid /= component_vol
    neigh_area = [area_left, area_right, area_front, area_back, area_top, area_bot]

    return is_valid_intersection, component_vol, component_centroid, neigh_area, particle_surface_area


def generate_gsp(file_type = 'sq', # single variable
                 sq_inputs = [(0.01,0.01,0.01,3,3)],
                 qs = [(1,0,0,0)],
                 centers = [(0,0,0)],
                 usr_input_file_names = [""],
                 discretizations = [[2,2,2]],
                 wall_extend_ratios = [0.001],
                 ngps = [1],
                 mesh_decimations=[False],
                 r_ratios=[0.1],
                 isRemoveTinyVolume = False,
                 critical_vol_fraction=0.001,
                 critical_totalvol_fraction=0.99,
                 isVTKs=[False],
                 isUsrVels = [True],
                 usrvels = [[0.0,0.0,0.0]],
                 rho_gps = [682],
                 units = [1],
                 isTp = False, # single variable
                 isRandomQ = False,
                 Tps = [273.15],
                 startPhase = None,
                 origin=[0.0, 0.0, 0.0],
                 end=[0.04, 0.04, 0.04],
                 spacing=[0.02, 0.02, 0.02],
                 sq_resolution = 30,
                 save_dir=os.path.dirname(os.path.abspath(__file__))):

    # below are hard coded parameters for initialization use
    #---------------------------------------
    tmp_number = 0                         #
    start_gid = 1                          #
    global_df_config = pd.DataFrame()      #
    global_df_particle = pd.DataFrame()    #
    df_gluepos = pd.DataFrame()            #
    df_glueposcopy = pd.DataFrame()        #
    phase_id = []                          #
    spheres_per_particle = []              #
    Niter = 1                              #
    nsphere_after_per_phase = 0            #
    #---------------------------------------
    if startPhase is not None:
        Nphase = startPhase - 1 # Will be incremented
    else:
        Nphase = 0
    # perform a valuefill to keep zip function working
    arrays = [sq_inputs, qs, centers, usr_input_file_names, discretizations, wall_extend_ratios, ngps,
              mesh_decimations, r_ratios, isVTKs,
              isUsrVels, usrvels, rho_gps, units, Tps]

    max_size = max(len(arr) for arr in arrays)

    for i, arr in enumerate(arrays):
        if len(arr) < max_size:
            arrays[i] = arr + [arr[0]] * (max_size - len(arr))

    for (sq_input, q, center, fname, discretization, wall_extend_ratio, ngp, mesh_decimation,
         r_ratio, isVTK, isUsrVel, usrvel,
         rho_gp, unit, Tp) in zip(*arrays):

        if stop_flag:
            return

        ret = seed_particles_gsp(
            file_type = file_type,
            sq_input = sq_input,
            q = q,
            center = center,
            fname = fname,
            discretization = discretization,
            wall_extend_ratio = wall_extend_ratio,
            mesh_decimation = mesh_decimation,
            r_ratio = r_ratio,
            isVTK = isVTK,
            isUsrVel = isUsrVel,
            usrvel = usrvel,
            rho_gp = rho_gp,
            unit = unit,
            sq_resolution = sq_resolution,
            save_dir=save_dir)
        if not ret:
            return
        nspheres, df_usr0, df_particle, equiv_diameter, bounding_diameter, area, volume, com, surface_devia, vol_devia = ret

        # start calling combined_seed
        Nphase += 1
        phase_id.append(Nphase)

        df_config_tmp, df_particle_tmp, nsphere_after_per_phase = combined_seed(
            df_usr0, df_particle, df_glueposcopy,
            bounding_diameter, Nphase, rho_gp,
            isUsrVel, usrvel, isTp, Tp, isRandomQ,
            isRemoveTinyVolume, critical_vol_fraction, critical_totalvol_fraction,
            save_dir=save_dir)

        spheres_per_particle.append(nsphere_after_per_phase)

        global_df_config = pd.concat([global_df_config, df_config_tmp], ignore_index=True)


    if stop_flag:
        return

    write_df(global_df_config, global_df_particle, phase_id,
             spheres_per_particle, isTp, save_dir=save_dir)

    return  nspheres, equiv_diameter, bounding_diameter, area, volume, com, surface_devia, vol_devia

def q_to_m(q):
    h = np.sqrt(np.dot(q,q))
    q = np.array(q) / h # normalize to unit quaternion
    M = np.zeros((3,3))
    M[0,0] = 2*(q[0]**2 + q[1]**2) - 1
    M[1,1] = 2*(q[0]**2 + q[2]**2) - 1
    M[2,2] = 2*(q[0]**2 + q[3]**2) - 1
    M[1,0] = 2*(q[1]*q[2] + q[0]*q[3])
    M[0,1] = 2*(q[1]*q[2] - q[0]*q[3])
    M[0,2] = 2*(q[1]*q[3] + q[0]*q[2])
    M[2,0] = 2*(q[1]*q[3] - q[0]*q[2])
    M[2,1] = 2*(q[2]*q[3] + q[0]*q[1])
    M[1,2] = 2*(q[2]*q[3] - q[0]*q[1])
    return M
