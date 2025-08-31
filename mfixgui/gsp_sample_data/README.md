# Glued sphere generator

## general procedure: given the input geometry, this script will firstly perform a iteratively intersecting, then convert the intersection into a volume equivalent sphere, this sphere is then a component sphere of the glued particle

Required third-party libraries:
vtk/9.2.6    numpy/1.26.0    pandas/2.1.3
trimesh/4.0.4

```python
import math
import os
import sys
import vtk
import numpy as np
import shutil
import time
import pandas as pd
import random
from vtk.util import numpy_support
import trimesh
```
## op level function

```python
def gsp_data_generator(file_type = 'stl', # string
                       sq_inputs = [[1.0, 2.0, 3.0, 2, 2]], # 2D vector
                       usr_input_file_names = ["r4h5dot3.stl"], # string vector
                       discretizations = [[2,2,2]], # 2D vector
                       wallextendratios = [0.001], # vector
                       ngps = [1], # vector
                       mesh_decimations=[False], # logical vector
                       r_ratios=[0.1], # vector
                       isVTKs=[False], # logical vector
                       isUsrVels = [True], # logical vector
                       usrvels = [[0.0,0.0,0.0]], # 2D vector
                       rho_gps = [682], # vector
                       units = [1e-3], # vector
                       isTp = False, # logical variable
                       Tps = [273.15], # vector
                       isAutoSeeding = False, # logical variable
                       isRandomQ = False, # logical variable
                       origin=[0.0, 0.0, 0.0], # vector
                       end=[0.04, 0.04, 0.04], # vector
                       spacing=[0.02, 0.02, 0.02], # vector
                       save_folder=os.path.dirname(os.path.abspath(__file__))): # string
```

## Detailed procedure to generate glued sphere data:

1. provide the necessary input variables:

    **Note: except isTp,origin,end and spacing, other variables should be a 2D array, their row numbers represent the number of solid phases in the system, and they should be equal to each other.**
    - file_type = 'stl': current the script support "stl", "step" and superquadratic input.

    - usr_input_file_names = ["gsp1.stl","gsp2.stl",...]: char vectors, each represents a input stl file, each stl file represents a solid phase, this is used when file_type = 'stl' or 'step'.

    - sq_inputs = [[a,b,c,m,n],[a1,b1,c1,m1,n1],...]: when file_type = 'sq', this sq_inputs should be used, providing a list of superquadratic shape parameters vector.

    - discretizations = [[2,2,2],[3,3,3],...]: the discretization in each direction of a glued sphere, the bounding box of the input geometry will be obtained first then divided by discretizations to generate a cube, this cube is then intersecting with input geometry and the intersection will be converted into a volume equivalent component sphere.

    - wallextendratios = [0.001,0.001,...]: extend the boundary box of the glued sphere to ensure the safe intersection, it is highly recommended keep the default value 0.001

    - ngps = [1,2,...]: number of glued spheres in each solid phase

    - mesh_decimations=[False,True,...]: flag to perform a mesh decimation to the input stl file, remesh the surface mesh of the input stl file and decrease the number of elements.

    - r_ratios=[0.1,0.5,...]: the ratio of number of surface elements before mesh decimation to after mesh decimation.

    - isVTKs=[False,True...]: flag to write the vtk of each intersections, use to check if intersecting properly, default to be False.

    - isUsrVels = [True,False,...]: flag to use user-defined sphere velocities, if false, scripts will assign some random velocity from 0~1 m/s to glued spheres in the same solid phase.

    - usrvels = [[0.0,0.0,0.0],[0.0,0.0,-0.1],...]: the value of user-defined sphere velocities.

    - rho_gps = [682,782,...]: the density of glued sphere in each solid phase

    - units = [1e-3,1e-6,...]: the origin unit system of the input geometries, they will be converted into SI.

    - isTp = False: a global flag to include temperature in output data files.

    - Tps = [273.15,293.15,...]: assign a usr-defined initial temperature to glued spheres in each solid phase.

    - isAutoSeeding = False: If true, MFiX will perform automatic seeding, otherwise the seeding locations should be specify.

    - isRandomQ = False: If true, q1 ~ q4 will be assigned random values, otherwise, they are [1,0,0,0]

    - origin=[0.0, 0.0, 0.0],end=[0.1, 0.1, 0.1],spacing=[0.01, 0.01, 0.01]: creat a uniformly distributed point set in space to manually seed the glued spheres.
    **(ideally, these are only used to create a point set to manually seed the glued spheres,  we should provide the freedom to user for providing a data file of point sets, which we can used to seed)**

    - save_folder=os.path.dirname(os.path.abspath(__file__)): if not explicitly assign a folder, all files will save to working directory

2. Run the python script from command line:
    - python gsp_generator.py
    - the script will generate the glued sphere configuration based on the order of input files or superquadractic vectors

3. Output the data file for MFIX use:
    - if autoseeding is true:
        -  only 'gsp_config.dat' will be saved
    - if autoseeding is false:
        - 'particle_input.dat' and 'gsp_config.dat'

4. Visualization for validation use:
    - the vtp files will be saved in the designated folder with the name 'phase_1.vtp','phase_2.vtp',...

## Limitations:
- for stl file who has a large number of surface elements, it will take a long time to process the intersection, so it's recommended using mesh decimation for this situation.
- the script can take superquadratic vector as input, however, it only works for convex shape now.
- in order to avoid the dependency of the scipy library, delaunay triangulation from vtk library is used to generate the geo surface from superquadratic vector, however, sometimes it will not perfectly resolve the desired geometry, user have to adjust the u_sample and v_sample to obtain better results.
