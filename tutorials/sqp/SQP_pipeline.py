from paraview.simple import *
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
from vtkmodules.vtkCommonCore import vtkFloatArray

# Step 1: Get the selected source and fetch its full data
my_sqp = GetActiveSource()
UpdatePipeline()
data = servermanager.Fetch(my_sqp)

# Step 2: Shift phase index so it can be used in the Glyph table
# If phase index is not present in vtp file, it should nmean there is only one solids phase,
# create a dummy value 0.
point_data = data.GetPointData()
calculator = Calculator(Input=my_sqp, registrationName="SQP_Calculator")
calculator.ResultArrayName = "Particle_Phase_ID_minus_1"
if "Particle_Phase_ID" in point_data.keys():
    calculator.Function = "Particle_Phase_ID - 1"
else:
    calculator.Function = "0"
UpdatePipeline()

# Step 3: Get superquadric parameters for each phase from Field Data and save them in a dictionary
field_data = data.GetFieldData()
sqp_dict = {}
keys = ['a', 'b', 'c', 'm', 'n']

# Loop through all field arrays
for i in range(field_data.GetNumberOfArrays()):
    array = field_data.GetArray(i)
    name = array.GetName()
    num_tuples = array.GetNumberOfTuples()
    num_components = array.GetNumberOfComponents()

    # Process only if array name starts with "SQP_"
    if name.startswith("SQP_") and num_tuples == 5 and num_components == 1:
        phase = int(name[-2:])
    # Extract 5 scalar values
        values = [array.GetComponent(t, 0) for t in range(5)]
    # Map them to a,b,c,m,n
        sqp_dict[name] = dict(zip(keys, values))

if not sqp_dict:
    raise RuntimeError("SQP_** field data not found.")

# Step 4: Create superquadric sources for each phase and group them
sqp_sources = []

for key, sqp_param in sqp_dict.items():
#    print(f"{key}: {sqp_param}")
    superquadric = Superquadric(registrationName=key)
    superquadric.Set(
        Scale=[sqp_param['a'], sqp_param['b'], sqp_param['c']],
        ThetaResolution=16,
        PhiResolution=16,
        ThetaRoundness=2.0/sqp_param['m'],
        PhiRoundness=2.0/sqp_param['n'],
        Size=1.0,
        Toroidal=0  )
    sqp_sources.append(superquadric)

groupDatasets = GroupDatasets(registrationName='SQP_GROUP', Input=sqp_sources)

# Step 5: Rotate all SQPs 90 degrees around x-axix
transform = Transform(registrationName='SQP_GROUP_RX90', Input=groupDatasets)
transform.Transform.Rotate = [90.0, 0.0, 0.0]
view = GetActiveViewOrCreate('RenderView')
display = Show(transform, view)
display.Representation = 'Outline'
Hide(transform, view)

# Step 6: Set Calculator representation to 3D Glyph, Orient, Scale and assign sqp sources based on shifted phase index
renderView1 = GetActiveViewOrCreate('RenderView')
Hide(my_sqp, renderView1)
SetActiveSource(calculator)
SQPDisplay = GetRepresentation(calculator, view=renderView1)
SQPDisplay.SetRepresentationType('3D Glyphs')
SQPDisplay.Orient = 1
SQPDisplay.OrientationMode = 'Quaternion'
SQPDisplay.SelectOrientationVectors = 'QUATERNION'
SQPDisplay.Scaling = 1
SQPDisplay.ScaleMode = 'Magnitude'
SQPDisplay.ScaleFactor = 1.0
SQPDisplay.SelectScaleArray = 'SQP_POLY_SCALE'
SQPDisplay.GlyphType = 'Pipeline Connection'
SQPDisplay.GlyphType.Input = transform
SQPDisplay.UseGlyphTable = 1
SQPDisplay.GlyphTableIndexArray = 'Particle_Phase_ID_minus_1'
SQPDisplay.UseCompositeGlyphTable = 1

Show(calculator, renderView1)
SetActiveSource(calculator)


