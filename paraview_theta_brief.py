# trace generated using paraview version 5.10.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
import glob

# Output Directory
filenames = glob.glob('/Users/user/Downloads/OUTPUT.13201065/ins_velo_*')
out_path = '/Users/user/Downloads/test/theta.png'

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


# import pdb; pdb.set_trace()
# create a new 'VisItTecplotBinaryReader'
ins_velo_000000 = VisItTecplotBinaryReader(registrationName='ins_velo_000000*', FileName=filenames)

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on ins_velo_000000
ins_velo_000000.PointArrayStatus = ['theta', 'u', 'v', 'w', 'x', 'y', 'z']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
ins_velo_000000Display = Show(ins_velo_000000, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
ins_velo_000000Display.Representation = 'Outline'

# reset view to fit data
renderView1.ResetCamera(False)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Iso Volume'
isoVolume1 = IsoVolume(registrationName='IsoVolume1', Input=ins_velo_000000)

# Properties modified on isoVolume1
isoVolume1.ThresholdRange = [0.01, 1.0]

# show data in view
isoVolume1Display = Show(isoVolume1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
isoVolume1Display.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(ins_velo_000000)

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=ins_velo_000000)

# toggle 3D widget visibility (only when running from the GUI)
Hide3DWidgets(proxy=clip1.ClipType)

# Properties modified on clip1
clip1.Invert = 0

# Properties modified on clip1.ClipType
clip1.ClipType.Normal = [0.0, 1.0, 0.0]

# show data in view
clip1Display = Show(clip1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
clip1Display.Representation = 'Surface'

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'u'))

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')
uLUT.RGBPoints = [1.1960824728012085, 0.231373, 0.298039, 0.752941, 12.022967994213104, 0.865003, 0.865003, 0.865003, 22.849853515625, 0.705882, 0.0156863, 0.14902]
uLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'u'
uPWF = GetOpacityTransferFunction('u')
uPWF.Points = [1.1960824728012085, 0.0, 0.5, 0.0, 22.849853515625, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# reset view to fit data
renderView1.ResetCamera(False)

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(972, 1153)

# current camera placement for renderView1
renderView1.CameraPosition = [-1.3404546376832744, -13.913813777467306, 4.206886251320443]
renderView1.CameraFocalPoint = [3.1170489788055415, 1.5585244894027708, 0.9741111993789672]
renderView1.CameraViewUp = [0.062002041424191974, 0.18698230745939765, 0.9804046937649739]
renderView1.CameraParallelScale = 4.292371598877116

# save animation
SaveAnimation(out_path, renderView1, ImageResolution=[972, 1153],
    TransparentBackground=1,
    FrameWindow=[0, 20])

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(972, 1153)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [-1.3404546376832744, -13.913813777467306, 4.206886251320443]
renderView1.CameraFocalPoint = [3.1170489788055415, 1.5585244894027708, 0.9741111993789672]
renderView1.CameraViewUp = [0.062002041424191974, 0.18698230745939765, 0.9804046937649739]
renderView1.CameraParallelScale = 4.292371598877116

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).