# trace generated using paraview version 5.10.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
import glob
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'VisItTecplotBinaryReader'
filenames = glob.glob('/scratch4/qwang4/ext-zyou6474/channel_flow_multi_scalar.14129029/OUTPUT/ins_velo_*')
ins_velo_000 = VisItTecplotBinaryReader(registrationName='ins_velo_000*', FileName=filenames)
ins_velo_000.MeshStatus = ['SOULTION']
ins_velo_000.PointArrayStatus = []

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on ins_velo_000
ins_velo_000.PointArrayStatus = ['theta1', 'theta2', 'theta3', 'u', 'v', 'w', 'x', 'y', 'z']

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
ins_velo_000Display = Show(ins_velo_000, renderView1, 'StructuredGridRepresentation')

# trace defaults for the display properties.
ins_velo_000Display.Representation = 'Outline'
ins_velo_000Display.ColorArrayName = [None, '']
ins_velo_000Display.SelectTCoordArray = 'None'
ins_velo_000Display.SelectNormalArray = 'None'
ins_velo_000Display.SelectTangentArray = 'None'
ins_velo_000Display.OSPRayScaleArray = 'theta1'
ins_velo_000Display.OSPRayScaleFunction = 'PiecewiseFunction'
ins_velo_000Display.SelectOrientationVectors = 'None'
ins_velo_000Display.ScaleFactor = 0.6258641719818115
ins_velo_000Display.SelectScaleArray = 'None'
ins_velo_000Display.GlyphType = 'Arrow'
ins_velo_000Display.GlyphTableIndexArray = 'None'
ins_velo_000Display.GaussianRadius = 0.03129320859909058
ins_velo_000Display.SetScaleArray = ['POINTS', 'theta1']
ins_velo_000Display.ScaleTransferFunction = 'PiecewiseFunction'
ins_velo_000Display.OpacityArray = ['POINTS', 'theta1']
ins_velo_000Display.OpacityTransferFunction = 'PiecewiseFunction'
ins_velo_000Display.DataAxesGrid = 'GridAxesRepresentation'
ins_velo_000Display.PolarAxes = 'PolarAxesRepresentation'
ins_velo_000Display.ScalarOpacityUnitDistance = 0.03496112989414026

# reset view to fit data
renderView1.ResetCamera(False)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=ins_velo_000)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'theta1']
clip1.Value = 0.5

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [3.1293208599090576, 1.5646604299545288, 0.4935269057750702]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [3.1293208599090576, 1.5646604299545288, 0.4935269057750702]

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
clip1Display.ColorArrayName = [None, '']
clip1Display.SelectTCoordArray = 'None'
clip1Display.SelectNormalArray = 'None'
clip1Display.SelectTangentArray = 'None'
clip1Display.OSPRayScaleArray = 'theta1'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 0.6258641719818115
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 0.03129320859909058
clip1Display.SetScaleArray = ['POINTS', 'theta1']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'theta1']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 0.04062730331089352
clip1Display.OpacityArrayName = ['POINTS', 'theta1']

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

# get opacity transfer function/opacity map for 'u'
uPWF = GetOpacityTransferFunction('u')

# reset view to fit data
renderView1.ResetCamera(False)

# reset view to fit data
renderView1.ResetCamera(True)

# set active source
SetActiveSource(ins_velo_000)

# create a new 'Iso Volume'
isoVolume1 = IsoVolume(registrationName='IsoVolume1', Input=ins_velo_000)
isoVolume1.InputScalars = ['POINTS', 'theta1']
isoVolume1.ThresholdRange = [0.0, 1.0]

# Properties modified on isoVolume1
isoVolume1.ThresholdRange = [0.01, 1.0]

# show data in view
isoVolume1Display = Show(isoVolume1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
isoVolume1Display.Representation = 'Surface'
isoVolume1Display.ColorArrayName = [None, '']
isoVolume1Display.SelectTCoordArray = 'None'
isoVolume1Display.SelectNormalArray = 'None'
isoVolume1Display.SelectTangentArray = 'None'
isoVolume1Display.OSPRayScaleArray = 'theta1'
isoVolume1Display.OSPRayScaleFunction = 'PiecewiseFunction'
isoVolume1Display.SelectOrientationVectors = 'None'
isoVolume1Display.ScaleFactor = 0.008796191215515137
isoVolume1Display.SelectScaleArray = 'None'
isoVolume1Display.GlyphType = 'Arrow'
isoVolume1Display.GlyphTableIndexArray = 'None'
isoVolume1Display.GaussianRadius = 0.0004398095607757568
isoVolume1Display.SetScaleArray = ['POINTS', 'theta1']
isoVolume1Display.ScaleTransferFunction = 'PiecewiseFunction'
isoVolume1Display.OpacityArray = ['POINTS', 'theta1']
isoVolume1Display.OpacityTransferFunction = 'PiecewiseFunction'
isoVolume1Display.DataAxesGrid = 'GridAxesRepresentation'
isoVolume1Display.PolarAxes = 'PolarAxesRepresentation'
isoVolume1Display.ScalarOpacityUnitDistance = 0.017000458568067745
isoVolume1Display.OpacityArrayName = ['POINTS', 'theta1']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
isoVolume1Display.ScaleTransferFunction.Points = [0.009999999776482582, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
isoVolume1Display.OpacityTransferFunction.Points = [0.009999999776482582, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(ins_velo_000)

# create a new 'Iso Volume'
isoVolume2 = IsoVolume(registrationName='IsoVolume2', Input=ins_velo_000)
isoVolume2.InputScalars = ['POINTS', 'theta1']
isoVolume2.ThresholdRange = [0.0, 1.0]

# Properties modified on isoVolume2
isoVolume2.InputScalars = ['POINTS', 'theta2']
isoVolume2.ThresholdRange = [0.01, 1.0]

# show data in view
isoVolume2Display = Show(isoVolume2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
isoVolume2Display.Representation = 'Surface'
isoVolume2Display.ColorArrayName = [None, '']
isoVolume2Display.SelectTCoordArray = 'None'
isoVolume2Display.SelectNormalArray = 'None'
isoVolume2Display.SelectTangentArray = 'None'
isoVolume2Display.OSPRayScaleArray = 'theta1'
isoVolume2Display.OSPRayScaleFunction = 'PiecewiseFunction'
isoVolume2Display.SelectOrientationVectors = 'None'
isoVolume2Display.ScaleFactor = 0.007325243949890137
isoVolume2Display.SelectScaleArray = 'None'
isoVolume2Display.GlyphType = 'Arrow'
isoVolume2Display.GlyphTableIndexArray = 'None'
isoVolume2Display.GaussianRadius = 0.00036626219749450687
isoVolume2Display.SetScaleArray = ['POINTS', 'theta1']
isoVolume2Display.ScaleTransferFunction = 'PiecewiseFunction'
isoVolume2Display.OpacityArray = ['POINTS', 'theta1']
isoVolume2Display.OpacityTransferFunction = 'PiecewiseFunction'
isoVolume2Display.DataAxesGrid = 'GridAxesRepresentation'
isoVolume2Display.PolarAxes = 'PolarAxesRepresentation'
isoVolume2Display.ScalarOpacityUnitDistance = 0.017316841437462226
isoVolume2Display.OpacityArrayName = ['POINTS', 'theta1']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
isoVolume2Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
isoVolume2Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(isoVolume2Display, ('POINTS', 'theta2'))

# rescale color and/or opacity maps used to include current data range
isoVolume2Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
isoVolume2Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'theta2'
theta2LUT = GetColorTransferFunction('theta2')

# get opacity transfer function/opacity map for 'theta2'
theta2PWF = GetOpacityTransferFunction('theta2')

# Properties modified on theta2LUT
theta2LUT.EnableOpacityMapping = 1

# Properties modified on theta2LUT
theta2LUT.EnableOpacityMapping = 0

# Properties modified on theta2LUT
theta2LUT.EnableOpacityMapping = 1

# Apply a preset using its name. Note this may not work as expected when presets have duplicate names.
theta2LUT.ApplyPreset('X Ray', True)

# Properties modified on theta2LUT
theta2LUT.RGBPoints = [0.009999999776482582, 1.0, 1.0, 1.0, 0.17396874725818634, 0.8352941176470589, 0.8352941176470589, 0.8352941176470589, 1.0, 0.0, 0.0, 0.0]

# Properties modified on theta2LUT
theta2LUT.RGBPoints = [0.009999999776482582, 1.0, 1.0, 1.0, 0.17396874725818634, 0.8352941176470589, 0.8352941176470589, 0.8352941176470589, 0.25440624356269836, 0.7529411764705882, 0.7529411764705882, 0.7529411764705882, 1.0, 0.0, 0.0, 0.0]

# Properties modified on theta2LUT
theta2LUT.RGBPoints = [0.009999999776482582, 1.0, 1.0, 1.0, 0.17396874725818634, 0.8352941176470589, 0.8352941176470589, 0.8352941176470589, 0.25440624356269836, 0.7529411764705882, 0.7529411764705882, 0.7529411764705882, 0.3905312418937683, 0.615686274509804, 0.615686274509804, 0.615686274509804, 1.0, 0.0, 0.0, 0.0]

# Properties modified on theta2LUT
theta2LUT.EnableOpacityMapping = 0

# set scalar coloring using an separate color/opacity maps
ColorBy(isoVolume2Display, ('POINTS', 'theta2'), True)

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(theta2LUT, renderView1)

# rescale color and/or opacity maps used to include current data range
isoVolume2Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
isoVolume2Display.SetScalarBarVisibility(renderView1, True)

# get separate color transfer function/color map for 'theta2'
separate_isoVolume2Display_theta2LUT = GetColorTransferFunction('theta2', isoVolume2Display, separate=True)

# get separate opacity transfer function/opacity map for 'theta2'
separate_isoVolume2Display_theta2PWF = GetOpacityTransferFunction('theta2', isoVolume2Display, separate=True)

# Properties modified on isoVolume2Display
isoVolume2Display.Opacity = 0.49

# set active source
SetActiveSource(ins_velo_000)

# create a new 'Iso Volume'
isoVolume3 = IsoVolume(registrationName='IsoVolume3', Input=ins_velo_000)
isoVolume3.InputScalars = ['POINTS', 'theta1']
isoVolume3.ThresholdRange = [0.0, 1.0]

# Properties modified on isoVolume3
isoVolume3.InputScalars = ['POINTS', 'theta3']
isoVolume3.ThresholdRange = [0.01, 1.0]

# show data in view
isoVolume3Display = Show(isoVolume3, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
isoVolume3Display.Representation = 'Surface'
isoVolume3Display.ColorArrayName = [None, '']
isoVolume3Display.SelectTCoordArray = 'None'
isoVolume3Display.SelectNormalArray = 'None'
isoVolume3Display.SelectTangentArray = 'None'
isoVolume3Display.OSPRayScaleArray = 'theta1'
isoVolume3Display.OSPRayScaleFunction = 'PiecewiseFunction'
isoVolume3Display.SelectOrientationVectors = 'None'
isoVolume3Display.ScaleFactor = 0.008796167373657227
isoVolume3Display.SelectScaleArray = 'None'
isoVolume3Display.GlyphType = 'Arrow'
isoVolume3Display.GlyphTableIndexArray = 'None'
isoVolume3Display.GaussianRadius = 0.00043980836868286134
isoVolume3Display.SetScaleArray = ['POINTS', 'theta1']
isoVolume3Display.ScaleTransferFunction = 'PiecewiseFunction'
isoVolume3Display.OpacityArray = ['POINTS', 'theta1']
isoVolume3Display.OpacityTransferFunction = 'PiecewiseFunction'
isoVolume3Display.DataAxesGrid = 'GridAxesRepresentation'
isoVolume3Display.PolarAxes = 'PolarAxesRepresentation'
isoVolume3Display.ScalarOpacityUnitDistance = 0.017000435844989618
isoVolume3Display.OpacityArrayName = ['POINTS', 'theta1']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
isoVolume3Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
isoVolume3Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

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
SaveAnimation('/scratch4/qwang4/ext-zyou6474/channel_flow_multi_scalar.14129029/OUTPUT/theta.png', renderView1, ImageResolution=[1368, 1354],
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
