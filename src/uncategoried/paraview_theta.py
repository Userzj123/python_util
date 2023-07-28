# trace generated using paraview version 5.10.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get active source.
ins_velo_000 = GetActiveSource()

# destroy ins_velo_000
Delete(ins_velo_000)
del ins_velo_000

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'VisItTecplotBinaryReader'
ins_velo_000 = VisItTecplotBinaryReader(registrationName='ins_velo_000*', FileName=['/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00000000.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00000500.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00001000.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00001500.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00002000.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00002500.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00003000.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00003500.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00004000.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00004500.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00005000.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00005500.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00006000.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00006500.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00007000.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00007500.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00008000.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00008500.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00009000.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00009500.plt', '/home/ext-zyou6474/Projects/lesgo_adjoint_tutorial_bundle/tests/2b_channel_flow_scalar_Qi/OUTPUT.13359817/ins_velo_00010000.plt'])
ins_velo_000.MeshStatus = ['SOULTION']
ins_velo_000.PointArrayStatus = []

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# Properties modified on ins_velo_000
ins_velo_000.PointArrayStatus = ['theta', 'u', 'v', 'w', 'x', 'y', 'z']

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
ins_velo_000Display.OSPRayScaleArray = 'theta'
ins_velo_000Display.OSPRayScaleFunction = 'PiecewiseFunction'
ins_velo_000Display.SelectOrientationVectors = 'None'
ins_velo_000Display.ScaleFactor = 0.6258641719818115
ins_velo_000Display.SelectScaleArray = 'None'
ins_velo_000Display.GlyphType = 'Arrow'
ins_velo_000Display.GlyphTableIndexArray = 'None'
ins_velo_000Display.GaussianRadius = 0.03129320859909058
ins_velo_000Display.SetScaleArray = ['POINTS', 'theta']
ins_velo_000Display.ScaleTransferFunction = 'PiecewiseFunction'
ins_velo_000Display.OpacityArray = ['POINTS', 'theta']
ins_velo_000Display.OpacityTransferFunction = 'PiecewiseFunction'
ins_velo_000Display.DataAxesGrid = 'GridAxesRepresentation'
ins_velo_000Display.PolarAxes = 'PolarAxesRepresentation'
ins_velo_000Display.ScalarOpacityUnitDistance = 0.03596971311433

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
ins_velo_000Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 48604.3203125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
ins_velo_000Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 48604.3203125, 1.0, 0.5, 0.0]

# reset view to fit data
renderView1.ResetCamera(False)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Iso Volume'
isoVolume1 = IsoVolume(registrationName='IsoVolume1', Input=ins_velo_000)
isoVolume1.InputScalars = ['POINTS', 'theta']
isoVolume1.ThresholdRange = [0.0, 48604.3203125]

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
isoVolume1Display.OSPRayScaleArray = 'theta'
isoVolume1Display.OSPRayScaleFunction = 'PiecewiseFunction'
isoVolume1Display.SelectOrientationVectors = 'None'
isoVolume1Display.ScaleFactor = 0.022763976454734804
isoVolume1Display.SelectScaleArray = 'None'
isoVolume1Display.GlyphType = 'Arrow'
isoVolume1Display.GlyphTableIndexArray = 'None'
isoVolume1Display.GaussianRadius = 0.00113819882273674
isoVolume1Display.SetScaleArray = ['POINTS', 'theta']
isoVolume1Display.ScaleTransferFunction = 'PiecewiseFunction'
isoVolume1Display.OpacityArray = ['POINTS', 'theta']
isoVolume1Display.OpacityTransferFunction = 'PiecewiseFunction'
isoVolume1Display.DataAxesGrid = 'GridAxesRepresentation'
isoVolume1Display.PolarAxes = 'PolarAxesRepresentation'
isoVolume1Display.ScalarOpacityUnitDistance = 0.02519202431106749
isoVolume1Display.OpacityArrayName = ['POINTS', 'theta']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
isoVolume1Display.ScaleTransferFunction.Points = [0.009999999776482582, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
isoVolume1Display.OpacityTransferFunction.Points = [0.009999999776482582, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(ins_velo_000)

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=ins_velo_000)
clip1.ClipType = 'Plane'
clip1.HyperTreeGridClipper = 'Plane'
clip1.Scalars = ['POINTS', 'theta']
clip1.Value = 24302.16015625

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [3.1293208599090576, 1.5646604299545288, 0.9870538115501404]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [3.1293208599090576, 1.5646604299545288, 0.9870538115501404]

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
clip1Display.OSPRayScaleArray = 'theta'
clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
clip1Display.SelectOrientationVectors = 'None'
clip1Display.ScaleFactor = 0.6258641719818115
clip1Display.SelectScaleArray = 'None'
clip1Display.GlyphType = 'Arrow'
clip1Display.GlyphTableIndexArray = 'None'
clip1Display.GaussianRadius = 0.03129320859909058
clip1Display.SetScaleArray = ['POINTS', 'theta']
clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
clip1Display.OpacityArray = ['POINTS', 'theta']
clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
clip1Display.DataAxesGrid = 'GridAxesRepresentation'
clip1Display.PolarAxes = 'PolarAxesRepresentation'
clip1Display.ScalarOpacityUnitDistance = 0.04199814252064949
clip1Display.OpacityArrayName = ['POINTS', 'theta']

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
clip1Display.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 48604.3203125, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
clip1Display.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 48604.3203125, 1.0, 0.5, 0.0]

# update the view to ensure updated data information
renderView1.Update()

# reset view to fit data
renderView1.ResetCamera(False)

# set scalar coloring
ColorBy(clip1Display, ('POINTS', 'u'))

# rescale color and/or opacity maps used to include current data range
clip1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
clip1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'u'
uLUT = GetColorTransferFunction('u')
uLUT.RGBPoints = [0.0770508348941803, 0.231373, 0.298039, 0.752941, 11.30783922970295, 0.865003, 0.865003, 0.865003, 22.53862762451172, 0.705882, 0.0156863, 0.14902]
uLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'u'
uPWF = GetOpacityTransferFunction('u')
uPWF.Points = [0.0770508348941803, 0.0, 0.5, 0.0, 22.53862762451172, 1.0, 0.5, 0.0]
uPWF.ScalarRangeInitialized = 1

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(1928, 1153)

# current camera placement for renderView1
renderView1.CameraPosition = [-2.5320411875824074, -8.995654025651632, 8.315507661076825]
renderView1.CameraFocalPoint = [3.129320859909057, 1.5646604299545286, 0.9870538115501405]
renderView1.CameraViewUp = [0.38167822560517645, 0.38007114553793725, 0.842536442196042]
renderView1.CameraParallelScale = 3.6352560476840026

# save animation
SaveAnimation('/home/zejiany/Documents/SynologyDrive/SDSU_Onedrive/research/lesgo_outputs/remote_test/theta.png', renderView1, ImageResolution=[1928, 1153],
    TransparentBackground=1,
    FrameWindow=[0, 20])

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(1928, 1153)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [-2.5320411875824074, -8.995654025651632, 8.315507661076825]
renderView1.CameraFocalPoint = [3.129320859909057, 1.5646604299545286, 0.9870538115501405]
renderView1.CameraViewUp = [0.38167822560517645, 0.38007114553793725, 0.842536442196042]
renderView1.CameraParallelScale = 3.6352560476840026

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).