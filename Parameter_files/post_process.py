try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

Custom_0028_vtk = GetActiveSource()
Threshold1 = Threshold()

Threshold1.Scalars = ['POINTS', 'pressure']
Threshold1.ThresholdRange = [-0.81475600000000004, 0.91933500000000001]

Threshold1.Scalars = ['POINTS', 'geom']
Threshold1.ThresholdRange = [0.17000000000000001, 1.0]

AnimationScene1 = GetAnimationScene()
RenderView1 = CreateRenderView()
RenderView1.HeadLightWarmth = 1.0
RenderView1.KeyLightIntensity = 1.0
RenderView1.UseLight = 1
RenderView1.FillLightWarmth = 1.0
RenderView1.CameraPosition = [7.4850001335144043, 0.73500001430511475, 29.058915388621521]
RenderView1.LightSwitch = 1
RenderView1.CameraClippingRange = [28.768326234735305, 29.494799119450846]
RenderView1.ViewTime = 0.0
RenderView1.LODThreshold = 5.0
RenderView1.BackLightWarmth = 1.0
RenderView1.CenterOfRotation = [7.4850001335144043, 0.73500001430511475, 0.0]
RenderView1.MaintainLuminance = 1
RenderView1.CameraFocalPoint = [7.4850001335144043, 0.73500001430511475, 0.0]
RenderView1.RenderInterruptsEnabled = 0
RenderView1.CameraParallelScale = 7.521000732597968
RenderView1.KeyLightWarmth = 1.0

a1_pressure_PVLookupTable = GetLookupTableForArray( "pressure", 1 )

DataRepresentation2 = Show()
DataRepresentation2.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation2.ColorAttributeType = 'POINT_DATA'
DataRepresentation2.ScalarOpacityFunction = []
DataRepresentation2.ColorArrayName = 'pressure'
DataRepresentation2.ScalarOpacityUnitDistance = 0.57183913878223691
DataRepresentation2.LookupTable = a1_pressure_PVLookupTable

AnimationScene1.ViewModules = [ a2DRenderView1, RenderView1 ]

Delete(a2DRenderView1)
DataRepresentation1 = GetDisplayProperties(Custom_0028_vtk)
Delete(DataRepresentation1)
AnimationScene1.ViewModules = RenderView1

SetActiveSource(Custom_0028_vtk)
Threshold2 = Threshold()

Threshold2.Scalars = ['POINTS', 'pressure']
Threshold2.ThresholdRange = [-0.81475600000000004, 0.91933500000000001]

Threshold2.Scalars = ['POINTS', 'geom']
Threshold2.ThresholdRange = [0.0, 0.84999999999999998]

DataRepresentation3 = Show()
DataRepresentation3.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation3.ColorAttributeType = 'POINT_DATA'
DataRepresentation3.ScalarOpacityFunction = []
DataRepresentation3.ColorArrayName = 'pressure'
DataRepresentation3.ScalarOpacityUnitDistance = 0.4135471980049702
DataRepresentation3.LookupTable = a1_pressure_PVLookupTable

DataRepresentation3.ColorArrayName = ''
DataRepresentation3.DiffuseColor = [0.66666666666666663, 0.0, 0.0]

SetActiveSource(Threshold1)
StreamTracer1 = StreamTracer( SeedType="Point Source" )

ScalarBarWidgetRepresentation1 = CreateScalarBar( Orientation='Horizontal', Title='pressure', Position2=[0.49999999999999972, 0.12999999999999989], Enabled=1, LabelFontSize=8, LookupTable=a1_pressure_PVLookupTable, TitleFontSize=10, Position=[0.25264106050305918, 0.15444121071012812] )
GetRenderView().Representations.append(ScalarBarWidgetRepresentation1)

a1_pressure_PVLookupTable.RGBPoints = [-0.092078699999999999, 0.0, 0.0, 1.0, 0.91933500000000001, 1.0, 0.0, 0.0]

StreamTracer1.SeedType.Center = [7.4850001335144043, 0.73500001430511475, 0.0]
StreamTracer1.SeedType.Radius = 1.4970000267028809
StreamTracer1.Vectors = ['POINTS', 'velocity']
StreamTracer1.SeedType = "Point Source"
StreamTracer1.MaximumStreamlineLength = 14.970000267028809

StreamTracer1.SeedType.NumberOfPoints = 100

StreamTracer1.MaximumStreamlineLength = 14.9700002670288
StreamTracer1.SeedType = "High Resolution Line Source"

DataRepresentation4 = Show()
DataRepresentation4.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation4.ColorAttributeType = 'POINT_DATA'
DataRepresentation4.ColorArrayName = 'pressure'
DataRepresentation4.Texture = []
DataRepresentation4.LookupTable = a1_pressure_PVLookupTable

DataRepresentation2.Visibility = 0

StreamTracer1.SeedType.Point2 = [7.5, 0.0, 0.0]
StreamTracer1.SeedType.Resolution = 50
StreamTracer1.SeedType.Point1 = [15.0, 1.5, 0.0]

DataRepresentation2.Visibility = 1

SetActiveSource(Threshold1)
StreamTracer2 = StreamTracer( SeedType="Point Source" )

DataRepresentation4.ColorArrayName = ''
DataRepresentation4.DiffuseColor = [0.0, 0.0, 0.0]

StreamTracer2.SeedType.Center = [7.4850001335144043, 0.73500001430511475, 0.0]
StreamTracer2.SeedType.Radius = 1.4970000267028809
StreamTracer2.Vectors = ['POINTS', 'velocity']
StreamTracer2.SeedType = "Point Source"
StreamTracer2.MaximumStreamlineLength = 14.970000267028809

StreamTracer2.SeedType.NumberOfPoints = 100

StreamTracer2.MaximumStreamlineLength = 14.9700002670288
StreamTracer2.SeedType = "High Resolution Line Source"

DataRepresentation5 = Show()
DataRepresentation5.EdgeColor = [0.0, 0.0, 0.50000762951094835]
DataRepresentation5.ColorAttributeType = 'POINT_DATA'
DataRepresentation5.ColorArrayName = 'pressure'
DataRepresentation5.Texture = []
DataRepresentation5.LookupTable = a1_pressure_PVLookupTable

DataRepresentation2.Visibility = 0

StreamTracer2.SeedType.Point2 = [15.0, 0.75, 0.0]
StreamTracer2.SeedType.Resolution = 25
StreamTracer2.SeedType.Point1 = [0.0, 1.5, 0.0]

StreamTracer1.SeedType.Resolution = 25

DataRepresentation2.Visibility = 1

StreamTracer2.SeedType.Point2 = [0.0, 0.75, 0.0]
StreamTracer2.SeedType.Point1 = [15.0, 1.5, 0.0]

DataRepresentation5.ColorArrayName = ''
DataRepresentation5.DiffuseColor = [0.0, 0.0, 0.0]

Render()
