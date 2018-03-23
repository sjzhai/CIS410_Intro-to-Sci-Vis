
def SetUpAnnotations():
  # Logging for SetAnnotationObjectOptions is not implemented yet.
  AnnotationAtts = AnnotationAttributes()
  AnnotationAtts.axes2D.visible = 1
  AnnotationAtts.axes2D.autoSetTicks = 1
  AnnotationAtts.axes2D.autoSetScaling = 1
  AnnotationAtts.axes2D.lineWidth = 0
  AnnotationAtts.axes2D.tickLocation = AnnotationAtts.axes2D.Outside  # Inside, Outside, Both
  AnnotationAtts.axes2D.tickAxes = AnnotationAtts.axes2D.BottomLeft  # Off, Bottom, Left, BottomLeft, All
  AnnotationAtts.axes2D.xAxis.title.visible = 1
  AnnotationAtts.axes2D.xAxis.title.font.font = AnnotationAtts.axes2D.xAxis.title.font.Courier  # Arial, Courier, Times
  AnnotationAtts.axes2D.xAxis.title.font.scale = 1
  AnnotationAtts.axes2D.xAxis.title.font.useForegroundColor = 1
  AnnotationAtts.axes2D.xAxis.title.font.color = (0, 0, 0, 255)
  AnnotationAtts.axes2D.xAxis.title.font.bold = 1
  AnnotationAtts.axes2D.xAxis.title.font.italic = 1
  AnnotationAtts.axes2D.xAxis.title.userTitle = 0
  AnnotationAtts.axes2D.xAxis.title.userUnits = 0
  AnnotationAtts.axes2D.xAxis.title.title = "X-Axis"
  AnnotationAtts.axes2D.xAxis.title.units = ""
  AnnotationAtts.axes2D.xAxis.label.visible = 1
  AnnotationAtts.axes2D.xAxis.label.font.font = AnnotationAtts.axes2D.xAxis.label.font.Courier  # Arial, Courier, Times
  AnnotationAtts.axes2D.xAxis.label.font.scale = 1
  AnnotationAtts.axes2D.xAxis.label.font.useForegroundColor = 1
  AnnotationAtts.axes2D.xAxis.label.font.color = (0, 0, 0, 255)
  AnnotationAtts.axes2D.xAxis.label.font.bold = 1
  AnnotationAtts.axes2D.xAxis.label.font.italic = 1
  AnnotationAtts.axes2D.xAxis.label.scaling = 0
  AnnotationAtts.axes2D.xAxis.tickMarks.visible = 1
  AnnotationAtts.axes2D.xAxis.tickMarks.majorMinimum = 0
  AnnotationAtts.axes2D.xAxis.tickMarks.majorMaximum = 1
  AnnotationAtts.axes2D.xAxis.tickMarks.minorSpacing = 0.02
  AnnotationAtts.axes2D.xAxis.tickMarks.majorSpacing = 0.2
  AnnotationAtts.axes2D.xAxis.grid = 0
  AnnotationAtts.axes2D.yAxis.title.visible = 1
  AnnotationAtts.axes2D.yAxis.title.font.font = AnnotationAtts.axes2D.yAxis.title.font.Courier  # Arial, Courier, Times
  AnnotationAtts.axes2D.yAxis.title.font.scale = 1
  AnnotationAtts.axes2D.yAxis.title.font.useForegroundColor = 1
  AnnotationAtts.axes2D.yAxis.title.font.color = (0, 0, 0, 255)
  AnnotationAtts.axes2D.yAxis.title.font.bold = 1
  AnnotationAtts.axes2D.yAxis.title.font.italic = 1
  AnnotationAtts.axes2D.yAxis.title.userTitle = 0
  AnnotationAtts.axes2D.yAxis.title.userUnits = 0
  AnnotationAtts.axes2D.yAxis.title.title = "Y-Axis"
  AnnotationAtts.axes2D.yAxis.title.units = ""
  AnnotationAtts.axes2D.yAxis.label.visible = 1
  AnnotationAtts.axes2D.yAxis.label.font.font = AnnotationAtts.axes2D.yAxis.label.font.Courier  # Arial, Courier, Times
  AnnotationAtts.axes2D.yAxis.label.font.scale = 1
  AnnotationAtts.axes2D.yAxis.label.font.useForegroundColor = 1
  AnnotationAtts.axes2D.yAxis.label.font.color = (0, 0, 0, 255)
  AnnotationAtts.axes2D.yAxis.label.font.bold = 1
  AnnotationAtts.axes2D.yAxis.label.font.italic = 1
  AnnotationAtts.axes2D.yAxis.label.scaling = 0
  AnnotationAtts.axes2D.yAxis.tickMarks.visible = 1
  AnnotationAtts.axes2D.yAxis.tickMarks.majorMinimum = 0
  AnnotationAtts.axes2D.yAxis.tickMarks.majorMaximum = 1
  AnnotationAtts.axes2D.yAxis.tickMarks.minorSpacing = 0.02
  AnnotationAtts.axes2D.yAxis.tickMarks.majorSpacing = 0.2
  AnnotationAtts.axes2D.yAxis.grid = 0
  AnnotationAtts.axes3D.visible = 1
  AnnotationAtts.axes3D.autoSetTicks = 1
  AnnotationAtts.axes3D.autoSetScaling = 1
  AnnotationAtts.axes3D.lineWidth = 0
  AnnotationAtts.axes3D.tickLocation = AnnotationAtts.axes3D.Inside  # Inside, Outside, Both
  AnnotationAtts.axes3D.axesType = AnnotationAtts.axes3D.ClosestTriad  # ClosestTriad, FurthestTriad, OutsideEdges, StaticTriad, StaticEdges
  AnnotationAtts.axes3D.triadFlag = 1
  AnnotationAtts.axes3D.bboxFlag = 1
  AnnotationAtts.axes3D.xAxis.title.visible = 1
  AnnotationAtts.axes3D.xAxis.title.font.font = AnnotationAtts.axes3D.xAxis.title.font.Arial  # Arial, Courier, Times
  AnnotationAtts.axes3D.xAxis.title.font.scale = 1
  AnnotationAtts.axes3D.xAxis.title.font.useForegroundColor = 1
  AnnotationAtts.axes3D.xAxis.title.font.color = (0, 0, 0, 255)
  AnnotationAtts.axes3D.xAxis.title.font.bold = 0
  AnnotationAtts.axes3D.xAxis.title.font.italic = 0
  AnnotationAtts.axes3D.xAxis.title.userTitle = 0
  AnnotationAtts.axes3D.xAxis.title.userUnits = 0
  AnnotationAtts.axes3D.xAxis.title.title = "X-Axis"
  AnnotationAtts.axes3D.xAxis.title.units = ""
  AnnotationAtts.axes3D.xAxis.label.visible = 1
  AnnotationAtts.axes3D.xAxis.label.font.font = AnnotationAtts.axes3D.xAxis.label.font.Arial  # Arial, Courier, Times
  AnnotationAtts.axes3D.xAxis.label.font.scale = 1
  AnnotationAtts.axes3D.xAxis.label.font.useForegroundColor = 1
  AnnotationAtts.axes3D.xAxis.label.font.color = (0, 0, 0, 255)
  AnnotationAtts.axes3D.xAxis.label.font.bold = 0
  AnnotationAtts.axes3D.xAxis.label.font.italic = 0
  AnnotationAtts.axes3D.xAxis.label.scaling = 0
  AnnotationAtts.axes3D.xAxis.tickMarks.visible = 1
  AnnotationAtts.axes3D.xAxis.tickMarks.majorMinimum = 0
  AnnotationAtts.axes3D.xAxis.tickMarks.majorMaximum = 1
  AnnotationAtts.axes3D.xAxis.tickMarks.minorSpacing = 0.02
  AnnotationAtts.axes3D.xAxis.tickMarks.majorSpacing = 0.2
  AnnotationAtts.axes3D.xAxis.grid = 0
  AnnotationAtts.axes3D.yAxis.title.visible = 1
  AnnotationAtts.axes3D.yAxis.title.font.font = AnnotationAtts.axes3D.yAxis.title.font.Arial  # Arial, Courier, Times
  AnnotationAtts.axes3D.yAxis.title.font.scale = 1
  AnnotationAtts.axes3D.yAxis.title.font.useForegroundColor = 1
  AnnotationAtts.axes3D.yAxis.title.font.color = (0, 0, 0, 255)
  AnnotationAtts.axes3D.yAxis.title.font.bold = 0
  AnnotationAtts.axes3D.yAxis.title.font.italic = 0
  AnnotationAtts.axes3D.yAxis.title.userTitle = 0
  AnnotationAtts.axes3D.yAxis.title.userUnits = 0
  AnnotationAtts.axes3D.yAxis.title.title = "Y-Axis"
  AnnotationAtts.axes3D.yAxis.title.units = ""
  AnnotationAtts.axes3D.yAxis.label.visible = 1
  AnnotationAtts.axes3D.yAxis.label.font.font = AnnotationAtts.axes3D.yAxis.label.font.Arial  # Arial, Courier, Times
  AnnotationAtts.axes3D.yAxis.label.font.scale = 1
  AnnotationAtts.axes3D.yAxis.label.font.useForegroundColor = 1
  AnnotationAtts.axes3D.yAxis.label.font.color = (0, 0, 0, 255)
  AnnotationAtts.axes3D.yAxis.label.font.bold = 0
  AnnotationAtts.axes3D.yAxis.label.font.italic = 0
  AnnotationAtts.axes3D.yAxis.label.scaling = 0
  AnnotationAtts.axes3D.yAxis.tickMarks.visible = 1
  AnnotationAtts.axes3D.yAxis.tickMarks.majorMinimum = 0
  AnnotationAtts.axes3D.yAxis.tickMarks.majorMaximum = 1
  AnnotationAtts.axes3D.yAxis.tickMarks.minorSpacing = 0.02
  AnnotationAtts.axes3D.yAxis.tickMarks.majorSpacing = 0.2
  AnnotationAtts.axes3D.yAxis.grid = 0
  AnnotationAtts.axes3D.zAxis.title.visible = 1
  AnnotationAtts.axes3D.zAxis.title.font.font = AnnotationAtts.axes3D.zAxis.title.font.Arial  # Arial, Courier, Times
  AnnotationAtts.axes3D.zAxis.title.font.scale = 1
  AnnotationAtts.axes3D.zAxis.title.font.useForegroundColor = 1
  AnnotationAtts.axes3D.zAxis.title.font.color = (0, 0, 0, 255)
  AnnotationAtts.axes3D.zAxis.title.font.bold = 0
  AnnotationAtts.axes3D.zAxis.title.font.italic = 0
  AnnotationAtts.axes3D.zAxis.title.userTitle = 0
  AnnotationAtts.axes3D.zAxis.title.userUnits = 0
  AnnotationAtts.axes3D.zAxis.title.title = "Z-Axis"
  AnnotationAtts.axes3D.zAxis.title.units = ""
  AnnotationAtts.axes3D.zAxis.label.visible = 1
  AnnotationAtts.axes3D.zAxis.label.font.font = AnnotationAtts.axes3D.zAxis.label.font.Arial  # Arial, Courier, Times
  AnnotationAtts.axes3D.zAxis.label.font.scale = 1
  AnnotationAtts.axes3D.zAxis.label.font.useForegroundColor = 1
  AnnotationAtts.axes3D.zAxis.label.font.color = (0, 0, 0, 255)
  AnnotationAtts.axes3D.zAxis.label.font.bold = 0
  AnnotationAtts.axes3D.zAxis.label.font.italic = 0
  AnnotationAtts.axes3D.zAxis.label.scaling = 0
  AnnotationAtts.axes3D.zAxis.tickMarks.visible = 1
  AnnotationAtts.axes3D.zAxis.tickMarks.majorMinimum = 0
  AnnotationAtts.axes3D.zAxis.tickMarks.majorMaximum = 1
  AnnotationAtts.axes3D.zAxis.tickMarks.minorSpacing = 0.02
  AnnotationAtts.axes3D.zAxis.tickMarks.majorSpacing = 0.2
  AnnotationAtts.axes3D.zAxis.grid = 0
  AnnotationAtts.axes3D.setBBoxLocation = 0
  AnnotationAtts.axes3D.bboxLocation = (0, 1, 0, 1, 0, 1)
  AnnotationAtts.userInfoFlag = 0
  AnnotationAtts.userInfoFont.font = AnnotationAtts.userInfoFont.Arial  # Arial, Courier, Times
  AnnotationAtts.userInfoFont.scale = 1
  AnnotationAtts.userInfoFont.useForegroundColor = 1
  AnnotationAtts.userInfoFont.color = (0, 0, 0, 255)
  AnnotationAtts.userInfoFont.bold = 0
  AnnotationAtts.userInfoFont.italic = 0
  AnnotationAtts.databaseInfoFlag = 0
  AnnotationAtts.timeInfoFlag = 1
  AnnotationAtts.databaseInfoFont.font = AnnotationAtts.databaseInfoFont.Arial  # Arial, Courier, Times
  AnnotationAtts.databaseInfoFont.scale = 1
  AnnotationAtts.databaseInfoFont.useForegroundColor = 1
  AnnotationAtts.databaseInfoFont.color = (0, 0, 0, 255)
  AnnotationAtts.databaseInfoFont.bold = 0
  AnnotationAtts.databaseInfoFont.italic = 0
  AnnotationAtts.databaseInfoExpansionMode = AnnotationAtts.File  # File, Directory, Full, Smart, SmartDirectory
  AnnotationAtts.databaseInfoTimeScale = 1
  AnnotationAtts.databaseInfoTimeOffset = 0
  AnnotationAtts.legendInfoFlag = 1
  AnnotationAtts.backgroundColor = (255, 255, 255, 255)
  AnnotationAtts.foregroundColor = (0, 0, 0, 255)
  AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
  AnnotationAtts.gradientColor1 = (0, 0, 255, 255)
  AnnotationAtts.gradientColor2 = (0, 0, 0, 255)
  AnnotationAtts.backgroundMode = AnnotationAtts.Solid  # Solid, Gradient, Image, ImageSphere
  AnnotationAtts.backgroundImage = ""
  AnnotationAtts.imageRepeatX = 1
  AnnotationAtts.imageRepeatY = 1
  AnnotationAtts.axesArray.visible = 1
  AnnotationAtts.axesArray.ticksVisible = 1
  AnnotationAtts.axesArray.autoSetTicks = 1
  AnnotationAtts.axesArray.autoSetScaling = 1
  AnnotationAtts.axesArray.lineWidth = 0
  AnnotationAtts.axesArray.axes.title.visible = 1
  AnnotationAtts.axesArray.axes.title.font.font = AnnotationAtts.axesArray.axes.title.font.Arial  # Arial, Courier, Times
  AnnotationAtts.axesArray.axes.title.font.scale = 1
  AnnotationAtts.axesArray.axes.title.font.useForegroundColor = 1
  AnnotationAtts.axesArray.axes.title.font.color = (0, 0, 0, 255)
  AnnotationAtts.axesArray.axes.title.font.bold = 0
  AnnotationAtts.axesArray.axes.title.font.italic = 0
  AnnotationAtts.axesArray.axes.title.userTitle = 0
  AnnotationAtts.axesArray.axes.title.userUnits = 0
  AnnotationAtts.axesArray.axes.title.title = ""
  AnnotationAtts.axesArray.axes.title.units = ""
  AnnotationAtts.axesArray.axes.label.visible = 1
  AnnotationAtts.axesArray.axes.label.font.font = AnnotationAtts.axesArray.axes.label.font.Arial  # Arial, Courier, Times
  AnnotationAtts.axesArray.axes.label.font.scale = 1
  AnnotationAtts.axesArray.axes.label.font.useForegroundColor = 1
  AnnotationAtts.axesArray.axes.label.font.color = (0, 0, 0, 255)
  AnnotationAtts.axesArray.axes.label.font.bold = 0
  AnnotationAtts.axesArray.axes.label.font.italic = 0
  AnnotationAtts.axesArray.axes.label.scaling = 0
  AnnotationAtts.axesArray.axes.tickMarks.visible = 1
  AnnotationAtts.axesArray.axes.tickMarks.majorMinimum = 0
  AnnotationAtts.axesArray.axes.tickMarks.majorMaximum = 1
  AnnotationAtts.axesArray.axes.tickMarks.minorSpacing = 0.02
  AnnotationAtts.axesArray.axes.tickMarks.majorSpacing = 0.2
  AnnotationAtts.axesArray.axes.grid = 0
  SetAnnotationAttributes(AnnotationAtts)
  
SetUpAnnotations()

actually_do_a_save=0
def MySaveWindow():
  global actually_do_a_save
  if (actually_do_a_save):
     SaveWindow()

a = CreateAnnotationObject("Text2D")
a.position = (0.45, 0.9)
a.text = "Y SLICE"

frame_idx=0

def PauseOneSecond():
  global frame_idx
  for i in range(30):
    sw.fileName="frames/m%02d" %(frame_idx)
    frame_idx = frame_idx + 1
    SetSaveWindowAttributes(sw)
    MySaveWindow()
   
OpenDatabase("~/410/proj7.vtk")
AddPlot("Pseudocolor", "hardyglobal")
AddOperator("Slice")
s = SliceAttributes()
DrawPlots()
sw = SaveWindowAttributes()
sw.family = 0
sw.width = 800
sw.height = 800
sw.format = sw.PPM

PauseOneSecond()

for i in range(60):
  t=i/59.0
  s.originIntercept = -10 + 20*t
  SetOperatorOptions(s)
  sw.fileName="frames/m%02d" %(frame_idx)
  frame_idx = frame_idx+1
  SetSaveWindowAttributes(sw)
  MySaveWindow()

PauseOneSecond()

s.axisType = s.XAxis
a.text = "X SLICE"
for i in range(60):
  t=i/59.0
  s.originIntercept = -10 + 20*t
  SetOperatorOptions(s)
  sw.fileName="frames/m%02d" %(frame_idx)
  frame_idx = frame_idx+1
  SetSaveWindowAttributes(sw)
  MySaveWindow()

PauseOneSecond()

s.axisType = s.ZAxis
a.text = "Z SLICE"
for i in range(60):
  t=i/59.0
  s.originIntercept = -10 + 20*t
  SetOperatorOptions(s)
  sw.fileName="frames/m%02d" %(frame_idx)
  frame_idx = frame_idx+1
  SetSaveWindowAttributes(sw)
  MySaveWindow()

PauseOneSecond()

import sys
#sys.exit()
