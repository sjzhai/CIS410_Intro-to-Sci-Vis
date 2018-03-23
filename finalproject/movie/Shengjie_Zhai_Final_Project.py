sw = SaveWindowAttributes()
sw.family = 0
sw.width = 800
sw.height = 800
sw.format = sw.PNG

OpenDatabase("localhost:/Users/shengjiezhai/Desktop/finalproject/movie/dbreak3d_data/dbreak3d_boundries.silo", 0)

##add subset domains
AddPlot("Subset", "domains", 1, 1)
SetActivePlots(0)
silr = SILRestriction()
silr.SuspendCorrectnessChecking()
silr.TurnOnAll()

##close walls for the subset domain
for silSet in (1,3,5,6):
    silr.TurnOffSet(silSet)
silr.EnableCorrectnessChecking()
SetPlotSILRestriction(silr ,1)
DrawPlots()

##change colors for subset domain
SubsetAtts = SubsetAttributes()
SubsetAtts.colorType = SubsetAtts.ColorBySingleColor
SubsetAtts.singleColor = (153, 204, 255, 255)
SubsetAtts.SetMultiColor(0, (255, 0, 0, 255))
SubsetAtts.SetMultiColor(1, (0, 255, 0, 255))
SubsetAtts.SetMultiColor(2, (0, 0, 255, 255))
SubsetAtts.SetMultiColor(3, (0, 255, 255, 255))
SubsetAtts.SetMultiColor(4, (255, 0, 255, 255))
SubsetAtts.SetMultiColor(5, (255, 255, 0, 255))
SetPlotOptions(SubsetAtts)

##add Pseudocolor
OpenDatabase("localhost:/Users/shengjiezhai/Desktop/finalproject/movie/dbreak3d_data/dbreak3d_fluid.visit", 0)
AddPlot("Pseudocolor", "alpha1", 1, 1)

##add operator--Isovolume
AddOperator("Isovolume", 0)
IsovolumeAtts = IsovolumeAttributes()
IsovolumeAtts.lbound = 0.5
IsovolumeAtts.ubound = 1e+37
IsovolumeAtts.variable = "default"
SetOperatorOptions(IsovolumeAtts, 0)
DrawPlots()

##change color for Isovolume
PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.colorTableName = "GnBu"
PseudocolorAtts.opacityType = PseudocolorAtts.Constant  # ColorTable, FullyOpaque, Constant, Ramp, VariableRange
PseudocolorAtts.opacityVariable = ""
PseudocolorAtts.opacity = 0.654902
PseudocolorAtts.opacityVarMin = 0
PseudocolorAtts.opacityVarMax = 1
PseudocolorAtts.opacityVarMinFlag = 0
PseudocolorAtts.opacityVarMaxFlag = 0
SetPlotOptions(PseudocolorAtts)

##adjust information
AnnotationAtts = AnnotationAttributes()
AnnotationAtts.axes3D.visible = 0
AnnotationAtts.axes3D.triadFlag = 0
AnnotationAtts.axes3D.bboxFlag = 1
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.timeInfoFlag = 1
AnnotationAtts.legendInfoFlag = 0
AnnotationAtts.axesArray.axes.title.visible = 1
AnnotationAtts.axesArray.axes.label.visible = 1
AnnotationAtts.axesArray.axes.tickMarks.visible = 1

AnnotationAtts.backgroundColor = (255, 255, 255, 255)
AnnotationAtts.foregroundColor = (153, 204, 255, 255)
AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
AnnotationAtts.gradientColor1 = (192, 192, 192, 255)
AnnotationAtts.gradientColor2 = (51, 51, 51, 255)
AnnotationAtts.backgroundMode = AnnotationAtts.Gradient
SetAnnotationAttributes(AnnotationAtts)

##modify and save windows images for first scene
i=0
m=0
offset=0.1
up = 0
angle = 20
left = 0
boolean = False
View3DAtts = View3DAttributes()

text0 = CreateAnnotationObject("Text2D")
text0.height = 0.02
text0.position = (0.45, 0.95)
text0.text = "Water flow"

while(m<=160):
    View3DAtts.viewNormal = (left, up, 1)
    View3DAtts.focus = (0.29246, 0.29146, 0.219095)
    View3DAtts.viewAngle = angle
    SetView3D(View3DAtts)

    for j in range(0, 4):
        SetTimeSliderState(m)
        sw.fileName="frames/m%02d" %(i+j)
        SetSaveWindowAttributes(sw)
        SaveWindow()
    i+=3
    m+=1
    if(m>32):
        boolean = True
    if (boolean == False):
        up+=offset
        angle+=10*offset
    if (boolean == True):
        if (m>64):
            left+=offset
        else:
            up-=offset
            angle-=10*offset
            left-=offset
text0.Delete()

#add plot vector
AddPlot("Vector", "U", 1, 0)
VectorAtts = VectorAttributes()
VectorAtts.glyphLocation = VectorAtts.AdaptsToMeshResolution  # AdaptsToMeshResolution, UniformInSpace
VectorAtts.useStride = 1
VectorAtts.stride = 70
VectorAtts.nVectors = 400
VectorAtts.lineStyle = VectorAtts.SOLID  # SOLID, DASH, DOT, DOTDASH
VectorAtts.scale = 0.75
VectorAtts.headSize = 0.25
VectorAtts.vectorColor = (0, 0, 0, 255)
VectorAtts.colorTableName = "RdYlBu"
VectorAtts.vectorOrigin = VectorAtts.Tail  # Head, Middle, Tail
VectorAtts.limitsMode = VectorAtts.OriginalData  # OriginalData, CurrentPlot
VectorAtts.lineStem = VectorAtts.Line  # Cylinder, Line
VectorAtts.geometryQuality = VectorAtts.High  # Fast, High
VectorAtts.glyphType = VectorAtts.Arrow  # Arrow, Ellipsoid
SetPlotOptions(VectorAtts)
DrawPlots()

##adjust information
AnnotationAtts = AnnotationAttributes()
AnnotationAtts.userInfoFlag = 0
AnnotationAtts.timeInfoFlag = 1
AnnotationAtts.axesArray.axes.tickMarks.visible = 0

AnnotationAtts.backgroundColor = (255, 255, 255, 255)
AnnotationAtts.foregroundColor = (0, 0, 0, 255)
AnnotationAtts.gradientBackgroundStyle = AnnotationAtts.Radial  # TopToBottom, BottomToTop, LeftToRight, RightToLeft, Radial
AnnotationAtts.gradientColor1 = (255, 255, 255, 255)
AnnotationAtts.gradientColor2 = (192, 192, 192, 255)
AnnotationAtts.backgroundMode = AnnotationAtts.Gradient
SetAnnotationAttributes(AnnotationAtts)

##modify and save windows images for second scene
n = 0
right = 0
left = 1
up = 0
angle = 20
boolean = False
View3DAtts = View3DAttributes()

#modify 2D text on the views
text1 = CreateAnnotationObject("Text2D")
text1.height = 0.02
text1.position = (0.35, 0.95)
text1.text = "Vector Plot on Water flow"
text1.visible = 1

text2 = CreateAnnotationObject("Text2D")
text2.height = 0.02
text2.position = (0.15, 0.95)
text2.text = "The actual extents are (0.0073, 0.5773, 0.0119997, 0.584, 0, 0.438)"
text2.visible = 0

text3 = CreateAnnotationObject("Text2D")
text3.height = 0.02
text3.position = (0.25, 0.1)
text3.text = "Vectors are in the same direction as water waves"
text3.visible = 0

text4 = CreateAnnotationObject("Text2D")
text4.height = 0.02
text4.position = (0.15, 0.94)
text4.text = "The directions of vectors in the water waves are going up"
text4.visible = 0

text5 = CreateAnnotationObject("Text2D")
text5.height = 0.02
text5.position = (0.15, 0.1)
text5.text = "The direction of vectors still follow the direction of the water waves"
text5.visible = 0

while(n<=100):
    View3DAtts.viewNormal = (right, up, left)
    View3DAtts.focus = (0.29246, 0.29146, 0.219095)
    View3DAtts.viewAngle = angle
    SetView3D(View3DAtts)
    # static for 5 frames
    if (n == 0):
        for s in range(5):
            SetView3D(View3DAtts)
            SetTimeSliderState(n)
            sw.fileName="frames/m%02d" %(i+s)
            SetSaveWindowAttributes(sw)
            SaveWindow()
        i+=5
    if (n == 20):
        text3.visible = 0
        for k in range(200):
            View3DAtts.viewNormal = (right, up, left)
            View3DAtts.viewAngle = angle
            if (k < 40):
                angle+=1
                right+=0.04 #0.077
                up+=0.05  #0.093
                left-=0.05   #0.0923
            elif (k >= 40 and k <= 100):
                text1.visible = 0
                text2.visible = 1
            elif (k > 100 and k <= 160):
                text2.position = (0.15, 0.97)
                text4.visible = 1
            elif (k > 160):
                text1.visible = 1
                text2.Delete()
                text4.Delete()
                angle-=1
                right-=0.04
                up-=0.05
                left+=0.05
            SetView3D(View3DAtts)
            SetTimeSliderState(n)
            sw.fileName="frames/m%02d" %(i+k)
            SetSaveWindowAttributes(sw)
            SaveWindow()
        i+=200
    elif (n != 20):
        if (n > 5 and n < 20):
            text3.visible = 1
        elif (n > 20 and n < 40):
            text5.visible = 1
        elif (n >= 50):
            text5.visible = 0
        for j in range(5):
            SetTimeSliderState(n)
            sw.fileName="frames/m%02d" %(i+j)
            SetSaveWindowAttributes(sw)
            SaveWindow()
        i+=5
    n+=1
    if(n>32):
        boolean = True
    if(n>96):
        boolean = False
    if (boolean == False):
        right+=offset
    if (boolean == True):
        right-=offset

import sys
sys.exit()
