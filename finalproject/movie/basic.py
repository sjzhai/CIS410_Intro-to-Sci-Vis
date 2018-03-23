OpenDatabase("/Users/shengjiezhai/Desktop/finalproject/movie/proj7.vtk")
AddPlot("Pseudocolor", "hardyglobal")
DrawPlots()
sw = SaveWindowAttributes()
sw.family = 0
sw.width = 800
sw.height = 800
sw.format = sw.PPM
for i in range(30):
  sw.fileName="frames/m%02d" %(i)
  SetSaveWindowAttributes(sw)
  SaveWindow()
import sys
sys.exit()
