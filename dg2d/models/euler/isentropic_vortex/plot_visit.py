#----------------------------------------------------------------------------
# To run this script, use following command
#       visit -cli -s scriptname.py
#---------------------------------------------------------------------------

# import visit utilitis in python
from visit_utils import *

OpenGUI()
database = "solution.xdmf"
OpenDatabase(database)

DefineScalarExpression("Mach", 
                "sqrt(XVelocity^2+YVelocity^2)/sqrt(1.4*Pressure/Density)")

# Step 2: Add plots. Generally this step is followed by adding operator and 
# setting Attributes for plot and operator.
AddPlot("Pseudocolor", "Density")
p = PseudocolorAttributes()
p.minFlag = 1                      
p.maxFlag = 1
p.min = 0.5                           # set minimum value for variable
p.max = 1.0                           # set maximum value of variable
p.colorTableName = "hot_desaturated"  # set color table
SetPlotOptions(p)

#formatting legend object for contours
plotName= "pl"
legend = GetAnnotationObject(GetPlotList().GetPlots(0).plotName)
legend.managePosition = 0
legend.position = (0.05, 0.75)
legend.fontFamily = 2  # Arial, Courier, Times
legend.fontBold = 1
legend.drawBoundingBox = 0
legend.numberFormat = "%# -9.2g"
legend.drawLabels = 1 # None, Values, Labels, Both
legend.fontHeight = 0.02
legend.xScale = 1.0
legend.yScale = 1.5
legend.drawTitle = 1
legend.drawMinMax = 1
legend.orientation = 1  # VerticalRight, VerticalLeft, HorizontalTop, HorizontalBottom
legend.controlTicks = 1
legend.numTicks = 5          # set number of ticks
legend.minMaxInclusive = 1

AddPlot("Contour", "Density")
c = ContourAttributes()
c.lineWidth = 1                    # Set line width of contours from here
c.colorType = 0                    # (0, 1): This selects color by (single, multiple) color/s: 
c.singleColor = (0, 0, 0, 255)     # This combination is for black color
c.contourNLevels = 20              # Set number of contour levels from here.
c.legendFlag = 0                   # (0, 1): Disable/enable display of legend for contours
SetPlotOptions(c)

# Controls annotations of figure area
a = AnnotationAttributes()
a.axes2D.xAxis.title.visible=0     # (0, 1): Disable/enable display of axis title
a.axes2D.yAxis.title.visible=0
a.userInfoFlag = 0                 # disables username display
a.databaseInfoFlag = 0             # diables database display
SetAnnotationAttributes(a)

# Object to create title of figure
t = CreateAnnotationObject("Text2D")
t.visible = 1
t.active = 1
t.position = (0.5, 0.95)        # set position on the window
t.height = 0.02                 # set height of title text, sets font size.
t.textColor = (0, 0, 0, 255)    # set color of title text
t.useForegroundForTextColor = 1
t.text = "Time = $time"         # string to be used for title
t.fontFamily = 2                # (0, 1, 2) Arial, Courier, Times
t.fontBold = 1

# Step 3: Draw the plots
DrawPlots()                 # The desired plots must be added before this command

exit()
