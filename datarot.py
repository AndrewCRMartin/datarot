#!/usr/bin/python3
import os
import sys
import math

SMALL = 0.001
m = 0.285
c = -32.477


def FindIntersection(m1, c1, m2, c2):
   """
      Line 1: y=m1*x + c1
      Line 2: y=m2*x + c2
      
      When y matches,
      m1*x + c1 = m2*x + c2
      So...
      m1*x - m2*x = c2 - c1
      x(m1-m2) = c2 - c1
      x = (c2-c1)/(m1-m2)
      
      and
      y = m1*x + c1
   """

   if(abs(m1-m2) < SMALL):
       return(0.0, 0.0, False)

   x = (c2-c1)/(m1-m2)
   y = m1*x + c1
   
   return(x, y, True)


def TranslateData(xData, yData, ndata, transX, transY):
    for i in range(ndata):
        xData[i] = xData[i] + transX
        yData[i] = yData[i] + transY
    return(xData, yData)


def FindAngle(m1, m2):
    """ $ tan(\theta) = (m1-m2)/(1 + m1*m2) $ """
    if(abs(m1-m2) < SMALL):
        return(0.0)

    tanTheta = (m1-m2)/(1 + m1*m2)
    return(math.atan(tanTheta))


def RotateData(xData, yData, nData, theta):
    for i in range(nData):
        x = xData[i] * math.cos(theta) - yData[i] * math.sin(theta)
        y = xData[i] * math.sin(theta) + yData[i] * math.cos(theta)
        xData[i] = x
        yData[i] = y
    return(xData, yData)

#if(__name__ == '__MAIN__'):
if(True):
    # Main program
    with open ("test/Everything_NR2_SklearnGBReg.csv") as file:
        ndata = 0
        xData = []
        xDataOrig = []
        yData = []
    
        for buffer in file:
            if(buffer[0] != '#'):
                fields = buffer.split(',')
                x = float(fields[1])
                y = float(fields[2])
                xData.append(x)
                xDataOrig.append(x)
                yData.append(y)
                ndata = ndata+1

        print("Original Data:")
      
        for i in range(ndata):
            print("%f,%f" % (xData[i], yData[i]))

        intersectX, intersectY, ok = FindIntersection(m, c, 1, 0)
        if(ok):
            xData, yData = TranslateData(xData, yData, ndata,
                                         -intersectX, -intersectY)
            angle = FindAngle(m, 1)
            xData, yData = RotateData(xData, yData, ndata, -angle)
            xData, yData = TranslateData(xData, yData, ndata,
                                         intersectX, intersectY)

        print("Int: %.3f %.3f" % (intersectX, intersectY))
        print("Ang: %.3f" % (180.0 * angle / math.pi))

        print("Rotated Data:")
        for i in range(ndata):
            print("%f,%f" % (xData[i], yData[i]))

        print("Rotated Y Data (X original):")
        for i in range(ndata):
            print("%f,%f" % (xDataOrig[i], yData[i]))


