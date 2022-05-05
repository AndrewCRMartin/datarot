#define TEST
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "svdfit.h"
#include "bioplib/array.h"


#define SMALL 0.001
#define MAXDATA 5000
#define MAXBUFF 160

int main(int argc, char **argv);
BOOL FindIntersection(REAL m1, REAL c1, REAL m2, REAL c2,
		      REAL *intersectX, REAL *intersectY);
void TranslateData(REAL *xData, REAL *yData, int ndata, REAL transX, REAL transY);
REAL FindAngle(REAL m1, REAL c1, REAL m2, REAL c2);
void RotateData(REAL *xData, REAL *yData, int ndata, REAL angle);








int main(int argc, char **argv)
{
   FILE *fp = NULL;
   
   if((fp=fopen("corrected.csv", "r"))!=NULL)
   {
      int ndata = 0, i;
      REAL xData[MAXDATA],
	yData[MAXDATA],
	m = 0.285,
	c = -32.477,
	intersectX,
	intersectY;

#ifdef TEST
      for(ndata=0; ndata<10; ndata++)
      {
         xData[ndata] = ndata;
         /*yData[ndata] = ndata;*/
	 yData[ndata] = 2 * ndata + 5;
      }
#else
      char buffer[MAXBUFF];
      while(fgets(buffer, MAXBUFF, fp))
      {
         REAL x, y;
         
         sscanf(buffer, "%lf,%lf", &x, &y);
         xData[ndata] = x;
         yData[ndata] = y;
         ndata++;
      }
#endif

      if(FindIntersection(m, c, 1, 0, &intersectX, &intersectY))
     {
       REAL angle;
       
       TranslateData(xData, yData, ndata, -intersectX, -intersectY);
       angle = FindAngle(m, 0, 1, 0);
       RotateData(xData, yData, ndata, angle);
       TranslateData(xData, yData, ndata, intersectX, intersectY);
     }
      for (i=0; i<ndata; i++)
      {
         printf("%f,%f\n", xData[i], yData[i]);
      }
   }
   return(0);
}


BOOL FindIntersection(REAL m1, REAL c1, REAL m2, REAL c2,
		      REAL *intersectX, REAL *intersectY)
{
  /* Line 1: y=m1*x + c1
     Line 2: y=m2*x + c2

     When y matches,
     m1*x + c1 = m2*x + c2
     So...
     m1*x - m2*x = c2 - c1
     x(m1-m2) = c2 - c1
     x = (c2-c1)/(m1-m2)

     and
     y = m1*x + c1;
  */
  REAL x, y;

  if(abs(m1-m2) < SMALL)
    {
      return(FALSE);
    }
  
  x = (c2-c1)/(m1-m2);
  y = m1*x + c1;

  *intersectX = x;
  *intersectY = y;
  return(TRUE);
}

void TranslateData(REAL *xData, REAL *yData, int ndata, REAL transX, REAL transY)
{
  int i;
  for(i=0; i<ndata; i++)
    {
      xData[i] += transX;
      yData[i] += transY;
    }
}

REAL FindAngle(REAL m1, REAL c1, REAL m2, REAL c2)
{
  return(0.0);
}

void RotateData(REAL *xData, REAL *yData, int ndata, REAL angle)
{
}

