#define DEBUG
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "svdfit.h"
#include "bioplib/array.h"
#include "bioplib/general.h"


#define SMALL 0.001
#define MAXDATA 5000
#define MAXBUFF 160

int main(int argc, char **argv);
BOOL FindIntersection(REAL m1, REAL c1, REAL m2, REAL c2,
		      REAL *intersectX, REAL *intersectY);
void TranslateData(REAL *xData, REAL *yData, int ndata, REAL transX, REAL transY);
REAL FindAngle(REAL m1, REAL m2);
void RotateData(REAL *xData, REAL *yData, int ndata, REAL angle);
void GetCSVFields(char *buffer,int field1, int field2, REAL *x, REAL *y);





int main(int argc, char **argv)
{
   FILE *fp = NULL;
   

   if((fp=fopen("test/Everything_NR2_SklearnGBReg.csv", "r"))!=NULL)
   {
      int ndata = 0, i;
      REAL xData[MAXDATA],
         yData[MAXDATA],
         xDataOrig[MAXDATA],
         m = 0.285,
         c = -32.477,
         intersectX,
         intersectY;
      
#ifdef TEST
      ndata = 0;
      
      for(i=-60; i<-30; i+=5)
      {
         xData[ndata] = i;
         xDataOrig[ndata] = i;
         /*yData[ndata] = i;*/
	 yData[ndata] = 2 * i + 5;
         ndata++;
      }
#else
      char buffer[MAXBUFF];
      while(fgets(buffer, MAXBUFF, fp))
      {
         if(buffer[0] != '#')
         {
            REAL x, y;
            GetCSVFields(buffer,1,2,&x,&y);
            xData[ndata] = x;
            xDataOrig[ndata] = x;
            yData[ndata] = y;
            ndata++;
         }
      }
#endif

#ifdef DEBUG
      printf("Original Data:\n");
      
      for (i=0; i<ndata; i++)
      {
         printf("%f,%f\n", xData[i], yData[i]);
      }
#endif
      
      if(FindIntersection(m, c, 1, 0, &intersectX, &intersectY))
      {
         REAL angle;
         
         TranslateData(xData, yData, ndata, -intersectX, -intersectY);
         angle = FindAngle(m, 1);
         RotateData(xData, yData, ndata, -angle);
         TranslateData(xData, yData, ndata, intersectX, intersectY);
#ifdef DEBUG
         printf("Int: %.3f %.3f\n", intersectX, intersectY);
         printf("Ang: %.3f\n", 180*angle/PI);
#endif
         
      }
      printf("Rotated Data:\n");
      for (i=0; i<ndata; i++)
      {
         printf("%f,%f\n", xData[i], yData[i]);
      }
      printf("Rotated Y Data (X original):\n");
      for (i=0; i<ndata; i++)
      {
         printf("%f,%f\n", xDataOrig[i], yData[i]);
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

   if(ABS(m1-m2) < SMALL)
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

REAL FindAngle(REAL m1, REAL m2)
{
   REAL tanTheta;
   
   /* $ tan(\theta) = (m1-m2)/(1 + m1*m2) $ */
   if(ABS(m1-m2) < SMALL)
   {
      return(0.0);
   }

   tanTheta = (m1-m2)/(1 + m1*m2);
   return(atan(tanTheta));
}

void RotateData(REAL *xData, REAL *yData, int nData, REAL theta)
{
   int i;
   for(i=0; i<nData; i++)
   {
      REAL x, y;
      x = xData[i] * cos(theta) - yData[i] * sin(theta);
      y = xData[i] * sin(theta) + yData[i] * cos(theta);
      xData[i] = x;
      yData[i] = y;
   }
}

void GetCSVFields(char *buffer,int field1, int field2, REAL *x, REAL *y)
{
   int field = 0;
   char *chp,
      word[MAXBUFF];
   
   for(chp=buffer; chp!=NULL; )
   {
      chp = blGetWord(chp, word, MAXBUFF);
      if(field == field1)
      {
         sscanf(word, "%lf", x);
      }
      else if(field == field2)
      {
         sscanf(word, "%lf", y);
      }
      
      field++;
   }
   
}

