#define TEST
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "svdfit.h"
#include "bioplib/array.h"


#define TOL 0.1
#define MAXDATA 5000
#define MAXBUFF 160



int main(int argc, char **argv)
{
   FILE *fp = NULL;
   
   if((fp=fopen("corrected.csv", "r"))!=NULL)
   {
      char buffer[MAXBUFF];
      int ndata = 0,
         ma=2;
      REAL xData[MAXDATA],
         yData[MAXDATA],
         sig[MAXDATA],
         a[MAXDATA],
         **u,
         **v,
         w[MAXDATA],
         chisq;


      u= (REAL **)blArray2D(sizeof(REAL), MAXDATA, MAXDATA);
      v= (REAL **)blArray2D(sizeof(REAL), MAXDATA, MAXDATA);
      

#ifdef TEST
      for(ndata=0; ndata<10; ndata++)
      {
         xData[ndata] = ndata;
         yData[ndata] = ndata;
#         yData[ndata] = 2 * ndata + 5;
         sig[ndata] = 1;
      }
#else
      while(fgets(buffer, MAXBUFF, fp))
      {
         REAL x, y;
         
         sscanf(buffer, "%lf,%lf", &x, &y);
         xData[ndata] = x;
         yData[ndata] = y;
         sig[ndata] = 1;
         ndata++;
      }
#endif

      svdfit(xData,yData,sig, ndata, TOL, a, ma, u, v, w, MAXDATA, MAXDATA,
             &chisq, &TheFunction);

      printf("Slope: %.2f\n", a[0]);
      printf("Inter: %.2f\n", a[1]);
      
      
   }
   return(0);
   
}
