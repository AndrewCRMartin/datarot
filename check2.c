#include<math.h>
#include<stdio.h>
#include "bioplib/MathType.h"
void regress(REAL *x, REAL *y, int ndata);

#ifdef SIMPLE
main(){
   int n,i;
   float x,y,m,c,d;
   float sumx=0,sumxsq=0,sumy=0,sumxy=0;
   printf("enter the number of values for n:");
   scanf("%d",&n);
   for(i=0;i<n;i++){
      printf("enter values of x and y");
      scanf("%f%f",&x,&y);
      sumx=sumx+x;
      sumxsq=sumxsq+(x*x);
      sumy=sumy+y;
      sumxy=sumxy+(x*y);
   }
   d=n*sumxsq-sumx*sumx;
   m=(n*sumxy-sumx*sumy)/d;
   c=(sumy*sumxsq-sumx*sumxy)/d;
   printf("M=%f\tC=%f\n",m,c);
}
#endif

#define TOL 0.1
#define MAXDATA 5000
#define MAXBUFF 160



int main(int argc, char **argv)
{
   FILE *fp = NULL;
   
   if((fp=fopen("corrected.csv", "r"))!=NULL)
   {
      int ndata = 0;
      REAL xData[MAXDATA],
         yData[MAXDATA];

#ifdef TEST
      for(ndata=0; ndata<10; ndata++)
      {
         xData[ndata] = ndata;
/*         yData[ndata] = ndata; */
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

      printf("NData: %d\n", ndata);
      
      regress(xData, yData, ndata);

      
   }
   return(0);
   
}


void regress(REAL *x, REAL *y, int ndata)
{
   int i;
   REAL m,c,d;
   REAL sumx=0,sumxsq=0,sumy=0,sumxy=0;
   for(i=0;i<ndata;i++)
   {
      sumx=sumx+x[i];
      sumxsq=sumxsq+(x[i]*x[i]);
      sumy=sumy+y[i];
      sumxy=sumxy+(x[i]*y[i]);
   }
   d=ndata*sumxsq-sumx*sumx;
   m=(ndata*sumxy-sumx*sumy)/d;
   c=(sumy*sumxsq-sumx*sumxy)/d;
   printf("M=%f\tC=%f\n",m,c);

}
