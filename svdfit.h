/*************************************************************************

   Program:    AMPlot2
   File:       svdfit.h
   
   Version:    V1.3
   Date:       27.07.94
   Function:   SVD Fitting routines
   
   Copyright:  SciTech Software 1992
   Author:     Dr. Andrew C. R. Martin
   Address:    SciTech Software
               23, Stag Leys,
               Ashtead,
               Surrey,
               KT21 2TD.
   Phone:      +44 (0372) 275775
   EMail:      UUCP: cbmuk!cbmuka!scitec!amartin
               JANET: andrew@uk.ac.ox.biop
               
**************************************************************************

   This program is copyright. Any copying without the permission of
   SciTech Software is illegal.

**************************************************************************

   Description:
   ============

**************************************************************************

   Usage:
   ======

**************************************************************************

   Revision History:
   =================

   V1.0  14.04.92 Tidied up from AMPlot V1.0.
   V1.1  27.04.92 Various fixes for dynamic memory allocation for big 
                  arrays.
   V1.2  01.05.92 Added routines to calculate variance and covariance 
                  matrices
   V1.3  27.07.94 Changed to REAL

*************************************************************************/
#ifndef _SVDFIT_H
#define _SVDFIT_H

#include <math.h>
#include <stdlib.h>

#include "bioplib/MathType.h"
#include "bioplib/SysDefs.h"

void TheFunction(REAL x,REAL *p,int np);
REAL expand(REAL *p,REAL *a,int n);
int svdfit(REAL *x,REAL *y,REAL *sig,int ndata,REAL tol,REAL *a,
           int ma,REAL **u,REAL **v,REAL *w,int mp,int np,REAL *chisq,
           void (*TheFunction)(REAL xf, REAL *pf, int nf));

#endif
