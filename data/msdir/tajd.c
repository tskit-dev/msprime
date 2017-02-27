/************************* tajima.c *************************************************************
 This program calculates Tajima's D when number of sequences, number of segregating sites,
   and average pairwise differences (pi) are known.  It also reports all the coefficients for Tajima's
   D (a1, a2, b1, b2, c1, c2, e1, e2). 
**************************************************************************************************/


#include <stdio.h>
#include <math.h>



double a1f(int);
double a2f(int);
double b1f(int);
double b2f(int);
double c1f(double, double);
double c2f(int, double, double, double);
double e1f(double, double);
double e2f(double, double, double);


	double
tajd(int nsam, int segsites, double sumk)
{

double  a1, a2, b1, b2, c1, c2, e1, e2; 
 
if( segsites == 0 ) return( 0.0) ;

a1 = a1f(nsam);
a2 = a2f(nsam);
b1 = b1f(nsam);
b2 = b2f(nsam);
c1 = c1f(a1, b1);
c2 = c2f(nsam, a1, a2, b2);
e1 = e1f(a1, c1);
e2 = e2f(a1, a2, c2);

return( (sumk - (segsites/a1))/sqrt((e1*segsites) + ((e2*segsites)*(segsites
-1))) ) ;

}

double a1f(int nsam)
{
double a1;
int i;
a1 = 0.0;
for (i=1; i<=nsam-1; i++) a1 += 1.0/i;
return (a1);
}


double a2f(int nsam) 
{
double a2;
int i;
a2 = 0.0;
for (i=1; i<=nsam-1; i++) a2 += 1.0/(i*i);
return (a2);
}


double b1f(int nsam)
{
double b1;
b1 = (nsam + 1.0)/(3.0*(nsam-1.0));
return (b1);
}


double b2f(int nsam)
{
double b2;
b2 = (2*(nsam*nsam + nsam + 3.0))/(9*nsam*(nsam - 1));
return (b2);
}


double e1f(double a1, double c1) 
{
double e1;
e1 = c1/a1;
return (e1);
}

double e2f(double a1, double a2, double c2)
{ 
double e2;
e2 = c2/((a1*a1)+a2);
return (e2);
}


double c1f(double a1, double b1)
{
double c1;
c1 = b1 - (1/a1);
return (c1);
}


double c2f(int nsam, double a1, double a2, double b2)
{
double c2;
c2 = b2 - ((nsam+2)/(a1*nsam)) + (a2/(a1 * a1));
return (c2);
}











