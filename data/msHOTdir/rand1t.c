/*  Link in this file for random number generation using drand48() */

#include <stdio.h>
#include <stdlib.h>

         double
ran1()
{
        double drand48();
        return( drand48() );
}               


	void seedit( char *flag )
{
	FILE *fopen(), *pfseed;
	unsigned short seedv[3], seedv2[3],  *seed48(), *pseed ;
	int i;


	if( flag[0] == 's' ){
	  seedv[0] =  (unsigned short)time( NULL ) ;
            seedv[1] = 27011; seedv[2] = 59243; 
          seed48( seedv );   

       printf("\n%d %d %d\n", seedv[0], seedv[1], seedv[2] );    
	}
}
