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

	if( flag[0] == 's' ) {
	   pfseed = fopen("seedms","r");
	   if( pfseed == NULL ) {
           seedv[0] = 3579 ; seedv[1] = 27011; seedv[2] = 59243; 
	   }
	   else {
	       seedv2[0] = 3579; seedv2[1] = 27011; seedv2[2] = 59243; 
           for(i=0;i<3;i++){ 
		       if(  fscanf(pfseed," %hd",seedv+i) < 1 )
		            seedv[i] = seedv2[i] ;
		   }
	       fclose( pfseed);
	   }
	   seed48( seedv );   

       // DG: suppress output
       //printf("\n%d %d %d\n", seedv[0], seedv[1], seedv[2] );
	}
	else {
	     pfseed = fopen("seedms","w");
         pseed = seed48(seedv);
         fprintf(pfseed,"%d %d %d\n",pseed[0], pseed[1],pseed[2]);     
	}
}

	void
commandlineseed( char **seeds)
{
	unsigned short seedv[3], *seed48();
	int i ;

	seedv[0] = atoi( seeds[0] );
	seedv[1] = atoi( seeds[1] );
	seedv[2] = atoi( seeds[2] );
    // DG: suppress output
	//printf("\n%d %d %d\n", seedv[0], seedv[1], seedv[2] );

	seed48(seedv);
}
