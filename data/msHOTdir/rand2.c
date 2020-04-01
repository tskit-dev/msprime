/*  Link in this file for random number generation with rand() */

#include <stdio.h>
#include <stdlib.h>

	double
ran1()
{
	int rand();
	return( rand()/(RAND_MAX+1.0)  );
}


	void seedit( const char *flag)
{
	FILE *fopen(), *pfseed;
	unsigned int seed2 ;

  if( flag[0] == 's' ) {
    pfseed = fopen("seedms","r");
        if( pfseed == NULL ) {
           seed2 = 59243; }
        else {
          fscanf(pfseed," %d",&seed2);
          fclose( pfseed);
          }
          srand( seed2) ;

        printf("\n%d\n", seed2 );    
	}
   else {
	pfseed = fopen("seedms","w");
        fprintf(pfseed,"%d \n",rand());  

	}
}
