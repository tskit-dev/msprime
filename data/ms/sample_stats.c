#include <stdio.h>
#include <stdlib.h>
#include <string.h>

	double nucdiv(int, int, char **);
	double tajd(int, int, double) ;
	double hfay(int, int, char **);
	double thetah(int, int, char **);

int maxsites = 1000 ;

main(argc,argv)
	int argc;
	char *argv[];
{
	int nsam, j ,nsites, i,  howmany  ;
	char **list, **cmatrix(), allele,na, line[1001], slashline[1001]  ;
	FILE *pf, *fopen(), *pfin ;
	double *posit   ;
	int   segsites, count  , nadv, probflag  ;
	double pi , h, th  ,prob ;
	char dum[20], astr[100] ;
	int  nsegsub, segsub( int nsam, int segsites, char **list ) ;

/* read in first two lines of output  (parameters and seed) */
  pfin = stdin ;
  fgets( line, 1000, pfin);
  sscanf(line," %s  %d %d", dum,  &nsam, &howmany);
  fgets( line, 1000, pfin);

	if( argc > 1 ) { 
	   nadv = atoi( argv[1] ) ; 
	}

  list = cmatrix(nsam,maxsites+1);
  posit = (double *)malloc( maxsites*sizeof( double ) ) ;

  count=0;
	probflag = 0 ;
while( howmany-count++ ) {

/* read in a sample */
  do {
     if( fgets( line, 1000, pfin) == NULL ){
	   exit(0);
	 }
	 if( line[0] == '/' )  strcpy(slashline,line+2);
  }while ( (line[0] != 's') && (line[0] != 'p' ) ) ;
 
  if( line[0] == 'p'){
      sscanf( line, "  prob: %lf", &prob );
	  probflag = 1 ;
	  if( fgets( line, 1000, pfin) == NULL ){
	    exit(0);
	  }
  }
  sscanf( line, "  segsites: %d", &segsites );
  if( segsites >= maxsites){
	maxsites = segsites + 10 ;
	posit = (double *)realloc( posit, maxsites*sizeof( double) ) ;
        biggerlist(nsam,maxsites, list) ;
        }
   if( segsites > 0) {
	fscanf(pfin," %s", astr);
	for( i=0; i<segsites ; i++) fscanf(pfin," %lf",posit+i) ;
	for( i=0; i<nsam;i++) fscanf(pfin," %s", list[i] );
	}
/* analyse sample ( do stuff with segsites and list) */
	if( argc > 1 ) nsegsub = segsub( nadv, segsites, list) ;
	pi = nucdiv(nsam, segsites, list) ;
	h = hfay(nsam, segsites, list) ;
	th = thetah(nsam, segsites, list) ;
	if( argc > 1 )
	printf("pi: %lf ss: %d  D: %lf H: %lf thetah: %lf segsub: %d \n", pi,segsites, tajd(nsam,segsites,pi) , h , th, nsegsub ) ;
	else if( probflag == 1 ) 
	  printf("pi:\t%lf\tss:\t%d\tD:\t%lf\tthetaH:\t%lf\tH:\t%lf\tprob:\t%g%s",
	          pi,segsites, tajd(nsam,segsites,pi) , th , h, prob , slashline ) ;
	else 
	  printf("pi:\t%lf\tss:\t%d\tD:\t%lf\tthetaH:\t%lf\tH:\t%lf%s", pi,segsites, tajd(nsam,segsites,pi) , th , h,slashline  ) ;
	

  }
}

	

/* allocates space for gametes (character strings) */
	char **
cmatrix(nsam,len)
	int nsam, len;
{
	int i;
	char **m;

	if( ! ( m = (char **) malloc( (unsigned)( nsam*sizeof( char* )) ) ) )
	   perror("alloc error in cmatrix") ;
	for( i=0; i<nsam; i++) {
		if( ! ( m[i] = (char *) malloc( (unsigned) (len*sizeof( char )) )))
			perror("alloc error in cmatric. 2");
		}
	return( m );
}

        int
biggerlist(nsam, nmax, list )
        int nsam ;
        unsigned nmax ;
        char ** list ;
{
        int i;

        maxsites = nmax  ;
        for( i=0; i<nsam; i++){
           list[i] = (char *)realloc( list[i],maxsites*sizeof(char) ) ;
           if( list[i] == NULL ) perror( "realloc error. bigger");
           }
}                        


	double
nucdiv( int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;

	nd = nsam;
	nnm1 = nd/(nd-1.0) ;
   	for( s = 0; s <segsites; s++){
		p1 = frequency('1', s,nsam,list)/nd ;
		pi += 2.0*p1*(1.0 -p1)*nnm1 ;
		}
	return( pi ) ;
}

/*   thetah - pi   */
	double
hfay( int nsam, int segsites, char **list)
{
	int s, frequency( char, int, int, char**);
	double pi, p1, nd, nnm1  ;

	pi = 0.0 ;

	nd = nsam;
	nnm1 = nd/(nd-1.0) ;
   	for( s = 0; s <segsites; s++){
		p1 = frequency('1', s,nsam,list)/nd ;
		pi += 2.0*p1*(2.*p1 - 1.0 )*nnm1 ;
		}
	return( -pi ) ;
}

/* Fay's theta_H  */
        double
thetah( int nsam, int segsites, char **list)
{
        int s, frequency( char, int, int, char**);
        double pi, p1, nd, nnm1  ;

        pi = 0.0 ;

        nd = nsam;
        nnm1 = nd/(nd-1.0) ;
        for( s = 0; s <segsites; s++){
                p1 = frequency('1', s,nsam,list) ;
                pi += p1*p1 ; 
                }
        return( pi*2.0/( nd*(nd-1.0) )  ) ;
}


        int
frequency( char allele,int site,int nsam,  char **list)
{
        int i, count=0;
        for( i=0; i<nsam; i++) count += ( list[i][site] == allele ? 1: 0 ) ;
        return( count);
}        

	int
segsub( int nsub, int segsites, char **list )
{
	int i, count = 0 , c1 ;
	int frequency( char, int, int, char**) ;

	for(i=0; i < segsites ; i++){
	  c1 = frequency('1',i,nsub, list);
	  if( ( c1 > 0 ) && ( c1 <nsub )  ) count++;
	  }
	return( count ) ;
}
	
