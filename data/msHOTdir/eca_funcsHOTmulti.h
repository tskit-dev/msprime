#ifdef UN_EXTERN
	#define GLOB
#else
	#define GLOB extern
#endif


/* function prototypes */

int re_v(int nsam);  /* function that performs recombination, taking into account 
						recombinational hotspots  */
int cinr_V(int nsam, int nsite);

int cleftr_V(int nsam);

int crightr_V(int nsam);

/*void fprint_coal_times(struct segl *seglist, int nsam, int nsegs, int nsites);
  void pwise_coal_times(struct node *ptree, double **A, int nsam);*/

/* Global variables */
                  /* recombination:*/
GLOB int *gHotSpotStart;
GLOB int *gHotSpotEnd;

GLOB int gNumHotSpots;  /* the number of recombination sites that have relative rates different
							than one */
GLOB double *gHSRates;
/* the number of TIMES the background rate the crossover hotspot should be */

                  /* gene conversion:*/
GLOB int *gHotSpotStartGC;
GLOB int *gHotSpotEndGC;

//GLOB int gHotSpotStartGC;
//GLOB int gHotSpotEndGC;

GLOB int gNumHotSpotsGC;  /* the number of recombination sites that have relative rates different
							than one */
//GLOB double gHSRatesGC;
GLOB double *gHSRatesGC;
/* the number of TIMES the background rate the GC hotspot should be */

