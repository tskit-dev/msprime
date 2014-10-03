struct devent {
	double time;
	int popi;
	int popj;
	double paramv;
	double **mat ;
	char detype ;
	struct devent *nextde;
	} ;
struct c_params {
	int npop;
	int nsam;
	int *config;
	double **mig_mat;
	double r;
	int nsites;
	double f;
	double track_len;
	double *size;
	double *alphag;
	struct devent *deventlist ;
	} ;
struct m_params {
	 double theta;
	int segsitesin;
	int treeflag;
	int timeflag;
	int mfreq;
	 } ;
struct params { 
	struct c_params cp;
	struct m_params mp;
	int commandlineseedflag ;
	int output_precision;
	};
	

