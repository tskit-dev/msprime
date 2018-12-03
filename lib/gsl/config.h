/* NO-OP */

#define HAVE_IEEE_COMPARISONS 1
#define RETURN_IF_NULL(x) if (!x) { return ; }
#define HAVE_EXTENDED_PRECISION_REGISTERS 1

#if HAVE_EXTENDED_PRECISION_REGISTERS
#define GSL_COERCE_DBL(x) (gsl_coerce_double(x))
#else
#define GSL_COERCE_DBL(x) (x)
#endif
