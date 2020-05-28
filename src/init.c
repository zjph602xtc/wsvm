
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


void
svmtrain (double *x, int *r, int *c,
          double *y,
          double *W,
          int    *rowindex, int *colindex,
          int    *svm_type,
          int    *kernel_type,
          int    *degree,
          double *gamma,
          double *coef0,
          double *cost,
          double *nu,
          int    *weightlabels,
          double *weights,
          int    *nweights,
          double *cache,
          double *tolerance,
          double *epsilon,
          int    *shrinking,
          int    *cross,
          int    *sparse,
          int    *probability,

          int    *nclasses,
          int    *nr,
          int    *index,
          int    *labels,
          int    *nSV,
          double *rho,
          double *coefs,
          double *sigma,
          double *probA,
          double *probB,

          double *cresults,
          double *ctotal1,
          double *ctotal2,
          char   **error);

void
svmpredict  (int    *decisionvalues,
    	     int    *probability,

             double *v, int *r, int *c,
	     int    *rowindex,
	     int    *colindex,
	     double *coefs,
	     double *rho,
	     int    *compprob,
	     double *probA,
	     double *probB,
	     int    *nclasses,
	     int    *totnSV,
	     int    *labels,
	     int    *nSV,
	     int    *sparsemodel,

	     int    *svm_type,
	     int    *kernel_type,
	     int    *degree,
	     double *gamma,
	     double *coef0,

	     double *x, int *xr,
	     int    *xrowindex,
	     int    *xcolindex,
	     int    *sparsex,

	     double *ret,
	     double *dec,
	     double *prob);

void
svmwrite (double *v, int *r, int *c,
	  int    *rowindex,
	  int    *colindex,
	  double *coefs,
	  double *rho,
	  int    *compprob,
          double *probA,
          double *probB,
	  int    *nclasses,
	  int    *totnSV,
	  int    *labels,
	  int    *nSV,
	  int    *sparsemodel,

	  int    *svm_type,
	  int    *kernel_type,
	  int    *degree,
	  double *gamma,
	  double *coef0,

	  char **filename);


static const R_CMethodDef CEntries[] = {
    {"wsvmpredict", (DL_FUNC) &svmpredict, 30},
    {"wsvmtrain", (DL_FUNC) &svmtrain, 38},
    {"wsvmwrite", (DL_FUNC) &svmwrite, 21},
    {NULL, NULL, 0}
};

void R_init_WeightSVM(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
