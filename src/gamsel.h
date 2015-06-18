#ifndef GAMSEL_H
#define GAMSEL_H

/** Print contents of vector */
void printVector(double *x, int lengthx);

/** Calculate z = x - y */
void vectorDifference(int *n, double *x, double *y, double *z);

/** Divergence fix */
double calculateThetaMin(double lam1, double lam2, double alpha, double beta, double A);

/** Calculate lambda_max */
double calculateLambdaMax(int *n, int *p, double *X, double *U, double *y, 
                          double *D, int *degrees, int *cum_degrees, int *numcolsU, 
                          int *family, double gamma);

/** Compute objective function at given values of the parameters. */
double calculateObjective(int *n, int *p, double *X, double *U, double *y, 
                          double *D, int *degrees, int *cum_degrees, int *numcolsU, 
                          double *lambdas_alpha, double *lambdas_beta, double *psis,
                          double *alpha0, double *alphas, double *betas,
                          int *family, double *fit, int *active_alpha, int *active_beta);
                          
/** Intercept update step */
void updateIntercept(double *alpha0, int *n, double *y, double *fit, int *family);

/** Update alpha value */
int updateAlpha(int j, int *n, double *y, double *x, double *fit, double *lambda, double *alphas,
                double *z, double *w, int *family);

/** Update beta value */
int updateBeta(int j, int *n, double *y, double *U, double *fit, double *lambdas, double *lambdas_alpha, double *psis, double *D, int *degrees, int *cum_degrees, double *betas, double *alphas, double *z, int *family);

/** Carry out line search  to determine ||theta|| */
double lineSearch(double m, double *D, double *Vr, double lambda);

#endif /** GAMSEL_H */
