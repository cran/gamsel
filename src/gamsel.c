#include <assert.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/BLAS.h>
#include <stdlib.h>
#include "gamsel.h"
#include <R_ext/Rdynload.h>

static const double one=1.0;
static const double zero=0.0;
static const int inc_one = 1;
static const double EPS = 1e-4;
// static const double SCALE_FACTOR = 0.9;
static const double MAX_DEV_EXPLAINED = 0.99;

#ifndef max
	#define max( a, b ) ( ((a) >= (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) <= (b)) ? (a) : (b) )
#endif

void printVector(double *x, int lengthx) {
  Rprintf("[");
  for(int i=0; i<lengthx; i++) {
    Rprintf("%f, ", x[i]);
  }
  Rprintf("]\n");
}

void printIntVector(int *x, int lengthx) {
  Rprintf("[");
  for(int i=0; i<lengthx; i++) {
    Rprintf("%d, ", x[i]);
  }
  Rprintf("]\n");
}

/* Computes difference between vectors x and y, both of length n
   and stores it in vector z */
void vectorDifference(int *n, double *x, double *y, double *z) {
  for(int i=0; i<*n; i++) {
    z[i] = x[i] - y[i];
  }
}

void vectorSum(int *n, double *x, double *y, double *z) {
  for(int i=0; i<*n; i++) {
    z[i] = x[i] + y[i];
  }
}

double calculateThetaMin(double lam1, double lam2, double alpha, double beta, double A) {
  double candidates[3];
  candidates[0] = -alpha;
  double pmterm = sqrt(fabs((pow(lam1,2)*A)/(pow(lam1,2) - pow(lam2, 2))));
  candidates[1] = beta - pmterm;
  candidates[2] = beta + pmterm;
  double cur_min_value = lam2 * sqrt(pow(beta - candidates[0],2) + A);
  int cur_min_index = 0;
  double fval;
  for(int i=1; i<3; i++) {
    fval = lam1*fabs(alpha + candidates[i]) + lam2 * sqrt(pow(beta - candidates[i],2) + A);
    if(fval < cur_min_value) {
      cur_min_value = fval;
      cur_min_index=i;
    }
  }
  return(candidates[cur_min_index]);
}

double implicitFunction(double c, double m, double *D, double *Vr, double lambda) {
  double out = -1;
  for(int i=0; i<m; i++) {
    out += pow(Vr[i]/(D[i]*c + lambda),2);
  }
  return out;
}

double lineSearch(double m, double *D, double *Vr, double lambda) {
  double tol = 0.0001;
  double start = 1;  
  double end = 2;
  //double scale = 2;
  double start_val = implicitFunction(start, m, D, Vr, lambda);
  double end_val = implicitFunction(end, m, D, Vr, lambda);
  while(sign(start_val) == sign(end_val)) {
    if(start_val > end_val && sign(start_val) > 0) {
      end*=2.0;
      end_val = implicitFunction(end, m, D, Vr, lambda);
    } else {
      start*=0.5;
      start_val = implicitFunction(start, m, D, Vr, lambda);
    }
  }
  double midpt = 0.5*(start + end);
  double midpt_val = implicitFunction(midpt, m, D, Vr, lambda);
  while(end - start > tol) {
    if(sign(midpt_val) == sign(start_val)) {
      start = midpt;
      start_val = midpt_val;
    } else {
      end = midpt;
      end_val = midpt_val;
    }
    midpt = 0.5*(start + end);
    midpt_val = implicitFunction(midpt, m, D, Vr, lambda);
  }
  return midpt; 
}

double calculateLambdaMax(int *n, int *p, double *X, double *U, double *y, 
                          double *D, int *degrees, int *cum_degrees, int *numcolsU, 
                          int *family, double gamma) {
  double curr_max = 0.0;
  double norm = 0.0;
  double trDinv;
  for(int j=0;j<*p;j++){
    trDinv = 0.0;
    double *Ujy = malloc(degrees[j]*sizeof(double));
    // Calculate alpha norm
    norm = fabs(F77_CALL(ddot)(n, X+(*n)*j, &inc_one, y, &inc_one))/gamma;
    curr_max = max(curr_max, norm);
    // Calculate beta norm
    F77_CALL(dgemv)("T",n,degrees+j,&one,U+(*n)*(cum_degrees[j]),n,y,
      &inc_one, &zero, Ujy, &inc_one);
    for(int i=0; i<degrees[j];i++) {
      trDinv += 1/D[cum_degrees[j] + i];
    }
    // Calculate norm of D^{-1/2}Ujy and scale
    free(Ujy);
  }
  return curr_max;
}


double calculateDeviance(int n, double *fit, double *y) {
  double *probs = malloc(n*sizeof(double));
  double loglik = 0.0;
  for(int i=0; i<n; i++) {
    probs[i] = 1/(1+exp(-fit[i]));
    loglik += log(1 - probs[i]) + y[i]*fit[i];
  }
  free(probs);
  return(-2*loglik);
}
/* Computes objective function at given values of the parameters. */
double calculateObjective(int *n, int *p, double *X, double *U, double *y, 
                          double *D, int *degrees, int *cum_degrees, int *numcolsU, 
                          double *lambdas_alpha, double *lambdas_beta, double *psis,
                          double *alpha0, double *alphas, double *betas,
                          int *family, double *fit,
                          int *active_alpha, int *active_beta) {
  double obj = 0.0;
  double lasso_penalty = 0.0;
  double group_penalty = 0.0;
  double spline_penalty = 0.0; 
  double group_norm_squared = 0.0;
  int group_start_index = 0;
  double *resid = Calloc(*n, double);
  
  memset(fit, 0.0, (*n)*sizeof(double));
  /** Calculate linear fit:    */
  for(int j=0; j<*p; j++) {
    if(active_alpha[j] == 1) {
      for(int i=0; i<*n; i++) {
        fit[i] += alphas[j]*X[(*n)*j + i];
      }
    }
    if(active_beta[j] == 1) {
      F77_CALL(dgemv)("N", n, degrees+j, &one, U+(*n)*(cum_degrees[j]), n, betas + cum_degrees[j], &inc_one, &one, fit, &inc_one);
    }
    // Increment fit by X*alpha
    // F77_CALL(dgemv)("N", n, p, &one, X, n, alphas, &inc_one, &one, fit, &inc_one);
    // Increment by U*beta
    // F77_CALL(dgemv)("N", n, numcolsU, &one, U, n, betas, &inc_one, &one, fit, &inc_one);
  }
  // Increment by intercept
  for(int i=0; i<*n; i++) {
    fit[i] += alpha0[0];
  }
  /** End linear fit calculation */
  
  if(*family == 0) {  // Linear model
    // Calculate residual
    vectorDifference(n, y, fit, resid);
    // Increment objective by RSS
    obj += F77_CALL(ddot)(n, resid, &inc_one, resid, &inc_one);
  } else if(*family == 1){  // Logistic model
    // Add negative log-likelihood to objective
    for(int i=0; i<*n; i++) {
      obj += -(y[i]*fit[i] - log(1 + exp(fit[i])));
    }
  }

  // Calculate lasso penalty
  for(int i=0; i<*p; i++) {
    if(active_alpha[i] == 1) {
      lasso_penalty += lambdas_alpha[i]*fabs(alphas[i]);
    }
  }
    
  double *Dbeta = Calloc(*numcolsU, double);
  for(int i=0; i<*numcolsU; i++) {
    Dbeta[i] = betas[i]*D[i];
  }
  
  // Calculate group penalty and spline penalty
  for(int i=0; i<*p; i++) {
    if(active_beta[i] == 1) {
      // group penalty
      group_norm_squared = F77_CALL(ddot)(degrees+i, Dbeta + group_start_index, &inc_one, betas + group_start_index, &inc_one);
      group_penalty += lambdas_beta[i]*sqrt(group_norm_squared);
      // spline penalty
      spline_penalty += psis[i]*(group_norm_squared- Dbeta[group_start_index]*betas[group_start_index]);
    }
    group_start_index += degrees[i];
  }

  // Calculate objective
  if(*family == 0) obj*=0.5;
  obj += lasso_penalty + group_penalty + 0.5*spline_penalty;
  // obj = obj/(*n);  // Apply 1/n factor
  Free(resid);
  Free(Dbeta);
  return obj;
}

// Intercept update
void updateIntercept(double *alpha0, int *n, double *y, double *fit, int *family) {
  if(*family == 0) {
    double diff = 0.0;
    for(int i=0; i<*n; i++) {
      diff += (y[i] - fit[i] + alpha0[0]);
    }
    alpha0[0] = diff/(*n);
  } else if (*family == 1) {
    // 1-d Newton-Raphson algorithm for logistic intercept
    double *fit_no_intercept = Calloc(*n, double);
    double tol = 1e-5;
    double old_alpha0 = alpha0[0];
    double new_alpha0 = alpha0[0] + 1;
    double deriv;
    double second_deriv;
    double exp_term;
    for(int i=0; i<*n; i++) {
      fit_no_intercept[i] = fit[i] - alpha0[0];
    }
    // Do Newton steps until convergence
    while(fabs(new_alpha0 - old_alpha0) > tol) {
      old_alpha0 = new_alpha0;
      deriv = 0.0; second_deriv=0.0;
      // Calculate first and second derivatives
      for(int i=0; i<*n; i++) {
        exp_term = exp(-(old_alpha0 + fit_no_intercept[i]));
        deriv += y[i] - 1/(1 + exp_term);
        second_deriv += -exp_term/pow(1+exp_term,2);
      }
      // Newton-Raphson step
      new_alpha0 = old_alpha0 - deriv/second_deriv;
    }
    // hard code:
    // alpha0[0] = 0.05;
    alpha0[0] = new_alpha0;
    Free(fit_no_intercept);
  } 
}

// Alpha vector update
// For logisitc regression, all weights taken to be 0.25
int updateAlpha(int j, int *n, double *y, double *x, double *fit, double *lambdas, double *alphas,           
                double *z, double *w, int *family) {
  int nonzero = 1;
  double dotprod = 0.0;
  double *resid = Calloc(*n, double);
  double alpha_old = alphas[j];
  if(*family == 0) {  // Linear regression
    vectorDifference(n, y, fit, resid);  // Calculate full residual
    // Adjust residual for contribution of alpha_j, and calculate inner product
    for(int i=0; i<*n; i++) {
      resid[i] += x[i]*alphas[j];
      dotprod += x[i]*resid[i];
    }
    // Apply soft-thresholding
    if(fabs(dotprod) < lambdas[j]) {
      alphas[j] = 0;
      nonzero = 0;
    } else if(dotprod > lambdas[j]) {
      alphas[j] = dotprod - lambdas[j]; 
    } else {
      alphas[j] = dotprod + lambdas[j];
    }
  } else if(*family == 1) {  // Logistic regression
    vectorDifference(n, z, fit, resid);  // Calculate full residual
    double w_sum = 0.0;
    // Adjust fit for contribution of alpha_j, and calculate inner product
    for(int i=0; i<*n; i++) {
      resid[i] += x[i]*alphas[j];
      dotprod += w[i]*x[i]*resid[i];
      w_sum += w[i]*pow(x[i],2);
    }
    if(fabs(dotprod) < lambdas[j]) {
      alphas[j] = 0; 
      nonzero=0;
    } else if(dotprod > lambdas[j]) {
      alphas[j] = (dotprod - lambdas[j])/w_sum;
    } else {
      alphas[j] = (dotprod + lambdas[j])/w_sum;
    }
  }
  // Update fitted value
  if(alphas[j] != alpha_old) {
    for(int i=0; i<*n; i++) {
      fit[i] += x[i]*(alphas[j] - alpha_old);
    }
  }
  Free(resid);
  return nonzero;
}



// Beta vector update
int updateBeta(int j, int *n, double *y, double *U, double *fit, double *lambdas, double *lambdas_alpha, double *psis, double *D, int *degrees, int *cum_degrees, double *betas, double *alphas, double *z, int *family) {
  double *resid = Calloc(*n, double);
  double *UjBj = Calloc(*n, double);
  double *DinvUjr = Calloc(degrees[j], double);
  double *old_beta = Calloc(degrees[j], double);
  int nonzero = 1;
  
  // Store old value of beta
  for(int i=0; i<degrees[j]; i++) {
    old_beta[i] = betas[cum_degrees[j] + i];
  }
  if(*family == 0) {  // Linear regression updates
    // Calculate residual
    vectorDifference(n,y,fit,resid);
    // Adjust resdiual for contribution of beta_j
    F77_CALL(dgemv)("N", n, degrees+j, &one, U+(*n)*(cum_degrees[j]), n, betas+cum_degrees[j], &inc_one, &zero, UjBj, &inc_one);
    for(int i=0; i<*n; i++) {
      resid[i] += UjBj[i];
    }
    // Compute D^{-1/2} * U_j^T * resid
    F77_CALL(dgemv)("T", n, degrees+j, &one, U+(*n)*(cum_degrees[j]), n, resid, 
      &inc_one, &zero, DinvUjr, &inc_one);
    for(int i=0; i<degrees[j]; i++) {
      DinvUjr[i] = DinvUjr[i] / sqrt(D[cum_degrees[j]+i]);
    }
    double norm_DinvUjr = sqrt(F77_CALL(ddot)(degrees+j, DinvUjr, &inc_one, DinvUjr, &inc_one));
    // Check if beta_j is 0
    if(norm_DinvUjr < lambdas[j]) {
      for(int i=0; i<degrees[j]; i++) {
        betas[cum_degrees[j] + i] = 0.0;
      }
      nonzero = 0;
    } else {
      double *Dstar = Calloc(degrees[j], double);
      // Form D*
      for(int i=0; i<degrees[j]; i++) {
        Dstar[i] = psis[j] + 1/D[cum_degrees[j]+i];
      }
      Dstar[0] -= psis[j];  // Adjust first entry
      // Line search to find c
      double c = lineSearch(degrees[j], Dstar, DinvUjr, lambdas[j]);
      // Update beta using equation (4.3)
      for(int i=0; i<degrees[j]; i++) {
        betas[cum_degrees[j]+i] = (DinvUjr[i]/sqrt(D[cum_degrees[j]+i]))/(Dstar[i] + lambdas[j]/c);
      }
      Free(Dstar);
    }
  } else if(*family == 1){
    // Calculate residual
    vectorDifference(n,z,fit,resid);
    // Adjust resdiual for contribution of beta_j
    F77_CALL(dgemv)("N", n, degrees+j, &one, U+(*n)*(cum_degrees[j]), n, betas+cum_degrees[j], &inc_one, &zero, UjBj, &inc_one);
    for(int i=0; i<*n; i++) {
      resid[i] += UjBj[i];
    }
    // Compute D^{-1/2} * U_j^T * resid
    F77_CALL(dgemv)("T", n, degrees+j, &one, U+(*n)*(cum_degrees[j]), n, resid, 
      &inc_one, &zero, DinvUjr, &inc_one);
    for(int i=0; i<degrees[j]; i++) {
      DinvUjr[i] = DinvUjr[i] / sqrt(D[cum_degrees[j]+i]);
    }
    double norm_DinvUjr = sqrt(F77_CALL(ddot)(degrees+j, DinvUjr, &inc_one, DinvUjr, &inc_one));
    // Check if beta_j is 0
    if(norm_DinvUjr < 4*lambdas[j]) {
      for(int i=0; i<degrees[j]; i++) {
        betas[cum_degrees[j] + i] = 0.0;
      }
      nonzero = 0;
    } else {
      double *Dstar = Calloc(degrees[j], double);
      // Form D*
      for(int i=0; i<degrees[j]; i++) {
        Dstar[i] = 4*psis[j] + 1/D[cum_degrees[j]+i];
      }
      Dstar[0] -= 4*psis[j];  // Adjust first entry
      // Line search to find c
      double c = lineSearch(degrees[j], Dstar, DinvUjr, 4*lambdas[j]);
      // Update beta using equation (4.3) and do
      // for(int i=0; i<degrees[j]; i++) {
      //   betas[cum_degrees[j]+i] = (DinvUjr[i]/sqrt(D[cum_degrees[j]+i]))/(Dstar[i] + 4*lambdas[j]/c);
      // }
  
      double *new_beta = malloc(degrees[j]*sizeof(double));
      for(int i=0; i < degrees[j]; i++) {
        new_beta[i] = (DinvUjr[i]/sqrt(D[cum_degrees[j]+i]))/(Dstar[i] + 4*lambdas[j]/c);
      }
      // Adjust beta and alpha to match
      // if(sign(alphas[j])*sign(new_beta[0]) < -0.1) {
      //   double A = 0.0;
      //   for(int i=1; i<degrees[j]; i++) {
      //     A += pow(new_beta[i], 2);
      //   }
      //   double mintheta = calculateThetaMin(lambdas_alpha[j], lambdas[j], alphas[j], 
      //     new_beta[0], A);
      //   alphas[j] += mintheta;
      //   new_beta[0] -= mintheta;
      // }
      
      // Backtracking step (backtrack until sign is not opposite of alpha[j])
      // int loop_counter = 0;
      // while(loop_counter < 20 && sign(alphas[j])*sign(new_beta[0]) < -0.1) {
      //   for(int i=0; i<degrees[j]; i++) {
      //     new_beta[i] = old_beta[i] + SCALE_FACTOR*(new_beta[i] - old_beta[i]);
      //   }
      //   loop_counter++;
      // }
      // if(sign(alphas[j])*sign(new_beta[0]) < -0.1 && fabs(new_beta[0]) > 1) {
      //   new_beta[0] += alphas[j];
      // }
      
      for(int i=0; i<degrees[j]; i++) {
        betas[cum_degrees[j]+i] = new_beta[i];
      }
      Free(Dstar);
      free(new_beta);
    }
  } else if(*family == 2) {  // Logistic regression updates
    double *gradient = Calloc(degrees[j], double);
    double *DinvTheta = Calloc(degrees[j], double);
    double *prob = Calloc(*n, double);
    double *grad_plus_dinv = Calloc(degrees[j], double);
    double norm_grad_plus_dinv;
    for(int i=0; i<degrees[j]; i++) {
      prob[i] = 1/(1 + exp(-fit[i]));
      if(prob[i] < EPS) {
        prob[i] = EPS;
      } else if(prob[i] > 1-EPS) {
        prob[i] = 1;
      }
    }
    vectorDifference(n, y, prob, resid);  // y - yhat
    // F77_CALL(dgemv)("N", n, degrees+j, &one, U+(*n)*(cum_degrees[j]), n, betas+cum_degrees[j], &inc_one, &zero, UjBj, &inc_one);
    // for(int i=0; i<*n; i++) {
    //   resid[i] += UjBj[i];
    //   // resid[i] *= 0.25;
    // }
    // Compute gradient = D^{-1/2} * U_j^T * resid
    F77_CALL(dgemv)("T", n, degrees+j, &one, U+(*n)*(cum_degrees[j]), n, resid, 
      &inc_one, &zero, gradient, &inc_one);
    for(int i=0; i<degrees[j]; i++) {
      gradient[i] = gradient[i] / sqrt(D[cum_degrees[j]+i]);
      // Compute D^{-1} theta = D^{-1/2} beta
      DinvTheta[i] = betas[cum_degrees[j]+i] / sqrt(D[cum_degrees[j]+i]);
      // DinvTheta[i] = betas[cum_degrees[j]+i] * sqrt(D[cum_degrees[j]+i]);
      grad_plus_dinv[i] = gradient[i] + DinvTheta[i]/4;
    }
    norm_grad_plus_dinv = sqrt(F77_CALL(ddot)(degrees+j, grad_plus_dinv, &inc_one, grad_plus_dinv, &inc_one));

    if(norm_grad_plus_dinv <= lambdas[j]) {
      for(int i=0; i<degrees[j]; i++) {
        betas[cum_degrees[j] + i] = 0.0;
      } 
    } else {
      double *Dstar = Calloc(degrees[j], double);
      // Form D*
      for(int i=0; i<degrees[j]; i++) {
        Dstar[i] = psis[j] + 0.25/D[cum_degrees[j]+i];
      }
      Dstar[0] -= psis[j];  // Adjust first entry
      // Line search to find c
      double c = lineSearch(degrees[j], Dstar, grad_plus_dinv, lambdas[j]);
      // Update beta using equation (4.3)
      for(int i=0; i<degrees[j]; i++) {
        betas[cum_degrees[j]+i] = (grad_plus_dinv[i]/sqrt(D[cum_degrees[j]+i]))/(Dstar[i] + lambdas[j]/c);
      }
      Free(Dstar);
    }
    Free(gradient);
    Free(DinvTheta);
    Free(prob);
  }
  /** Fitted value update step */
  double *beta_diff = Calloc(degrees[j], double);
  int equal_Flag = 0;
  for(int i=0; i<degrees[j]; i++) {
    beta_diff[i] = betas[cum_degrees[j]+i] - old_beta[i];
    if(old_beta[i] != betas[cum_degrees[j]+i]) {
      equal_Flag = 1;
    }
  }
  if(equal_Flag != 0) {
    double *Ujbeta_diff = Calloc(*n, double);
    F77_CALL(dgemv)("N", n, degrees+j, &one, U+(*n)*(cum_degrees[j]), n, beta_diff, 
      &inc_one, &zero, Ujbeta_diff, &inc_one);
    for(int i=0; i<*n; i++) {
      fit[i] += Ujbeta_diff[i];
    }
    Free(Ujbeta_diff);
  }
  
  Free(resid);
  Free(UjBj);
  Free(old_beta);
  Free(DinvUjr);
  Free(beta_diff);
  return nonzero;
}

SEXP gamselFit(SEXP R_y, SEXP R_X, SEXP R_U, SEXP R_tol,
               SEXP R_degrees, SEXP R_D, SEXP R_gamma, SEXP R_psis,
               SEXP R_family, SEXP R_MAXITER, SEXP R_lambda_seq, SEXP R_num_lambda,
               SEXP R_traceit) {
  PROTECT(R_y = coerceVector(R_y, REALSXP));
  PROTECT(R_X = coerceVector(R_X, REALSXP));
  PROTECT(R_U = coerceVector(R_U, REALSXP));
  PROTECT(R_degrees = coerceVector(R_degrees, INTSXP));
  PROTECT(R_D = coerceVector(R_D, REALSXP));
  
  double tol = *REAL(R_tol);
  
  SEXP dims =  getAttrib(R_X, R_DimSymbol);
  SEXP dimsU = getAttrib(R_U, R_DimSymbol);
  int n = INTEGER(dims)[0];             // Number of observations
  int p = INTEGER(dims)[1];             // Number of parameters
  int numcolsU = INTEGER(dimsU)[1];     // Number of columns in U
  double *X = REAL(R_X);                // Pointer to design matrix
  double *U = REAL(R_U);                // Pointer to basis functions matrix
  double *y = REAL(R_y);                // Pointer to response y
  int *degrees = INTEGER(R_degrees);    // Pointer to degrees
  double *D = REAL(R_D);                // Pointer to diagonal entries of D
  double *psis_real = REAL(R_psis);          // Pointer to psi values
  int family = *INTEGER(R_family);      // 0 = Gaussian; 1 = Logistic
  double gamma = *REAL(R_gamma);        // Value of gamma
  // double *lambda_seq = REAL(R_lambda_seq);  // User-specified lambda values
  int num_lambda = *INTEGER(R_num_lambda);            // Number of lambdas
  int traceit = *INTEGER(R_traceit);                  // Print flag
  int *cum_degrees = malloc(p*sizeof(int));  // Cumulative degrees for beta
  double *psis = malloc(p*sizeof(double));
  double *lambda_seq = malloc(num_lambda*sizeof(double));
  double *residual_deviance = Calloc(num_lambda, double);
  double *frac_deviance_explained = Calloc(num_lambda, double);
  // Set cum_degrees
  int curr_sum = 0;
  for(int i=0; i<p; i++) {
    cum_degrees[i] = curr_sum;
    curr_sum += degrees[i];
  }
  
  // If lambda sequence is not provided, calculate it here
  if(REAL(R_lambda_seq)[0] < 0) {
    // Set starting value of lambda at lambda_max
    double lambda_max = calculateLambdaMax(&n, &p, X, U, y, D, degrees, cum_degrees, &numcolsU, &family, gamma);  
    double lambda_factor = family == 0 ? 0.01 : 0.05;
    // lambda_factor *= (n > p ? 0.1 : 1);
    double lambda_min = lambda_factor*lambda_max;   // Smallest lambda value
    // unused double lambda = lambda_max;
    // Construct lambda sequence
    double t;
    for(int i=0; i<num_lambda; i++) {
      t = ((double)i)/(num_lambda - 1);
      // logarithmic sequence
      lambda_seq[num_lambda-i-1] = exp(log(lambda_min) + t*log(lambda_max/lambda_min));
    }
  } else {
    for(int i=0; i<num_lambda; i++) {
      lambda_seq[i] = REAL(R_lambda_seq)[i];
    }
  }
  // if(traceit == 1) {
  //   printVector(lambda_seq, num_lambda); 
  // }
  
  
  // Penalty parameter values
  double *lambdas_alpha = malloc(p*sizeof(double));    
  double *lambdas_beta = malloc(p*sizeof(double)); 
  // Fitted values
  double *fit = Calloc(n, double);
  // Working response (logistic regression)
  double *z = Calloc(n, double);
  // Probabilities (logistic regression)
  double *prob = Calloc(n, double);
  // Weights (logistic regression)
  double *w = Calloc(n, double);
  
  // Allocate output
  SEXP R_all_alpha0 = PROTECT(allocVector(REALSXP, num_lambda));
  SEXP R_all_alphas = PROTECT(allocMatrix(REALSXP, p, num_lambda));
  SEXP R_all_betas =  PROTECT(allocMatrix(REALSXP, numcolsU, num_lambda));
  
  double *alpha0 = Calloc(1, double);
  double *alphas = Calloc(p, double);
  double *betas = Calloc(numcolsU, double);
  
  double old_obj = 0.0;
  double new_obj = 1.0;
  double pre_update_obj;
  double post_update_obj;
  int changedFlag = 1;     // Flag to determine if active set needs to be updated
  int *active_alpha = Calloc(p, int);
  int *active_beta = Calloc(p, int);
  int is_nonzero = 0;
  for(int i=0; i<p; i++) {
    active_alpha[i] = 0;
    active_beta[i] = 0;
  }
  
  // Calculate null deviance for logistic model
  double null_deviance = 0 /* -Wall */;
  if(family == 1) {
    double *null_fit = malloc(n*sizeof(double));
    double y_sum = 0.0;
    for(int i = 0; i<n; i++) {
      y_sum += y[i];
    }
    double fitted_val = log((y_sum/n)/(1-(y_sum/n)));
    for(int i=0; i<n; i++) {
      null_fit[i] = fitted_val;
    }
    null_deviance = calculateDeviance(n, null_fit, y);
    free(null_fit);
    // Rprintf("Null deviance = %f \n", null_deviance);
  }
  
  int iter=1; // Count number of passes through the variables
  int num_quadratic_updates = 0;
  int num_outer_loops = 0;
  double dotprod;
  double *resid = malloc(n*sizeof(double));
  double lambda;
  for(int idx = 0; idx < num_lambda; idx++) {
    changedFlag = 1;
    // Set lambdas
    lambda = lambda_seq[idx];
    // lambda = lambda_max;
    int curr_index = 0;
    for(int i=0; i<p; i++) {
      double trDinv = 0.0;
      for(int j=0; j<degrees[i]; j++) {
        trDinv += 1/D[curr_index];
        curr_index++;
      }
      if(degrees[i] > 1) {
        lambdas_alpha[i] = lambda*gamma;
      } else {
        lambdas_alpha[i] = lambda*min(gamma, 1-gamma);
      }
      lambdas_beta[i] = (1-gamma)*lambda*sqrt(trDinv);
      // psis[i] = psis_real[i]/(1 + lambda);
      psis[i] = psis_real[i];
    }
    
    /** Determine active set using strong rules*/
    if(family == 0 && idx > 0) {
      // Calculate residual
      vectorDifference(&n, y, fit, resid);
      for(int j=0; j<p; j++) {
        dotprod = F77_CALL(ddot)(&n, X+(n*j), &inc_one, resid, &inc_one);
        if(fabs(dotprod) < 2*lambdas_alpha[j] - lambda_seq[idx-1]*gamma) {
          active_alpha[j] = 0;
        } else {
          active_alpha[j] = 1;
        }
      }
    }
    /** End strong rules calculation */
    
    while(changedFlag > 0) {
      changedFlag = 0;
      num_quadratic_updates = 0.0;
      
      /** Update working response and weights */
      if(family == 1) {
        for(int i=0; i<n; i++) {
          prob[i] = 1/(1+exp(-fit[i]));
          w[i] = prob[i]*(1-prob[i]);
          if(prob[i] <= EPS) {
            prob[i] = 0;
            w[i] = EPS;
          } else if(prob[i] >= 1-EPS) {
            prob[1] = 1;
            w[i] = EPS;
          } else {
            w[i] = prob[i]*(1-prob[i]);
          }
          w[i] = 0.25;
          z[i] = fit[i] + (y[i] - prob[i])/w[i];
        }
      }
      
      updateIntercept(alpha0, &n, y, fit, &family);
      
      /** Determine active set by checking which are non-0*/
      for(int j=0; j<p; j++) {
        if(family == 1) {
          is_nonzero = updateAlpha(j, &n, y, X+(n*j), fit, lambdas_alpha, alphas, z, w, &family);
          if(is_nonzero > active_alpha[j]) {
            active_alpha[j] = 1;
            changedFlag = 1;
          }
        }
        if(degrees[j] > 1) {
          is_nonzero = updateBeta(j, &n, y, U, fit, lambdas_beta, lambdas_alpha, psis, D, degrees, cum_degrees, betas, alphas, z, &family);
          if(is_nonzero > active_beta[j]) {
            active_beta[j] = 1;
            changedFlag = 1;
          }
        } 
      }
      
      /** End of active set determination */
      
      pre_update_obj = calculateObjective(&n, &p, X, U, y, D, degrees, cum_degrees, &numcolsU,
                        lambdas_alpha, lambdas_beta, psis, alpha0, alphas, betas, &family, fit,
                        active_alpha, active_beta);
      post_update_obj = -1;
      while(pre_update_obj - post_update_obj > 1e-3) {
        if(post_update_obj > - 1) {  // After first cycle, update pre_object_obj
          pre_update_obj = post_update_obj;
        }
        /** Begin quadratic approximation update */
        if(family == 1) {
          for(int i=0; i<n; i++) {
            prob[i] = 1/(1+exp(-fit[i]));
            if(prob[i] <= EPS) {
              prob[i] = 0;
              w[i] = EPS;
            } else if(prob[i] >= 1-EPS) {
              prob[1] = 1;
              w[i] = EPS;
            } else {
              w[i] = prob[i]*(1-prob[i]);
              // w[i] = 0.25;
            }
            w[i] = 0.25;
            z[i] = fit[i] + (y[i] - prob[i])/w[i];
          }
          num_quadratic_updates++;
        }
        /** End quadratic approximation update */
  
        /** Loop over active set until convergence */
        old_obj = calculateObjective(&n, &p, X, U, y, D, degrees, cum_degrees, &numcolsU,
                          lambdas_alpha, lambdas_beta, psis, alpha0, alphas, betas, &family, fit,
                          active_alpha, active_beta);
        new_obj = old_obj + 1;
        while(iter<*INTEGER(R_MAXITER)) {
          /** End quadratic update */
          updateIntercept(alpha0, &n, y, fit, &family);
          for(int j=0; j<p; j++){
            // if(sign(alphas[j]*betas[cum_degrees[j]])<0) {
            //   alphas[j] = alphas[j] + betas[cum_degrees[j]];
            //   betas[cum_degrees[j]] = 0.0;
            // }
            if(active_alpha[j] == 1) {
              updateAlpha(j, &n, y, X+(n*j), fit, lambdas_alpha, alphas, z, w, &family);
            }  
            if(degrees[j] > 1 && active_beta[j] == 1) {
              updateBeta(j, &n, y, U, fit, lambdas_beta, lambdas_alpha, psis, D, degrees, cum_degrees, betas, alphas, z,  &family);
            }    
          }
          new_obj = calculateObjective(&n, &p, X, U, y, D, degrees, cum_degrees, &numcolsU,
                          lambdas_alpha, lambdas_beta, psis, alpha0, alphas, betas, &family, fit,
                          active_alpha, active_beta);
          // Rprintf("%f \n", new_obj);
          if(fabs(old_obj - new_obj) < tol) {
            break;
          }
          old_obj = new_obj;
          iter++;
          // Rprintf("Completed iteration %d \n", iter);
        }
        post_update_obj = new_obj;
      } /** End loop over quadratic approximation */
      // if(family == 1 && traceit == 1) {
      //   Rprintf("Number of quadratic updates = %d \n", num_quadratic_updates);
      // }    
      num_outer_loops++;
    } /** End loop over active set */
    // Print final objective 
    if(traceit == 1) {
      // Rprintf("Final objective value: %f \n", new_obj);
    }
    // Print output
    // Rprintf("Final estimates: \n");
    // printVector(alphas, p);
    // printVector(betas, numcolsU);
    /** Prepare output for current lambda value */
    REAL(R_all_alpha0)[idx] = alpha0[0];
    for(int i=0; i<p; i++) {
      REAL(R_all_alphas)[p*idx + i] = alphas[i];
    }
    for(int i=0; i<numcolsU; i++) {
      REAL(R_all_betas)[numcolsU*idx + i] = betas[i];
    }
    // Residual deviance calculation
    if(family == 1) {
      residual_deviance[idx] = calculateDeviance(n, fit, y);
      frac_deviance_explained[idx] = (null_deviance - residual_deviance[idx])/null_deviance;
      if(traceit == 1) {
        Rprintf("Percent deviance explained = %f \n", frac_deviance_explained[idx]);
      }
      // Terminate if deviance explained gets too high
      if(frac_deviance_explained[idx] > MAX_DEV_EXPLAINED) {
        Rprintf("Did not reach end of path.  Deviance explained exceeded %f \n", MAX_DEV_EXPLAINED);
        break;
      }
    }
    /** End output */
  }
  if(traceit == 1) {
    Rprintf("Number of iterations: %d \n", iter+1);
    // if(p <= 100) {
    //   Rprintf("Final active set: \n");
    //   printIntVector(active_alpha, p);
    //   printIntVector(active_beta, p);
    // }
    Rprintf("Number of outer loops: %d \n", num_outer_loops);
  }
  
  
  
  // Rprintf("Final quadratic approximation: \n");
  // printVector(z, n);
  // printVector(w, n);
  // printVector(prob, n);
  
  
  
  
  // Prepare output 
  SEXP R_lambda_vals = PROTECT(allocVector(REALSXP, num_lambda));
  SEXP R_dev_explained = PROTECT(allocVector(REALSXP, num_lambda));
  for(int i=0; i<num_lambda; i++) {
    REAL(R_lambda_vals)[i] = lambda_seq[i];
    REAL(R_dev_explained)[i] = frac_deviance_explained[i];
  }
  SEXP result = PROTECT(allocVector(VECSXP,5));
  SET_VECTOR_ELT(result, 0, R_all_alpha0);
  SET_VECTOR_ELT(result, 1, R_all_alphas);
  SET_VECTOR_ELT(result, 2, R_all_betas);
  SET_VECTOR_ELT(result, 3, R_lambda_vals);
  SET_VECTOR_ELT(result, 4, R_dev_explained);
  const char *result_names[5] = {"intercept", "alphas", "betas", "lambdas", "dev.explained"};
  SEXP sNames = PROTECT(allocVector(STRSXP, 5));
  for(int i=0; i<5; i++) {
    SET_STRING_ELT(sNames, i, mkChar(result_names[i]));
  }
  setAttrib(result, R_NamesSymbol, sNames);
    
  // Unprotect memory and return
  UNPROTECT(12);
  Free(active_alpha);
  Free(active_beta);
  Free(alpha0);
  Free(alphas);
  Free(betas);
  Free(fit);
  Free(frac_deviance_explained);
  Free(prob);
  Free(residual_deviance);
  Free(w);
  Free(z);
  free(cum_degrees);
  free(lambda_seq);
  free(lambdas_alpha);
  free(lambdas_beta);
  free(psis);
  free(resid);
  return result;
} 


extern SEXP gamselFit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] =
  {
   {"gamselFit", (DL_FUNC) &gamselFit, 13},
   {NULL, NULL, 0}
  };

void R_init_gamsel(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
