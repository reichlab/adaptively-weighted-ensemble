/*
 ============================================================================
 Name        : utility.c
 Author      : Evan Ray
 Version     :
 Copyright   : 
 Description : utility functions for numeric calculations and interacting with R
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>

#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include <math.h>
#include <float.h>

#include "utility.h"

SEXP get_dbl_max() {
	SEXP retval;
	retval = PROTECT(allocVector(REALSXP, 1));
	
	*(REAL(retval)) = DBL_MAX;
	
	UNPROTECT(1);
	return retval;
}

SEXP getListElement(SEXP list, const char *str) {
	// Return the element with name matching str from the list provided.
	// Based on http://cran.r-project.org/doc/manuals/R-exts.html#Handling-lists
	//
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

	for (R_len_t i = 0; i < length(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}

	return elmt;
}


SEXP logspace_add_C(SEXP log_x, SEXP log_y) {
	// An interface to R's C function logspace_add
	// Computes log(exp(log_x) + exp(log_y))
	
	SEXP retval;
	retval = PROTECT(allocVector(REALSXP, 1));
	
	double lx = *(REAL(log_x)), ly = *(REAL(log_y));

	if(lx == R_NegInf && ly == R_NegInf) {
//		Rprintf("Both = -Infty\n");
		*(REAL(retval)) = R_NegInf;
	} else {
		*(REAL(retval)) = logspace_add(lx, ly);
	}
	
	UNPROTECT(1);
	return retval;
}


double logspace_add_safe(double log_x, double log_y) {
	// An interface to R's C function logspace_add
	// Computes log(exp(log_x) + exp(log_y))
	
	if(log_x == -INFINITY && log_y == -INFINITY) {
		return(-INFINITY);
	} else {
		return(logspace_add(log_x, log_y));
	}
}


SEXP logspace_sub_C(SEXP log_x, SEXP log_y) {
	// An interface to R's C function logspace_sub
	// Computes log(exp(log_x) - exp(log_y))
	
	SEXP retval;
	retval = PROTECT(allocVector(REALSXP, 1));
	
	double lx = *(REAL(log_x)), ly = *(REAL(log_y));

	*(REAL(retval)) = logspace_sub(lx, ly);
	
	UNPROTECT(1);
	return retval;
}

SEXP logspace_sum_matrix_rows_C(SEXP Xp, SEXP N_rowp, SEXP N_colp) {
	int i, j, n_row = *INTEGER(N_rowp), n_col = *INTEGER(N_colp);
	SEXP retval = PROTECT(allocVector(REALSXP, n_row));
	double *dblptr = REAL(retval), *X = REAL(Xp);
	
	for(i = 0; i < n_row; i++) {
		*(dblptr + i) = *(X + i);
	}

	for(j = 1; j < n_col; j++) {
		for(i = 0; i < n_row; i++) {
			if(!(*(dblptr + i) == R_NegInf && *(X + i + j*n_row) == R_NegInf))
				*(dblptr + i) = logspace_add_safe(*(dblptr + i), *(X + i + j*n_row));
		}
	}

	UNPROTECT(1);
	return retval;
}

SEXP logspace_sub_matrix_rows_C(SEXP Xp, SEXP N_rowp) {
	int i, n_row = *INTEGER(N_rowp);
	SEXP retval = PROTECT(allocVector(REALSXP, n_row));
	double *dblptr = REAL(retval), *X = REAL(Xp);
	
	for(i = 0; i < n_row; i++) {
		*(dblptr + i) = logspace_sub(*(X + i), *(X + i + n_row));
	}

	UNPROTECT(1);
	return retval;
}


R_CallMethodDef callMethods[] =
{
    {"logspace_add_C", (DL_FUNC)&logspace_sub_C, 2},
    {"logspace_sum_matrix_rows_C", (DL_FUNC)&logspace_sum_matrix_rows_C, 3},
    {"logspace_sub_C", (DL_FUNC)&logspace_sub_C, 2},
    {"logspace_sub_matrix_rows_C", (DL_FUNC)&logspace_sub_matrix_rows_C, 2},
	{NULL,NULL, 0}
};

void R_init_awes(DllInfo *dll)
{
   R_registerRoutines(dll,NULL,callMethods,NULL,NULL);
}
