/* 
 ============================================================================
 Name        : utility.h
 Author      : Evan Ray
 Version     :
 Copyright   : 
 Description : utility functions for numeric calculations and interacting with R
 ============================================================================
 */

#include <Rdefines.h>


SEXP get_dbl_max();

SEXP getListElement(SEXP list, const char *str);

double logspace_add_safe(double log_x, double log_y);

SEXP logspace_add_C(SEXP x, SEXP y);

SEXP logspace_sub_C(SEXP x, SEXP y);

SEXP logspace_sum_matrix_rows(SEXP Xp, SEXP N_rowp, SEXP N_colp);

SEXP logspace_sub_matrix_rows(SEXP Xp, SEXP N_rowp);
