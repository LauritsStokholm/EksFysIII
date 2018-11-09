#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sort_vector.h>

#include "myheaders.h"
#include "singlegauss.c"

int main (int argc, char **argv){
  /* .......................................................
   * Given the number of peaks (int) and a datafile (string)
   * this function will return the filtered peaks (gaussians).
   * ....................................................... */

  if (argc != 3){
    printf("Expected 3 arguments\n\
        0) ./main\n\
        1) (int) number of peaks,\n\
        2) (string) datafile\n\
    But received %i arguments.", argc);
  };

  int numberofpeaks = atoi(argv[1]);
  FILE *fp = fopen(argv[2], "r");

  if (fp == NULL) {printf("Could not open datafile");};

  singlegauss(numberofpeaks, fp);

  return 0;
}

