#include <stdio.h>
#include <gsl/gsl_math.h>

// For doubles
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

int main (void){
  /* .......................................................
   * Here we fit Gaussians
   * ....................................................... */

  // Al 000 000
  //405764
  len = strlen(name);
  dirp = opendir(".");
  while ((dp = readdir(dirp)) != NULL)
  if (dp->d_namlen == len && !strcmp(dp->d_name, name)) {
  (void)closedir(dirp);
  return FOUND;
           }
   (void)closedir(dirp);
   return NOT_FOUND;

  return 0;
}

