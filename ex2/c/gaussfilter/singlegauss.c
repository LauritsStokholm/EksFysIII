/*int singlegauss(int numberofpeaks, char *fname, FILE *fp){ */
int singlegauss(int numberofpeaks, FILE *fp){
  /* .......................................................
   * Given the integer number of peaks to find and a dataset 
   * (string) this function filters the single gaussians of 
   * greatest amplitude. It is divided into following the following;
   *
   * 1) Initialize & Load data.
   * 2) Greatest peak VAL and IDX is found.
   * 3) Set half max val (VAL/2) and find all points near this point.
   * 4) Subtract from original data and square (now all is positive)
   * 5) Minima (zero points) will be on the line of half max, these are
   * weighted, by distance from centroid (closer to centroid is better)
   * and by most importantly, the val (best to be at zero corresponding to
   * VAL/2).
   * 6) Now maximas are found (corresponding to real zeros, as they are maximal
   * away from the VAL/2. The maxima within 2 FWHM is the first zero, which is
   * used.
   * 7) Full interval of gauss is twice the distance from centroid to index of
   * max (or zero point), and the interval is offset the index of zeropoint.
   * 8) Vector values (Not saved to file) within gauss interval is set to zero,
   * and next peak is found. Itterates through given amount of peaks.
   * ....................................................... */


  /* This is to open a file, which stores original stdout (we use stdout to save
   * only relevant gauss data (gnuplot graphics) */
  FILE* stdout2 = fopen("stdout2.txt", "w");

  /* initializing variables */
  double x, y;

  /* Obtain size of file */
  int size_of_file = 0;
  while(fscanf(fp, "%lg, %lg", &x, &y)!=EOF){
    size_of_file++;
  }
  rewind(fp);

  fprintf(stdout2,"%i\n", size_of_file);
  gsl_vector *idx = gsl_vector_alloc(size_of_file);
  gsl_vector *val = gsl_vector_alloc(size_of_file);

  /* Obtain data in vector */
  int i = 0;
  while(fscanf(fp, "%lg, %lg", &x, &y)!=EOF){
    gsl_vector_set(idx, i, x); /* Channel numbers */
    gsl_vector_set(val, i, y); /* Count numbers */
    i++;
  }
  rewind(fp);

  /* =========================================================================
  NOW WE HAVE: Data (IDX and VAL) loaded into vectors
  ============================================================================*/


  /* IMPORTANT ITTERATOR (ALL GAUSSIANS DEPEND ON THIS) */
  int J;
  for (J=0; J<numberofpeaks; J++){
  printf("# <<START%i>>\n", J);
  fprintf(stdout2, "<<START%i>>\n", J);


  /* Find data peak and corresponding index.*/
  int max_val = gsl_vector_max(val);
  int max_val_idx = gsl_vector_max_index(val);

  /* Save half maximum value in variable */
  int half_max_val = max_val / 2;

  /* Check values */
  fprintf(stdout2, "The greatest countnumber is %d at channelnumber %d, and the halfmax is %d\n", max_val, max_val_idx, half_max_val);

  /* =========================================================================
  NOW WE HAVE: Data (MAXIDX and MAXVAL)
  NOW WE WANT: Subtract halfmax from data, square (data>0), and find FWHM
  (Other irrelevant data points might have y = VAL/2, and these are to be
  filtered. Thus a weight of both (delta_val = 0) and a small distance to
  gaussian centroid is important!)
  ============================================================================*/

  /* Initialize new vector */
  gsl_vector *delta_val = gsl_vector_alloc(size_of_file);

  /* Copy elements of vector src val into vector dest new_val */
  gsl_vector_memcpy(delta_val, val);
  gsl_vector_add_constant(delta_val, -half_max_val);
  gsl_vector_mul(delta_val, delta_val);

  /* If there is a zero, it must be the half max */
  /* Store all data, and sort them by closest to Half max and idx near centroid
   * */

  /* Not many values will be at half max (compared to size of file). We guess
   * (many more than needed) about 100 */

  int n = 10;

  double delta_val_small[n]; /* Smallest values of delta_val (ascending) */
  size_t delta_val_small_idx[n]; /* Indices of smallest values */

 /* delta_val is very important (square), such we are at Half Max value.
  * Index is also very important to make sure, data is at right gaussian.*/

  /* Weights vals and indices*/
  gsl_vector *weights = gsl_vector_alloc(n);
  size_t weights_idx[n];

  /* Find n smallest vals in delta_val (VAL/2) and save values/indices */
  gsl_sort_vector_smallest(delta_val_small, n, delta_val);
  gsl_sort_vector_smallest_index(delta_val_small_idx, n, delta_val);

  /* Calculate weights (points near halfmaxval and idx near centroid). */
  int k; double a; double b;
  for (k=0; k<n; k++){
    /* Deviate from halfmaxval (squared) */
    a = delta_val_small[k] * delta_val_small[k];

    /* Deviate from centroid (squared) */
    b = abs((delta_val_small_idx[k] - max_val_idx)); //
    /* Calculate and set weight */

    gsl_vector_set(weights, k, a + b);
    fprintf(stdout2, "%lg \t %li \t %lg \n",
        delta_val_small[k],
        delta_val_small_idx[k],
        gsl_vector_get(weights, k));
  }

  /* Point of smallest weight is best fitting as FWHM point */
  gsl_sort_vector_smallest_index(weights_idx, n, weights);

  fprintf(stdout2, "Point of FWHM: (%li, %lg) equivalent to (%ti, %lg)\n",
      weights_idx[0],
      delta_val_small[weights_idx[0]],
      delta_val_small_idx[weights_idx[0]],
      gsl_vector_get(val, delta_val_small_idx[weights_idx[0]]));

  /* Testing abnormalities */
  if(delta_val_small[weights_idx[0]] == delta_val_small[0]){
    fprintf(stdout2, "Everything as usual with the FWHM\n");
  }
  else{
    fprintf(stdout2, "Something is off! Take a look at FWHM (graph?)\n");
  }

  double 
  FWHM = abs(2 * (max_val_idx - (int) delta_val_small_idx[weights_idx[0]]));

  fprintf(stdout2, "The full width half maxima is %lg\n", FWHM);

  /* =========================================================================
  NOW WE HAVE: FWHM of Gauss.
  NOW WE WANT: Find the lowest points (zeros) of gauss.
  ============================================================================*/

  /*
   When gauss reaches the roots, it will be at a maximal distance from the
   half_max_val (as non negative vals stored in delta_val),
   so now we find the index of the closest maxima from point of FWHM
   (away from max_val_idx) to obtain index of gauss zero
  */


  /* Prepare vector for gaussian tail */

  /* As calculated, weights_idx[0] is index ranging from 0 to n-1
   * of FWHM point in p array. Translated to original indices: */
  int offset = delta_val_small_idx[weights_idx[0]];

  /* We look after 0 points (originally) by finding peak (in delta_val)
   * and we guess this is "l" away from FWHM */
  size_t l = 2*FWHM;
  gsl_vector * tail_vec = gsl_vector_alloc(l);

  /* FWHM might be right or left to centroid:
   * if left then lower index
   * If right then increase index */

  if (max_val_idx - offset > 0){
    fprintf(stdout2, "The FWHM point is Left to center\n");
    for (i = 0; i < (int) l ; i++){
      gsl_vector_set(tail_vec, i, gsl_vector_get(delta_val, offset - i));
    }
  };
  if (max_val_idx - offset < 0){
    fprintf(stdout2, "The FWHM point is Right to center\n");
    for (i = 0; i < (int) l ; i++){
      gsl_vector_set(tail_vec, i, gsl_vector_get(delta_val, offset + i));
    }
  };

  /* The first zero should do alright*/
  double endpoint[l]; /* To hold the (many possible) end point vals of gauss tail */
  size_t endpoint_idx[l]; /* To hold the indices of values */

  /* Tail_vec is sorted in descending order */
  gsl_sort_vector_largest(endpoint, l, tail_vec);
  gsl_sort_vector_largest_index(endpoint_idx, l, tail_vec);


  int root_index;
  if (max_val_idx - offset > 0){
    root_index = offset -  (int) endpoint_idx[0];
  }

  if (max_val_idx - offset < 0){
    root_index = offset + (int) endpoint_idx[0];
  }

  double end_val = gsl_vector_get(val, root_index);
  double end_idx = gsl_vector_get(idx, root_index);

  /* Point of smallest weight is best fitting as FWHM point */
//  gsl_sort_vector_smallest_index(weights_tail_idx, l, weights_tail);

  fprintf(stdout2, "The endpoint index is (%ti, %lg) correspondingly (%lg, %lg)\n",
      endpoint_idx[0], endpoint[0], end_idx, end_val);

  /* NOW the whole interval is obtained */

  int size_of_gauss;
  size_of_gauss = 2 * abs(max_val_idx - end_idx);
  fprintf(stdout2, "Size of gaussian is %i\n", size_of_gauss);

  /* Saving gaussian data */
  gsl_matrix *gauss_data = gsl_matrix_alloc(size_of_gauss, 2);
  int index;
  int idx_low = abs(max_val_idx - (size_of_gauss / 2));
  int idx_hgh = abs(max_val_idx + (size_of_gauss / 2));

  i = 0;
  for (index = idx_low ; index < idx_hgh; index++){
    gsl_matrix_set(gauss_data, i, 0, gsl_vector_get(idx, index));
    gsl_matrix_set(gauss_data, i, 1, gsl_vector_get(val, index));
    gsl_vector_set(val, index, 0);
    i++;
  };

  for (i = 0; i<size_of_gauss; i++){

    printf("%lg, %lg\n",
        gsl_matrix_get(gauss_data, i, 0),
        gsl_matrix_get(gauss_data, i, 1));
//    printf("%lg\n", gsl_matrix_get(gauss_data, i, 0));
  };

  printf("<<END%i>>\n", J);
  /* This should be freed before next itterator */
  gsl_vector_free(delta_val);
  gsl_vector_free(weights);
  gsl_vector_free(tail_vec);
  gsl_matrix_free(gauss_data);
  };

  /* This is closed when done with all peaks */
  fclose(stdout2);
  fclose(fp);
  gsl_vector_free(val);
  gsl_vector_free(idx);
  return 0;
}
