void data_import(char[] filename){
  /* .......................................................
   * Here we import a data file and save data in arrays
   * ....................................................... */

  // File size (down)
  int count = 0
  while(fscanf(filename, "%lg", &x) != EOF){count++}
}
  // Initialize matrix (memory allocation)
//  gsl_matrix * M1 = gsl_matrix_alloc(n, m);
//
//  // Initialize stream (read filename)
//  FILE * mystream = fopen("matrixvals.txt", "r");
//
//  // File scan into Matrix
//  gsl_matrix_fscanf(mystream, M1);
//  fclose(mystream);
//
//
//  // Initialize read data (domain of functions)
//  double x;
//
//  // Function pointers
//  double (* f)(double);
//  double (* g)(double);
//
//  // Functions
// f = &gsl_sf_airy_Ai;
// g = &gsl_sf_airy_Bi;
//
//  // Reading input
//  FILE * mystream1 = fopen("input_airy.txt", "r");
//  FILE * mystream2 = fopen("output_airy.txt", "w");
//  while(fscanf(mystream1, "%lg", &x) != EOF){
//    fprintf(mystream2, "%lg \t %lg \t %lg\n", x,f(x), g(x));
//  }
//}

//void task2(void){
  /* .......................................................
   * Task 2: Solve a system of linear equation in the matrix
   * form.
   * ( 6.13  -2.90   5.86)(x0)   (6.23)
   * ( 8.08  -6.31  -3.89)(x1) = (5.37)
   * (-4.36   1.00   0.19)(x2)   (2.29)
   * ....................................................... */

  // Type 1: Matrix from input stream
  // Initialize matrix
 // size_t n = 3; // # of rows
 // size_t m = 3; // # of columns

 // // Initialize matrix (memory allocation)
 // gsl_matrix * M1 = gsl_matrix_alloc(n, m);

 // // Initialize stream (read filename)
 // FILE * mystream = fopen("matrixvals.txt", "r");

 // // File scan into Matrix
 // gsl_matrix_fscanf(mystream, M1);
 // fclose(mystream);


  // Type 2: Matrix from array
  // Initialize matrix
//  double matrix_data[] = {6.13, -2.90,  5.86,
//                          8.08, -6.31, -3.89,
//                         -4.36,  1.00,  0.19};
//
//  gsl_matrix_view M2 = gsl_matrix_view_array(matrix_data, n, m);
//
//  for(int i=0; i<n; i++){
//  for(int j=0; j<m; j++){
//    printf("M1(%i, %i) = %lg \t M2(%i, %i) = %lg \n", i, j, gsl_matrix_get(&M2.matrix, i, j), i, j, gsl_matrix_get(M1,i,j));
//  }};
//
//  // Copy of M2 into M3 for solving in next part (need 3 identical matrices)
//  // as each method destroys the given matrix
//  gsl_matrix * M3 = gsl_matrix_alloc(n, m);
//  M3 = gsl_matrix_memcpy(M3, M2);
//
//
//  // Initialize vector (solution)
//  // Mx = b
//  // Initialize b (vector)
//  double vector_data[] = {6.23, 5.37, 2.29};
//  // double * base (here array) and size_t n (size of vector = 3)
//  gsl_vector_view b2 = gsl_vector_view_array(vector_data, n);
//
//  //Solutions:
//  //gsl_linalg_HH_solve
//  //
//  // gsl_linalg_LU_decomp
//  // gsl_linalg_LU_solve 
//  //
//  // gsl_linalg_QR_decomp
//  // gsl_linalg_QR_solve
//
//
//
//  // As for every (M)alloc: remember to free memory
//  gsl_matrix_free(M1);
//}
