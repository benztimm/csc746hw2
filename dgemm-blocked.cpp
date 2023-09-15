const char *dgemm_desc = "Blocked dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm_blocked(int n, int block_size, double *A, double *B, double *C)
{
   // These outer loops iterate over the blocks of the matrices.
   // si, sj, and sk represent the starting row, column, and depth indices of the current blocks.
   for (int si = 0; si < n; si += block_size)
   {
      for (int sj = 0; sj < n; sj += block_size)
      {
         for (int sk = 0; sk < n; sk += block_size)
         {
            // Multiply the individual blocks
            // These inner loops iterate over the individual elements inside the current blocks.
            // i and j iterate over the rows and columns of the current block in matrix C.
            for (int i = si; i < si + block_size; ++i)
            {
               for (int j = sj; j < sj + block_size; ++j)
               {
                  // Initialize the accumulator with the current value in matrix C.
                  double cij = C[i + j * n]; // C is in column-major format
                  // k loop computes the dot product for the blocks of the i-th row of A and the j-th column of B.
                  for (int k = sk; k < sk + block_size; ++k)
                     cij += A[i + k * n] * B[k + j * n]; // A and B are in column-major format
                  C[i + j * n] = cij; // Store the accumulated result back in matrix C.

               }
            }
         }
      }
   }
}
