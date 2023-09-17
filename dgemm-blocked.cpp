const char *dgemm_desc = "Blocked dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are n-by-n matrices stored in column-major format.
 * On exit, A and B maintain their input values. */
void square_dgemm_blocked(int n, int block_size, double *A, double *B, double *C)
{
    double A_block[block_size * block_size];
    double B_block[block_size * block_size];

    for (int si = 0; si < n; si += block_size)
    {
        for (int sj = 0; sj < n; sj += block_size)
        {
            for (int sk = 0; sk < n; sk += block_size)
            {
                // Copy blocks of A and B into local storage
                for (int i = 0; i < block_size; ++i)
                {
                    for (int j = 0; j < block_size; ++j)
                    {
                        A_block[i + j * block_size] = A[(si + i) + (sk + j) * n];
                        B_block[i + j * block_size] = B[(sk + i) + (sj + j) * n];
                    }
                }

                // Multiply the individual blocks using local storage
                for (int i = 0; i < block_size; ++i)
                {
                    for (int j = 0; j < block_size; ++j)
                    {
                        double cij = C[(si + i) + (sj + j) * n];
                        for (int k = 0; k < block_size; ++k)
                        {
                            cij += A_block[i + k * block_size] * B_block[k + j * block_size];
                        }
                        C[(si + i) + (sj + j) * n] = cij;
                    }
                }
            }
        }
    }
}
