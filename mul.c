/*
 * mul.c - Matrix multiply C = A * B
 *
 * Each multiply function must have a prototype of the form:
 * void mul(int M, int N, int A[N][M], int B[M][N], int C[N][N]);
 *
 * A multiply function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */
#include <stdio.h>
#include "cachelab.h"

int is_mul(int M, int N, int A[N][M], int B[M][N], int C[N][N]);

/*
 * mul_submit - This is the solution multiply function that you will be
 *     graded on for Part C of the assignment. Do not change the
 *     description string "Multiply submission", as the driver searches
 *     for that string to identify the multiply function to be graded.
 *     C starts out holding garbage, so every element of C must be written.
 *
 * The cache is 32 sets x 1 line x 32 bytes (s = 5, E = 1, b = 5).  With
 * everything packed at a 32-int stride, A row i and C row i both live in
 * sets 4i..4i+3, so the straight i-k-j traversal thrashes on them.  Each
 * tile of 3 rows x 4 columns is accumulated in twelve registers: every k
 * loads one block of B (B[k][j..j+3]) and one block per A row, and the C
 * tile is written exactly once.  The limit is sixteen int locals; twelve
 * accumulators plus the four B values use all sixteen.
 */
char mul_submit_desc[] = "Multiply submission";
void mul_submit(int M, int N, int A[N][M], int B[M][N], int C[N][N])
{
    int acc[12];
    int i, j, k, y, x;
    int b0, b1, b2, b3;

    for (i = 0; i < N; i += 3) {
        for (j = 0; j + 4 <= N; j += 4) {
            for (x = 0; x < 12; x++)
                acc[x] = 0;
            for (k = 0; k < M; k++) {
                b0 = B[k][j + 0];
                b1 = B[k][j + 1];
                b2 = B[k][j + 2];
                b3 = B[k][j + 3];
                for (y = 0; y < 3 && i + y < N; y++) {
                    acc[y * 4 + 0] += A[i + y][k] * b0;
                    acc[y * 4 + 1] += A[i + y][k] * b1;
                    acc[y * 4 + 2] += A[i + y][k] * b2;
                    acc[y * 4 + 3] += A[i + y][k] * b3;
                }
            }
            for (y = 0; y < 3 && i + y < N; y++)
                for (x = 0; x < 4; x++)
                    C[i + y][j + x] = acc[y * 4 + x];
        }
        /* Columns that do not fill a 4-wide tile: plain dot product. */
        for (; j < N; j++) {
            for (y = 0; y < 3 && i + y < N; y++) {
                acc[0] = 0;
                for (k = 0; k < M; k++)
                    acc[0] += A[i + y][k] * B[k][j];
                C[i + y][j] = acc[0];
            }
        }
    }
}

/*
 * mul - A simple baseline multiply function, not optimized for the cache.
 */
char mul_desc[] = "Simple i-j-k multiply";
void mul(int M, int N, int A[N][M], int B[M][N], int C[N][N])
{
    int i, j, k, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            tmp = 0;
            for (k = 0; k < M; k++)
                tmp += A[i][k] * B[k][j];
            C[i][j] = tmp;
        }
    }
}

/*
 * registerFunctions - This function registers your multiply functions
 *     with the driver. At runtime, the driver will evaluate each of the
 *     registered functions and summarize their performance. This is a
 *     handy way to experiment with different multiply strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerMulFunction(mul_submit, mul_submit_desc);

    /* Register any additional multiply functions */
    registerMulFunction(mul, mul_desc);
}

/*
 * is_mul - This helper function checks if C is the product of A and B.
 *     You can check the correctness of your multiply by calling
 *     it before returning from the multiply function.
 */
int is_mul(int M, int N, int A[N][M], int B[M][N], int C[N][N])
{
    int i, j, k, num;

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            num = 0;
            for (k = 0; k < M; k++)
                num += A[i][k] * B[k][j];
            if (num != C[i][j])
                return 0;
        }
    }
    return 1;
}
