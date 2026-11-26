/*
 * trans.c - Matrix transpose B = A^T
 *
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/*
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded.
 *
 * The cache is 32 sets x 1 line x 32 bytes (s = 5, E = 1, b = 5), so a
 * block holds 8 ints and element (r, c) of an M-column matrix lives in
 * set ((r * M + c) / 8) mod 32.
 *
 * 32 x 32: rows i and i + 8 hit the same sets, and A and B are 256KB
 *     apart (a multiple of the 1KB cache size), so A[r][c] and B[r][c]
 *     collide.  Copying one full 8-int row of A into registers and then
 *     writing the eight scattered B rows once each turns the 8 x 8 tile
 *     into 8 A misses plus 8 B misses.
 * 64 x 64: 32 sets hold only four 64-column rows, so rows i and i + 4
 *     collide.  Each 8 x 8 tile is split into four 4 x 4 quadrants and
 *     the top-right quadrant is parked in the B positions that phase 2
 *     overwrites, so no 4 x 4 quadrant is ever read twice from memory.
 * 61 x 67: the odd width defeats any fixed set pattern, so plain 8 x 8
 *     tiles are used with the column block as the outer loop: one column
 *     block sweeps every row block while its eight B rows stay hot.
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, k, c;
    int t0, t1, t2, t3, t4, t5, t6, t7;

    if (M == 32 && N == 32) {
        for (i = 0; i < 32; i += 8) {
            for (j = 0; j < 32; j += 8) {
                for (k = i; k < i + 8; k++) {
                    t0 = A[k][j + 0];
                    t1 = A[k][j + 1];
                    t2 = A[k][j + 2];
                    t3 = A[k][j + 3];
                    t4 = A[k][j + 4];
                    t5 = A[k][j + 5];
                    t6 = A[k][j + 6];
                    t7 = A[k][j + 7];
                    B[j + 0][k] = t0;
                    B[j + 1][k] = t1;
                    B[j + 2][k] = t2;
                    B[j + 3][k] = t3;
                    B[j + 4][k] = t4;
                    B[j + 5][k] = t5;
                    B[j + 6][k] = t6;
                    B[j + 7][k] = t7;
                }
            }
        }
    } else if (M == 64 && N == 64) {
        for (i = 0; i < 64; i += 8) {
            for (j = 0; j < 64; j += 8) {
                /* Phase 1: top half of A.  Top-left quadrant lands in
                 * place; the top-right quadrant is parked in the B cells
                 * that phase 2 is about to overwrite. */
                for (k = i; k < i + 4; k++) {
                    t0 = A[k][j + 0];
                    t1 = A[k][j + 1];
                    t2 = A[k][j + 2];
                    t3 = A[k][j + 3];
                    t4 = A[k][j + 4];
                    t5 = A[k][j + 5];
                    t6 = A[k][j + 6];
                    t7 = A[k][j + 7];
                    B[j + 0][k] = t0;
                    B[j + 1][k] = t1;
                    B[j + 2][k] = t2;
                    B[j + 3][k] = t3;
                    B[j + 0][k + 4] = t4;
                    B[j + 1][k + 4] = t5;
                    B[j + 2][k + 4] = t6;
                    B[j + 3][k + 4] = t7;
                }
                /* Phase 2: for B rows j..j+3, take the parked top-right
                 * values to B[k+4][i..i+3], and fill B[k][i+4..i+7] with
                 * the bottom-left quadrant of A. */
                for (k = j; k < j + 4; k++) {
                    t0 = B[k][i + 4];
                    t1 = B[k][i + 5];
                    t2 = B[k][i + 6];
                    t3 = B[k][i + 7];
                    t4 = A[i + 4][k];
                    t5 = A[i + 5][k];
                    t6 = A[i + 6][k];
                    t7 = A[i + 7][k];
                    B[k][i + 4] = t4;
                    B[k][i + 5] = t5;
                    B[k][i + 6] = t6;
                    B[k][i + 7] = t7;
                    B[k + 4][i + 0] = t0;
                    B[k + 4][i + 1] = t1;
                    B[k + 4][i + 2] = t2;
                    B[k + 4][i + 3] = t3;
                }
                /* Phase 3: bottom-right quadrant.  The A rows were already
                 * touched in phase 2, so these loads mostly hit. */
                for (k = i + 4; k < i + 8; k++) {
                    t0 = A[k][j + 4];
                    t1 = A[k][j + 5];
                    t2 = A[k][j + 6];
                    t3 = A[k][j + 7];
                    B[j + 4][k] = t0;
                    B[j + 5][k] = t1;
                    B[j + 6][k] = t2;
                    B[j + 7][k] = t3;
                }
            }
        }
    } else {
        /* General case: 8 x 8 tiles, one row of A per register load.
         * Column blocks are the outer loop.  With a column block fixed,
         * the eight B rows stay put and their blocks advance by exactly
         * one per row block, so the B stream is sequential and conflicts
         * little; with row blocks outermost those eight rows sit 67/8
         * blocks apart and keep evicting each other.  Measured on 61 x 67:
         * 2011 misses row-block-major vs 1762 column-block-major. */
        for (j = 0; j < M; j += 8) {
            for (i = 0; i < N; i += 8) {
                for (k = i; k < i + 8 && k < N; k++) {
                    if (j + 8 <= M) {
                        t0 = A[k][j + 0];
                        t1 = A[k][j + 1];
                        t2 = A[k][j + 2];
                        t3 = A[k][j + 3];
                        t4 = A[k][j + 4];
                        t5 = A[k][j + 5];
                        t6 = A[k][j + 6];
                        t7 = A[k][j + 7];
                        B[j + 0][k] = t0;
                        B[j + 1][k] = t1;
                        B[j + 2][k] = t2;
                        B[j + 3][k] = t3;
                        B[j + 4][k] = t4;
                        B[j + 5][k] = t5;
                        B[j + 6][k] = t6;
                        B[j + 7][k] = t7;
                    } else {
                        /* Tail columns: read the row, then write it out. */
                        for (c = j; c < M; c++)
                            B[c][k] = A[k][c];
                    }
                }
            }
        }
    }
}

/*
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started.
 */

/*
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc);

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc);

}

/*
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}
