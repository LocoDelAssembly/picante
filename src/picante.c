#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <time.h>

#define BITS_PER_WORD (sizeof(int) * 8)
#define WORD_OFFSET(b) ((b) / BITS_PER_WORD)
#define BIT_OFFSET(b) ((b) % BITS_PER_WORD)
#define BITMASK(b) (1 << BIT_OFFSET(b))
#define GET_BIT(a, b) ((a[WORD_OFFSET(b)] & BITMASK(b)) != 0)
#define SET_BIT(a, b) (a[WORD_OFFSET(b)] |= BITMASK(b))
#define CLEAR_BIT(a, b) (a[WORD_OFFSET(b)] &= ~BITMASK(b))
#define ASSIGN_BIT(a, b, value) (value ? SET_BIT(a, b) : CLEAR_BIT(a, b))

/*This takes an integer n and returns one random integer*/
int intrand(int n)
{
    double u;
    u = unif_rand();
    return ((int)(u * n)); /*Cast the double as an integer*/
}

int count_ones(double *v, int row, int column)
{
    int count = 0;

    for (int i = 0; i < row * column; i++)
        count += (v[i] == 1.0);
    return count;
}

/* inefficient but probably not important;
 * allocate and copy in vector to matrix */
double **vectomat(double *v, int row, int column)
{
    int i, j;
    double **m;

    m = (double **)R_alloc(row, sizeof(double *));
    for (i = 0; i < row; i++)
    {
        m[i] = (double *)R_alloc(column, sizeof(double));
        for (j = 0; j < column; j++)
        {
            m[i][j] = v[row * j + i]; /* R uses column-first ordering */
        }
    }
    return (m);
}

int *vectomat_pam(double *v, int row, int column)
{
    int i, j;

    int *bit_vector = (int *)R_alloc(((row * column + BITS_PER_WORD - 1) / BITS_PER_WORD), sizeof(int));
    printf("Allocated %lu bytes for bit vector\n", ((row * column + BITS_PER_WORD - 1) / BITS_PER_WORD) * sizeof(int));
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < column; j++)
        {
            if (v[row * j + i] > 0.0)
                SET_BIT(bit_vector, i * column + j);
            else
                CLEAR_BIT(bit_vector, i * column + j);
        }
    }
    return bit_vector;
}

/* copy matrix back into vector */
void mattovec(double *v, double **m, int row, int column)
{
    int i, j, k;
    k = 0;
    for (j = 0; j < column; j++)
    {
        for (i = 0; i < row; i++)
        {
            v[k++] = m[i][j];
        }
    }
}

/* copy matrix back into vector */
void mattovec_pam(double *v, int *bit_vector, int row, int column)
{
    int i, j, k;
    k = 0;
    for (j = 0; j < column; j++)
    {
        for (i = 0; i < row; i++)
        {
            v[k++] = (double)GET_BIT(bit_vector, i * column + j);
        }
    }
}

void trialswap(double *v, int *pintervals, int *prow, int *pcolumn)
{
    long int trial;
    int i, j, k, l;
    int row, column;
    int intervals;
    double tmp;
    double **m;

    row = *prow;
    column = *pcolumn;
    intervals = *pintervals;

    m = vectomat(v, row, column);

    GetRNGstate();
    for (trial = 0; trial < intervals; trial++)
    {
        i = intrand(row); // Choose a random row
        while ((j = intrand(row)) == i)
            ;                // make sure that you do not randomly choose the same row as before
        k = intrand(column); // Choose a random column
        while ((l = intrand(column)) == k)
            ; // make sure that you do not randomly choose the same column as before
        if ((m[i][k] > 0.0 && m[j][l] > 0.0 && m[i][l] + m[j][k] == 0.0) || (m[i][k] + m[j][l] == 0.0 && m[i][l] > 0.0 && m[j][k] > 0.0))
        {
            // currently swaps abundances within columns (=species)
            // should have a switch to swap abundances within rows, columns, or random
            tmp = m[i][k];
            m[i][k] = m[j][k];
            m[j][k] = tmp;
            tmp = m[i][l];
            m[i][l] = m[j][l];
            m[j][l] = tmp;
        }
    }
    mattovec(v, m, row, column);
    PutRNGstate();
}

void richness(double *v, int *prow, int *pcolumn)
{

    int i, j, k;
    int row, column;
    double tmp;
    double **m;

    row = *prow;
    column = *pcolumn;

    m = vectomat(v, row, column);

    GetRNGstate();
    for (i = 0; i < row; i++)
    {
        for (j = 0; j < column; j++)
        {
            k = intrand(column); // choose another column (species) at random
            tmp = m[i][j];
            m[i][j] = m[i][k];
            m[i][k] = tmp;
        }
    }
    mattovec(v, m, row, column);
    PutRNGstate();
}

void frequency(double *v, int *prow, int *pcolumn)
{

    int i, j, k;
    int row, column;
    double tmp;
    double **m;

    row = *prow;
    column = *pcolumn;

    m = vectomat(v, row, column);

    GetRNGstate();
    for (i = 0; i < column; i++)
    {
        for (j = 0; j < row; j++)
        {
            k = intrand(row); // choose another row (sample) at random
            tmp = m[j][i];
            m[j][i] = m[k][i];
            m[k][i] = tmp;
        }
    }
    mattovec(v, m, row, column);
    PutRNGstate();
}

void independentswap_pam(double *v, int *pintervals, int *prow, int *pcolumn)
{
    long int swap;
    int swapped;
    int i, j, k, l;
    int row, column;
    int intervals;
    double tmp;
    int *bit_vector;

    row = *prow;
    column = *pcolumn;
    intervals = *pintervals;

    bit_vector = vectomat_pam(v, row, column);

    GetRNGstate();
    time_t start_time = time(NULL);
    long int iter_count = 0;
    for (swap = 0; swap < intervals; swap++)
    {
        swapped = 0;
        while (swapped == 0)
        {
            iter_count++;
            if (iter_count % (2 << 24) == 0)
            {
                time_t current_time = time(NULL);
                double seconds = difftime(current_time, start_time);
                if (seconds > 0)
                {
                    fprintf(stderr, "\rIterations per second: %.2f. Swaps: %ld/%d(%.2f%%). Total iterations: %ld|", iter_count / seconds, swap, intervals, swap * 100.0 / intervals, iter_count);
                }
            }
            i = intrand(row); // Choose a random row
            while ((j = intrand(row)) == i)
                ;                // make sure that you do not randomly choose the same row as before
            k = intrand(column); // Choose a random column
            while ((l = intrand(column)) == k)
                ; // make sure that you do not randomly choose the same column as before
            if ((GET_BIT(bit_vector, i * column + k) && GET_BIT(bit_vector, j * column + l) &&
                 !GET_BIT(bit_vector, i * column + l) && !GET_BIT(bit_vector, j * column + k)) ||
                (!GET_BIT(bit_vector, i * column + k) && !GET_BIT(bit_vector, j * column + l) &&
                 GET_BIT(bit_vector, i * column + l) && GET_BIT(bit_vector, j * column + k)))
            {
                // currently swaps abundances within columns (=species)
                // should have a switch to swap abundances within rows, columns, or random
                // Perform the bit swaps
                int tmp = GET_BIT(bit_vector, i * column + k);
                ASSIGN_BIT(bit_vector, i * column + k, GET_BIT(bit_vector, j * column + k));
                ASSIGN_BIT(bit_vector, j * column + k, tmp);

                tmp = GET_BIT(bit_vector, i * column + l);
                ASSIGN_BIT(bit_vector, i * column + l, GET_BIT(bit_vector, j * column + l));
                ASSIGN_BIT(bit_vector, j * column + l, tmp);

                swapped = 1;
            }
        }
    }
    mattovec_pam(v, bit_vector, row, column);
    PutRNGstate();
}

void independentswap(double *v, int *pintervals, int *prow, int *pcolumn)
{
    printf("\n");
    int old_count = count_ones(v, *prow, *pcolumn);
    printf("Number of ones: %d\n", old_count);
    independentswap_pam(v, pintervals, prow, pcolumn);
    if (old_count != count_ones(v, *prow, *pcolumn))
    {
        fprintf(stderr, "Error: count_ones() failed\n");
        exit(1);
    }
    printf("\n");
    return;
    long int swap;
    int swapped;
    int i, j, k, l;
    int row, column;
    int intervals;
    double tmp;
    double **m;

    row = *prow;
    column = *pcolumn;
    intervals = *pintervals;

    m = vectomat(v, row, column);

    GetRNGstate();
    time_t start_time = time(NULL);
    long int iter_count = 0;
    for (swap = 0; swap < intervals; swap++)
    {
        swapped = 0;
        while (swapped == 0)
        {
            iter_count++;
            if (iter_count % (2 << 24) == 0)
            {
                time_t current_time = time(NULL);
                double seconds = difftime(current_time, start_time);
                if (seconds > 0)
                {
                    fprintf(stderr, "\rIterations per second: %.2f. Swaps: %ld/%d(%.2f%%). Total iterations: %ld|", iter_count / seconds, swap, intervals, swap * 100.0 / intervals, iter_count);
                }
            }
            i = intrand(row); // Choose a random row
            while ((j = intrand(row)) == i)
                ;                // make sure that you do not randomly choose the same row as before
            k = intrand(column); // Choose a random column
            while ((l = intrand(column)) == k)
                ; // make sure that you do not randomly choose the same column as before
            if ((m[i][k] > 0.0 && m[j][l] > 0.0 && m[i][l] + m[j][k] == 0.0) || (m[i][k] + m[j][l] == 0.0 && m[i][l] > 0.0 && m[j][k] > 0.0))
            {
                // currently swaps abundances within columns (=species)
                // should have a switch to swap abundances within rows, columns, or random
                tmp = m[i][k];
                m[i][k] = m[j][k];
                m[j][k] = tmp;
                tmp = m[i][l];
                m[i][l] = m[j][l];
                m[j][l] = tmp;
                swapped = 1;
            }
        }
    }
    mattovec(v, m, row, column);
    PutRNGstate();
}
