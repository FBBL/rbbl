/*  This file is part of RBBL (RAM-Based BKW for LWE).
 *
 *  RBBL is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RBBL is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Nome-Programma.  If not, see <http://www.gnu.org/licenses/>
 */

#include "lwe_instance.h"
#include "utils.h"
#include "random_utils.h"

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define PI 3.14159265359

static int roundInt(double d)
{
    return d > 0.0 ? d + 0.5 : d - 0.5;
}

/* adapted malb */
int chi(double sigma, rand_ctx *rnd)
{
#if 1
    const double aa = randomUtilDouble(rnd);
    const double a = aa == 0 ? 0 : sqrt(-2 * log(aa));
    const double b = 2 * PI * randomUtilDouble(rnd);
    const double x = sigma * a * cos(b);
    //const double y = sigma * a * sin(b);
    int ret = roundInt(x);
    return ret;
#else
    const double aa = rand() / (double)RAND_MAX;
    const double a = aa == 0 ? 0 : sqrt(-2 * log(aa));
    const double b = 2 * PI * rand() / (double)RAND_MAX;
    const double x = sigma * a * cos(b);
    //const double y = sigma * a * sin(b);
    int ret = roundInt(x);
    return ret;
#endif
}

/* Initialize LWE struct, store secret (distributed as the noise) */
void lwe_init(lweInstance *lwe, u16 n, u16 q, double alpha){

    lwe->n = n;
    lwe->q = q;
    lwe->alpha = alpha;
    lwe->sigma = alpha*q;

    srand(time(NULL));
    randomUtilRandomize();
    randomUtilInit(&lwe->ctx);

    for (u16 i = 0; i < n; i++)
        lwe->s[i] = (u16)(randomUtil64(&lwe->ctx));

    unsigned int i, j;

    for (i = 0; i < n; i++)
    {
        lwe->s[i] = (chi(lwe->sigma, NULL) + q) % q;
        printf("%d ", lwe->s[i]);
    }printf("\n");
}

/* Create and allocate n_samples in list */
void create_lwe_samples(unsortedSamplesList *Samples, lweInstance *lwe, u64 n_samples){

    u16 n = lwe->n;
    u16 q = lwe->q;

    Samples->n_samples = n_samples;
    Samples->max_samples = n_samples;
    Samples->a_list = calloc(n_samples*n, sizeof(u16));
    Samples->z_list = calloc(n_samples, sizeof(u16));
    Samples->e_list = calloc(n_samples, sizeof(u16));

    int seed = get_seed(), in;
    u16 error;
    srand((unsigned) seed);

    for (int i = 0; i < n_samples; i++)
    {
        in = i*n;
        Samples->z_list[i] = 0;

        for (int j = 0; j < n; j++)
        {
            Samples->a_list[in+j] = (u16)randomUtilInt(&lwe->ctx, q);
            Samples->z_list[i] = (Samples->z_list[i] + Samples->a_list[in+j]*lwe->s[j]) % q;
        }

        // sample error
        Samples->e_list[i] = (chi(lwe->sigma, &lwe->ctx) + q) % q;
 
        // z = a*s + e mod q
        Samples->z_list[i] = (Samples->z_list[i] + Samples->e_list[i]) % q;
    }

}


/* allocate memory for sorted samples. n_samples is the nuber of total sampels in input before sorting/bkwstep */
void allocate_unsorted_samples_list(unsortedSamplesList *Samples, lweInstance *lwe, u64 n_samples){

    Samples->n_samples = 0;
    Samples->a_list = malloc(n_samples*lwe->n * sizeof(u16));
    Samples->z_list = malloc(n_samples * sizeof(u16));
    Samples->e_list = malloc(n_samples * sizeof(u16));

    ASSERT(Samples->a_list != NULL, "Failed allocation");
    ASSERT(Samples->z_list != NULL, "Failed allocation");
    ASSERT(Samples->e_list != NULL, "Failed allocation");

    Samples->max_samples = n_samples;
}


/* free samples */
void free_samples(unsortedSamplesList *Samples){

    free(Samples->a_list);
    free(Samples->z_list);
    free(Samples->e_list);

    Samples->n_samples = 0;
}

/* free sortedSamplesList */
void free_sorted_samples(sortedSamplesList *Samples, u64 max_categories){

    free(Samples->a_list);
    free(Samples->z_list);
    free(Samples->e_list);
    free(Samples->n_in_categories);

    Samples->n_samples = 0;
    Samples->max_samples = 0;
    Samples->n_categories = 0;

}

void clean_sorted_samples(sortedSamplesList *Samples){

    memset(Samples->n_in_categories, 0, sizeof(u8)*Samples->n_categories);

    Samples->n_samples = 0;
    Samples->max_samples = 0;
    Samples->n_categories = 0;
}


/*************** BKW STEP PARAMETERS ******************/

#define MIN(a,b) (((a)<(b))?(a):(b))

u64 num_categories(lweInstance *lwe, bkwStepParameters *bkwStepPar)
{
    u64 numCategories = 0;
    int p = bkwStepPar->p;
    int p1 = bkwStepPar->p1;
    int p2 = bkwStepPar->p2;
    ASSERT(lwe->q & 1, "Modulo power of 2 not handled!");
    int q_ = lwe->q%2 == 1 ? (lwe->q+1)/2 : lwe->q/2;
    u64 c = ((2*q_-1) % p) == 0 ? ((2*q_-1) / p) : ((2*q_-1) / p) + 1;
    u64 c1 = ((2*q_-1) % p1) == 0 ? ((2*q_-1) / p1) : ((2*q_-1) / p1) + 1;
    if (bkwStepPar->prev_p1 == -1)   // first step
    {
        numCategories = 1;
        int lastPosition = MIN(bkwStepPar->numPositions, bkwStepPar->numPositions + 1);
        for (int i=0; i<lastPosition; i++){
            numCategories *= c;
        }
        numCategories *= c1;
    }
    else if (bkwStepPar->startIndex + bkwStepPar->numPositions == lwe->n)     // last step
    {
        q_ = bkwStepPar->prev_p1;
        u64 c2 = ((2*q_-1) % p2) == 0 ? ((2*q_-1) / p2) : ((2*q_-1) / p2) + 1;
        numCategories = c2;
        for (int i=1; i<bkwStepPar->numPositions; i++)
            numCategories *= c;
    }
    else     // other steps
    {
        q_ = bkwStepPar->prev_p1;
        u64 c2 = ((2*q_-1) % p2) == 0 ? ((2*q_-1) / p2) : ((2*q_-1) / p2) + 1;
        numCategories = c2;
        int lastPosition = MIN(bkwStepPar->numPositions, bkwStepPar->numPositions + 1);
        for (int i=1; i<lastPosition; i++)
            numCategories *= c;
        numCategories *= c1;
    }
    return numCategories; 
}


/* allocate memory for sorted samples. n_samples is the nuber of total sampels in input before sorting/bkwstep */
void allocate_sorted_samples_list(sortedSamplesList *Samples, lweInstance *lwe, bkwStepParameters *bkwStepPar, u64 n_samples, u64 max_categories){

    // careful with the following setting... could be modified
    Samples->n_categories = num_categories(lwe, bkwStepPar);
    Samples->max_samples = MIN(n_samples + ceil(SAMPLES_INCREASE_FACTOR*n_samples), max_categories*6); // allow some more samples at each step

    Samples->a_list = malloc(lwe->n*max_categories*6 * sizeof(u16));
    Samples->z_list = malloc(max_categories*6 * sizeof(u16));
    Samples->e_list = malloc(max_categories*6 * sizeof(u16));
    ASSERT(Samples->a_list != NULL, "Failed allocation");
    ASSERT(Samples->z_list != NULL, "Failed allocation");
    ASSERT(Samples->e_list != NULL, "Failed allocation");

    Samples->n_in_categories = calloc(max_categories, sizeof(u8));
    Samples->n_samples = 0;
}

/* allocate memory for sorted samples. n_samples is the nuber of total sampels in input before sorting/bkwstep */
void set_sorted_samples_list(sortedSamplesList *Samples, lweInstance *lwe, bkwStepParameters *bkwStepPar, u64 n_samples, u64 max_categories){

    // careful with the following setting... could be modified
    Samples->n_categories = num_categories(lwe, bkwStepPar);
    Samples->max_samples = n_samples + ceil(SAMPLES_INCREASE_FACTOR*n_samples); // allow some more samples at each step

    memset(Samples->n_in_categories, 0, sizeof(u8)*max_categories);
    Samples->n_samples = 0;
}








