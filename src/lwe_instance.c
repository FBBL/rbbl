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

#include <stdio.h>
#include <math.h>
#include <time.h>

#define LEN_MAX_CDF 200
#define SAMPLES_PER_CATEGORY 6

static u32 CDF_TABLE[LEN_MAX_CDF];

static int CDF_TABLE_LEN = 0;

/* Gaussian Distribution PDF */
const double pi = 3.141592653589793;
const double e = 2.718281828459045;
double normal_PDF(int x, double mean, double sigma)
{
    return (1/(sigma*sqrt(2*pi)))*exp(-0.5*(((x-mean)/sigma)*((x-mean)/sigma)));
}

void precompute_cdf_table(double sigma)
{

    CDF_TABLE_LEN = 4*sigma+2;
    if (CDF_TABLE_LEN>LEN_MAX_CDF)
    {
        printf("ERROR: CDF_TABLE_LEN>LEN_MAX_CDF\n");
        exit(0);
    }

    int len_x = 16;
    double sum;

    CDF_TABLE[0] = (1<<(len_x-1))*normal_PDF(0, 0, sigma)-1;
    for (int z = 1; z < CDF_TABLE_LEN; z++)
    {
        sum = 0;
        for (int i = 1; i <= z; i++)
            sum += normal_PDF(i, 0, sigma);
        CDF_TABLE[z] = CDF_TABLE[0]+(1<<len_x)*sum;
    }

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

    if(!CDF_TABLE_LEN){
        printf("ERROR: one must first precompute CDF table\n");
        exit(0);
    }

    for (i = 0; i < n; i++)
    {
        u16 sample = 0;
        u16 prnd = lwe->s[i] >> 1;    // Drop the least significant bit
        u16 sign = lwe->s[i] & 0x1;    // Pick the least significant bit

        // No need to compare with the last value.
        for (j = 0; j < (unsigned int)(CDF_TABLE_LEN - 1); j++)
        {
            // Constant time comparison: 1 if CDF_TABLE[j] < s, 0 otherwise. Uses the fact that CDF_TABLE[j] and s fit in 15 bits.
            sample += (u16)(CDF_TABLE[j] - prnd) >> 15;
        }
        // Assuming that sign is either 0 or 1, flips sample iff sign = 1. WARNING: NOT CONSTANT TIME
        lwe->s[i] = sign ? sample : (q-sample) % q;
    }

}

/* Create and allocate n_samples in list */
void create_lwe_samples(samplesList *Samples, lweInstance *lwe, int n_samples){

    Samples->n_samples = n_samples;
    Samples->list = calloc(n_samples, sizeof(sample));

    u16 n = lwe->n;
    u16 q = lwe->q;

    int seed = get_seed();
    srand((unsigned) seed);

    for (int i = 0; i < n_samples; i++)
    {
        Samples->list[i].a = calloc(lwe->n, sizeof(u16));
        Samples->list[i].z = 0;

        for (int j = 0; j < n; j++)
        {
            Samples->list[i].a[j] = (u16)randomUtilInt(&lwe->ctx, q);
            Samples->list[i].z = (Samples->list[i].z + Samples->list[i].a[j]*lwe->s[j]) % q;
        }

        // sample error
        Samples->list[i].error = (u16)(rand());
        u16 sample = 0;
        u16 prnd = Samples->list[i].error >> 1;    // Drop the least significant bit
        u16 sign = Samples->list[i].error & 0x1;    // Pick the least significant bit

        // No need to compare with the last value.
        for (int j = 0; j < (unsigned int)(CDF_TABLE_LEN - 1); j++)
        {
            // Constant time comparison: 1 if CDF_TABLE[j] < s, 0 otherwise. Uses the fact that CDF_TABLE[j] and s fit in 15 bits.
            sample += (u16)(CDF_TABLE[j] - prnd) >> 15;
        }
        // Assuming that sign is either 0 or 1, flips sample iff sign = 1. WARNING: NOT CONSTANT TIME
        Samples->list[i].error = sign ? sample : (q-sample) % q;

        // z = a*s + e mod q
        Samples->list[i].z = (Samples->list[i].z + Samples->list[i].error) % q;
    }

}


/* allocate memory for sorted samples. n_samples is the nuber of total sampels in input before sorting/bkwstep */
void allocate_samples_list(samplesList *Samples, lweInstance *lwe, int n_samples){

    Samples->n_samples = n_samples;
    Samples->list = calloc(n_samples, sizeof(sample));

    u16 n = lwe->n;
    u16 q = lwe->q;
}


/* free samples */
void free_samples(samplesList *Samples){

    for (int i = 0; i < Samples->n_samples; i++)
        free(Samples->list[i].a);

    free(Samples->list);

    Samples->n_samples = 0;
}

/* free sortedSamplesList */
void free_sorted_samples(sortedSamplesList *Samples, u64 max_categories){

    for (int i = 0; i < max_categories; i++)
    {
        for (int j = 0; j < SAMPLES_PER_CATEGORY; j++)
            free(Samples->list_categories[i].list[j].a);
        free(Samples->list_categories[i].list);
        Samples->list_categories[i].n_samples = 0;
    }
    free(Samples->list_categories);

    Samples->n_samples = 0;
    Samples->n_samples_per_category = 0;
    Samples->max_samples = 0;
    Samples->n_categories = 0;

}

void clean_sorted_samples(sortedSamplesList *Samples){

    for (int i = 0; i < Samples->n_categories; i++)
        Samples->list_categories[i].n_samples = 0;

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
    Samples->n_samples_per_category = SAMPLES_PER_CATEGORY;//5*ceil((double)n_samples/Samples->n_categories)+1; // no need to store too many samples for each category
    ASSERT(Samples->n_samples_per_category >= 2, "Not enough samples");

    Samples->max_samples = n_samples + ceil(0.05*n_samples); // allow some more samples at each step

    Samples->list_categories = calloc(max_categories, sizeof(category));
    ASSERT(Samples->list_categories != NULL, "Failed allocation");
    for (int i = 0; i < max_categories; ++i){
        Samples->list_categories[i].list = calloc(Samples->n_samples_per_category, sizeof(sample));
        ASSERT(Samples->list_categories[i].list != NULL, "Failed allocation");
        for (int j = 0; j < Samples->n_samples_per_category; j++)
            Samples->list_categories[i].list[j].a = malloc(lwe->n*sizeof(u16));
    }
    Samples->n_samples = 0;
}

/* allocate memory for sorted samples. n_samples is the nuber of total sampels in input before sorting/bkwstep */
void set_sorted_samples_list(sortedSamplesList *Samples, lweInstance *lwe, bkwStepParameters *bkwStepPar, int n_samples){

    // careful with the following setting... could be modified
    Samples->n_categories = num_categories(lwe, bkwStepPar);
    Samples->n_samples_per_category = SAMPLES_PER_CATEGORY; // no need to store too many samples for each category
    Samples->max_samples = n_samples + ceil(0.05*n_samples); // allow some more samples at each step
    Samples->n_samples = 0;
}








