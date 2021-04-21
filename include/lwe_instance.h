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

#ifndef LWE_INSTANCE_H
#define LWE_INSTANCE_H

#include "utils.h"
#include "config.h"
#include "random_utils.h"

/* max number of rows */
#define MAX_N 50

typedef struct
{
    u16 n;
    u16 q;
    double alpha;
    double sigma;
    short s[MAX_N];

    // random
    rand_ctx ctx;

} lweInstance;

// struct for a simple LWE sample
typedef struct
{
    u16 *a; // (a_1...a_n)
    u16 z;  // sum(a_i * s_i) + error
    u16 error;
} sample;

typedef struct
{
    sample *list;
    u32 n_samples;
} samplesList;

typedef struct
{
    sample *list;
    u16 n_samples;
} category;

typedef struct
{
    category *list_categories;
    u64 n_categories;
    u8 n_samples_per_category;
    u64 n_samples;
    u64 max_samples;
} sortedSamplesList;

/******************* BKW STEP PARAMETERS ***********************/

typedef struct
{
    int startIndex;
    int numPositions; /* = Ni, in case of codedBKW */
    short p; /* general reduction factor */
    short p1; /* additional reduction factor for position Ni+1 */
    short p2; /* additional reduction factor for first position after the first step */
    short prev_p1; /* reduction factor for position Ni+1 in the previous step - used on all but the first steps */
    short un_selection;
} bkwStepParameters;

/******************* BKW STEP PARAMETERS ***********************/

// return the number of categories
u64 num_categories(lweInstance *lwe, bkwStepParameters *bkwStepPar);

// precompute CDF tables - must be run before lwe_init
void precompute_cdf_table(double sigma);

// sample s at random
void lwe_init(lweInstance *lwe, u16 n, u16 q, double alpha);

// create and allocate lwe samples
void create_lwe_samples(samplesList *Samples, lweInstance *lwe, u64 n_samples);

// allocate struct for unsorted lwe samples
void allocate_samples_list(samplesList *Samples, lweInstance *lwe, u64 n_samples);

// allocate memory for sorted samples
void allocate_sorted_samples_list(sortedSamplesList *Samples, lweInstance *lwe, bkwStepParameters *bkwStepPar, u64 n_samples, u64 max_categories);

// free samples
void free_samples(samplesList *Samples);

// free sorted samples
void free_sorted_samples(sortedSamplesList *Samples, u64 max_categories);

void set_sorted_samples_list(sortedSamplesList *Samples, lweInstance *lwe, bkwStepParameters *bkwStepPar, u64 n_samples);

void clean_sorted_samples(sortedSamplesList *Samples);


#endif
