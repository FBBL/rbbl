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


#include "config.h"
#include "lwe_instance.h"
#include "utils.h"
#include "transition_times2_modq.h"
#include "transition_bkw_step_smooth_lms.h"
#include "transition_bkw_step_final.h"
#include "lookup_tables.h"
#include "solve_fwht.h"

#include <stdio.h>
#include <stdlib.h>

#define NUM_REDUCTION_STEPS 13
#define BRUTE_FORCE_POSITIONS 0
#define ZERO_POSITIONS 16

int main()
{
    u64 n_samples = 45000000;
    u64 samples_for_guessing = 45000000;

    lweInstance lwe;
    int n = 40;
    int q = 1601;
    double alpha = 0.005;

    time_stamp("LWE parameters: n: %d, q: %d, sigma: %lf*q. Initial samples: %llu", n, q, alpha, n_samples);

    time_stamp("Precomputation");
    precompute_cdf_table(alpha*q);
    lwe_init(&lwe, n, q, alpha);
    if (createSumAndDiffTables(lwe.q)){
        printf("ERROR precomputing Sums and Diffs tables\n");
        return 0;
    }

//                                                    0      1     2    3    4    5    6    7    8    9   10    11    12
    // BRUTEFORCE = 3
   // int start_index[NUM_REDUCTION_STEPS] =            {0,     2,    4,   6,   9,  11,  13,  15,  18,  21,  24,   28,   32};
   // int len_step[NUM_REDUCTION_STEPS] =               {2,     2,    2,   2,   2,   2,   2,   2,   3,   3,   4,    4,    5};
   // int p_step[NUM_REDUCTION_STEPS] =                 {1,     1,    1,   1,   1,   1,   1,   1,   9,  14,  19,   25,   35};
   // int p1_step[NUM_REDUCTION_STEPS] =                {165,  30,    6,   1, 165,  30,   6,   1, 270,  25, 100,  120,  1601};
   // int prev_p1_step[NUM_REDUCTION_STEPS] =           {-1,  165,   30,   6,  -1, 165,  30,   6,  -1, 270,  25,  100,  120};

//                                                    0      1     2    3    4    5    6    7    8      9    10     11     12
   // int start_index[NUM_REDUCTION_STEPS] =            {0,     2,    4,   6,   9,  11,  13,  15,   18,   21,   25,    29,    33};
   // int len_step[NUM_REDUCTION_STEPS] =               {2,     2,    2,   2,   2,   2,   2,   2,    3,    4,    4,     4,     5};
   // int p_step[NUM_REDUCTION_STEPS] =                 {1,     1,    1,   1,   1,   1,   1,   1,    7,   19,   23,    27,    35};
   // int p1_step[NUM_REDUCTION_STEPS] =                {165,  30,    6,   1, 165,  30,   6,   1,  346,  800,  800,   268,  1601};
   // int prev_p1_step[NUM_REDUCTION_STEPS] =           {-1,  165,   30,   6,  -1, 165,  30,   6,   -1,  346,  800,   800,   268};
   // int un_selection[NUM_REDUCTION_STEPS] =           {0,     0,    0,   0,   0,   0,   0,   0,    0,    0,    0,   130,   130};

//                                                     0      1     2     3    4     5     6     7     8     9    10    11    12
    int start_index[NUM_REDUCTION_STEPS] =            {0,     2,    4,    7,   9,   11,   14,   16,   19,   22,   26,   30,   35};
    int len_step[NUM_REDUCTION_STEPS] =               {2,     2,    2,    2,   2,    2,    2,    3,    3,    4,    4,    5,    5};
    int p_step[NUM_REDUCTION_STEPS] =                 {1,     1,    1,    1,   1,    1,    1,    5,    9,   13,   22,   38,   43};
    int p1_step[NUM_REDUCTION_STEPS] =                {90,   10,    1,   90,  10,    1,   75,   81,   14,  140,  115,  800, 1601};
    int prev_p1_step[NUM_REDUCTION_STEPS] =           {-1,   90,   10,   -1,  90,   10,   -1,   75,   81,   14,  140,  115,  800};
    int un_selection[NUM_REDUCTION_STEPS] =           {0,     0,    0,    0,   0,    0,    0,    0,    0,    0,   60,   60,   60};

    bkwStepParameters bkwStepPar[NUM_REDUCTION_STEPS];
    /* Set steps: smooth LMS */
    for (int i=0; i<NUM_REDUCTION_STEPS; i++)
    {
        bkwStepPar[i].startIndex = start_index[i];// i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = len_step[i];
        bkwStepPar[i].p = p_step[i];
        bkwStepPar[i].p1 = p1_step[i];
        bkwStepPar[i].p2 = bkwStepPar[i].p;
        bkwStepPar[i].prev_p1 = prev_p1_step[i];//i == 0 ? -1 : bkwStepPar[i-1].p1;
        bkwStepPar[i].un_selection = un_selection[i];
        ASSERT(bkwStepPar[i].p2 != 0, "smooth-LMS p2 parameter not valid");
        printf("step %02d categories %llu\n", i, num_categories(&lwe, &bkwStepPar[i]));
    }

    // exit(0);

    int bf_positions = BRUTE_FORCE_POSITIONS;
    int zero_positions = ZERO_POSITIONS;
    int fwht_positions = lwe.n - zero_positions - bf_positions;

    u8 binary_solution[fwht_positions];
    short bf_solution[bf_positions];

    time_stamp("LWE instance created");

    samplesList Samples;
    create_lwe_samples(&Samples, &lwe, n_samples);

    time_stamp("Samples allocated");


    time_stamp("Start reduction phase");

    sortedSamplesList sortedSamples1;
    sortedSamplesList sortedSamples2;

    sortedSamplesList *srcSamples, *dstSamples, *tmpSamples;

    /* multiply times 2 mod q and sort (unsorted) samples */
    time_stamp("Multiply samples times 2 modulo q");
    int ret = transition_times2_modq(&lwe, &bkwStepPar[0], &sortedSamples1, &Samples);
    time_stamp("Number of samples: %d - %d samples per category", sortedSamples1.n_samples, sortedSamples1.n_samples_per_category);

    // free original samples - save up memory
    free_samples(&Samples);

    srcSamples = &sortedSamples1;
    dstSamples = &sortedSamples2;

    // perform smooth LMS steps
    int numReductionSteps = NUM_REDUCTION_STEPS;
    for (int i=0; i<numReductionSteps-1; i++){

    	time_stamp("Perform smooth LMS reduction step %d/%d", i+1, numReductionSteps);
        ret = transition_bkw_step_smooth_lms(&lwe, &bkwStepPar[i], &bkwStepPar[i+1], srcSamples, dstSamples);

        // clean past list
        tmpSamples = srcSamples;
        free_sorted_samples(tmpSamples);
        srcSamples = dstSamples;
        dstSamples = tmpSamples;
        time_stamp("Number of samples: %d - %d samples per category", srcSamples->n_samples, srcSamples->n_samples_per_category);
    }

    /* perform last reduction step */
    int i = numReductionSteps-1;
    time_stamp("Perform last smooth LMS reduction step %d/%d", numReductionSteps, numReductionSteps);
    ret = transition_bkw_step_final(&lwe, &bkwStepPar[i], srcSamples, &Samples, samples_for_guessing);

    // clean past list
    tmpSamples = srcSamples;
    free_sorted_samples(tmpSamples);

    freeSumAndDiffTables();

    time_stamp("Number of samples: %d", Samples.n_samples);

    /* compute binary secret */
    u8 original_binary_secret[lwe.n];
    // printf("(");
    for (int i = 0; i < lwe.n; ++i)
    {
        // printf("%d ", lwe.s[i]);
        if (lwe.s[i] < q/2)
            original_binary_secret[i] = lwe.s[i] % 2;
        else
            original_binary_secret[i] = (lwe.s[i]+1) % 2;
    }
    // printf(")\n");

    /* Solving phase - using Fast Walsh Hadamard Tranform */

    time_stamp("Apply Fast Walsh Hadamard Tranform");
    ret = solve_fwht_search(binary_solution, zero_positions, fwht_positions, &Samples, &lwe);
    // ret = solve_fwht_search_bruteforce(binary_solution, bf_solution, zero_positions, bf_positions, fwht_positions, &Samples, &lwe);
    if(ret)
    {
        printf("error %d in solve_fwht_search_hybrid\n", ret);
        exit(-1);
    }

    printf("EXAMPLE\n");
    for (int i = 0; i < 3; ++i)
    {
        printf("(");
        for (int j = 0; j < n; ++j)
        {
            printf("%d  ", Samples.list[i].a[j]);
        }printf(")\n");
    }

    free_samples(&Samples);


    printf("\nFound Solution   \n");
    for(int i = 0; i< fwht_positions; i++)
        printf("%d ", binary_solution[i]);
    printf("- ");
    for(int i = 0; i < bf_positions; i++)
        printf("%d ", bf_solution[i]);

    printf("\nOriginal Solution\n");
    for(int i = zero_positions; i < zero_positions+fwht_positions; i++)
        printf("%d ",original_binary_secret[i]);
    printf("- ");
    for(int i = zero_positions+fwht_positions; i < n; i++)
        printf("%d ", lwe.s[i]);
    printf("\n\n");

    time_stamp("Terminate program");

    return 0;
}
