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

#define _POSIX_C_SOURCE 199309L

#include "config.h"
#include "lwe_instance.h"
#include "utils.h"
#include "transition_times2_modq.h"
#include "transition_bkw_step_smooth_lms.h"
#include "transition_bkw_step_final.h"
#include "solve_fwht.h"
#include "error_rate.h"

#include <stdio.h>
#include <stdlib.h>

#define NUM_REDUCTION_STEPS 13
#define BRUTE_FORCE_POSITIONS 1
#define ZERO_POSITIONS 17

int main()
{
    u64 n_samples = 18000000;
    u64 samples_for_guessing = 3000000;

    lweInstance lwe;
    int n = 40;
    int q = 1601;
    double alpha = 0.005;

    time_stamp("LWE parameters: n: %d, q: %d, sigma: %lf*q. Initial samples: %llu, n_cores: %d", n, q, alpha, n_samples, (int)NUM_THREADS);

    // initialize random
    // time_t start = time(NULL);
    srand(time(NULL));
    randomUtilRandomize();

    time_stamp("Create LWE instance");
    lwe_init(&lwe, n, q, alpha);

//   n = 40, alpha = 0.005, 8 smoothplainBKW + 5 smoothLMS      0      1     2    3    4     5     6    7    8    9   10    11    12
    int start_index[NUM_REDUCTION_STEPS] =                     {0,     2,    4,   6,   8,   11,   13,  15,  18,  21,  25,   29,   34};
    int len_step[NUM_REDUCTION_STEPS] =                        {2,     2,    2,   2,   2,    2,    2,   2,   3,   4,   4,    4,    5};
    int p_step[NUM_REDUCTION_STEPS] =                          {1,     1,    1,   1,   1,    1,    1,   1,  13,  24,  33,   48,   66};
    int p1_step[NUM_REDUCTION_STEPS] =                         {280,  80,   20,   5,   1,  178,   45,  12, 180, 479, 336,   48, 1601};
    int prev_p1_step[NUM_REDUCTION_STEPS] =                    {-1,  280,   80,  20,   5,   -1,  178,  45,  -1, 180, 479,  336, 1601};
    int un_selection[NUM_REDUCTION_STEPS] =                    {0,     0,    0,   0,   0,    0,    0,   0,  16,  23,  33,   40,   40};

    float increase_factor[NUM_REDUCTION_STEPS] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    bkwStepParameters bkwStepPar[NUM_REDUCTION_STEPS];
    /* Set steps: smooth LMS */
    u64 max_categories = 0, tmp_categories[NUM_REDUCTION_STEPS];

    for (int i=0; i<NUM_REDUCTION_STEPS; i++)
    {
        bkwStepPar[i].startIndex = start_index[i];// i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = len_step[i];//2;
        bkwStepPar[i].p = p_step[i];//3; // test
        bkwStepPar[i].p1 = p1_step[i]; //19; // test
        bkwStepPar[i].p2 = bkwStepPar[i].p;
        bkwStepPar[i].prev_p1 = prev_p1_step[i];// i ==  0 ? -1 : bkwStepPar[i-1].p1;
        bkwStepPar[i].un_selection = un_selection[i];
        ASSERT(bkwStepPar[i].p2 != 0, "smooth-LMS p2 parameter not valid");
        tmp_categories[i] = num_categories(&lwe, &bkwStepPar[i]);
        printf("step %d categories %lu\n", i, tmp_categories[i]);
        if (tmp_categories[i] > max_categories)
            max_categories = tmp_categories[i];
    }

    int bf_positions = BRUTE_FORCE_POSITIONS;
    int zero_positions = ZERO_POSITIONS;
    int fwht_positions = lwe.n - zero_positions - bf_positions;

    u8 binary_solution[fwht_positions];
    short bf_solution[bf_positions];

    time_stamp("Generate %llu samples", n_samples);
    unsortedSamplesList Samples;
    create_lwe_samples(&Samples, &lwe, n_samples);

    time_stamp("Start reduction phase");

    sortedSamplesList sortedSamples1;
    sortedSamplesList sortedSamples2;

    sortedSamplesList *srcSamples, *dstSamples, *tmpSamples;

    allocate_sorted_samples_list(&sortedSamples1, &lwe, &bkwStepPar[0], tmp_categories[0]*increase_factor[0], max_categories);
    set_sorted_samples_list(&sortedSamples1, &lwe, &bkwStepPar[0], tmp_categories[0]*increase_factor[0], max_categories);

    /* multiply times 2 mod q and sort (unsorted) samples */
    time_stamp("Multiply samples times 2 modulo q");
    int ret = transition_times2_modq(&lwe, &bkwStepPar[0], &sortedSamples1, &Samples);
    time_stamp("Number of samples: %d", sortedSamples1.n_samples);

    // free original samples - save up memory
    free_samples(&Samples);

    srcSamples = &sortedSamples1;
    dstSamples = &sortedSamples2;

    allocate_sorted_samples_list(dstSamples, &lwe, &bkwStepPar[1], tmp_categories[0]*increase_factor[0], max_categories);

    struct timespec begin, end;
    clock_gettime(CLOCK_REALTIME, &begin);

    // perform smooth LMS steps
    int numReductionSteps = NUM_REDUCTION_STEPS;
    for (int i=0; i<numReductionSteps-1; i++)
    {

        time_stamp("Perform smooth LMS reduction step %d/%d", i+1, numReductionSteps);
        set_sorted_samples_list(dstSamples, &lwe, &bkwStepPar[i+1], tmp_categories[i+1]*increase_factor[i+1], max_categories);

        ret = transition_bkw_step_smooth_lms(&lwe, &bkwStepPar[i+1], srcSamples, dstSamples);

        if (i != numReductionSteps-2)
        {
            // clean past list
            tmpSamples = srcSamples;
            clean_sorted_samples(tmpSamples);
            srcSamples = dstSamples;
            dstSamples = tmpSamples;
        }
        else
        {
            free_sorted_samples(srcSamples, max_categories);
            srcSamples = dstSamples;
        }

        time_stamp("Number of samples: %d", srcSamples->n_samples);
    }

    /* perform last reduction step */
    int i = numReductionSteps-1;
    time_stamp("Perform last smooth LMS reduction step %d/%d", numReductionSteps, numReductionSteps);

    allocate_unsorted_samples_list(&Samples, &lwe, samples_for_guessing); // actually one could have more or less samples
    ret = transition_bkw_step_final(&lwe, &bkwStepPar[i], srcSamples, &Samples, samples_for_guessing);

    time_stamp("Number of samples: %d", Samples.n_samples);

    // clean past list
    free_sorted_samples(srcSamples, max_categories);

    /* compute binary secret */
    u8 original_binary_secret[lwe.n];
    for (int i = 0; i < lwe.n; ++i)
    {
        if (lwe.s[i] < q/2)
            original_binary_secret[i] = lwe.s[i] % 2;
        else
            original_binary_secret[i] = (lwe.s[i]+1) % 2;
    }

    if (Samples.n_samples == 0)
    {
        time_stamp("ERROR: not enough samples. Terminate the program.");
        free_samples(&Samples);
        exit(1);
    }

    /* Solving phase - using Fast Walsh Hadamard Tranform */
    time_stamp("Apply Fast Walsh Hadamard Tranform");
    ret = solve_fwht_search_bruteforce(binary_solution, bf_solution, zero_positions, bf_positions, fwht_positions, &Samples, &lwe);
    if(ret)
    {
        printf("error %d in solve_fwht_search_hybrid\n", ret);
        exit(-1);
    }
    free_samples(&Samples);

    clock_gettime(CLOCK_REALTIME, &end);
    long seconds = end.tv_sec - begin.tv_sec;
    long nanoseconds = end.tv_nsec - begin.tv_nsec;
    double elapsed = seconds + nanoseconds*1e-9;

    printf("\nFound Solution   \n");
    for(int i = 0; i<fwht_positions; i++)
        printf("%d ",binary_solution[i]);
    for(int i = 0; i<bf_positions; i++)
        printf("%d ",bf_solution[i]);
    printf("\n");

    printf("\nOriginal Solution\n");
    for(int i = zero_positions; i<zero_positions+fwht_positions; i++)
        printf("%d ",original_binary_secret[i]);
    for(int i = zero_positions+fwht_positions; i<lwe.n; i++)
        printf("%d ",lwe.s[i]);
    printf("\n");

    printf("Time measured: %.3f seconds.\n", elapsed);

    for (int i = 0; i < fwht_positions; ++i)
    {
        if (binary_solution[i] != original_binary_secret[zero_positions+i])
        {
            printf("Test Failed");
            return 1;
        }
    }
    for (int i = 0; i<bf_positions; i++)
        if (bf_solution[i] != lwe.s[zero_positions+fwht_positions+i])
        {
            printf("Test Failed");
            return 1;
        }

    printf("Test Passed");

    return 0;
}
