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

#define NUM_REDUCTION_STEPS 5
#define BRUTE_FORCE_POSITIONS 0

int main()
{
    u64 n_samples = 10000;

    lweInstance lwe;
    int n = 10;
    int q = 101;
    double alpha = 0.005;

    time_stamp("Precomputation");
    precompute_cdf_table(alpha*q);
    lwe_init(&lwe, n, q, alpha);
    if (createSumAndDiffTables(lwe.q)){
        printf("ERROR precomputing Sums and Diffs tables\n");
        return 0;
    }


    time_stamp("LWE instance created");

    samplesList Samples;
    create_lwe_samples(&Samples, &lwe, n_samples);

    time_stamp("Samples allocated");

    bkwStepParameters bkwStepPar[NUM_REDUCTION_STEPS];

    /* Set steps: smooth LMS */
    for (int i=0; i<NUM_REDUCTION_STEPS; i++)
    {
        bkwStepPar[i].startIndex = i == 0 ? 0 : bkwStepPar[i-1].startIndex + bkwStepPar[i-1].numPositions;
        bkwStepPar[i].numPositions = 2;
        bkwStepPar[i].p = 11; // test
        bkwStepPar[i].p1 = 19; // test
        bkwStepPar[i].p2 = bkwStepPar[i].p;
        bkwStepPar[i].prev_p1 = i == 0 ? -1 : bkwStepPar[i-1].p1;
        ASSERT(bkwStepPar[i].p2 != 0, "smooth-LMS p2 parameter not valid");
    }

    int bruteForcePositions = BRUTE_FORCE_POSITIONS;
    int fwht_positions = lwe.n - bruteForcePositions;
    int zero_positions = 0;

    u8 binary_solution[fwht_positions];
    short bf_solution[bruteForcePositions];

    time_stamp("Start reduction phase");

    sortedSamplesList sortedSamples1;
    sortedSamplesList sortedSamples2;

    sortedSamplesList *srcSamples, *dstSamples, *tmpSamples;

    /* multiply times 2 mod q and sort (unsorted) samples */
    time_stamp("Multiply samples times 2 modulo q");
    int ret = transition_times2_modq(&lwe, &bkwStepPar[0], &sortedSamples1, &Samples);

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
        time_stamp("Number of samples: %d", srcSamples->n_samples);
    }

    /* perform last reduction step */
    int i = numReductionSteps-1;
    time_stamp("Perform last smooth LMS reduction step %d/%d", numReductionSteps, numReductionSteps);
    ret = transition_bkw_step_final(&lwe, &bkwStepPar[i], &bkwStepPar[i+1], srcSamples, &Samples, srcSamples->n_samples);

    time_stamp("Number of samples: %d", Samples.n_samples);

    // clean past list
    tmpSamples = srcSamples;
    free_sorted_samples(tmpSamples);

    /* compute binary secret */
    u8 real_binary_secret[lwe.n];
    // printf("(");
    for (int i = 0; i < lwe.n; ++i)
    {
        // printf("%d ", lwe.s[i]);
        if (lwe.s[i] < q/2)
            real_binary_secret[i] = lwe.s[i] % 2;
        else
            real_binary_secret[i] = (lwe.s[i]+1) % 2;
    }
    // printf(")\n");

    freeSumAndDiffTables();

    /* Solving phase - using Fast Walsh Hadamard Tranform */

    time_stamp("Apply Fast Walsh Hadamard Tranform");
    ret = solve_fwht_search(binary_solution, zero_positions, fwht_positions, &Samples, &lwe);
    if(ret)
    {
        printf("error %d in solve_fwht_search_hybrid\n", ret);
        exit(-1);
    }
    free_samples(&Samples);


    printf("Binary Solution Found (");
    for(int i = 0; i<lwe.n-bruteForcePositions; i++)
        printf("%hhu ",binary_solution[i]);
    printf("- ");
    for(int i = 0; i<bruteForcePositions; i++)
        printf("%hi ",bf_solution[i]);
    printf(")\n");

    printf("Real Binary Solution  (");
    for(int i = 0; i<lwe.n-bruteForcePositions; i++)
        printf("%hhu ",real_binary_secret[i]);
    printf("- ");
    for(int i = 0; i<bruteForcePositions; i++)
        printf("%hi ",lwe.s[i+zero_positions+fwht_positions]);
    printf(")\n");

    for(int i = 0; i<lwe.n-bruteForcePositions; i++)
    {
        if (binary_solution[i] != real_binary_secret[i])
        {
            printf("WRONG retrieved solution!\n");
            return 1;
        }
    }
    for(int i = 0; i<bruteForcePositions; i++)
    {
        if (bf_solution[i] != lwe.s[i+zero_positions+fwht_positions])
        {
            printf("WRONG retrieved solution!\n");
            return 1;
        }
    }

    time_stamp("Terminate program");

    return 0;
}
