/*  This file is part of FBBL (File-Based BKW for LWE).
 *
 *  FBBL is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  FBBL is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Nome-Programma.  If not, see <http://www.gnu.org/licenses/>
 */

#include "transition_times2_modq.h"

#include <string.h>

/* perform multiplication times 2 mod q of srcSample and store in dstSample. return the category */
int sample_times2_modq(sample *dstSample, sample *srcSample, lweInstance *lwe, bkwStepParameters *bkwStepPar)
{
    int n = lwe->n;
    int q = lwe->q;
    for (int i=0; i<n; i++)
        dstSample->a[i] = (2*srcSample->a[i]) % q;
    dstSample->z = (2*srcSample->z) % q;
    dstSample->error = (2*srcSample->error) % q;

    int index = position_values_2_category_index(lwe, bkwStepPar, dstSample->a);

    return index;
}

/* multiply each sample times 2 mod q in sortedSamples. Then store the result in dstSortedSamplesList according to its category */
int transition_times2_modq(lweInstance *lwe, bkwStepParameters *bkwStepPar, sortedSamplesList *sortedSamples, samplesList* unsortedSamples)
{

    allocate_sorted_samples_list(sortedSamples, lwe, bkwStepPar, unsortedSamples->n_samples);

    sample tmpSample;
    tmpSample.a = calloc(lwe->n, sizeof(u16));

    int count = 0;
    int category;
    int n_samples_in_category;
    while(count < unsortedSamples->n_samples){
    
        category = sample_times2_modq(&tmpSample, &unsortedSamples->list[count], lwe, bkwStepPar);
        n_samples_in_category = sortedSamples->list_categories[category].n_samples;

        if (n_samples_in_category < sortedSamples->n_samples_per_category && sortedSamples->n_samples < sortedSamples->max_samples)
        {
            if (category > sortedSamples->n_categories)
            {
                printf("ERROR: category %d tot categories %d \n", category, sortedSamples->n_categories );
                exit(0);
            }
            
            sortedSamples->list_categories[category].list[n_samples_in_category].a = calloc(lwe->n, sizeof(u16));
            memcpy(sortedSamples->list_categories[category].list[n_samples_in_category].a, tmpSample.a, lwe->n*sizeof(u16));
            sortedSamples->list_categories[category].list[n_samples_in_category].z = tmpSample.z;
            sortedSamples->list_categories[category].list[n_samples_in_category].error = tmpSample.error;
            sortedSamples->list_categories[category].n_samples++;
            sortedSamples->n_samples++;
        }
        count++;
    }

    // for (int i = 0; i < sortedSamples->n_categories; ++i)
    //     printf("%d ", sortedSamples->list_categories[i].n_samples);
    // printf("\n");

    free(tmpSample.a);

    return 1;

}

