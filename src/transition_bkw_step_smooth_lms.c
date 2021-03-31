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

#include "transition_bkw_step_smooth_lms.h"

#include <math.h>
#include <string.h>

#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

static u64 subtractSamples(lweInstance *lwe, sample *outSample, sample *sample1, sample *sample2, bkwStepParameters *dstBkwStepPar, bkwStepParameters *srcBkwStepPar)
{
    int q = lwe->q;

    for (int i=0; i < lwe->n; i++)
        outSample->a[i] = (sample1->a[i] - sample2->a[i] + q) %q;//     diffTable(sample1->a[i], sample2->a[i]);
    outSample->z = (sample1->z - sample2->z + q) %q;
    // outSample->error = diffTable(sample1->error, sample2->error);

    u64 index = position_values_2_category_index(lwe, dstBkwStepPar, outSample->a + dstBkwStepPar->startIndex);

    return index;
}

static u64 addSamples(lweInstance *lwe, sample *outSample, sample *sample1, sample *sample2, bkwStepParameters *dstBkwStepPar, bkwStepParameters *srcBkwStepPar)
{
    int n = lwe->n;
    int q = lwe->q;

    for (int i=0; i < n; i++)
        outSample->a[i] = (sample1->a[i] + sample2->a[i]) % q;//sumTable(sample1->a[i], sample2->a[i]);
    outSample->z = (sample1->z + sample2->z) % q;//sumTable(sample1->z, sample2->z);
    // outSample->error = sumTable(sample1->error, sample2->error);

    u64 index = position_values_2_category_index(lwe, dstBkwStepPar, outSample->a + dstBkwStepPar->startIndex);

    return index;
}

int transition_bkw_step_smooth_lms(lweInstance *lwe, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, sortedSamplesList *srcSamples, sortedSamplesList *dstSamples)
{

    set_sorted_samples_list(dstSamples, lwe, dstBkwStepPar, srcSamples->n_samples);

    sample tmpSample;
    tmpSample.a = calloc(lwe->n, sizeof(u16));

    u64 index1, index2, category;
    int n_samples_in_category;

    u64 discarded = 0;

    if (srcSamples->n_categories & 1)
    {
        // process single category
        for (int i = 0; i < srcSamples->list_categories[0].n_samples; i++)
        {
            for (int j=i+1; j < srcSamples->list_categories[0].n_samples; j++)
            {

                category = subtractSamples(lwe, &tmpSample, &srcSamples->list_categories[0].list[i], &srcSamples->list_categories[0].list[j], dstBkwStepPar, srcBkwStepPar);
                
                if (category > dstSamples->n_categories || category < 0)
                {
                    printf("ERROR: category %llu tot categories %llu \n", category, dstSamples->n_categories );
                    exit(0);
                }

                // add it to the new list
                n_samples_in_category = dstSamples->list_categories[category].n_samples;
                if (n_samples_in_category < dstSamples->n_samples_per_category)
                {
                	if (!checkzero((char*)tmpSample.a, sizeof(u16)*lwe->n))
                	{   
	                    memcpy(dstSamples->list_categories[category].list[n_samples_in_category].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list_categories[category].list[n_samples_in_category].z = tmpSample.z;
	                    // dstSamples->list_categories[category].list[n_samples_in_category].error = tmpSample.error;
	                    dstSamples->list_categories[category].n_samples++;
	                    dstSamples->n_samples++;
	                    if (dstSamples->n_samples == dstSamples->max_samples)
		                    goto exit;
	                }
                }
                else
                    discarded++;
            }
        }
        
        index1 = 1;
        index2 = 2;
    }
    else
    {
        index1 = 0;
        index2 = 1;
    }

    /* process samples with LF2 method */
    while (index2 < srcSamples->n_categories && dstSamples->n_samples < dstSamples->max_samples)
    {
        // process single category
        for (int i = 0; i < srcSamples->list_categories[index1].n_samples; i++)
        {
            for (int j=i+1; j < srcSamples->list_categories[index1].n_samples; j++)
            {
                category = subtractSamples(lwe, &tmpSample, &srcSamples->list_categories[index1].list[i], &srcSamples->list_categories[index1].list[j], dstBkwStepPar, srcBkwStepPar);
                
                if (category > dstSamples->n_categories || category < 0)
                {
                    printf("ERROR: category %llu tot categories %llu \n", category, dstSamples->n_categories );
                    exit(0);
                }

                // add it to the new list
                n_samples_in_category = dstSamples->list_categories[category].n_samples;
                if (n_samples_in_category < dstSamples->n_samples_per_category)
                {
                	if (!checkzero((char*)tmpSample.a, sizeof(u16)*lwe->n))
                	{
	                    memcpy(dstSamples->list_categories[category].list[n_samples_in_category].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list_categories[category].list[n_samples_in_category].z = tmpSample.z;
	                    // dstSamples->list_categories[category].list[n_samples_in_category].error = tmpSample.error;
	                    dstSamples->list_categories[category].n_samples++;
	                    dstSamples->n_samples++;

	                    if (dstSamples->n_samples == dstSamples->max_samples)
		                    goto exit;
	                }
                }
                else
                    discarded++;
            }
        }

        // process single category
        for (int i = 0; i < srcSamples->list_categories[index2].n_samples; i++)
        {
            for (int j=i+1; j < srcSamples->list_categories[index2].n_samples; j++)
            {
                category = subtractSamples(lwe, &tmpSample, &srcSamples->list_categories[index2].list[i], &srcSamples->list_categories[index2].list[j], dstBkwStepPar, srcBkwStepPar);

                if (category > dstSamples->n_categories || category < 0)
                {
                    printf("ERROR: category %llu tot categories %llu \n", category, dstSamples->n_categories );
                    exit(0);
                }
         
                // add it to the new list
                n_samples_in_category = dstSamples->list_categories[category].n_samples;
                if (n_samples_in_category < dstSamples->n_samples_per_category)
                {
                	if (!checkzero((char*)tmpSample.a, sizeof(u16)*lwe->n))
                	{
	                    memcpy(dstSamples->list_categories[category].list[n_samples_in_category].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list_categories[category].list[n_samples_in_category].z = tmpSample.z;
	                    // dstSamples->list_categories[category].list[n_samples_in_category].error = tmpSample.error;
	                    dstSamples->list_categories[category].n_samples++;
	                    dstSamples->n_samples++;

	                    if (dstSamples->n_samples == dstSamples->max_samples)
		                    goto exit;	                    
	                }
                }
                else
                    discarded++;
            }
        }

        // process two adjacent category
        for (int i=0; i<srcSamples->list_categories[index1].n_samples; i++)
        {
            for (int j=0; j<srcSamples->list_categories[index2].n_samples; j++)
            {
                category = addSamples(lwe, &tmpSample, &srcSamples->list_categories[index1].list[i], &srcSamples->list_categories[index2].list[j], dstBkwStepPar, srcBkwStepPar);

                if (category > dstSamples->n_categories || category < 0)
                {
                    printf("ERROR: category %llu tot categories %llu \n", category, dstSamples->n_categories );
                    exit(0);
                }
                
                // add it to the new list
	                n_samples_in_category = dstSamples->list_categories[category].n_samples;
                if (n_samples_in_category < dstSamples->n_samples_per_category)
                {
                	if (!checkzero((char*)tmpSample.a, sizeof(u16)*lwe->n))
                	{	                    
	                    memcpy(dstSamples->list_categories[category].list[n_samples_in_category].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list_categories[category].list[n_samples_in_category].z = tmpSample.z;
	                    // dstSamples->list_categories[category].list[n_samples_in_category].error = tmpSample.error;
	                    dstSamples->list_categories[category].n_samples++;
	                    dstSamples->n_samples++;

	                    if (dstSamples->n_samples == dstSamples->max_samples)
		                    goto exit;
	                }
                }
                else
                    discarded++;
            }
        }

        index1 += 2;
        index2 += 2;
    }

exit:

    time_stamp("discarded samples: %ld", discarded);
    free(tmpSample.a);

    return 0;
}

