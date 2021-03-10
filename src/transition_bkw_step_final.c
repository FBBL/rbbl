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

static u64 subtractSamples(lweInstance *lwe, sample *outSample, sample *sample1, sample *sample2, bkwStepParameters *dstBkwStepPar)
{
    int n = lwe->n;
    int q = lwe->q;

    int startIndex = dstBkwStepPar->startIndex;
    int numPositions = dstBkwStepPar->numPositions;

    for (int i=0; i < n; i++)
        outSample->a[i] = diffTable(sample1->a[i], sample2->a[i]);
    outSample->z = diffTable(sample1->z, sample2->z);
    outSample->error = diffTable(sample1->error, sample2->error);

    return 0;
}

static int addSamples(lweInstance *lwe, sample *outSample, sample *sample1, sample *sample2, bkwStepParameters *dstBkwStepPar)
{
    int n = lwe->n;
    int q = lwe->q;

    int startIndex = dstBkwStepPar->startIndex;
    int numPositions = dstBkwStepPar->numPositions;

    for (int i=0; i < n; i++)
        outSample->a[i] = sumTable(sample1->a[i], sample2->a[i]);
    outSample->z = sumTable(sample1->z, sample2->z);
    outSample->error = sumTable(sample1->error, sample2->error);

    return 0;
}

int transition_bkw_step_final(lweInstance *lwe, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, sortedSamplesList *srcSamples, samplesList *dstSamples, int tot_final_samples)
{

    allocate_samples_list(dstSamples, lwe, tot_final_samples); // actually one could have more or less samples

    sample tmpSample;
    tmpSample.a = calloc(lwe->n, sizeof(u16));

    u64 index1, index2, count = 0;

    dstSamples->n_samples = 0;

    if (srcSamples->n_categories & 1){

        // process single category
        for (int i = 0; i < srcSamples->list_categories[0].n_samples; i++)
        {
            for (int j=i+1; j < srcSamples->list_categories[0].n_samples; j++)
            {
                subtractSamples(lwe, &tmpSample, &srcSamples->list_categories[0].list[i], &srcSamples->list_categories[0].list[j], dstBkwStepPar);
                
                if (!checkzero((char*)tmpSample.a, sizeof(u16)*lwe->n))
                {
	                // add it to the new list
	                if (count < tot_final_samples)
	                {
	                    dstSamples->list[count].a = calloc(lwe->n, sizeof(u16));
	                    memcpy(dstSamples->list[count].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list[count].z = tmpSample.z;
	                    dstSamples->list[count].error = tmpSample.error;
	                    dstSamples->n_samples++;
                        count++;
	                }
	                else
	                    goto exit;
	            }
            }
        }

        index1 = 1;
        index2 = 2;
    }
    else{
        index1 = 0;
        index2 = 1;
    }

    /* process samples with LF2 method */
    while (index2 < srcSamples->n_categories && count < tot_final_samples)
    {

        // process single category
        for (int i = 0; i < srcSamples->list_categories[index1].n_samples; i++)
        {
            for (int j=i+1; j < srcSamples->list_categories[index1].n_samples; j++)
            {
                subtractSamples(lwe, &tmpSample, &srcSamples->list_categories[index1].list[i], &srcSamples->list_categories[index1].list[j], dstBkwStepPar);
                
                if (!checkzero((char*)tmpSample.a, sizeof(u16)*lwe->n))
                {
	                // add it to the new list
	                if (count < tot_final_samples)
	                {
	                    dstSamples->list[count].a = malloc(lwe->n*sizeof(u16));
	                    memcpy(dstSamples->list[count].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list[count].z = tmpSample.z;
	                    dstSamples->list[count].error = tmpSample.error;
                        dstSamples->n_samples++;
	                    count++;
	                }
	                else
	                    goto exit;
				}
            }
        }

        // process single category
        for (int i = 0; i < srcSamples->list_categories[index2].n_samples; i++)
        {
            for (int j=i+1; j < srcSamples->list_categories[index2].n_samples; j++)
            {
                subtractSamples(lwe, &tmpSample, &srcSamples->list_categories[index2].list[i], &srcSamples->list_categories[index2].list[j], dstBkwStepPar);
                
                if (!checkzero((char*)tmpSample.a, sizeof(u16)*lwe->n))
                {

	                // add it to the new list
	                if (count < tot_final_samples)
	                {
	                    dstSamples->list[count].a = malloc(lwe->n*sizeof(u16));
	                    memcpy(dstSamples->list[count].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list[count].z = tmpSample.z;
	                    dstSamples->list[count].error = tmpSample.error;
                        dstSamples->n_samples++;
	                    count++;
	                }
	                else
	                    goto exit;
				}
            }
        }

        // process two adjacent category
        for (int i=0; i<srcSamples->list_categories[index1].n_samples; i++)
        {
            for (int j=0; j<srcSamples->list_categories[index2].n_samples; j++)
            {
                addSamples(lwe, &tmpSample, &srcSamples->list_categories[index1].list[i], &srcSamples->list_categories[index2].list[j], dstBkwStepPar);
                
                if (!checkzero((char*)tmpSample.a, sizeof(u16)*lwe->n))
                {
	                // add it to the new list
	                if (count < tot_final_samples)
	                {
	                    dstSamples->list[count].a = malloc(lwe->n*sizeof(u16));
	                    memcpy(dstSamples->list[count].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list[count].z = tmpSample.z;
	                    dstSamples->list[count].error = tmpSample.error;
                        dstSamples->n_samples++;
	                    count++;
	                }
	                else
	                    goto exit;
				}
            }
        }

        index1 += 2;
        index2 += 2;
    }

exit:

    free(tmpSample.a);

    return 0;
}

