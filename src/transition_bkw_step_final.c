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

static u64 subtractSamples(lweInstance *lwe, sample *outSample, sample *sample1, sample *sample2)
{
    int n = lwe->n;
    int q = lwe->q;

    for (int i=0; i < n; i++)
        outSample->a[i] = diffTable(sample1->a[i], sample2->a[i]);
    outSample->z = diffTable(sample1->z, sample2->z);
    // outSample->error = diffTable(sample1->error, sample2->error);

    return 0;
}

static int addSamples(lweInstance *lwe, sample *outSample, sample *sample1, sample *sample2)
{
    int n = lwe->n;
    int q = lwe->q;

    for (int i=0; i < n; i++)
        outSample->a[i] = sumTable(sample1->a[i], sample2->a[i]);
    outSample->z = sumTable(sample1->z, sample2->z);
    // outSample->error = sumTable(sample1->error, sample2->error);

    return 0;
}

/* Return 1 if passes unnatural selection, 0 otherwise */
static int unnaturalSelection(lweInstance *lwe, sample *Sample, bkwStepParameters *srcBkwStepPar)
{
    int n = lwe->n;
    int q = lwe->q;

    int tmp;

    for (int i=0; i < srcBkwStepPar->startIndex + srcBkwStepPar->numPositions; i++)
    {
        tmp = Sample->a[i] < q/2 ? Sample->a[i] : (int)Sample->a[i] -q;
        if (tmp > srcBkwStepPar->un_selection || tmp < -srcBkwStepPar->un_selection)
            return 0;
    }

    return 1;
}

int transition_bkw_step_final(lweInstance *lwe, bkwStepParameters *srcBkwStepPar, sortedSamplesList *srcSamples, samplesList *dstSamples, int tot_final_samples)
{

    allocate_samples_list(dstSamples, lwe, tot_final_samples); // actually one could have more or less samples

    sample tmpSample;
    tmpSample.a = calloc(lwe->n, sizeof(u16));

    u64 index1, index2, discarded = 0;

    dstSamples->n_samples = 0;

    if (srcSamples->n_categories & 1){

        // process single category
        for (int i = 0; i < srcSamples->list_categories[0].n_samples; i++)
        {
            for (int j=i+1; j < srcSamples->list_categories[0].n_samples; j++)
            {
                subtractSamples(lwe, &tmpSample, &srcSamples->list_categories[0].list[i], &srcSamples->list_categories[0].list[j]);
                
                if (!checkzero((char*)tmpSample.a, sizeof(u16)*lwe->n) && unnaturalSelection(lwe, &tmpSample, srcBkwStepPar))
                {
	                // add it to the new list
	                if (dstSamples->n_samples < tot_final_samples)
	                {
	                    dstSamples->list[dstSamples->n_samples].a = calloc(lwe->n, sizeof(u16));
	                    memcpy(dstSamples->list[dstSamples->n_samples].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list[dstSamples->n_samples].z = tmpSample.z;
	                    // dstSamples->list[dstSamples->n_samples].error = tmpSample.error;
	                    dstSamples->n_samples++;
	                }
	                else
	                    goto exit;
	            }
	            else
	            	discarded++;
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
    while (index2 < srcSamples->n_categories && dstSamples->n_samples < tot_final_samples)
    {

        // process single category
        for (int i = 0; i < srcSamples->list_categories[index1].n_samples; i++)
        {
            for (int j=i+1; j < srcSamples->list_categories[index1].n_samples; j++)
            {
                subtractSamples(lwe, &tmpSample, &srcSamples->list_categories[index1].list[i], &srcSamples->list_categories[index1].list[j]);
                
                if (!checkzero((char*)tmpSample.a, sizeof(u16)*lwe->n) && unnaturalSelection(lwe, &tmpSample, srcBkwStepPar))
                {
	                // add it to the new list
	                if (dstSamples->n_samples < tot_final_samples)
	                {
	                    dstSamples->list[dstSamples->n_samples].a = malloc(lwe->n*sizeof(u16));
	                    memcpy(dstSamples->list[dstSamples->n_samples].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list[dstSamples->n_samples].z = tmpSample.z;
	                    // dstSamples->list[dstSamples->n_samples].error = tmpSample.error;
                        dstSamples->n_samples++;
	                }
	                else
	                    goto exit;
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
                subtractSamples(lwe, &tmpSample, &srcSamples->list_categories[index2].list[i], &srcSamples->list_categories[index2].list[j]);
                
                if (!checkzero((char*)tmpSample.a, sizeof(u16)*lwe->n) && unnaturalSelection(lwe, &tmpSample, srcBkwStepPar))
                {

	                // add it to the new list
	                if (dstSamples->n_samples < tot_final_samples)
	                {
	                    dstSamples->list[dstSamples->n_samples].a = malloc(lwe->n*sizeof(u16));
	                    memcpy(dstSamples->list[dstSamples->n_samples].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list[dstSamples->n_samples].z = tmpSample.z;
	                    // dstSamples->list[dstSamples->n_samples].error = tmpSample.error;
                        dstSamples->n_samples++;
	                }
	                else
	                    goto exit;
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
                addSamples(lwe, &tmpSample, &srcSamples->list_categories[index1].list[i], &srcSamples->list_categories[index2].list[j]);
                
                if (!checkzero((char*)tmpSample.a, sizeof(u16)*lwe->n) && unnaturalSelection(lwe, &tmpSample, srcBkwStepPar))
                {
	                // add it to the new list
	                if (dstSamples->n_samples < tot_final_samples)
	                {
	                    dstSamples->list[dstSamples->n_samples].a = malloc(lwe->n*sizeof(u16));
	                    memcpy(dstSamples->list[dstSamples->n_samples].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list[dstSamples->n_samples].z = tmpSample.z;
	                    // dstSamples->list[dstSamples->n_samples].error = tmpSample.error;
                        dstSamples->n_samples++;
	                }
	                else
	                    goto exit;
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

