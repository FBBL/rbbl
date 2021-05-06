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
#include <pthread.h>

#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

/* define mutexes to protect common resources from concurrent access */
static pthread_mutex_t screen_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t tot_count_mutex = PTHREAD_MUTEX_INITIALIZER;
static u64 discarded = 0;

typedef struct {
    lweInstance *lwe;
    bkwStepParameters *bkwStepPar;
    sortedSamplesList *srcSamples;
    samplesList* dstSamples;
    u64 minIndex1;
    u64 maxIndex1;
    u64 maxTotSamples;
} Params;

static u64 subtractSamples(lweInstance *lwe, sample *outSample, sample *sample1, sample *sample2)
{
    int q = lwe->q;

    for (int i=0; i < lwe->n; i++)
        outSample->a[i] = (sample1->a[i] - sample2->a[i] + q) %q;//     diffTable(sample1->a[i], sample2->a[i]);
    outSample->z = (sample1->z - sample2->z + q) %q;
    // outSample->error = diffTable(sample1->error, sample2->error);

    return 0;
}

static int addSamples(lweInstance *lwe, sample *outSample, sample *sample1, sample *sample2)
{
    int q = lwe->q;

    for (int i=0; i < lwe->n; i++)
        outSample->a[i] = (sample1->a[i] + sample2->a[i]) % q;//sumTable(sample1->a[i], sample2->a[i]);
    outSample->z = (sample1->z + sample2->z) % q;//sumTable(sample1->z, sample2->z);
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

void *single_thread_final_lf2_work(void *params){

    Params *p = (Params*)params;

    sample tmpSample;
    u16 sample_a[p->lwe->n];
    tmpSample.a = sample_a;

    u64 index1, index2, category, sample_index;
    int n_samples_in_category, mutex_index;

    /* process samples with LF2 method */
    for (index1 = p->minIndex1; index1 < p->maxIndex1; index1 = index1+2)
    {
        index2 = index1+1;

        if (!(p->dstSamples->n_samples < p->maxTotSamples)){
            return NULL;
        }

        // process single category
        for (int i = 0; i < p->srcSamples->list_categories[index1].n_samples; i++)
        {
            for (int j=i+1; j < p->srcSamples->list_categories[index1].n_samples; j++)
            {
                subtractSamples(p->lwe, &tmpSample, &p->srcSamples->list_categories[index1].list[i], &p->srcSamples->list_categories[index1].list[j]);

                // add it to the new list
            	if (!checkzero((char*)tmpSample.a, sizeof(u16)*p->lwe->n) && unnaturalSelection(p->lwe, &tmpSample, p->bkwStepPar))
                {
                	pthread_mutex_lock(&tot_count_mutex);
                    sample_index = p->dstSamples->n_samples;
                    if (!(p->dstSamples->n_samples < p->maxTotSamples))
                    {
                    	pthread_mutex_unlock(&tot_count_mutex);
	            		return NULL;
                    }
                    p->dstSamples->n_samples++;
                    pthread_mutex_unlock(&tot_count_mutex);                    
                    p->dstSamples->list[sample_index].a = malloc(p->lwe->n*sizeof(u16));
                    memcpy(p->dstSamples->list[sample_index].a, tmpSample.a, p->lwe->n*sizeof(u16));
                    p->dstSamples->list[sample_index].z = tmpSample.z;
                    // dstSamples->list[sample_index].error = tmpSample.error;
                }
                else
                {
                  	discarded++;
                }
            }
        }

        // process single category
        for (int i = 0; i < p->srcSamples->list_categories[index2].n_samples; i++)
        {
            for (int j=i+1; j < p->srcSamples->list_categories[index2].n_samples; j++)
            {
                subtractSamples(p->lwe, &tmpSample, &p->srcSamples->list_categories[index2].list[i], &p->srcSamples->list_categories[index2].list[j]);
                
                // add it to the new list
            	if (!checkzero((char*)tmpSample.a, sizeof(u16)*p->lwe->n) && unnaturalSelection(p->lwe, &tmpSample, p->bkwStepPar))
                {
                	pthread_mutex_lock(&tot_count_mutex);
                    sample_index = p->dstSamples->n_samples;
                    if (!(p->dstSamples->n_samples < p->maxTotSamples))
                    {
                    	pthread_mutex_unlock(&tot_count_mutex);
	            		return NULL;
                    }
                    p->dstSamples->n_samples++;
                    pthread_mutex_unlock(&tot_count_mutex);                   
                    p->dstSamples->list[sample_index].a = malloc(p->lwe->n*sizeof(u16));
                    memcpy(p->dstSamples->list[sample_index].a, tmpSample.a, p->lwe->n*sizeof(u16));
                    p->dstSamples->list[sample_index].z = tmpSample.z;
                    // dstSamples->list[sample_index].error = tmpSample.error;
                }
                else
                {
                    discarded++;
                }
            }
        }

        // process two adjacent category
        for (int i=0; i<p->srcSamples->list_categories[index1].n_samples; i++)
        {
            for (int j=0; j<p->srcSamples->list_categories[index2].n_samples; j++)
            {
                addSamples(p->lwe, &tmpSample, &p->srcSamples->list_categories[index1].list[i], &p->srcSamples->list_categories[index2].list[j]);
                
                // add it to the new list
            	if (!checkzero((char*)tmpSample.a, sizeof(u16)*p->lwe->n) && unnaturalSelection(p->lwe, &tmpSample, p->bkwStepPar))
                {
                	pthread_mutex_lock(&tot_count_mutex);
                    sample_index = p->dstSamples->n_samples;
                    if (!(p->dstSamples->n_samples < p->maxTotSamples))
                    {
                    	pthread_mutex_unlock(&tot_count_mutex);
	            		return NULL;
                    }
                    p->dstSamples->n_samples++;
                    pthread_mutex_unlock(&tot_count_mutex);
                    p->dstSamples->list[sample_index].a = malloc(p->lwe->n*sizeof(u16));
                    memcpy(p->dstSamples->list[sample_index].a, tmpSample.a, p->lwe->n*sizeof(u16));
                    p->dstSamples->list[sample_index].z = tmpSample.z;
                    // dstSamples->list[sample_index].error = tmpSample.error; 
                }
                else
                {
                  	discarded++;
                }
            }
        }
    }
}

/*
    VERY IMPORTANT: p->maxTotSamples must be < than the expected number of samples that can be generated!
*/

int transition_bkw_step_final(lweInstance *lwe, bkwStepParameters *srcBkwStepPar, sortedSamplesList *srcSamples, samplesList *dstSamples, u64 maxSamples)
{

    ASSERT(NUM_THREADS >= 1, "Unexpected number of threads!");

    sample tmpSample;
    u16 sample_a[lwe->n];
    tmpSample.a = sample_a;

    u64 index1, index2, minc;
    dstSamples->n_samples = 0;
    discarded = 0;

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
	                if (dstSamples->n_samples < maxSamples)
	                {
	                    dstSamples->list[dstSamples->n_samples].a = calloc(lwe->n, sizeof(u16));
	                    memcpy(dstSamples->list[dstSamples->n_samples].a, tmpSample.a, lwe->n*sizeof(u16));
	                    dstSamples->list[dstSamples->n_samples].z = tmpSample.z;
	                    // dstSamples->list[dstSamples->n_samples].error = tmpSample.error;
	                    dstSamples->n_samples++;
	                }
	                // else
	                //     goto exit;
	            }
	            else
	            	discarded++;
            }
        }

        minc = 1;
    }
    else{
        minc = 0;
    }

    pthread_t thread[NUM_THREADS];
    Params param[NUM_THREADS]; /* one set of in-/output paramaters per thread, so no need to lock these */

    u64 cat_per_thread = srcSamples->n_categories/NUM_THREADS;
    if(cat_per_thread & 1) // make it even
        cat_per_thread--;

    /* load input parameters */
    for (int i=0; i<NUM_THREADS; i++) {
        param[i].lwe = lwe; /* set input parameter to thread number */
        param[i].bkwStepPar = srcBkwStepPar;
        param[i].srcSamples = srcSamples;
        param[i].dstSamples = dstSamples;
        param[i].minIndex1 = i == 0 ? minc : param[i-1].minIndex1 + cat_per_thread;
        param[i].maxIndex1 = param[i].minIndex1 + cat_per_thread;
        param[i].maxTotSamples = maxSamples;
    }
    param[NUM_THREADS-1].maxIndex1 = srcSamples->n_categories-2;

    /* process samples with LF2 method */
    /* start threads */
    for (int i = 0; i < NUM_THREADS; ++i)
    {
        if (!pthread_create(&thread[i], NULL, single_thread_final_lf2_work, (void*)&param[i])) {
            // pthread_mutex_lock(&screen_mutex);
            // printf("Thread %d created!\n", i+1);
            // pthread_mutex_unlock(&screen_mutex);
        } else {
            // pthread_mutex_lock(&screen_mutex);
            // printf("Error creating thread %d!\n", i+1);
            // pthread_mutex_unlock(&screen_mutex);
        }
    }

    /* wait until all threads have completed */
    for (int i = 0; i < NUM_THREADS; i++) {
        if (!pthread_join(thread[i], NULL)) {
            // pthread_mutex_lock(&screen_mutex);
            // printf("Thread %d joined!\n", i+1);
            // pthread_mutex_unlock(&screen_mutex);
        } else {
            // pthread_mutex_lock(&screen_mutex);
            // printf("Error joining thread %d!\n", i+1);
            // pthread_mutex_unlock(&screen_mutex);
        }
    }

    time_stamp("discarded samples: %ld", discarded);

    return 0;
}

