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
static pthread_mutex_t *storage_mutex;
static int n_storage_mutex;
static u64 discarded = 0;

typedef struct {
    lweInstance *lwe;
    bkwStepParameters *bkwStepPar;
    sortedSamplesList *srcSamples;
    sortedSamplesList* dstSamples;
    u64 minIndex1;
    u64 maxIndex1;
} Params;

static u64 subtractSamples(lweInstance *lwe, sample *outSample, sample *sample1, sample *sample2, bkwStepParameters *bkwStepPar)
{
    int q = lwe->q;

    for (int i=0; i < lwe->n; i++)
        outSample->a[i] = (sample1->a[i] - sample2->a[i] + q) %q;//     diffTable(sample1->a[i], sample2->a[i]);
    outSample->z = (sample1->z - sample2->z + q) %q;
    // outSample->error = diffTable(sample1->error, sample2->error);

    u64 index = position_values_2_category_index(lwe, bkwStepPar, outSample->a + bkwStepPar->startIndex);

    return index;
}

static u64 addSamples(lweInstance *lwe, sample *outSample, sample *sample1, sample *sample2, bkwStepParameters *bkwStepPar)
{
    int n = lwe->n;
    int q = lwe->q;

    for (int i=0; i < n; i++)
        outSample->a[i] = (sample1->a[i] + sample2->a[i]) % q;//sumTable(sample1->a[i], sample2->a[i]);
    outSample->z = (sample1->z + sample2->z) % q;//sumTable(sample1->z, sample2->z);
    // outSample->error = sumTable(sample1->error, sample2->error);

    u64 index = position_values_2_category_index(lwe, bkwStepPar, outSample->a + bkwStepPar->startIndex);

    return index;
}


void *single_thread_lf2_work(void *params){

	Params *p = (Params*)params;

	u64 index1, index2, category;
	int n_samples_in_category, mutex_index;

    // pthread_mutex_lock(&screen_mutex);
    // printf("Min %lu Max %lu\n", p->minIndex1, p->maxIndex1);
    // pthread_mutex_unlock(&screen_mutex);

    sample tmpSample;
    u16 sample_a[p->lwe->n];
    tmpSample.a = sample_a;

	for (index1 = p->minIndex1; index1 < p->maxIndex1; index1 = index1+2)
	{
		index2 = index1+1;

        if (!(p->dstSamples->n_samples < p->dstSamples->max_samples)){        	
            return NULL;
        }

	    // process single category
	    for (int i = 0; i < p->srcSamples->list_categories[index1].n_samples; i++)
	    {
	        for (int j=i+1; j < p->srcSamples->list_categories[index1].n_samples; j++)
	        {
	            category = subtractSamples(p->lwe, &tmpSample, &p->srcSamples->list_categories[index1].list[i], &p->srcSamples->list_categories[index1].list[j], p->bkwStepPar);
	            mutex_index = category / n_storage_mutex;

	            if (category > p->dstSamples->n_categories || category < 0)
	            {
	            	pthread_mutex_lock(&screen_mutex);
	                printf("ERROR subtractSamples: category %lu tot categories %lu \n", category, p->dstSamples->n_categories );
	                pthread_mutex_unlock(&screen_mutex);
	                exit(0);
	            }

	            pthread_mutex_lock(&storage_mutex[mutex_index]);
	            // add it to the new list
	            n_samples_in_category = p->dstSamples->list_categories[category].n_samples;
	            if (n_samples_in_category < p->dstSamples->n_samples_per_category)
	            {
	            	if (!checkzero((char*)tmpSample.a, sizeof(u16)*p->lwe->n))
	            	{
	                    memcpy(p->dstSamples->list_categories[category].list[n_samples_in_category].a, tmpSample.a, p->lwe->n*sizeof(u16));
	                    p->dstSamples->list_categories[category].list[n_samples_in_category].z = tmpSample.z;
	                    // dstSamples->list_categories[category].list[n_samples_in_category].error = tmpSample.error;
	                    p->dstSamples->list_categories[category].n_samples++;
	                    p->dstSamples->n_samples++;

	                    if (p->dstSamples->n_samples >= p->dstSamples->max_samples){
	                    	pthread_mutex_unlock(&storage_mutex[mutex_index]);
		                    return NULL;
	                    }
	                }
	            }
	            else
	                discarded++;
	            pthread_mutex_unlock(&storage_mutex[mutex_index]);
	        }
	    }

	    // process single category
	    for (int i = 0; i < p->srcSamples->list_categories[index2].n_samples; i++)
	    {
	        for (int j=i+1; j < p->srcSamples->list_categories[index2].n_samples; j++)
	        {
	            category = subtractSamples(p->lwe, &tmpSample, &p->srcSamples->list_categories[index2].list[i], &p->srcSamples->list_categories[index2].list[j], p->bkwStepPar);
	            mutex_index = category / n_storage_mutex;

	            if (category > p->dstSamples->n_categories || category < 0)
	            {
	            	pthread_mutex_lock(&screen_mutex);
	                printf("ERROR subtractSamples: category %lu tot categories %lu \n", category, p->dstSamples->n_categories );
	                pthread_mutex_unlock(&screen_mutex);
	                exit(0);
	            }
	     		
	     		pthread_mutex_lock(&storage_mutex[mutex_index]);
	            // add it to the new list
	            n_samples_in_category = p->dstSamples->list_categories[category].n_samples;
	            if (n_samples_in_category < p->dstSamples->n_samples_per_category)
	            {
	            	if (!checkzero((char*)tmpSample.a, sizeof(u16)*p->lwe->n))
	            	{
	                    memcpy(p->dstSamples->list_categories[category].list[n_samples_in_category].a, tmpSample.a, p->lwe->n*sizeof(u16));
	                    p->dstSamples->list_categories[category].list[n_samples_in_category].z = tmpSample.z;
	                    // dstSamples->list_categories[category].list[n_samples_in_category].error = tmpSample.error;
	                    p->dstSamples->list_categories[category].n_samples++;
	                    p->dstSamples->n_samples++;

	                    if (p->dstSamples->n_samples >= p->dstSamples->max_samples){
	                    	pthread_mutex_unlock(&storage_mutex[mutex_index]);
		                    return NULL;	                    
	                    }
	                }
	            }
	            else
	                discarded++;
	            pthread_mutex_unlock(&storage_mutex[mutex_index]);
	        }
	    }

	    // process two adjacent category
	    for (int i=0; i<p->srcSamples->list_categories[index1].n_samples; i++)
	    {
	        for (int j=0; j<p->srcSamples->list_categories[index2].n_samples; j++)
	        {
	            category = addSamples(p->lwe, &tmpSample, &p->srcSamples->list_categories[index1].list[i], &p->srcSamples->list_categories[index2].list[j], p->bkwStepPar);
	            mutex_index = category / n_storage_mutex;

	            if (category > p->dstSamples->n_categories || category < 0)
	            {
	            	pthread_mutex_lock(&screen_mutex);
	                printf("ERROR addSamples: category %lu tot categories %lu \n", category, p->dstSamples->n_categories );
	                pthread_mutex_unlock(&screen_mutex);
	                exit(0);
	            }
	            
	            pthread_mutex_lock(&storage_mutex[mutex_index]);
	            // add it to the new list
	                n_samples_in_category = p->dstSamples->list_categories[category].n_samples;
	            if (n_samples_in_category < p->dstSamples->n_samples_per_category)
	            {
	            	if (!checkzero((char*)tmpSample.a, sizeof(u16)*p->lwe->n))
	            	{	                    
	                    memcpy(p->dstSamples->list_categories[category].list[n_samples_in_category].a, tmpSample.a, p->lwe->n*sizeof(u16));
	                    p->dstSamples->list_categories[category].list[n_samples_in_category].z = tmpSample.z;
	                    // dstSamples->list_categories[category].list[n_samples_in_category].error = tmpSample.error;
	                    p->dstSamples->list_categories[category].n_samples++;
	                    p->dstSamples->n_samples++;

	                    if (p->dstSamples->n_samples >= p->dstSamples->max_samples){
	                    	pthread_mutex_unlock(&storage_mutex[mutex_index]);
		                    return NULL;
	                    }
	                }
	            }
	            else
	                discarded++;
	            pthread_mutex_unlock(&storage_mutex[mutex_index]);
	        }
	    }
	}
}

int transition_bkw_step_smooth_lms(lweInstance *lwe, bkwStepParameters *bkwStepPar, sortedSamplesList *srcSamples, sortedSamplesList *dstSamples)
{

	ASSERT(NUM_THREADS >= 1, "Unexpected number of threads!");

    set_sorted_samples_list(dstSamples, lwe, bkwStepPar, srcSamples->n_samples);

    sample tmpSample;
    u16 sample_a[lwe->n];
    tmpSample.a = sample_a;

    u64 category, minc;
    int n_samples_in_category;

    discarded = 0;

    if (srcSamples->n_categories & 1)
    {
        // process single category
        for (int i = 0; i < srcSamples->list_categories[0].n_samples; i++)
        {
            for (int j=i+1; j < srcSamples->list_categories[0].n_samples; j++)
            {

                category = subtractSamples(lwe, &tmpSample, &srcSamples->list_categories[0].list[i], &srcSamples->list_categories[0].list[j], bkwStepPar);
                
                if (category > dstSamples->n_categories || category < 0)
                {
                    printf("ERROR subtractSamples: category %lu tot categories %lu \n", category, dstSamples->n_categories );
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
	                    // if (dstSamples->n_samples == dstSamples->max_samples)
		                   //  goto exit;
	                }
                }
                else
                    discarded++;
            }
        }
        minc = 1;
    }
    else
        minc = 0;

    pthread_t thread[NUM_THREADS];
    Params param[NUM_THREADS]; /* one set of in-/output paramaters per thread, so no need to lock these */

    n_storage_mutex = dstSamples->n_categories < MAX_NUM_STORAGE_MUTEXES ? dstSamples->n_categories : MAX_NUM_STORAGE_MUTEXES;
    storage_mutex = malloc(n_storage_mutex * sizeof(pthread_mutex_t));
    for (int i=0; i<n_storage_mutex; i++) { pthread_mutex_init(&storage_mutex[i], NULL); }

    u64 cat_per_thread = srcSamples->n_categories/NUM_THREADS;
    if(cat_per_thread & 1) // make it even
    	cat_per_thread--;

    /* load input parameters */
    for (int i=0; i<NUM_THREADS; i++) {
        param[i].lwe = lwe; /* set input parameter to thread number */
        param[i].bkwStepPar = bkwStepPar;
        param[i].srcSamples = srcSamples;
        param[i].dstSamples = dstSamples;
        param[i].minIndex1 = i == 0 ? minc : param[i-1].minIndex1 + cat_per_thread;
        param[i].maxIndex1 = param[i].minIndex1 + cat_per_thread;
    }
    param[NUM_THREADS-1].maxIndex1 = srcSamples->n_categories-2;

    /* process samples with LF2 method */
    /* start threads */
    for (int i = 0; i < NUM_THREADS; ++i)
    {
        if (!pthread_create(&thread[i], NULL, single_thread_lf2_work, (void*)&param[i])) {
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

    for (int i=0; i<n_storage_mutex; i++) { pthread_mutex_destroy(&storage_mutex[i]); }

    time_stamp("discarded samples: %ld", discarded);
    free(storage_mutex);

    return 0;
}

