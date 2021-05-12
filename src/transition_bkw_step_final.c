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
    unsortedSamplesList* dstSamples;
    u64 minIndex1;
    u64 maxIndex1;
    u64 maxTotSamples;
} Params;

static u64 subtractSamples(lweInstance *lwe, u16 *dst_a, u16 *dst_z, u16 *src1_a, u16 src1_z, u16 *src2_a, u16 src2_z)
{

    for (int i=0; i < lwe->n; i++){
        dst_a[i] = (src1_a[i] - src2_a[i] + lwe->q);
        dst_a[i] = dst_a[i] >= lwe->q ? dst_a[i] -lwe->q : dst_a[i];
    }
    *dst_z = (src1_z - src2_z + lwe->q);
    *dst_z = *dst_z >= lwe->q ? *dst_z -lwe->q : *dst_z;

    return 0;
}

static u64 addSamples(lweInstance *lwe, u16 *dst_a, u16 *dst_z, u16 *src1_a, u16 src1_z, u16 *src2_a, u16 src2_z)
{

    for (int i=0; i < lwe->n; i++){
        dst_a[i] = (src1_a[i] + src2_a[i]);
        dst_a[i] = dst_a[i] >= lwe->q ? dst_a[i] -lwe->q : dst_a[i];
    }
    *dst_z = (src1_z + src2_z);
    *dst_z = *dst_z >= lwe->q ? *dst_z -lwe->q : *dst_z;

    return 0;
}


/* Return 1 if passes unnatural selection, 0 otherwise */
static int unnaturalSelection(lweInstance *lwe, u16 *dst_a, bkwStepParameters *srcBkwStepPar)
{
    int n = lwe->n;
    int q = lwe->q;

    int tmp;

    for (int i=0; i < srcBkwStepPar->startIndex + srcBkwStepPar->numPositions; i++)
    {
        tmp = dst_a[i] < q/2 ? dst_a[i] : (int)dst_a[i] -q;
        if (tmp > srcBkwStepPar->un_selection || tmp < -srcBkwStepPar->un_selection)
            return 0;
    }

    return 1;
}

void *single_thread_final_lf2_work(void *params){

    Params *p = (Params*)params;

    u16 tmp_a[p->lwe->n];
    u16 tmp_z = 0;
    u16 tmp_e;

    u16 in, jn;

    u64 index1, index2, category, sample_index;
    int n_samples_in_category, mutex_index;

    int block_a = p->lwe->n*SAMPLES_PER_CATEGORY;
    int block_z = SAMPLES_PER_CATEGORY;

    /* process samples with LF2 method */
    for (index1 = p->minIndex1; index1 < p->maxIndex1; index1 = index1+2)
    {
        index2 = index1+1;

        if (!(p->dstSamples->n_samples < p->maxTotSamples)){
            return NULL;
        }

        // process single category
        for (int i = 0; i < p->srcSamples->n_in_categories[index1]; i++)
        {
            in = i*p->lwe->n;
            for (int j=i+1; j < p->srcSamples->n_in_categories[index1]; j++)
            {
                jn = j*p->lwe->n;
                subtractSamples(p->lwe, tmp_a, &tmp_z, &p->srcSamples->a_list[index1*block_a +in], p->srcSamples->z_list[index1*block_z +i], &p->srcSamples->a_list[index1*block_a+jn], p->srcSamples->z_list[index1*block_z +j]);

                // error - DEBUG
                tmp_e = (p->srcSamples->e_list[index1*block_z +i] - p->srcSamples->e_list[index1*block_z +j] +p->lwe->q) %p->lwe->q;

                // add it to the new list
            	if (!checkzero((char*)tmp_a, sizeof(u16)*p->lwe->n) && unnaturalSelection(p->lwe, tmp_a, p->bkwStepPar))
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
                    memcpy(&p->dstSamples->a_list[sample_index*p->lwe->n], tmp_a, p->lwe->n*sizeof(u16));
                    p->dstSamples->z_list[sample_index] = tmp_z;
                    p->dstSamples->e_list[sample_index] = tmp_e;
                }
                else
                {
                  	discarded++;
                }
            }
        }

        // process single category
        for (int i = 0; i < p->srcSamples->n_in_categories[index2]; i++)
        {
            in = i*p->lwe->n;
            for (int j=i+1; j < p->srcSamples->n_in_categories[index2]; j++)
            {
                jn = j*p->lwe->n;
                subtractSamples(p->lwe, tmp_a, &tmp_z, &p->srcSamples->a_list[index2*block_a +in], p->srcSamples->z_list[index2*block_z +i], &p->srcSamples->a_list[index2*block_a+jn], p->srcSamples->z_list[index2*block_z +j]);

                // error - DEBUG
                tmp_e = (p->srcSamples->e_list[index2*block_z +i] - p->srcSamples->e_list[index2*block_z +j] +p->lwe->q) %p->lwe->q;

                // add it to the new list
            	if (!checkzero((char*)tmp_a, sizeof(u16)*p->lwe->n) && unnaturalSelection(p->lwe, tmp_a, p->bkwStepPar))
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
                    memcpy(&p->dstSamples->a_list[sample_index*p->lwe->n], tmp_a, p->lwe->n*sizeof(u16));
                    p->dstSamples->z_list[sample_index] = tmp_z;
                    p->dstSamples->e_list[sample_index] = tmp_e;
                }
                else
                {
                    discarded++;
                }
            }
        }

        // process two adjacent category
        for (int i=0; i < p->srcSamples->n_in_categories[index1]; i++)
        {
            in = i*p->lwe->n;
            for (int j=0; j < p->srcSamples->n_in_categories[index2]; j++)
            {
                jn = j*p->lwe->n;
                addSamples(p->lwe, tmp_a, &tmp_z, &p->srcSamples->a_list[index1*block_a+in], p->srcSamples->z_list[index1*block_z +i], &p->srcSamples->a_list[index2*block_a +jn], p->srcSamples->z_list[index2*block_z +j]);

                // error - DEBUG
                tmp_e = (p->srcSamples->e_list[index1*block_z +i] + p->srcSamples->e_list[index2*block_z +j] +p->lwe->q) %p->lwe->q;
                
                // add it to the new list
            	if (!checkzero((char*)tmp_a, sizeof(u16)*p->lwe->n) && unnaturalSelection(p->lwe, tmp_a, p->bkwStepPar))
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
                    memcpy(&p->dstSamples->a_list[sample_index*p->lwe->n], tmp_a, p->lwe->n*sizeof(u16));
                    p->dstSamples->z_list[sample_index] = tmp_z;

                    p->dstSamples->e_list[sample_index] = tmp_e;
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

int transition_bkw_step_final(lweInstance *lwe, bkwStepParameters *srcBkwStepPar, sortedSamplesList *srcSamples, unsortedSamplesList *dstSamples, u64 maxSamples)
{

    ASSERT(NUM_THREADS >= 1, "Unexpected number of threads!");

    u16 tmp_a[lwe->n];
    u16 tmp_z = 0;
    u16 tmp_e;

    u64 index1, index2, minc, sample_index;
    dstSamples->n_samples = 0;
    discarded = 0;

    if (srcSamples->n_categories & 1){

        // process single category
        for (int i = 0; i < srcSamples->n_in_categories[0]; i++)
        {
            for (int j=i+1; j < srcSamples->n_in_categories[0]; j++)
            {
                subtractSamples(lwe, tmp_a, &tmp_z, &srcSamples->a_list[0+i*lwe->n], srcSamples->z_list[0+i], &srcSamples->a_list[0+j*lwe->n], srcSamples->z_list[0+j]);
                
                // error - DEBUG
                tmp_e = (srcSamples->z_list[0+i] - srcSamples->z_list[j] + lwe->q) %lwe->q;

                if (!checkzero((char*)tmp_a, sizeof(u16)*lwe->n) && unnaturalSelection(lwe, tmp_a, srcBkwStepPar))
                {
	                // add it to the new list
	                if (dstSamples->n_samples < maxSamples)
	                {
                        memcpy(&dstSamples->a_list[sample_index*lwe->n], tmp_a, lwe->n*sizeof(u16));
                        dstSamples->z_list[sample_index] = tmp_z;
                        dstSamples->e_list[sample_index] = tmp_e;
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

