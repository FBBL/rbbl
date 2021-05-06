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

#include "transition_times2_modq.h"
#include <string.h>
#include <pthread.h>

// global variables
typedef struct {
    lweInstance *lwe;
    bkwStepParameters *bkwStepPar;
    sortedSamplesList *sortedSamples;
    samplesList* unsortedSamples;
    u64 min;
    u64 max;
} Params;

/* define mutexes to protect common resources from concurrent access */
static pthread_mutex_t screen_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t *storage_mutex;
static int n_storage_mutex;

/* perform multiplication times 2 mod q of srcSample and store in dstSample. return the category */
int sample_times2_modq(sample *dstSample, sample *srcSample, lweInstance *lwe, bkwStepParameters *bkwStepPar)
{
    int n = lwe->n;
    int q = lwe->q;
    for (int i=0; i<n; i++)
        dstSample->a[i] = (srcSample->a[i] << 1) % q;
    dstSample->z = (srcSample->z << 1) % q;
    // dstSample->error = (2*srcSample->error) % q;

    int index = position_values_2_category_index(lwe, bkwStepPar, dstSample->a);

    return index;
}

void *single_thread_work(void *params){

    Params *p = (Params*)params;

    sample tmpSample;
    u64 count = 0, category;
    int n_samples_in_category;
    int mutex_index;

    u16 sample_a[p->lwe->n];
    tmpSample.a = sample_a;

    for(count = p->min; count < p->max; count++)
    {
        category = sample_times2_modq(&tmpSample, &p->unsortedSamples->list[count], p->lwe, p->bkwStepPar);

        n_samples_in_category = p->sortedSamples->list_categories[category].n_samples;

        if (n_samples_in_category < p->sortedSamples->n_samples_per_category && p->sortedSamples->n_samples < p->sortedSamples->max_samples)
        {
            mutex_index = category / n_storage_mutex;
            pthread_mutex_lock(&storage_mutex[mutex_index]);
            if (!checkzero((char*)tmpSample.a, sizeof(u16)*(p->lwe->n)))
            {
                if (category > p->sortedSamples->n_categories)
                {
                    pthread_mutex_lock(&screen_mutex);
                    printf("ERROR: category %lu tot categories %lu \n", category, p->sortedSamples->n_categories );
                    pthread_mutex_unlock(&screen_mutex);
                    pthread_mutex_unlock(&storage_mutex[mutex_index]);
                    exit(0);
                }
                memcpy(p->sortedSamples->list_categories[category].list[n_samples_in_category].a, tmpSample.a, p->lwe->n*sizeof(u16));
                p->sortedSamples->list_categories[category].list[n_samples_in_category].z = tmpSample.z;
                // sortedSamples->list_categories[category].list[n_samples_in_category].error = tmpSample.error;
                p->sortedSamples->list_categories[category].n_samples++;
                p->sortedSamples->n_samples++;
            }
            pthread_mutex_unlock(&storage_mutex[mutex_index]);
        }
    }
}


/* multiply each sample times 2 mod q in sortedSamples. Then store the result in dstSortedSamplesList according to its category */
int transition_times2_modq(lweInstance *lwe, bkwStepParameters *bkwStepPar, sortedSamplesList *sortedSamples, samplesList* unsortedSamples)
{

    ASSERT(NUM_THREADS >= 1, "Unexpected number of threads!");

    pthread_t thread[NUM_THREADS];
    Params param[NUM_THREADS]; /* one set of in-/output paramaters per thread, so no need to lock these */

    /* allocate memory for an optimal number of mutex */
    n_storage_mutex = sortedSamples->n_categories < MAX_NUM_STORAGE_MUTEXES ? sortedSamples->n_categories : MAX_NUM_STORAGE_MUTEXES;
    storage_mutex = malloc(n_storage_mutex * sizeof(pthread_mutex_t));
    for (int i=0; i<n_storage_mutex; i++) { pthread_mutex_init(&storage_mutex[i], NULL); }

    /* load input parameters */
    for (int i=0; i<NUM_THREADS; i++) {
        param[i].lwe = lwe; /* set input parameter to thread number */
        param[i].bkwStepPar = bkwStepPar;
        param[i].sortedSamples = sortedSamples;
        param[i].unsortedSamples = unsortedSamples;
        param[i].min = i*(unsortedSamples->n_samples/NUM_THREADS);
        param[i].max = (i+1)*(unsortedSamples->n_samples/NUM_THREADS);
    }
    param[NUM_THREADS-1].max = unsortedSamples->n_samples;


    /* start threads */
    for (int i = 0; i < NUM_THREADS; ++i)
    {
        if (!pthread_create(&thread[i], NULL, single_thread_work, (void*)&param[i])) {
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

    free(storage_mutex);

}

