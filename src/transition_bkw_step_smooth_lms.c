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
#include "error_rate.h"

#include <math.h>
#include <string.h>
#include <pthread.h>

#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

/* define mutexes to protect common resources from concurrent access */
static pthread_mutex_t screen_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t *storage_mutex;
static int n_storage_mutex;
static u64 max_num_categories_per_mutex;
static u64 discarded = 0;

typedef struct {
    lweInstance *lwe;
    bkwStepParameters *bkwStepPar;
    sortedSamplesList *srcSamples;
    sortedSamplesList* dstSamples;
    u64 minIndex1;
    u64 maxIndex1;
} Params;

static u64 subtractSamples(lweInstance *lwe, u16 *dst_a, u16 *dst_z, u16 *src1_a, u16 src1_z, u16 *src2_a, u16 src2_z,  bkwStepParameters *bkwStepPar)
{

    for (int i=0; i < lwe->n; i++){
        dst_a[i] = (src1_a[i] - src2_a[i] + lwe->q);
        dst_a[i] = dst_a[i] >= lwe->q ? dst_a[i] -lwe->q : dst_a[i];
    }
    *dst_z = (lwe->q + src1_z - src2_z) ;
    *dst_z = *dst_z >= lwe->q ? (*dst_z) -lwe->q : *dst_z;

    u64 index = position_values_2_category_index(lwe, bkwStepPar, dst_a + bkwStepPar->startIndex);

    return index;
}

static u64 addSamples(lweInstance *lwe, u16 *dst_a, u16 *dst_z, u16 *src1_a, u16 src1_z, u16 *src2_a, u16 src2_z,  bkwStepParameters *bkwStepPar)
{

    for (int i=0; i < lwe->n; i++){
        dst_a[i] = (src1_a[i] + src2_a[i]);
        dst_a[i] = dst_a[i] >= lwe->q ? dst_a[i] -lwe->q : dst_a[i];
    }
    *dst_z = (src1_z + src2_z);
    *dst_z = *dst_z >= lwe->q ? *dst_z -lwe->q : *dst_z;

    u64 index = position_values_2_category_index(lwe, bkwStepPar, dst_a + bkwStepPar->startIndex);

    return index;
}

#ifdef DEBUG
// return 1 if sample is wrong
int verify_sample(u16 *a, u16 z, u16 e, lweInstance *lwe){

	int sum = 0;
	for (int k = 0; k < lwe->n; ++k)
    {
        sum = (sum + a[k]*lwe->s[k]) % lwe->q;
    }
    sum = (sum + e) % lwe->q;
    if (sum != z)
    {
        printf("TADAAA in sample verification\n");
        return 1;
    }
    return 0;
}
#endif

void *single_thread_lf2_work(void *params){

	Params *p = (Params*)params;

	u64 index1, index2, category;
	int n_in_categories, mutex_index;

    // pthread_mutex_lock(&screen_mutex);
    // printf("Min %lu Max %lu\n", p->minIndex1, p->maxIndex1);
    // pthread_mutex_unlock(&screen_mutex);

    u16 tmp_a[p->lwe->n];
    u16 tmp_z = 0;
   	u16 in, jn;

#ifdef DEBUG
    u16 tmp_e;
#endif

    int block_a = p->lwe->n*SAMPLES_PER_CATEGORY;
    int block_z = SAMPLES_PER_CATEGORY;

	for (index1 = p->minIndex1; index1 < p->maxIndex1; index1 = index1+2)
	{
		index2 = index1+1;

        if (!(p->dstSamples->n_samples < p->dstSamples->max_samples)){        	
            return NULL;
        }

	    // process single category
	    for (int i = 0; i < p->srcSamples->n_in_categories[index1]; i++)
	    {
	    	in = i*p->lwe->n;
	        for (int j=i+1; j < p->srcSamples->n_in_categories[index1]; j++)
	        {
	        	jn = j*p->lwe->n;

	            category = subtractSamples(p->lwe, tmp_a, &tmp_z, &p->srcSamples->a_list[index1*block_a +in], p->srcSamples->z_list[index1*block_z +i], &p->srcSamples->a_list[index1*block_a +jn], p->srcSamples->z_list[index1*block_z +j], p->bkwStepPar);

#ifdef DEBUG
	            tmp_e = (p->srcSamples->e_list[index1*block_z +i] - p->srcSamples->e_list[index1*block_z +j] +p->lwe->q) %p->lwe->q;
	            if(verify_sample(tmp_a, tmp_z, tmp_e, p->lwe)){

	            	printf("Subtract 1\n");
	            	printf("%d %d < %d \n", i, j, p->srcSamples->n_in_categories[index1] );

		            for (int k = 0; k < p->lwe->n; k++)
		            {
		            	printf("%d - %d = %d\n", p->srcSamples->a_list[index1*block_a +in + k], p->srcSamples->a_list[index1*block_a +jn + k], tmp_a[k]);           	
		            }
		            printf("%d - %d = %d \n", p->srcSamples->z_list[index1*block_z +i], p->srcSamples->z_list[index1*block_z +j], tmp_z);
		            printf("%d - %d = %d \n", p->srcSamples->e_list[index1*block_z +i], p->srcSamples->e_list[index1*block_z +j], tmp_e);
		            exit(0);
	            }
#endif
	            mutex_index = category / max_num_categories_per_mutex;

	            if (category > p->dstSamples->n_categories || category < 0)
	            {
	            	pthread_mutex_lock(&screen_mutex);
	                printf("ERROR subtractSamples: category %lu tot categories %lu \n", category, p->dstSamples->n_categories );
	                pthread_mutex_unlock(&screen_mutex);
	                exit(0);
	            }

	            pthread_mutex_lock(&storage_mutex[mutex_index]);
	            // add it to the new list
	            n_in_categories = p->dstSamples->n_in_categories[category];
	            if (n_in_categories < SAMPLES_PER_CATEGORY)
	            {
	            	if (!checkzero((char*)tmp_a, sizeof(u16)*p->lwe->n))
	            	{
	                    memcpy(&p->dstSamples->a_list[category*block_a + n_in_categories*p->lwe->n], tmp_a, p->lwe->n*sizeof(u16));
	                    p->dstSamples->z_list[category*block_z +n_in_categories] = tmp_z;
#ifdef DEBUG
	                    p->dstSamples->e_list[category*block_z +n_in_categories] = tmp_e;
#endif
	                    p->dstSamples->n_in_categories[category]++;
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
	    for (int i = 0; i < p->srcSamples->n_in_categories[index2]; i++)
	    {
	    	in = i*p->lwe->n;
	        for (int j=i+1; j < p->srcSamples->n_in_categories[index2]; j++)
	        {
	        	jn = j*p->lwe->n;
	            category = subtractSamples(p->lwe, tmp_a, &tmp_z, &p->srcSamples->a_list[index2*block_a +in], p->srcSamples->z_list[index2*block_z +i], &p->srcSamples->a_list[index2*block_a +jn], p->srcSamples->z_list[index2*block_z +j], p->bkwStepPar);

#ifdef DEBUG
	            tmp_e = (p->srcSamples->e_list[index2*block_z +i] - p->srcSamples->e_list[index2*block_z +j] +p->lwe->q) %p->lwe->q;
	            if(verify_sample(tmp_a, tmp_z, tmp_e, p->lwe)){

	            	printf("Subtract 2\n");
	            	printf("%d %d < %d \n", i, j, p->srcSamples->n_in_categories[index2] );

		            for (int k = 0; k < p->lwe->n; k++)
		            {
		            	printf("%d - %d = %d\n", p->srcSamples->a_list[index2*block_a +in + k], p->srcSamples->a_list[index2*block_a +jn + k], tmp_a[k]);           	
		            }
		            printf("%d - %d = %d \n", p->srcSamples->z_list[index2*block_z +i], p->srcSamples->z_list[index2*block_z +j], tmp_z);
		            printf("%d - %d = %d \n", p->srcSamples->e_list[index2*block_z +i], p->srcSamples->e_list[index2*block_z +j], tmp_e);
		            exit(0);
	            }
#endif
	            mutex_index = category / max_num_categories_per_mutex;

	            if (category > p->dstSamples->n_categories || category < 0)
	            {
	            	pthread_mutex_lock(&screen_mutex);
	                printf("ERROR subtractSamples: category %lu tot categories %lu \n", category, p->dstSamples->n_categories );
	                pthread_mutex_unlock(&screen_mutex);
	                exit(0);
	            }
	     		
	     		pthread_mutex_lock(&storage_mutex[mutex_index]);
	            // add it to the new list
	            n_in_categories = p->dstSamples->n_in_categories[category];
	            if (n_in_categories < SAMPLES_PER_CATEGORY)
	            {
	            	if (!checkzero((char*)tmp_a, sizeof(u16)*p->lwe->n))
	            	{
	                    memcpy(&p->dstSamples->a_list[category*block_a +n_in_categories*p->lwe->n], tmp_a, p->lwe->n*sizeof(u16));
	                    p->dstSamples->z_list[category*block_z +n_in_categories] = tmp_z;
#ifdef DEBUG
	                    p->dstSamples->e_list[category*block_z +n_in_categories] = tmp_e;
#endif
	                    p->dstSamples->n_in_categories[category]++;
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
	    for (int i=0; i<p->srcSamples->n_in_categories[index1]; i++)
	    {
	    	in = i*p->lwe->n;
	        for (int j=0; j<p->srcSamples->n_in_categories[index2]; j++)
	        {
	        	jn = j*p->lwe->n;
	            category = addSamples(p->lwe, tmp_a, &tmp_z, &p->srcSamples->a_list[index1*block_a +in], p->srcSamples->z_list[index1*block_z +i], &p->srcSamples->a_list[index2*block_a +jn], p->srcSamples->z_list[index2*block_z +j], p->bkwStepPar);
#ifdef DEBUG
	            tmp_e = (p->srcSamples->e_list[index1*block_z +i] + p->srcSamples->e_list[index2*block_z +j] +p->lwe->q) %p->lwe->q;
	            if(verify_sample(tmp_a, tmp_z, tmp_e, p->lwe)){

	            	printf("Addition \n");
	            	printf("%d: %d < %d \n", index1, i, p->srcSamples->n_in_categories[index1] );
	            	printf("%d: %d < %d \n", index2, j, p->srcSamples->n_in_categories[index2] );

		            for (int k = 0; k < p->lwe->n; k++)
		            {
		            	printf("%d - %d = %d\n", p->srcSamples->a_list[index1*block_a +in + k], p->srcSamples->a_list[index2*block_a +jn + k], tmp_a[k]);           	
		            }
		            printf("%d - %d = %d \n", p->srcSamples->z_list[index1*block_z +i], p->srcSamples->z_list[index2*block_z +j], tmp_z);
		            printf("%d - %d = %d \n", p->srcSamples->e_list[index1*block_z +i], p->srcSamples->e_list[index2*block_z +j], tmp_e);
		            exit(0);
	            }
#endif
	            mutex_index = category / max_num_categories_per_mutex;

	            if (category > p->dstSamples->n_categories || category < 0)
	            {
	            	pthread_mutex_lock(&screen_mutex);
	                printf("ERROR addSamples: category %lu tot categories %lu \n", category, p->dstSamples->n_categories );
	                pthread_mutex_unlock(&screen_mutex);
	                exit(0);
	            }
	            
	            pthread_mutex_lock(&storage_mutex[mutex_index]);
	            // add it to the new list
	                n_in_categories = p->dstSamples->n_in_categories[category];
	            if (n_in_categories < SAMPLES_PER_CATEGORY)
	            {
	            	if (!checkzero((char*)tmp_a, sizeof(u16)*p->lwe->n))
	            	{	                    
	                    memcpy(&p->dstSamples->a_list[category*block_a +n_in_categories*p->lwe->n], tmp_a, p->lwe->n*sizeof(u16));
	                    p->dstSamples->z_list[category*block_z +n_in_categories] = tmp_z;
#ifdef DEBUG
	                    p->dstSamples->e_list[category*block_z +n_in_categories] = tmp_e;
#endif
	                    p->dstSamples->n_in_categories[category]++;
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

    u16 tmp_a[lwe->n];
    u16 tmp_z = 0;
#ifdef DEBUG
    u16 tmp_e;
#endif
    u64 category, minc;
    int n_in_categories, in, jn;

    discarded = 0;

    int block_a = lwe->n*SAMPLES_PER_CATEGORY;
    int block_z = SAMPLES_PER_CATEGORY;

    if (srcSamples->n_categories & 1)
    {
        // process single category
        for (int i = 0; i < srcSamples->n_in_categories[0]; i++)
        {
        	in = i*lwe->n;
            for (int j=i+1; j < srcSamples->n_in_categories[0]; j++)
            {
            	jn = j*lwe->n;
                category = subtractSamples(lwe, tmp_a, &tmp_z, &srcSamples->a_list[0+i*lwe->n], srcSamples->z_list[0+i], &srcSamples->a_list[0+j*lwe->n], srcSamples->z_list[0+j], bkwStepPar);

#ifdef DEBUG
                tmp_e = (srcSamples->e_list[i] - srcSamples->e_list[j] + lwe->q) %lwe->q;
				if(verify_sample(tmp_a, tmp_z, tmp_e, lwe)){
		            for (int k = 0; k < lwe->n; k++)
		            {
		            	printf("%d - %d = %d\n", srcSamples->a_list[in + k], srcSamples->a_list[jn + k], tmp_a[k]);           	
		            }
		            printf("z %d - %d = %d \n", srcSamples->z_list[i], srcSamples->z_list[j], tmp_z);
		            printf("error %d - %d = %d \n", srcSamples->e_list[i], srcSamples->e_list[j], tmp_e);
		            exit(0);
	            }
#endif

                if (category > dstSamples->n_categories || category < 0)
                {
                    printf("ERROR subtractSamples: category %lu tot categories %lu \n", category, dstSamples->n_categories );
                    exit(0);
                }

                // add it to the new list
                n_in_categories = dstSamples->n_in_categories[category];
                if (n_in_categories < SAMPLES_PER_CATEGORY)
                {
                	if (!checkzero((char*)tmp_a, sizeof(u16)*lwe->n))
                	{   
                memcpy(&srcSamples->a_list[block_a*category+n_in_categories*lwe->n], tmp_a, lwe->n*sizeof(u16));
                srcSamples->z_list[block_z*category+n_in_categories] = tmp_z;
#ifdef DEBUG
                srcSamples->e_list[block_z*category+n_in_categories] = tmp_e;
#endif
                srcSamples->n_in_categories[category]++;
                srcSamples->n_samples++;
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
    max_num_categories_per_mutex = (dstSamples->n_categories + n_storage_mutex - 1) / n_storage_mutex;
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

