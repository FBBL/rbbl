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

#include "solve_fwht.h"
#include "lwe_instance.h"

#include <stdint.h>
#include <string.h>
#include <math.h>
#include <pthread.h>

// #define PRINT_INTERMEDIATE

/*
 * integer to binary sequence
 */
void int_to_bin(u64 input, u8 *binary, int binlen)
{
    // zero all entries
    for (int i = 0; i<binlen; i++)
        binary[i] = 0;

    // counter for binary array
    int i = 0;
    while (input > 0)
    {
        // storing remainder in binary array
        binary[i] = input % 2;
        input = input / 2;
        i++;
    }
}

/*
 * Convert sample to binary sequence, then to integer (u64)
 */
u64 sample_to_int(u16 *input, int len, int q)
{

    u64 output = 0;
    short a;

    for (int i = 0; i<len; i++)
    {
        a = input[i] <= q/2 ? input[i] : abs(input[i]-q);
        if (a % 2 != 0)
            output += ((u64)1)<<i;
    }

    return output;
}

void FWHT (long* data, int size)
{
    int n = log2(size);
    long tmp;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < (1 << n); j += 1 << (i+1))
        {
            for (int k = 0; k < (1<<i); ++k)
            {
                int a = j + k;
                int b = j + k + (1<<i);

                tmp = data[a];
                data[a] += data[b];
                data[b] = tmp - data[b];
            }
        }
    }
}


/* Retrieve binary secret using Fast Walsh Hadamard Transform */
int solve_fwht_search(u8 *binary_solution, int zero_positions, int fwht_positions, unsortedSamplesList *reducedSamples, lweInstance *lwe)
{

    int n = lwe->n;
    int q = lwe->q;

    ASSERT(fwht_positions <= MAX_FWHT, "The number of positions for fwht is not supported in this implementation!\n");

    if(zero_positions == n)
    {
        printf("Coordinates all solved\n");
        return 1;
    }

    /* create initial list */
    u64 N = (u64)1<<fwht_positions; // N = 2^fwht_positions
    long* list = calloc(N, sizeof(long));
    if (!list)
    {
        printf("*** solve_fwht_search: failed to allocate memory for initial list\n");
        exit(-1);
    }

    short z, lsb_z;
    u64 intsample;

    // process all samples
    for (u64 i = 0; i < reducedSamples->n_samples; ++i)
    {

        intsample = sample_to_int(&reducedSamples->a_list[i*n + zero_positions], fwht_positions, q);
        z = reducedSamples->z_list[i] > (q-1)/2 ? (reducedSamples->z_list[i] -q) : (reducedSamples->z_list[i]);
        lsb_z = z%2 == 0 ? 0 : 1;
        if (lsb_z == 0)
            list[intsample] += 1;
        else
            list[intsample] -= 1;
    }

    /* Apply Fast Walsh Hadamard Tranform */
    FWHT(list, N);

    // find maximum
    u64 max_pos = -1;
    double max = 0;
    double tot = 0;
    for (u64 i = 0; i<N; i++)
    {
        tot += labs(list[i]);
        if (max < labs(list[i]))
        {
            max = labs(list[i]);
            max_pos = i;
        }
    }

    // Convert solution into binary
    int_to_bin(max_pos, binary_solution, fwht_positions);

    free(list);
    return 0;
}



/* THREADED HYBRID SEARCH */


/* min and max refer to the first bruteforce position, 
  return 0 if guesses are terminated, otherwise 1 */
int nextBruteForceGuess(int min, int max, int ratio, int *BFguess, int lenght)
{

    for(int i = lenght-1; i>=1; i--)
    {
        if(BFguess[i] < ratio)
        {
            BFguess[i]++;
            for(int j = lenght-1; j>i; j--)
                BFguess[j] = -ratio;
            return 1;
        }
    }

	if (BFguess[0] < max) {
		BFguess[0]++;
		for(int j = lenght-1; j>0; j--)
                BFguess[j] = -ratio;
		return 1;
	}
	else
		return 0;

}

// global variables
typedef struct {
    lweInstance *lwe;
    unsortedSamplesList *Samples;
    int bf_positions;
    int fwht_positions;
    int zero_positions;
    int min;
    int max;
    int ratio;
} Params;

static double global_max;
static u8 *p_binary_solution;
static short *p_bf_solution;

static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t screen_mutex = PTHREAD_MUTEX_INITIALIZER;

void *single_work_bruteforce_guess(void *params){

	Params *p = (Params*)params;
	int q = p->lwe->q;

	/* create initial list */
    u64 N = (u64)1 << p->fwht_positions; // N = 2^fwht_positions
    long* list = calloc(N,sizeof(long));
    if (!list)
    {
        printf("*** solve_fwht_search: failed to allocate memory for initial list\n");
        exit(-1);
    }

    long max_pos = -1;
    double max = 0;
    u64 intsample;
    short z, lsb_z;

    u8 bin_guess[p->fwht_positions];

    int BFguess[p->bf_positions];
    int tmp_si;

    // initialize brute force guess
    BFguess[0] = p->min;
    for(int i = 1; i<p->bf_positions; i++)
        BFguess[i] = -p->ratio;

	do{

       memset(list, 0, N*sizeof(long));

        // process all samples
        for (u64 i = 0; i < p->Samples->n_samples; ++i)
        {

            intsample = sample_to_int(&p->Samples->a_list[i*p->lwe->n+ p->zero_positions], p->fwht_positions, q);
            z = p->Samples->z_list[i];
            // update z with the bruteforce-guessed positions
            for(int j = 0; j<p->bf_positions; j++)
            {
                if(BFguess[j] < 0)
                    tmp_si = BFguess[j]+q;
                else
                    tmp_si = BFguess[j];
                z = (z - (p->Samples->a_list[i*p->lwe->n+p->zero_positions+p->fwht_positions+j]*tmp_si) % q );
                if (z < 0)
                    z += q;
            }
            z = z > (q-1)/2 ? z -q : z;
            intsample = sample_to_int(&p->Samples->a_list[i*p->lwe->n + p->zero_positions], p->fwht_positions, q);
            lsb_z = z%2 == 0 ? 0 : 1;

            if (lsb_z == 0)
                list[intsample] += 1;
            else
                list[intsample] -= 1;
        }

        /* Apply Fast Walsh Hadamard Tranform */
        FWHT(list, N);

        // find maximum
        max_pos = -1;
        max = 0;
        for (u64 i = 0; i<N; i++)
        {
            if (max < labs(list[i]))
            {
                max = labs(list[i]);
                max_pos = i;
            }
        }

        // Convert solution into binary
        int_to_bin(max_pos, bin_guess, p->fwht_positions);
#ifdef PRINT_INTERMEDIATE
        pthread_mutex_lock(&screen_mutex);
        printf("Index found %ld - max %f \n(", max_pos, max);
        for(int j = 0; j<p->bf_positions; j++)
            printf("%d ", BFguess[j]);
        printf(") - (");
        for(int j = 0; j<p->fwht_positions; j++)
            printf("%hu ", bin_guess[j]);
        printf(")\n");
        pthread_mutex_unlock(&screen_mutex);
#endif
        if(max > global_max)
        {
        	pthread_mutex_lock(&global_mutex);
            global_max = max;
            for(int j = 0; j<p->fwht_positions; j++)
                p_binary_solution[j] = bin_guess[j];
            for(int j = 0; j<p->bf_positions; j++)
                p_bf_solution[j] = BFguess[j] >= 0 ? BFguess[j] : BFguess[j] +q;
            pthread_mutex_unlock(&global_mutex);
        }

	} while(nextBruteForceGuess(p->min, p->max, p->ratio, BFguess, p->bf_positions));

	free(list);

}


/* Hybrid solver that uses brute-force for bf_positions number of positions and Fast Walsh Hadamard Transform for fwht_positions number of positions
 */
int solve_fwht_search_bruteforce(u8 *binary_solution, short *bf_solution, int zero_positions, int bf_positions, int fwht_positions, unsortedSamplesList *reducedSamples, lweInstance *lwe)
{
    int q = lwe->q;
    int n = lwe->n;

    ASSERT(1 <= fwht_positions && fwht_positions <= MAX_FWHT, "The number of positions for fwht is not supported in this implementation!\n");
    ASSERT(1 <= bf_positions && bf_positions <= MAX_BRUTE_FORCE, "The number of positions for bruteforce guessing is not supported in this implementation!\n");
    ASSERT(fwht_positions + bf_positions + zero_positions == n, "The number of positions for bruteforce and fwht is => n!\n");
    ASSERT(NUM_THREADS >= 1, "Unexpected number of threads!");

    time_stamp("Start FWHT: brute force %d positions, guess %d positions\n", bf_positions, fwht_positions);

    int ratio = round(lwe->sigma*4); // 2*3*standard_deviation is the interval length where to search

	pthread_t thread[NUM_THREADS];
    Params param[NUM_THREADS]; /* one set of in-/output paramaters per thread, so no need to lock these */


    /* load input parameters */
    for (int i=0; i<NUM_THREADS; i++) {
        param[i].lwe = lwe; /* set input parameter to thread number */
        param[i].Samples = reducedSamples;
        param[i].bf_positions = bf_positions;
        param[i].fwht_positions = fwht_positions;
        param[i].zero_positions = zero_positions;
        param[i].min = ((2*ratio+1)/NUM_THREADS)*i-ratio;
        param[i].max = ((2*ratio+1)/NUM_THREADS)*(i+1)-ratio;
        param[i].ratio = ratio;
    }
    param[NUM_THREADS-1].max = ratio+1;

    global_max = 0;
    p_bf_solution = bf_solution;
    p_binary_solution = binary_solution;

 
    /* start threads */
    for (int i = 0; i < NUM_THREADS; ++i)
    {
        if (!pthread_create(&thread[i], NULL, single_work_bruteforce_guess, (void*)&param[i])) {
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

    return 0;
}