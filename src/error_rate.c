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


/* Compute Error rate when all positions have been reduced and no bruteforce is used */
int error_rate(int zero_positions, unsortedSamplesList *Samples, lweInstance *lwe)
{

    int n = lwe->n;
    int q = lwe->q;

    u64 errors = 0;

    u8 bin_secret[n], bina, binz, sum;

    // make secret binary
    for (int i = 0; i < n; i++)
    {
        if (lwe->s[i] < q/2)
            bin_secret[i] = lwe->s[i] % 2;
        else
            bin_secret[i] = (lwe->s[i] +1) % 2;
    }

    // process all samples
    for (u64 i = 0; i < Samples->n_samples; i++)
    {
        sum = 0;
        for (int j = zero_positions; j < n; j++)
        {
            bina = Samples->a_list[i*n+j] < q/2 ? Samples->a_list[i*n+j] % 2 : (Samples->a_list[i*n+j] +1) % 2;
            sum = (sum + bina*bin_secret[j]) % 2;
        }
        binz = Samples->z_list[i] < q/2 ? Samples->z_list[i] % 2 : (Samples->z_list[i] +1) % 2;

        // if (i == 1){
        //         for (int k = 0; k < n; k++)
        //         {
        //             printf("%d ", Samples->a_list[k]);
        //         }printf(" - ");
        //     printf("%d %d \n", sum, binz );
        //     printf("\n");
        // }

        // if(i % 1000 == 0)
        // {
        //     int sum = 0;
        //     for (int k = 0; k < n; ++k)
        //     {
        //         sum = (sum + Samples->a_list[i*n+k]*lwe->s[k]) % q;
        //     }
        //     sum = (sum + Samples->e_list[i]) % q;
        //     printf("%d - %d\n", sum, Samples->z_list[i]);
        // }


        errors += sum == binz ? 0 : 1;
    }

    time_stamp("Error rate: %lf", (float)errors/(float)Samples->n_samples );

    return 0;
}



int verify_samples(int zero_positions, sortedSamplesList *Samples, lweInstance *lwe){

    int n = lwe->n;
    int q = lwe->q;

    for (int j = 0; j < Samples->n_categories; j++)
    {
        for (int i = 0; i < Samples->n_in_categories[j]; ++i)
        {
            int sum = 0;
            for (int k = 0; k < n; ++k)
            {
                sum = (sum + Samples->a_list[SAMPLES_PER_CATEGORY*n*j + i*n + k]*lwe->s[k]) % q;
            }
            sum = (sum + Samples->e_list[SAMPLES_PER_CATEGORY*j +i]) % q;
            if (sum != Samples->z_list[SAMPLES_PER_CATEGORY*j +i])
            {
                printf("ERROR in sample verification\n");
                exit(0);
            }
        }
        
    }

    return 0;

}


int verify_unsorted_samples(int zero_positions, unsortedSamplesList *Samples, lweInstance *lwe){

    int n = lwe->n;
    int q = lwe->q;

    for (int i = 0; i < Samples->n_samples; ++i)
    {
        int sum = 0;
        for (int k = 0; k < n; ++k)
        {
            sum = (sum + Samples->a_list[i*n + k]*lwe->s[k]) % q;
        }
        sum = (sum + Samples->e_list[i]) % q;
        if (sum != Samples->z_list[i])
        {
            printf("ERROR in unsorted sample verification\n");
            exit(0);
        }
    }

}



