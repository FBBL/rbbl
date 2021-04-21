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
int error_rate(int zero_positions, samplesList *Samples, lweInstance *lwe)
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
        for (int j = zero_positions; j < n; j++){
            bina = Samples->list[i].a[j] < q/2 ? Samples->list[i].a[j] % 2 : (Samples->list[i].a[j] +1) % 2;
            sum = (sum + bina*bin_secret[j]) % 2;
        }
        binz = Samples->list[i].z < q/2 ? Samples->list[i].z % 2 : (Samples->list[i].z +1) % 2;
        errors += sum == binz ? 0 : 1;
    }

    time_stamp("Error rate: %lf", (float)errors/(float)Samples->n_samples );

    return 0;
}

