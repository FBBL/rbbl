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

#include "position_values_2_category_index.h"

#include <math.h>

/* Smooth LMS */

int is_smooth_lms_singleton(u64 categoryIndex, u64 numCategories)
{
    if(numCategories & 1) // whether it is odd
    {
        return (categoryIndex == 0);
    }
    return 0;
}

// n is the length, t is a list containing the indices, p is a list containing the reduction factor for each position.
// efficiency may be improved with a look-up table
static u64 positionValuesToCategoryGeneralized(int n, int *t, int *c)
{

    u64 index;

    if (c[n-1] & 1)   // c odd case
    {

        if (n == 1)   // single position case
        {
            if (t[0] == 0)
                return 0;
            else if (2*t[0] < c[0])
                return 2*t[0] -1;
            else
                return 2*(c[0]-t[0]);

        }
        else     // other recursive cases
        {
            if (t[n-1] == 0)
            {
                return positionValuesToCategoryGeneralized(n-1, t, c);
            }
            else if ((t[n-1] << 1) < c[n-1])
            {
                index = ((t[n-1] << 1)-1);
                for (int i = 0; i<n-1; i++)
                    index *= c[i];
                index += (positionValuesToCategoryGeneralized(n-1, t, c) << 1);
            }
            else
            {
                index = ((c[n-1]-t[n-1]) << 1) -1;
                int nt[n-1];
                for (int i = 0; i<n-1; i++)
                {
                    if (c[i] & 1)   /* Category i is odd */
                    {
                        nt[i] = (c[i] - t[i])%c[i];
                    }
                    else     /* Category i is even */
                    {
                        nt[i] = (c[i] - t[i]-1);
                    }
                    index *= c[i];
                }
                index += 1 + (positionValuesToCategoryGeneralized(n-1, nt, c) << 1);
            }
        }

    }
    else     // c even  case
    {

        if (n == 1)    // single position case
        {
            if ((t[0] << 1) < c[0])
            {
                return (t[0] << 1);
            }
            else
            {
                return ((c[0]-t[0]) << 1)-1;
            }
        }
        else     // other recursive cases
        {
            if ((t[n-1] << 1) < (c[n-1]))
            {
                index = t[n-1] << 1;
                for (int i = 0; i<n-1; i++)
                    index *= c[i];
                index += (positionValuesToCategoryGeneralized(n-1, t, c) << 1);
            }
            else
            {
                index = (c[n-1]-t[n-1]-1) << 1;
                int nt[n-1];
                for (int i = 0; i<n-1; i++)
                {
                    if (c[i] & 1)   /* Category i is odd */
                    {
                        nt[i] = (c[i] - t[i])%c[i];
                    }
                    else     /* Category i is even */
                    {
                        nt[i] = (c[i] - t[i]-1);
                    }
                    index *= c[i];
                }
                index += 1 + (positionValuesToCategoryGeneralized(n-1, nt, c) << 1);
            }
        }

    }
    return index;
}

/* makes an LMS mapping of pi using q and p as moduli - version for smooth LMS */
/* here q_ must be ceil(lwe.q/2) if the position is not already reduced, otherwise it must be 2*p from previous step */
short positionSmoothLMSMap(u16 pn, u16 q, u16 q_, u16 p, int c)
{

    int Delta;

    if (c & 1)   // if c is odd
    {
        Delta = p*((c/2) +1) -q_;
        if (pn < q_)
            return (pn + Delta) / p;
        else
            return (c -(((q - pn) + Delta) / p)) % c;
    }
    else     // if c is even
    {
        Delta = p*(c/2)-q_; // shift
        if (pn < q_)
            return (pn+Delta)/p;
        else
            return c -1 -(((q - pn) + Delta) / p);
    }
}

u64 position_values_2_category_index(lweInstance *lwe, bkwStepParameters *dstBkwStepPar, u16 *pn)
{
    u16 q = lwe->q;
    u16 Ni = dstBkwStepPar->numPositions;
    u16 p = dstBkwStepPar->p;
    u16 p1 = dstBkwStepPar->p1;
    u16 p2;
    u64 index_cat = 0;
    u16 t2q_1;

    int t[MAX_SMOOTH_LMS_POSITIONS+1];
    int c[MAX_SMOOTH_LMS_POSITIONS+1];
    u16 q_;

    /* Differentiate the first step for general steps */
    if (dstBkwStepPar->prev_p1 == -1)   // first step
    {
        q_ = q%2 == 1 ? (q+1)/2 : q/2;
        t2q_1 = 2*q_-1;
        for(int i = 0; i < Ni; i++)
        {
            c[i] = ((t2q_1) % p) == 0 ? ((t2q_1) / p) : ((t2q_1) / p) + 1;
            t[i] = positionSmoothLMSMap(pn[i], q, q_, p, c[i]);
        }
        c[Ni] = ((t2q_1) % p1) == 0 ? ((t2q_1) / p1) : ((t2q_1) / p1) + 1;
        t[Ni] = positionSmoothLMSMap(pn[Ni], q, q_, p1, c[Ni]);
        index_cat = positionValuesToCategoryGeneralized(Ni+1, t, c);
    }
    else if (dstBkwStepPar->startIndex + Ni == lwe->n)     // last step
    {
        p2 = dstBkwStepPar->p2;
        q_= dstBkwStepPar->prev_p1;
        t2q_1 = 2*q_-1;
        c[0] = ((t2q_1) % p2) == 0 ? ((t2q_1) / p2) : ((t2q_1) / p2) + 1;
        t[0] = positionSmoothLMSMap(pn[0], q, q_, p2, c[0]);
        q_ = q%2 == 1 ? (q+1)/2 : q/2;
        t2q_1 = 2*q_-1;
        for(int i = 1; i < Ni; i++)
        {
            c[i] = ((t2q_1) % p) == 0 ? ((t2q_1) / p) : ((t2q_1) / p) + 1;
            t[i] = positionSmoothLMSMap(pn[i], q, q_, p, c[i]);
        }
        index_cat = positionValuesToCategoryGeneralized(Ni, t, c);
    }
    else      // middle steps
    {
        p2 = dstBkwStepPar->p2;
        q_= dstBkwStepPar->prev_p1;
        t2q_1 = 2*q_-1;
        c[0] = ((t2q_1) % p2) == 0 ? ((t2q_1) / p2) : ((t2q_1) / p2) + 1;
        t[0] = positionSmoothLMSMap(pn[0], q, q_, p2, c[0]);
        q_ = q%2 == 1 ? (q+1)/2 : q/2;
        t2q_1 = 2*q_-1;
        for(int i = 1; i < Ni; i++)
        {
            c[i] = ((t2q_1) % p) == 0 ? ((t2q_1) / p) : ((t2q_1) / p) + 1;
            t[i] = positionSmoothLMSMap(pn[i], q, q_, p, c[i]);
        }
        c[Ni] = ((t2q_1) % p1) == 0 ? ((t2q_1) / p1) : ((t2q_1) / p1) + 1;
        t[Ni] = positionSmoothLMSMap(pn[Ni], q, q_, p1, c[Ni]);
        index_cat = positionValuesToCategoryGeneralized(Ni+1, t, c);
    }
    return index_cat;
}
