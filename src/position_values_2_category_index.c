/*  This file is part of FBBL (File-Based BKW for LWE).
 *
 *  FBBL is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  FBBL is distributed in the hope that it will be useful,
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
static u64 positionValuesToCategoryGeneralized(int n, short *t, int *c)
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
            else if (2*t[n-1] < c[n-1])
            {
                index = (2*t[n-1]-1);
                for (int i = 0; i<n-1; i++)
                    index *= c[i];
                index += 2*positionValuesToCategoryGeneralized(n-1, t, c);
            }
            else
            {
                index = 2*(c[n-1]-t[n-1])-1;
                short nt[n-1];
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
                index += 1 + 2*positionValuesToCategoryGeneralized(n-1, nt, c);
            }
        }

    }
    else     // c even  case
    {

        if (n == 1)    // single position case
        {
            if (2*t[0] < c[0])
            {
                return 2*t[0];
            }
            else
            {
                return 2*(c[0]-t[0])-1;
            }
        }
        else     // other recursive cases
        {
            if (t[n-1] < (c[n-1]/2))
            {
                index = 2*t[n-1];
                for (int i = 0; i<n-1; i++)
                    index *= c[i];
                index += 2*positionValuesToCategoryGeneralized(n-1, t, c);
            }
            else
            {
                index = 2*(c[n-1]-t[n-1]-1);
                short nt[n-1];
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
                index += 1 + 2*positionValuesToCategoryGeneralized(n-1, nt, c);
            }
        }

    }
    return index;
}

/* makes an LMS mapping of pi using q and p as moduli - version for smooth LMS */
/* here q_ must be ceil(lwe.q/2) if the position is not already reduced, otherwise it must be 2*p from previous step */
short positionSmoothLMSMap(short pn, int q, int q_, int p, int c)
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
    int q = lwe->q;
    int Ni = dstBkwStepPar->numPositions;
    short p = dstBkwStepPar->p;
    short p1 = dstBkwStepPar->p1;
    short p2;
    u64 index_cat = 0;

    short t[MAX_SMOOTH_LMS_POSITIONS+1];
    int c[MAX_SMOOTH_LMS_POSITIONS+1];
    int q_;

    /* Differentiate the first step for general steps */
    if (dstBkwStepPar->prev_p1 == -1)   // first step
    {
        q_ = q%2 == 1 ? (q+1)/2 : q/2;
        for(int i = 0; i < Ni; i++)
        {
            c[i] = ((2*q_-1) % p) == 0 ? ((2*q_-1) / p) : ((2*q_-1) / p) + 1;
            t[i] = positionSmoothLMSMap(pn[i], q, q_, p, c[i]);
        }
        c[Ni] = ((2*q_-1) % p1) == 0 ? ((2*q_-1) / p1) : ((2*q_-1) / p1) + 1;
        t[Ni] = positionSmoothLMSMap(pn[Ni], q, q_, p1, c[Ni]);
            // printf("pn %d %d %d\n", pn[0], pn[1], pn[2] );
            // printf("t %d %d %d\n", t[0], t[1], t[2] );
            // printf("c %d %d %d\n", c[0], c[1], c[2] );


        index_cat = positionValuesToCategoryGeneralized(Ni+1, t, c);
    }
    else if (dstBkwStepPar->startIndex + Ni == lwe->n)     // last step
    {
        p2 = dstBkwStepPar->p2;
        q_= dstBkwStepPar->prev_p1;
        c[0] = ((2*q_-1) % p2) == 0 ? ((2*q_-1) / p2) : ((2*q_-1) / p2) + 1;
        t[0] = positionSmoothLMSMap(pn[0], q, q_, p2, c[0]);
        q_ = q%2 == 1 ? (q+1)/2 : q/2;
        for(int i = 1; i < Ni; i++)
        {
            c[i] = ((2*q_-1) % p) == 0 ? ((2*q_-1) / p) : ((2*q_-1) / p) + 1;
            t[i] = positionSmoothLMSMap(pn[i], q, q_, p, c[i]);
        }
        index_cat = positionValuesToCategoryGeneralized(Ni, t, c);
    }
    else      // middle steps
    {
        p2 = dstBkwStepPar->p2;
        q_= dstBkwStepPar->prev_p1;
        c[0] = ((2*q_-1) % p2) == 0 ? ((2*q_-1) / p2) : ((2*q_-1) / p2) + 1;
        t[0] = positionSmoothLMSMap(pn[0], q, q_, p2, c[0]);
        q_ = q%2 == 1 ? (q+1)/2 : q/2;
        for(int i = 1; i < Ni; i++)
        {
            c[i] = ((2*q_-1) % p) == 0 ? ((2*q_-1) / p) : ((2*q_-1) / p) + 1;
            t[i] = positionSmoothLMSMap(pn[i], q, q_, p, c[i]);
        }
        c[Ni] = ((2*q_-1) % p1) == 0 ? ((2*q_-1) / p1) : ((2*q_-1) / p1) + 1;
        t[Ni] = positionSmoothLMSMap(pn[Ni], q, q_, p1, c[Ni]);
        index_cat = positionValuesToCategoryGeneralized(Ni+1, t, c);
    }
    return index_cat;
}
