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

#ifndef SOLVE_FWHT_H
#define SOLVE_FWHT_H

#include "lwe_instance.h"

#define MAX_FWHT 35

#define MAX_BRUTE_FORCE 60

int solve_fwht_search(u8 *binary_solution, int zeroPositions, int fwht_positions, unsortedSamplesList *reducedSamples, lweInstance *lwe);

int solve_fwht_search_bruteforce(u8 *binary_solution, short *bf_solution, int zero_positions, int bf_positions, int fwht_positions, unsortedSamplesList *reducedSamples, lweInstance *lwe);

#endif /* SOLVE_FWHT_H */
