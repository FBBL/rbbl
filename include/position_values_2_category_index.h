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

#ifndef POSITION_VALUES_2_CATEGORY_INDEX
#define POSITION_VALUES_2_CATEGORY_INDEX

#include "lwe_instance.h"

#define MAX_SMOOTH_LMS_POSITIONS 10

/* smooth LMS */
u64 position_values_2_category_index(lweInstance *lwe, bkwStepParameters *dstBkwStepPar, u16 *pn);

#endif

