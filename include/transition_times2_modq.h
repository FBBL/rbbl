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

#ifndef SRC_FILE_BASED_LWE_MOD2ATTACK_UTILS_H_
#define SRC_FILE_BASED_LWE_MOD2ATTACK_UTILS_H_

#include "position_values_2_category_index.h"
#include "utils.h"
#include "lwe_instance.h"
#include <time.h>

int transition_times2_modq(lweInstance *lwe, bkwStepParameters *bkwStepPar, sortedSamplesList *sortedSamplesList, samplesList* unsortedSampleList);

#endif /* SRC_FILE_BASED_LWE_MOD2ATTACK_UTILS_H_ */
