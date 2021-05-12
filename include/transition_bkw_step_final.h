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

#ifndef TRANSITION_BKW_STEP_FINAL_H
#define TRANSITION_BKW_STEP_FINAL_H

#include "lwe_instance.h"
#include "position_values_2_category_index.h"
#include "utils.h"

int transition_bkw_step_final(lweInstance *lwe, bkwStepParameters *srcBkwStepPar, sortedSamplesList *srcSamples, unsortedSamplesList *dstSamples, u64 maxSamples);

#endif /* TRANSITION_BKW_STEP_FINAL_H */
