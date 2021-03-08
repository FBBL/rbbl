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

#ifndef TRANSITION_BKW_STEP_SMOOTH_LMS_H
#define TRANSITION_BKW_STEP_SMOOTH_LMS_H

#include "lwe_instance.h"
#include "position_values_2_category_index.h"
#include "utils.h"
#include "lookup_tables.h"

int transition_bkw_step_smooth_lms(lweInstance *lwe, bkwStepParameters *srcBkwStepPar, bkwStepParameters *dstBkwStepPar, sortedSamplesList *srcSamples, sortedSamplesList *dstSamples);

#endif /* TRANSITION_BKW_STEP_SMOOTH_LMS_H */
