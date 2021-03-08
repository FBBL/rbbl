/**
	utils.h
	
	Authors:
		Alessandro Budroni, May 2019
*/

#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>

#include "config.h"

void ASSERT(int condition, const char* string, ...);

void time_stamp (const char* string, ...);

int get_seed();

#endif