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

#ifndef UTILS_H
#define UTILS_H

#include "utils.h"
#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdarg.h>
#include <unistd.h>
#include <fcntl.h>

/* if condition is false, then print error message and exit */
void ASSERT(int condition, const char* string, ...)
{

    if (!condition)
    {
        va_list args;
        va_start(args, string);
        char err_msg[200];
        vsprintf(err_msg, string, args);
        printf("ERROR: %s\n", err_msg);
        va_end(args);
        exit(-1);
    }
}

/*
	time_stamp: make a time stamp
	INPUT
	- string: string to print
*/
void time_stamp (const char* string, ...)
{
    time_t my_time = time(NULL);
    va_list args;
    va_start(args, string);
    char msg[256];
    vsprintf(msg, string, args);
    printf("%.24s - %s\n", ctime(&my_time), msg);
    va_end(args);
}


/* get seed for RNG from /dev/urandom */
int get_seed()
{
    // get seed for RNG
    int randomData = open("/dev/urandom", O_RDONLY);
    ASSERT(randomData  > 0, "could not open /dev/urandom in get_seed");
    int seed, ret;
    ret = read(randomData, &seed, sizeof(int));
    ASSERT(ret > 0, "could not read from /dev/urandom in get_seed");
    close(randomData);

    return seed;
}

// /* return 1 if the given array is all zero */
// int checkzero(char *string, int length) {
//     int is_zero;
//     __asm__ (
//         "cld\n"
//         "xorb %%al, %%al\n"
//         "repz scasb\n"
//         : "=c" (is_zero)
//         : "c" (length), "D" (string)
//         : "eax", "cc"
//     );
//     return !is_zero;
// }

/* return 1 if the given array is all zero */
int checkzero(char *string, int length) {
    for (int i = 0; i < length; ++i)
    {
        if (string!= 0)
            return 0;
    }
    return 1;
}


#endif