/**
	utils.h
	
	Authors:
		Alessandro Budroni, May 2019
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

#endif