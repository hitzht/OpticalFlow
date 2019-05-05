#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include "Matchcommon.h"
using namespace std;

double size;
//
//  timer
//


double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}
//
//  keep density constant
//


void init_events( int n, int *g )
{
    srand48( time( NULL ) );

    for (int i = 0; i <  n; i++) 
        for (int j = 0; j < n; j++) 
            g[i*N+j] = drand48();
}

//
//  interact two particles
//

//  
//  I/O routines
//


//
//  command line option processing
//
