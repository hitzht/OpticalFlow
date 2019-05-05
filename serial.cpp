#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <vector>
#include <time.h>
#include <sys/time.h>
#include "common.h"
using namespace std;

double binSize, gridSize;
int binNum;
double size;

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
double set_size( int n )
{
    size = sqrt( density * n );
    return size;
}
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;

    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;

    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
//void apply_force( particle_t &particle, particle_t &neighbor )
{
// dmin: the distance of the nearest particle 
// dave: average distance of all the particles
// navg: total number of particles that are averaged
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
    if (r2 != 0)
        {
       if (r2/(cutoff*cutoff) < *dmin * (*dmin))
          *dmin = sqrt(r2)/cutoff;
           (*davg) += sqrt(r2)/cutoff;
           (*navg) ++;
        }
        
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );
    
    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option ) // compare input to the option that 
// we give in the program and return a value greater or equal to 1 when match and a value equal to
// -1 when fails to match
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

 int read_int( int argc, char **argv, const char *option, int default_value ) // 
{
    int iplace = find_option( argc, argv, option ); 
    if( iplace >= 0 && iplace < argc-1 ) // 
        return atoi( argv[iplace+1] ); // Ascii to integer conversion
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}


void buildBins(vector<vector<particle_t> >& bins, particle_t* particles, int n, double gridSize)
{
    // gridSize = sqrt(n*density);
    int gs = int(gridSize/cutoff) +1;
    binSize = cutoff * 2;  
    binNum = int(gridSize / binSize)+1; // Should be around sqrt(N/2)

    // printf("Grid Size: %.4lf\n",gridSize);
    // printf("Bin Size: %.2lf\n",binSize);
    // printf("number of grids: %d\n",gs);
    // printf("Number of Bins: %d*%d\n",binNum,binNum);
    // // Increase\Decrease binNum to be something like 2^k?

    bins.resize(binNum * binNum);

    for (int i = 0; i < n; i++)
    {
        int x = int(particles[i].x / binSize);
        int y = int(particles[i].y / binSize);
        // printf("x is %d while particles[i].x is %f\n",x,particles[i].x);
        // printf("y is %d while particles[i].y is %f\n",y,particles[i].y);
        //this adds a new element to a vector
        bins[x*binNum + y].push_back(particles[i]);
        // printf("size of bin: %lx\n",bins.size());
    }
}

 // benchmarking program

int main( int argc, char **argv )
{    
    int navg,nabsavg=0;
    double davg, dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int( argc, argv, "-n", 10 ); // the number of particles
    char *savename = read_string( argc, argv, "-o", NULL ); 
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL; // fopen binds a name resource 
    // specified by a filename to a stream and 
    //(w) place the file pointer at the (beginning) of the file and truncate the file to zero length
    //  If the file does not exist, attempt to create it. 
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;
    //(a) place the file pointer at the (end) of the file. If the file does not exist, 
    // attempt to create it. writes are always appended

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    double gridSize = set_size(n);

    init_particles( n, particles );


    vector<vector<particle_t> > particle_bins;
    vector<particle_t> temp;
    buildBins(particle_bins, particles, n, gridSize);


//  Initialize the particle positions and velocities
//
//
//  simulate a number of time steps
//
   double simulation_time = read_timer( );

//////////////////////////////////////////////////////////////////////////////
////////////// The Heart of the optimized serial code is situated here ///////
//////////////////////////////////////////////////////////////////////////////
    for (int step = 0; step < NSTEPS; step++ )
    {
        // printf("particles[5].x is %f, and particles[5].y is %f\n", particles[5].x, particles[5].y);

        navg = 0;
        davg = 0.0;
        dmin = 1.0;
        for (int i = 0; i < binNum; i++)
        {
            for (int j = 0; j < binNum; j++)
            {
                vector<particle_t>& vec = particle_bins[i*binNum+j];
                // printf("vec.size() is: %lx, and size of particle_bins[i*binNum+j] is:%lx\n",vec.size(), sizeof(particle_bins[i*binNum+j]));
                for (int k = 0; k < vec.size(); k++)
                    vec[k].ax = vec[k].ay = 0;
                for (int dx = -1; dx <= 1; dx++)   //Search over nearby 8 bins and itself
                {
                    for (int dy = -1; dy <= 1; dy++)
                    {
                        if (i + dx >= 0 && i + dx < binNum && j + dy >= 0 && j + dy < binNum)
                        {
                            vector<particle_t> &vec2 = particle_bins[(i+dx) * binNum + j + dy];
                            // Comparing two Adjacent Vectors
                            for (int k = 0; k < vec.size(); k++)
                                for (int l = 0; l < vec2.size(); l++)
                                    apply_force(vec[k], vec2[l], &dmin, &davg, &navg);
                        }
                    }
                }
            }
        }

        for (int i = 0; i < binNum; i++)
        {
            for(int j = 0; j < binNum; j++)
            {
                vector<particle_t>& vec = particle_bins[i * binNum + j];
                int tail = vec.size(), k = 0;
                for(; k < tail; )
                {
                    move( vec[k] );
                    int x = int(vec[k].x / binSize);  //Check the position
                    // printf("vec[k].x is %f, binSize is %d, x=%d\n",vec[k].x,binSize,x);
                    int y = int(vec[k].y / binSize);
                    if (x == i && y == j)  // Still inside original bin
                        k++;
                    else
                    {
                        temp.push_back(vec[k]);  // Store paricles that have changed bin. 
                        vec[k] = vec[--tail]; //Remove it from the current bin.
                     // printf("k is %d, tail is %d\n",k,tail);
                    }
                }
                vec.resize(k);
            }
        }

        for (int i = 0; i < temp.size(); i++)  // Put them into the new bin 
        {
            int x = int(temp[i].x / binSize);
            int y = int(temp[i].y / binSize);
            particle_bins[x*binNum+y].push_back(temp[i]);
        }
        temp.clear();
        
        if( find_option( argc, argv, "-no" ) == -1 )
        {
            if (navg) 
            {
                absavg +=  davg/navg;
                nabsavg++;
            }
              
            if (dmin < absmin) 
                absmin = dmin;

            if( fsave && (step%SAVEFREQ) == 0 )
                save( fsave, n, particles );
        }
        // printf("particles[5].x is %f, and particles[5].y is %f\n", particles[5].x, particles[5].y);
    }

//////////////////////////////////////////////////////////////////////////////
////////////// The Heart of the optimized serial code is situated here ///////
//////////////////////////////////////////////////////////////////////////////
   
    // for( int step = 0; step < NSTEPS; step++ )
    // {
    // navg = 0;
    // davg = 0.0;
    // dmin = 1.0;
    //     //
    //     //  compute forces
    //     //
    //     for( int i = 0; i < n; i++ )
    //     {
    //         // Reset acceleration
    //         particles[i].ax = particles[i].ay = 0;
    //          for (int j = 0; j < n; j++ )
    //          apply_force( particles[i], particles[j],&dmin,&davg,&navg);
    //     }
    //     //
    //     //  move particles
    //     //
    //     for( int i = 0; i < n; i++ )
    //         move( particles[i] );        

    //     if( find_option( argc, argv, "-no" ) == -1 )
    //     {
    //       //
    //       // Computing statistical data
    //       //
    //       if (navg) {
    //         absavg +=  davg/navg;
    //         nabsavg++;
    //       }
    //       if (dmin < absmin) absmin = dmin;
    //       //
    //       //  save if necessary
    //       //
    //       if( fsave && (step%SAVEFREQ) == 0 )
    //           save( fsave, n, particles );
    //     }
    // }
//////////////////////////////////////////////////////////////////////
///////////////////////////   The End ////////////////////////////////
//////////////////////////////////////////////////////////////////////

    simulation_time = read_timer( ) - simulation_time;
    
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
    if (nabsavg) absavg /= nabsavg;
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    return 0;
}