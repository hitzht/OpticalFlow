#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>
#include <algorithm>
#include "common.h"
#include <time.h>
#include <sys/time.h>

#define NUM_THREADS 256
//#define _cutoff 0.01    //Value copied from common.cpp
//#define _density 0.0005
#define MAXITEM 4 //Assume at most MAXITEM particles in one bin. Change depends on binSize
#define CUTOFF_SCALE 2 

double gridsize;

double binSize;
int binNum;

//
//  benchmarking program
//
double size;

//
//  tuned constants
//

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
double set_size( int n )
{
    size = sqrt( density * n );
    return size;
}

//
//  Initialize the particle positions and velocities
//
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
        //cd
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

void apply_force2( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
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
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor )
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
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
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}

__device__ void apply_force_gpu(particle_t &particle, particle_t &neighbor)
{
  double dx = neighbor.x - particle.x;
  double dy = neighbor.y - particle.y;
  double r2 = dx * dx + dy * dy;
  if( r2 > cutoff*cutoff )
      return;
  //r2 = fmax( r2, min_r*min_r );
  r2 = (r2 > min_r*min_r) ? r2 : min_r*min_r;
  double r = sqrt( r2 );

  //
  //  very simple short-range repulsive force
  //
  double coef = ( 1 - cutoff / r ) / r2 / mass;
  particle.ax += coef * dx;
  particle.ay += coef * dy;

}

__constant__ const int dir[9][2]={{0,0},{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1},{-1,0}};

__global__ void compute_forces_gpu(particle_t * particles, int*cnt,int n,double binSize,int binNum)
{
  // Get thread (particle) ID
  
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = gridDim.x*blockDim.x;
    for(int ii=tid;ii<n;ii+=offset)
    {
        particle_t p = particles[ii];
        p.ax=p.ay=0;
        int i = int(p.x / binSize);
        int j = int(p.y / binSize);
        for(int t=0;t<9;t++)
        {
            int x = i + dir[t][0];
            int y = j + dir[t][1];
            if (x >= 0 && x < binNum && y >= 0 && y < binNum)
            {
                int id = x*binNum+y;
                int start = cnt[id-1],end = cnt[id];
                for (int k = start; k < end; k++)
                    apply_force_gpu(p, particles[k]);
            }
        }
        particles[ii].ax = p.ax;
        particles[ii].ay = p.ay;
    }
}

__global__ void move_gpu (particle_t * __restrict__ particles, int n, double gridsize)
{

  // Get thread (particle) ID
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = gridDim.x*blockDim.x;
    for(int i=tid;i<n;i+=offset)
    {

        particle_t * p = &particles[i];
        //
        //  slightly simplified Velocity Verlet integration
        //  conserves energy better than explicit Euler method
        //
        p->vx += p->ax * dt;
        p->vy += p->ay * dt;
        p->x  += p->vx * dt;
        p->y  += p->vy * dt;
    
        //
        //  bounce from walls
        //
        while( p->x < 0 || p->x > gridsize )
        {
            p->x  = p->x < 0 ? -(p->x) : 2*gridsize-p->x;
            p->vx = -(p->vx);
        }
        while( p->y < 0 || p->y > gridsize )
        {
            p->y  = p->y < 0 ? -(p->y) : 2*gridsize-p->y;
            p->vy = -(p->vy);
        }
    }

}

// __global__ void move_gpu2 (particle_t * particles, int n, double gridsize)
// {

//   // Get thread (particle) ID
//   int tid = threadIdx.x + blockIdx.x * blockDim.x;
//   if(tid >= n) return;

//   //particle_t  p = particles[tid];
//     //
//     //  slightly simplified Velocity Verlet integration
//     //  conserves energy better than explicit Euler method
//     //
//     double vx = particles[tid].vx;
//     double ax = particles[tid].ax;
//     double vy = particles[tid].vy;
//     double ay = particles[tid].ay;
//     vx += ax * dt;
//     vy += ay * dt;
//     double x = particles[tid].x;
//     double y = particles[tid].y;
//     x  += vx * dt;
//     y  += vy * dt;
//     //
//     //  bounce from walls
//     //
//     while( x < 0 || x > size )
//     {
//         x  = x < 0 ? -(x) : 2*size-x;
//         vx = -(vx);
//     }
//     while( y < 0 || y > size )
//     {
//         y  = y < 0 ? -(y) : 2*size-y;
//         vy = -(vy);
//     }
//     /*particles[tid].x = x;
//     particles[tid].y = y;
//     particles[tid].vx = vx;
//     particles[tid].vy = vy;*/
// }


__global__ void getCount(particle_t* particles, int* count,int n,double binSize,int binNum)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = gridDim.x*blockDim.x;
    for(int i=tid;i<n;i+=offset)
    {
        int x = int(particles[i].x / binSize);
        int y = int(particles[i].y / binSize);
        atomicAdd(count+x*binNum+y,1);
    }
}

__global__ void buildBins(particle_t* particles,particle_t* tmp,int* count,int n,double binSize,int binNum)
{
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    int offset = gridDim.x*blockDim.x;
    for(int i=tid;i<n;i+=offset)
    {
        int x = int(particles[i].x / binSize);
        int y = int(particles[i].y / binSize);
        int id = atomicSub(count+x*binNum+y,1);
        tmp[id-1] = particles[i];
    }
}

int main( int argc, char **argv )
{    
    // This takes a few seconds to initialize the runtime
    cudaDeviceSynchronize(); 

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        return 0;
    }
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    // GPU particle data structure
    particle_t * d_particles,*tmp;
    cudaMalloc((void **) &d_particles, n * sizeof(particle_t));
    cudaMalloc((void **) &tmp, n * sizeof(particle_t));
    
    double gridsize = set_size( n );
    binSize = cutoff*CUTOFF_SCALE;  
    binNum = int(gridsize / binSize)+1; // Should be around sqrt(N/2)
    // printf("Grid Size: %.4lf\n",size);
    // printf("Number of Bins: %d*%d\n",binNum,binNum);
    // printf("Bin Size: %.2lf\n",binSize);
    
    
    int* cnt;
    cudaMalloc((void **) &cnt, (binNum*binNum+1) * sizeof(int));
    cudaMemset(cnt,0,(binNum*binNum+1)*sizeof(int));
    cnt+=1; //Add one therefore cnt[-1]==0
    int* count = (int*) malloc(binNum*binNum * sizeof(int));
      
    init_particles( n, particles );

    cudaDeviceSynchronize();
    double copy_time = read_timer( );

    // Copy the particles to the GPU
    cudaMemcpy(d_particles, particles, n * sizeof(particle_t), cudaMemcpyHostToDevice);

    cudaDeviceSynchronize();
    copy_time = read_timer( ) - copy_time;
    
    //
    //  simulate a number of time steps
    //
    cudaDeviceSynchronize();
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        //
        //  compute forces
        //
        int threadNum = NUM_THREADS;
        int blks = min(1024,(n + NUM_THREADS - 1) / NUM_THREADS);
      int blockNum = blks;//min(512,(n+threadNum-1)/threadNum);

        //Old methods that don't assume maxitems for each bin
        cudaMemset(cnt,0,binNum*binNum*sizeof(int));
        getCount<<<blockNum,threadNum>>>(d_particles,cnt,n,binSize,binNum);
        
        cudaMemcpy(count, cnt, binNum*binNum * sizeof(int), cudaMemcpyDeviceToHost);
        for(int i=1;i<binNum*binNum;i++)  //Prefix sum  could be faster using parallel....
            count[i]+=count[i-1];
        cudaMemcpy(cnt, count, binNum*binNum * sizeof(int), cudaMemcpyHostToDevice);
        buildBins<<<blockNum,threadNum>>>(d_particles,tmp,cnt,n,binSize,binNum);
        std::swap(d_particles,tmp);
        cudaMemcpy(cnt, count, binNum*binNum * sizeof(int), cudaMemcpyHostToDevice);
        
        //compute_forces_gpu <<< blks, NUM_THREADS >>> (d_particles, n);
      compute_forces_gpu<<<blks, NUM_THREADS>>> (d_particles,cnt,n,binSize,binNum);
        
        //
        //  move particles
        //
      
      move_gpu <<< blks, NUM_THREADS >>> (d_particles, n, gridsize);
        
        //
        //  save if necessary
        //
        if( fsave && (step%SAVEFREQ) == 0 ) {
      // Copy the particles back to the CPU
            cudaMemcpy(particles, d_particles, n * sizeof(particle_t), cudaMemcpyDeviceToHost);
            save( fsave, n, particles);
      }
    }
    cudaDeviceSynchronize();
    simulation_time = read_timer( ) - simulation_time;
    
    printf( "CPU-GPU copy time = %g seconds\n", copy_time);
    printf( "n = %d, simulation time = %g seconds\n", n, simulation_time );
    
    free( particles );
    cudaFree(d_particles);
    if( fsave )
        fclose( fsave );
    
    return 0;
}