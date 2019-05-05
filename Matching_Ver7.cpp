#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <sys/time.h>
#include <vector>
#include "mpi.h"

// #include "Matchcommon.h"
#define Radius 1
#define NPixel 9
#define Deltat 0.00001


int find_option( int argc, char **argv, const char *option )
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

//Initializing events
void init_events( int n, int *g )
{
    srand48( time( NULL ) );

    for (int i = 0; i <  n; i++) 
        for (int j = 0; j < n; j++) 
            g[i*n+j] = drand48()+.5;
}

void HammingDistance(int *tmd, int *tm2d, int tmdaddr, int tm2daddr, int *temp, int N )
{
         for (int k = 0; k < NPixel; k++)
                for (int l = 0; l < NPixel; l++)
                    temp += int (tmd[tmdaddr+(k*N+l)] ^ tm2d[tm2daddr+(k*N+l)]);
}

void OpticalFlow(int *tmd, int *tm2d, int tmdaddr, int tm2daddr, int *velocity, int N )
{
         for (int k = 0; k < NPixel; k++)
                for (int l = 0; l < NPixel; l++)
                    velocity[tmdaddr+(k*N+l)] = int (tmd[tmdaddr+(k*N+l)] - tm2d[tm2daddr+(k*N+l)]);
}


int main(int argc, char **argv)
{
int N = read_int( argc, argv, "-n", 81 ); // the number of particles

int *tmd = (int*) malloc( N * N * sizeof(int) );
int *tm2d = (int*) malloc( N * N * sizeof(int) );
int *Targetind = (int*) malloc( (N/NPixel) * (N/NPixel) * sizeof(int) );
int *targetvalue = (int*) malloc( (N/NPixel) * (N/NPixel) * sizeof(int) );
int *velocity = (int*) malloc( N * N * sizeof(int) );

    //
    //  set up MPI
    //
int n_proc, rank;
MPI_Init( &argc, &argv );
MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // //
    // //  allocate generic resources
    // //
    // FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    // FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;

// int *tmd = new int[N * N];
// int *tm2d = new int[N * N];
// int *Targetind = new int[(N/NPixel) * (N/NPixel)];
// int *targetvalue = new int[(N/NPixel) * (N/NPixel)];
// int *velocity = new int[N * N];

// MPI_Datatype int;
// MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
// MPI_Type_commit( &PARTICLE );

if (rank==0){
    init_events(N, tmd);
    init_events(N, tm2d);
}
    MPI_Bcast(tmd, N * N, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(tm2d, N * N, MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(Targetind, (N/NPixel) * (N/NPixel), MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(targetvalue, (N/NPixel) * (N/NPixel), MPI_INT, 0, MPI_COMM_WORLD);
    // MPI_Bcast(velocity, N * N, MPI_INT, 0, MPI_COMM_WORLD);

    // delete[] tmd;
    // delete[] tm2d;
    // delete[] Targetind;
    // delete[] targetvalue;
    // delete[] velocity;
    
    tmd = NULL;
    tm2d = NULL;

double simulation_time = read_timer( );
////////////////////////////////////////////////////////////////////////////////     
////////////////////////////  Target Defining      ////////////////////////////
////////////////////////////////////////////////////////////////////////////////

        for (int i = rank; i < (N/NPixel); i+=n_proc)
        {
            for (int j = 0; j < (N/NPixel); j++)
            {
            // count = ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel;
            Targetind[i*(N/NPixel) + j] = i*NPixel*N + j*NPixel; // Assign block in t-2d as target
            // printf("count= %d, ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel= %d\n", count, ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel);
            HammingDistance(tmd, tm2d, i*NPixel*N + j*NPixel, i*NPixel*N + j*NPixel, &targetvalue[i*(N/NPixel) + j], N );            
            }
        }

////////////////////////////////////////////////////////////////////////////////     
////////////////////////////  Optical Matching      ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
         for (int i = rank; i < (N/NPixel); i+=n_proc)
        {
            for (int j = 0; j < (N/NPixel); j++)
            {
                // count = ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel;
                // int *temp = (int*) calloc((2*Radius +1)*(2*Radius +1), sizeof(int));
                // Targetind[i*NBlocks+j] = (i * NBlocks + j);
                for (int dx = -Radius; dx <= Radius; dx++)   //Search over nearby 8 blocks and the target block 
                {
                    for (int dy = -Radius; dy <= Radius; dy++)
                    {
                        if (i*NPixel + dx*NPixel >= 0 && i*NPixel + dx*NPixel < N && j*NPixel + dy*NPixel >= 0 && j*NPixel + dy*NPixel < N)
                        {
                            int temp = 0;
                            // Calculating the hamming distance
                            HammingDistance(tmd, tm2d, i*NPixel*N + j*NPixel, (i*NPixel+dx*Radius)*N + (j*NPixel+dy*Radius), &temp, N);
                            // for (int k = 0; k < NPixel; k++)
                            //     for (int l = 0; l < NPixel; l++)
                            //         temp = temp + int (tmd[(i+k)*N+(j+l)] xor tm2d[(i+k+dx*Radius)*N+(j+l+dy*Radius)]);
                            // printf("temp is %d\n",temp);
                            if (temp < targetvalue[i*(N/NPixel) + j])
                            {
                                targetvalue[i*(N/NPixel) + j] = temp;
                                Targetind[i*(N/NPixel) + j] = (i*NPixel + dx*NPixel)*N + (j*NPixel + dy*NPixel);
                                // int Targetind[count] = 0;
                            }
                        }
                    }
                }
            }
        }
///////////////////////////////////////////////////////////////////////////////////     
////////////////////////////////////  Optical Flow Calculation ////////////////////
///////////////////////////////////////////////////////////////////////////////////
        for (int i = rank; i < (N/NPixel); i+=n_proc)
        {
            for (int j = 0; j < (N/NPixel); j++)
            {
                // count = ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel;
// Calculating the Optical Flow
                OpticalFlow(tmd, tm2d, i*NPixel*N + j*NPixel, Targetind[i*(N/NPixel) + j], velocity, N);                
                // for (int k = 0; k < NPixel; k++)
                //     for (int l = 0; l < NPixel; l++)
                //         OpticalFlow[(i+k)*N+(j+l)] = tmd[(i+k)*N+(j+l)] - tm2d[Targetind[count]+(k*N+l)];
            }
        }

   simulation_time = read_timer( ) - simulation_time;
   printf( "N = %d, OpenMP simulation time = %g seconds\n", N, simulation_time);


// int i, j, count = 0; 
// int *B[N];
//     for (i=0; i<M; i++)
//          B[i] = (int *)malloc(N * sizeof(int));

    // for (int i = 0; i <  N; i++) 
    //   for (int j = 0; j < N; j++)
    //      printf("tmd in position %d*%d is:%d\n",i, j, tmd[i*N+j]);

free( tmd );
free( tm2d );
free( velocity );
MPI_Finalize();

return 0;
}
//run using 
//mpicc Matching_Ver7.cpp -o Matching_Ver7 -lm
//mpic++ Matching_Ver7.cpp -o Matching_Ver7
//mpirun -np 100 ./Matching_Ver7
