// SHELL='bash'
// #
// # Bridges - PSC
// #
// # Intel Compilers are loaded by default
// # You will need to specifically switch to GNU Modules
// # With with `modules.sh`
// #

// CC = g++
// MPCC = mpic++
// OPENMP = -fopenmp
// CFLAGS = -O3
// LIBS =

// cuCC = nvcc
// cuCFLAGS = -std=c++11 -O3 -arch=compute_30 -code=sm_30
// cuNVCCFLAGS = -std=c++11 -O3 -arch=compute_30 -code=sm_30

// TARGETS = Matching_Ver1 Matching_Ver2 Matching_Ver3 Matching_Ver4 Matching_Ver5 Matching_Ver6 Matching_Ver7


// Matching_Ver6: Matching_Ver6.cpp
//     $(CC) -o $@ $(LIBS) $(OPENMP) Matching_Ver6.cpp

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include <omp.h>
#include <chrono>

// #include "Matchcommon.h"

#define Deltat 0.00001
#define N 128
#define NPixel 8
#define Radius 1


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
void init_events( int n, int *g1 , int *g2)
{
    srand48( time( NULL ) );

    for (int i = 0; i <  n; i++) {
        for (int j = 0; j < n; j++){ 
            g1[i*n+j] = drand48()+.5;
            g2[i*n+j] = drand48()+.5;
        }
    }
}

// void HammingDistance(int *tmd, int *tm2d, int tmdaddr, int tm2daddr, int *temp, int N, int NPixel )
// {
//          // temp=0;
//          // printf("temp is %d\n",temp);
//          for (int k = 0; k < NPixel; k++){
//                 for (int l = 0; l < NPixel; l++){
//                     // temp += int (tmd[tmdaddr+(k*N+l)] ^ tm2d[tm2daddr+(k*N+l)]);
//                     *temp += int (tmd[tmdaddr+(k*N+l)] ^ tm2d[tm2daddr+(k*N+l)]);
//                     // printf("tm2daddr is %d+ (k*N+l) is %d and tm2daddr+(k*N+l) is %d and temp is %d\n", tm2daddr, (k*N+l),tm2daddr+(k*N+l),temp);
//                 }
//          }
// }

// void OpticalFlow(int *tmd, int *tm2d, int tmdaddr, int tm2daddr, int *velocity, int N, int NPixel )
// {
//          for (int k = 0; k < NPixel; k++)
//                 for (int l = 0; l < NPixel; l++)
//                     velocity[tmdaddr+(k*N+l)] = int (tmd[tmdaddr+(k*N+l)] - tm2d[tm2daddr+(k*N+l)]);
// }


int main(int argc, char **argv)
{
int NBlocks = N/NPixel;
// int numthreads, id;
int *tmd = (int*) malloc( N * N * sizeof(int) );
int *tm2d = (int*) malloc( N * N * sizeof(int) );
int *velocity = (int*) malloc( N * N * sizeof(int) );

int *Targetind = (int*) calloc( NBlocks * NBlocks, sizeof(int) );
int *targetvalue = (int*) calloc( NBlocks * NBlocks, sizeof(int) );

int temp;

init_events(N, tmd,tm2d);

// /////////////////////////////////////////////////
// ///////////////// just for test /////////////////
// /////////////////////////////////////////////////
//         printf("tmd is a 4*4 matrix in 2*2 block:\n");
//         printf("|%d  %d | %d  %d|:\n",tmd[0],tmd[1],tmd[2],tmd[3]);
//         printf("|%d  %d | %d  %d|:\n",tmd[4],tmd[5],tmd[6],tmd[7]);
//         printf("---------------\n");
//         printf("|%d  %d | %d  %d|:\n",tmd[8],tmd[9],tmd[10],tmd[11]);
//         printf("|%d  %d | %d  %d|:\n",tmd[12],tmd[13],tmd[14],tmd[15]);

//         printf("tm2d is a 4*4 matrix in 2*2 block:\n");
//         printf("|%d  %d | %d  %d|:\n",tm2d[0],tm2d[1],tm2d[2],tm2d[3]);
//         printf("|%d  %d | %d  %d|:\n",tm2d[4],tm2d[5],tm2d[6],tm2d[7]);
//         printf("---------------\n");
//         printf("|%d  %d | %d  %d|:\n",tm2d[8],tm2d[9],tm2d[10],tm2d[11]);
//         printf("|%d  %d | %d  %d|:\n",tm2d[12],tm2d[13],tm2d[14],tm2d[15]);
// /////////////////////////////////////////////////
// /////////////////////////////////////////////////
// /////////////////////////////////////////////////

// auto begin_sim = std::chrono::high_resolution_clock::now();
////////////////////////////////////////////////////////////////////////////////     
////////////////////////////  Target Defining      ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
  #pragma omp parallel
     {
        // #pragma omp master
        // {
        //     numthreads = omp_get_num_threads();
        //     // temp.resize(numthreads);
        //     // printf("number if threads is: %d\n", numthreads);
        // }
        // // int id = omp_get_thread_num();
        temp = 0;
//         #pragma omp critical
//         {
//             numthreads = omp_get_num_threads();
// }
        #pragma omp for reduction (+:temp)
        for (int i = 0; i < NBlocks; i++)
        {
            for (int j = 0; j < NBlocks; j++)
            {
            // count = ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel;
            // Targetind[i*NBlocks + j] = i*NPixel*N + j*NPixel; // Assign block in t-2d as target
            Targetind[i*NBlocks + j] = i*NPixel*N + j*NPixel; // Assign block in t-2d as target
            // printf("Targetind[%d*NBlocks + %d] is %d\n",i,j, i*NBlocks + j);
            // printf("targetvalue[%d*NBlocks + %d] is %d\n",i,j, targetvalue[i*NBlocks + j]);
            temp=0;
            // printf("count= %d, ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel= %d\n", count, ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel);
            // HammingDistance(tmd, tm2d, i*NPixel*N + j*NPixel, i*NPixel*N + j*NPixel, &targetvalue[i*NBlocks + j], N , NPixel);
            // HammingDistance(tmd, tm2d, i*N*NPixel + j*NPixel, i*N*NPixel + j*NPixel, targetvalue[i*NBlocks + j], N );
         // printf("targetvalue[i*NBlocks + j] is %d\n",targetvalue[i*NBlocks + j]);
         for (int k = 0; k < NPixel; k++){
                for (int l = 0; l < NPixel; l++){
                    // temp += int (tmd[tmdaddr+(k*N+l)] ^ tm2d[tm2daddr+(k*N+l)]);
                    targetvalue[i*NBlocks + j] += int (tmd[i*N*NPixel + j*NPixel+(k*N+l)] ^ tm2d[i*N*NPixel + j*NPixel+(k*N+l)]);
                    // printf("tmdaddr is %d+ (k*N+l) is %d and tmdaddr+(k*NPixel+l) is %d and targetvalue[i*NBlocks + j] is %d\n", i*N*NPixel + j*NPixel, (k*N+l),i*N*NPixel + j*NPixel+(k*N+l),targetvalue[i*NBlocks + j]);
                }
         }

            // targetvalue[i*NBlocks + j] = temp;
            }
        }

// /////////////////////////////////////////////////
// ///////////////// just for test /////////////////
// /////////////////////////////////////////////////
//         #pragma omp critical
//         {
//         printf("Target Indices for 2*2 blocks:\n");
//         printf("|%d | %d|:\n",Targetind[0],Targetind[1]);
//         printf("|%d | %d|:\n",Targetind[2],Targetind[3]);
//         printf("Target Values for 2*2 blocks:\n");
//         printf("|%d | %d|:\n",targetvalue[0],targetvalue[1]);
//         printf("|%d | %d|:\n",targetvalue[2],targetvalue[3]);
//         }
// /////////////////////////////////////////////////
// /////////////////////////////////////////////////
// /////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////     
////////////////////////////  Optical Matching      ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
         #pragma omp for reduction (+:temp)
         for (int i = 0; i < NBlocks; i++)
        {
            for (int j = 0; j < NBlocks; j++)
            {
                for (int dx = -Radius; dx <= Radius; dx++) 
                {
                    for (int dy = -Radius; dy <= Radius; dy++)
                    {
                        if (((i+dx)*NPixel >= 0) && ((i+dx)*NPixel < N) && ((j+dy)*NPixel) >= 0 && ((j+dy)*NPixel) < N)
                        {
                            temp = 0;
                            // HammingDistance(tmd, tm2d, i*NPixel*N + j*NPixel, ((i+dx)*N*NPixel) + ((j+dy)*NPixel), &temp, N, NPixel);
                            for (int k = 0; k < NPixel; k++)
                                for (int l = 0; l < NPixel; l++)
                                    temp = temp + int (tmd[(i+k)*N+(j+l)] xor tm2d[(i+k+dx*Radius)*N+(j+l+dy*Radius)]);
                                // printf("size is %d\n",sizeof(tmd[(i+k)*N+(j+l)]));
                            if (temp < targetvalue[i*NBlocks + j])
                            {
                                targetvalue[i*NBlocks + j] = temp;
                                Targetind[i*NBlocks + j] = (i+dx)*NPixel*N + (j+dy)*NPixel;
                            }
                        }
                    }
                }
            }
        }

// /////////////////////////////////////////////////
// ///////////////// just for test /////////////////
// /////////////////////////////////////////////////
//         #pragma omp critical
//         {
//         printf("New Target Indices for 2*2 blocks:\n");
//         printf("|%d | %d|:\n",Targetind[0],Targetind[1]);
//         printf("|%d | %d|:\n",Targetind[2],Targetind[3]);
//         printf("New Target Values for 2*2 blocks:\n");
//         printf("|%d | %d|:\n",targetvalue[0],targetvalue[1]);
//         printf("|%d | %d|:\n",targetvalue[2],targetvalue[3]);
//         }
// /////////////////////////////////////////////////
// /////////////////////////////////////////////////
// /////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////     
////////////////////////////////////  Optical Flow Calculation ////////////////////
///////////////////////////////////////////////////////////////////////////////////
         #pragma omp for
        for (int i = 0; i < NBlocks; i++)
        {
            for (int j = 0; j < NBlocks; j++)
            {
                // count = ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel;
                // Calculating the Optical Flow
                int tmdaddr= i*NPixel*N + j*NPixel;
                // int tm2daddri= (Targetind[i*NBlocks + j]*NPixel)/N;
                // int tm2daddrj= Targetind[i*NBlocks + j]-tm2daddri*(N/NPixel);
                // int tm2daddr= tm2daddri*NPixel*N+tm2daddrj*NPixel;
                // int tm2daddr= Targetind[i*NBlocks + j];
                // OpticalFlow(tmd, tm2d, tmdaddr, Targetind[i*NBlocks + j], velocity, N, NPixel);
                for (int k = 0; k < NPixel; k++){
                    for (int l = 0; l < NPixel; l++){
                        velocity[tmdaddr+(k*N+l)] = int (tmd[tmdaddr+(k*N+l)] - tm2d[Targetind[i*NBlocks + j]+(k*N+l)]);
                        // printf("tmd addr:%d, tm2d addr:%d, velocity in %d is %d\n",tmdaddr+(k*N+l), tm2daddr+(k*N+l), tmdaddr+(k*N+l), velocity[tmdaddr+(k*N+l)]);
                    }
                }
            }
        }
        #pragma omp barrier
}

// /////////////////////////////////////////////////
// ///////////////// just for test /////////////////
// /////////////////////////////////////////////////
//         #pragma omp critical
//         {
//         printf("Diff for 4*4 matrix and 2*2 blocks:\n");
//         printf("|%d  %d | %d  %d|:\n",velocity[0],velocity[1],velocity[2],velocity[3]);
//         printf("|%d  %d | %d  %d|:\n",velocity[4],velocity[5],velocity[6],velocity[7]);
//         printf("---------------\n");
//         printf("|%d  %d | %d  %d|:\n",velocity[8],velocity[9],velocity[10],velocity[11]);
//         printf("|%d  %d | %d  %d|:\n",velocity[12],velocity[13],velocity[14],velocity[15]);
//         }
// /////////////////////////////////////////////////
// /////////////////////////////////////////////////
// /////////////////////////////////////////////////

  // auto end_sim = std::chrono::high_resolution_clock::now();
  // double simdur = std::chrono::duration <double> (end_sim - begin_sim).count();
  //  printf( "N: %d, threads: %d, NPixels: %d, NBlocks: %d, Radius: %d, OpenMP ST: %g\n", N, numthreads, NPixel, N/NPixel, Radius, simdur);


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
return 0;
}

         //        int tmdaddr=i*NPixel*N + j*NPixel;
         //        int tm2daddr=Targetind3[i*NBlocks + j];
         // for (int k = 0; k < NPixel; k++){
         //        for (int l = 0; l < NPixel; l++){
         //            velocity[tmdaddr+(k*N+l)] = int (tmd[tmdaddr+(k*N+l)] - tm2d[tm2daddr+(k*N+l)]);
         //            // upcxx::future<> fput5 = upcxx::rput(velocity[(i*NPixel*N + j*NPixel)+(k*N+l)], velocitypv[upcxx::rank_me()] + j);
         //            // fut_all5 = upcxx::when_all(fut_all5, fput5);
         //        }
         //    }
