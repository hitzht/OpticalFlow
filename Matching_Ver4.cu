#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include <chrono>

// #include "Matchcommon.h"
#define Radius 1
// #define NPixel 8
#define Deltat 0.00001
// int NBlocks;
// int Blocks;
// int binNum;

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

__device__ void HammingDistance(int *tmd, int *tm2d, int tmdaddr, int tm2daddr, int *temp, int N, int NPixel)
{
            int xtid = threadIdx.x + blockIdx.x * blockDim.x;
            int ytid = threadIdx.y + blockIdx.y * blockDim.y;
            int xoffset = gridDim.x*blockDim.x;
            int yoffset = gridDim.y*blockDim.y;
         for (int k = xtid; k < NPixel; k+=xoffset)
                for (int l = ytid; l < NPixel; l+=yoffset)
                    *temp += int (tmd[tmdaddr+(k*N+l)] ^ tm2d[tm2daddr+(k*N+l)]);
}

__device__ void OpticalFlow(int *tmd, int *tm2d, int tmdaddr, int tm2daddr, int *velocity, int N, int NPixel)
{
        int xtid = threadIdx.x + blockIdx.x * blockDim.x;
        int ytid = threadIdx.y + blockIdx.y * blockDim.y;
        int xoffset = gridDim.x*blockDim.x;
        int yoffset = gridDim.y*blockDim.y;
         for (int k = xtid; k < NPixel; k+=xoffset)
                for (int l = ytid; l < NPixel; l+=yoffset)                    
                    velocity[tmdaddr+(k*N+l)] = int (tmd[tmdaddr+(k*N+l)] - tm2d[tm2daddr+(k*N+l)]);
}

__global__ void TargetDefining( int* tmd, int* tm2d, int* Targetind, int* targetvalue, int N, int NPixel)
{
    // // NBlocks = int (N/NPixel);
        int xtid = threadIdx.x + blockIdx.x * blockDim.x;
        int ytid = threadIdx.y + blockIdx.y * blockDim.y;
        int xoffset = gridDim.x*blockDim.x;
        int yoffset = gridDim.y*blockDim.y;
        for (int i = xtid; i < (N/NPixel); i+=xoffset){
            for (int j = ytid; j < (N/NPixel); j+=yoffset){
            Targetind[i*(N/NPixel) + j] = i*NPixel*N + j*NPixel; // Assign block in t-2d as target
            HammingDistance(tmd, tm2d, i*NPixel*N + j*NPixel, i*NPixel*N + j*NPixel, &targetvalue[i*(N/NPixel) + j], N,NPixel );
            }
        }
}

__constant__ const int dir[8][2]={{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1},{-1,0}};
__global__ void OpticalMatching(int* __restrict__ tmd, int* __restrict__ tm2d, int* __restrict__ Targetind, int* __restrict__ targetvalue, int N, int NPixel)
{
        int xtid = threadIdx.x + blockIdx.x * blockDim.x;
        int ytid = threadIdx.y + blockIdx.y * blockDim.y;
        int xoffset = gridDim.x*blockDim.x;
        int yoffset = gridDim.y*blockDim.y;
        for (int i = xtid; i <  N/NPixel; i+=xoffset){
            for (int j = ytid; j <  N/NPixel; j+=yoffset){
                for(int t=0;t<8;t++){
                              int x = (i + dir[t][0]);
                              int y = (j + dir[t][1]);
                        if (x*NPixel >= 0 && x*NPixel < N && y*NPixel >= 0 && y*NPixel < N)
                        {
                            int temp = 0;
                            // Calculating the hamming distance
                            HammingDistance(tmd, tm2d, i*NPixel*N + j*NPixel, (x*NPixel*N) + (y*NPixel), &temp, N, NPixel);
                            // printf("temp for block %d where dx is %d and dy is %d is %d\n",i* N/NPixel+j,dx,dy,temp);
                            // for (int k = 0; k < NPixel; k++)
                            //     for (int l = 0; l < NPixel; l++)
                            //         temp = temp + int (tmd[(i+k)*N+(j+l)] xor tm2d[(i+k+dx*Radius)*N+(j+l+dy*Radius)]);
                            // printf("temp is %d\n",temp);
                            if (temp < targetvalue[i* N/NPixel + j])
                            {
                                targetvalue[i* N/NPixel + j] = temp;
                                // printf("targetvalue[%d] is %d\n",i* N/NPixel+j,targetvalue[i* N/NPixel + j]);
                                Targetind[i* N/NPixel + j] = x*NPixel*N + y*NPixel;
                            }
                        }
                }
            }
        }
}

__global__ void OpticalFlowCalculation(int *tmd, int *tm2d, int *Targetind, int *velocity, int N, int NPixel){
        int xtid = threadIdx.x + blockIdx.x * blockDim.x;
        int ytid = threadIdx.y + blockIdx.y * blockDim.y;
        int xoffset = gridDim.x*blockDim.x;
        int yoffset = gridDim.y*blockDim.y;
        for (int i = xtid; i < (N/NPixel); i+=xoffset){
            for (int j = ytid; j < (N/NPixel); j+=yoffset){
                OpticalFlow(tmd, tm2d, i*NPixel*N + j*NPixel, Targetind[i*(N/NPixel)+j], velocity, N, NPixel);
            }
        }
}

int main(int argc, char **argv)
{
int N = read_int( argc, argv, "-n", 64 ); // the number of particles
int NUM_THREADS = read_int( argc, argv, "-t", 256 ); // the number of particles
int blks = read_int( argc, argv, "-b", 1024 ); // the number of particles
int NPixel = read_int( argc, argv, "-p", 8 ); // the number of particles


int NBlocks = int (N/NPixel);
// binNum = int(N / NPixel); // Should be around sqrt(N/2)

int *tmd = (int*) malloc( N * N * sizeof(int) );
int *tm2d = (int*) malloc( N * N * sizeof(int) );
int *velocity = (int*) malloc( N * N * sizeof(int) );

int *Targetind = (int*) malloc( NBlocks * NBlocks * sizeof(int) );
int *targetvalue = (int*) malloc( NBlocks * NBlocks * sizeof(int) );

int * d_tmd,*d_tm2d, *d_velocity, *d_Targetind, *d_targetvalue;
cudaMalloc((void **) &d_tmd, N * N * sizeof(int));
cudaMalloc((void **) &d_tm2d, N * N * sizeof(int));
cudaMalloc((void **) &d_velocity, N * N * sizeof(int));
cudaMalloc((void **) &d_Targetind, NBlocks * NBlocks * sizeof(int));
cudaMalloc((void **) &d_targetvalue, NBlocks * NBlocks * sizeof(int));

init_events(N, tmd);
init_events(N, tm2d);

cudaDeviceSynchronize();
double copy_time = read_timer( );
cudaMemcpy(d_tmd, tmd, N * N * sizeof(int), cudaMemcpyHostToDevice);
cudaMemcpy(d_tm2d, tm2d, N * N * sizeof(int), cudaMemcpyHostToDevice);
cudaDeviceSynchronize();
copy_time = read_timer( ) - copy_time;

auto begin_sim = std::chrono::high_resolution_clock::now();
////////////////////////////////////////////////////////////////////////////////     
////////////////////////////  Target Defining      ////////////////////////////
////////////////////////////////////////////////////////////////////////////////
        int threadNum = NUM_THREADS;
        // int blks = min(1024,(N*N + NUM_THREADS - 1) / NUM_THREADS);
        int blockNum = blks;//min(512,(n+threadNum-1)/threadNum);

        cudaMemset(d_Targetind, 0, NBlocks * NBlocks * sizeof(int));
        cudaMemset(d_targetvalue, 0, NBlocks * NBlocks * sizeof(int));
        TargetDefining<<<blockNum,threadNum>>>(d_tmd, d_tm2d, d_Targetind, d_targetvalue, N, NPixel );
        // cudaMemcpy(Targetind, d_Targetind, NBlocks * NBlocks * sizeof(int), cudaMemcpyDeviceToHost);
        // cudaMemcpy(targetvalue, d_targetvalue, NBlocks * NBlocks * sizeof(int), cudaMemcpyDeviceToHost);

        // for (int i = 0; i < N; i+=NPixel)
        // {
        //     for (int j = 0; j < N; j+=NPixel)
        //     {
        //     count = ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel;
        //     Targetind[count] = i*N + j; // Assign block in t-2d as target
        //     // printf("count= %d, ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel= %d\n", count, ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel);
        //     HammingDistance(NPixel, N, i*N + j, i*N + j, tmd, tm2d, &targetvalue[count]);
        //     }
        // }

////////////////////////////////////////////////////////////////////////////////     
////////////////////////////  Optical Matching      ////////////////////////////
////////////////////////////////////////////////////////////////////////////////

        // cudaMemcpy(d_Targetind, Targetind, NBlocks * NBlocks * sizeof(int), cudaMemcpyHostToDevice);
        // cudaMemcpy(d_targetvalue, targetvalue, NBlocks * NBlocks * sizeof(int), cudaMemcpyHostToDevice);
        OpticalMatching<<<blockNum,threadNum>>>(d_tmd, d_tm2d, d_Targetind, d_targetvalue, N, NPixel);
        // cudaMemcpy(Targetind, d_Targetind, NBlocks * NBlocks * sizeof(int), cudaMemcpyDeviceToHost);

        //  for (int i = 0; i < N; i+=NPixel)
        // {
        //     for (int j = 0; j < N; j+=NPixel)
        //     {
        //         count = ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel;
        //         // int *temp = (int*) calloc((2*Radius +1)*(2*Radius +1), sizeof(int));
        //         // Targetind[i*NBlocks+j] = (i * NBlocks + j);
        //         for (int dx = -Radius; dx <= Radius; dx++)   //Search over nearby 8 blocks and the target block 
        //         {
        //             for (int dy = -Radius; dy <= Radius; dy++)
        //             {
        //                 if (i + dx*NPixel >= 0 && i + dx*NPixel < N && j + dy*NPixel >= 0 && j + dy*NPixel < N)
        //                 {
        //                     temp = 0;
        //                     // Calculating the hamming distance
        //                     HammingDistance(NPixel, N, i*N + j, (i+dx*Radius)*N + (j+dy*Radius), tmd, tm2d, &temp);
        //                     // for (int k = 0; k < NPixel; k++)
        //                     //     for (int l = 0; l < NPixel; l++)
        //                     //         temp = temp + int (tmd[(i+k)*N+(j+l)] xor tm2d[(i+k+dx*Radius)*N+(j+l+dy*Radius)]);
        //                     // printf("temp is %d\n",temp);
        //                     if (temp < targetvalue[count])
        //                     {
        //                     	targetvalue[count] = temp;
        //                         Targetind[count] = (i + dx*NPixel)*N + (j + dy*NPixel);
        //                         // int Targetind[count] = 0;
        //                     }
        //                 }
        //             }
        //         }
        //     }
        // }
///////////////////////////////////////////////////////////////////////////////////     
////////////////////////////////////  Optical Flow Calculation ////////////////////
///////////////////////////////////////////////////////////////////////////////////
        // cudaMemcpy(d_Targetind, Targetind, NBlocks * NBlocks * sizeof(int), cudaMemcpyHostToDevice);
        cudaMemset(d_velocity, 0, N * N * sizeof(int) );
        OpticalFlowCalculation<<<blockNum,threadNum>>>(d_tmd, d_tm2d, d_Targetind, d_velocity, N, NPixel);
        cudaMemcpy(velocity, d_velocity, N * N * sizeof(int), cudaMemcpyDeviceToHost);

//         for (int i = 0; i < N; i+=NPixel)
//         {
//             for (int j = 0; j < N; j+=NPixel)
//             {
//                 count = ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel;
// // Calculating the Optical Flow
//                 OpticalFlow(NPixel, N, i*N + j, Targetind[count], tmd, tm2d, velocity);                
//                 // for (int k = 0; k < NPixel; k++)
//                 //     for (int l = 0; l < NPixel; l++)
//                 //         OpticalFlow[(i+k)*N+(j+l)] = tmd[(i+k)*N+(j+l)] - tm2d[Targetind[count]+(k*N+l)];
//             }
//         }

   cudaDeviceSynchronize();
  auto end_sim = std::chrono::high_resolution_clock::now();
  double simdur = std::chrono::duration <double> (end_sim - begin_sim).count();
  printf( "N: %d, NPixels: %d, NBlocks: %d, Radius: %d, Threads: %d, Blocks: %d, GPGPU ST: %g, CT: %g\n", N, NPixel, NBlocks, Radius, NUM_THREADS, blockNum, simdur, copy_time);


// int i, j, count = 0; 
// int *B[N];
//     for (i=0; i<M; i++)
//          B[i] = (int *)malloc(N * sizeof(int));

    // for (int i = 0; i <  N; i++) 
    //   for (int j = 0; j < N; j++)
    //      printf("tmd in position %d*%d is:%d\n",i, j, tmd[i*N+j]);
cudaFree(d_tmd);
cudaFree(d_tm2d);
cudaFree(d_velocity);

free( tmd );
free( tm2d );
free( velocity );
free( Targetind );
free( targetvalue );
return 0;
}
// //module load cuda
// //salloc -N 1 -t 01:30:00 -p gpu

// #include <stdlib.h>
// #include <stdio.h>
// #include <assert.h>
// #include <float.h>
// #include <string.h>
// #include <math.h>
// #include <time.h>
// #include <sys/time.h>
// #include <vector>
// // #include "Matchcommon.h"
// #define NUM_THREADS 256
// #define Radius 1
// #define NPixel 8
// #define Deltat 0.00001
// int NBlocks;
// int Blocks;
// int binNum;

// int find_option( int argc, char **argv, const char *option )
// {
//     for( int i = 1; i < argc; i++ )
//         if( strcmp( argv[i], option ) == 0 )
//             return i;
//     return -1;
// }

// int read_int( int argc, char **argv, const char *option, int default_value ) // 
// {
//     int iplace = find_option( argc, argv, option ); 
//     if( iplace >= 0 && iplace < argc-1 ) // 
//         return atoi( argv[iplace+1] ); // Ascii to integer conversion
//     return default_value;
// }


// double read_timer( )
// {
//     static bool initialized = false;
//     static struct timeval start;
//     struct timeval end;
//     if( !initialized )
//     {
//         gettimeofday( &start, NULL );
//         initialized = true;
//     }
//     gettimeofday( &end, NULL );
//     return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
// }

// void init_events( int n, int *g1 , int *g2)
// {
//     srand48( time( NULL ) );

//     for (int i = 0; i <  n; i++) {
//         for (int j = 0; j < n; j++){ 
//             g1[i*n+j] = drand48()+.5;
//             g2[i*n+j] = drand48()+.5;
//         }
//     }
// }

// __device__ void HammingDistance(int *tmd, int *tm2d, int tmdaddr, int tm2daddr, int *temp, int N )
// {
//          //    int xtid = threadIdx.x + blockIdx.x * blockDim.x;
//          //    int ytid = threadIdx.y + blockIdx.y * blockDim.y;
//          //    int xoffset = gridDim.x*blockDim.x;
//          //    int yoffset = gridDim.y*blockDim.y;
//          // for (int k = xtid; k < NPixel; k+=xoffset){
//          //        for (int l = ytid; l < NPixel; l+=yoffset){
//          for (int k = 0; k < NPixel; k++){
//                 for (int l = 0; l < NPixel; l++){
//                     *temp += int (tmd[tmdaddr+(k*N+l)] ^ tm2d[tm2daddr+(k*N+l)]);
//                     // atomicAdd(temp,int (tmd[tmdaddr+(k*N+l)] ^ tm2d[tm2daddr+(k*N+l)]));
//         // printf("tm2daddr is %d+ (k*N+l) is %d and tm2daddr+(k*N+l) is %d and temp is %d\n", tm2daddr, (k*N+l),tm2daddr+(k*N+l),*temp);
//                 }
//          }
// }

// __device__ void OpticalFlow(int *tmd, int *tm2d, int tmdaddr, int tm2daddr, int *velocity, int N )
// {
//         // int xtid = threadIdx.x + blockIdx.x * blockDim.x;
//         // int ytid = threadIdx.y + blockIdx.y * blockDim.y;
//         // int xoffset = gridDim.x*blockDim.x;
//         // int yoffset = gridDim.y*blockDim.y;
//         //  for (int k = xtid; k < NPixel; k+=xoffset)
//         //         for (int l = ytid; l < NPixel; l+=yoffset)                    
//          for (int k = 0; k < NPixel; k++){
//                 for (int l = 0; l < NPixel; l++){
//                     velocity[tmdaddr+(k*N+l)] = int (tmd[tmdaddr+(k*N+l)] - tm2d[tm2daddr+(k*N+l)]);
//                 }
//             }
// }

// __global__ void TargetDefining( int* tmd, int* tm2d, int* Targetind, int* targetvalue, int N )
// {
//     // // NBlocks = int (N/NPixel);
//         int xtid = threadIdx.x + blockIdx.x * blockDim.x;
//         int ytid = threadIdx.y + blockIdx.y * blockDim.y;
//         int xoffset = gridDim.x*blockDim.x;
//         int yoffset = gridDim.y*blockDim.y;
//         for (int i = xtid; i < (N/NPixel); i+=xoffset){
//             for (int j = ytid; j < (N/NPixel); j+=yoffset){
//             Targetind[i*(N/NPixel) + j] = i*(N/NPixel) + j; // Assign block in t-2d as target
//             HammingDistance(tmd, tm2d, i*NPixel*N + j*NPixel, i*NPixel*N + j*NPixel, &targetvalue[i*(N/NPixel) + j], N );
//             // printf("targetvalue[i*NBlocks + j] is %d\n",targetvalue[i*(N/NPixel) + j]);
//             }
//         }
// }

// __constant__ const int dir[8][2]={{-1,-1},{0,-1},{1,-1},{1,0},{1,1},{0,1},{-1,1},{-1,0}};

// __global__ void OpticalMatching(int* __restrict__ tmd, int* __restrict__ tm2d, int* __restrict__ Targetind, int* __restrict__ targetvalue, int N){
//         int xtid = threadIdx.x + blockIdx.x * blockDim.x;
//         int ytid = threadIdx.y + blockIdx.y * blockDim.y;
//         int xoffset = gridDim.x*blockDim.x;
//         int yoffset = gridDim.y*blockDim.y;
//         for (int i = xtid; i <  N/NPixel; i+=xoffset){
//             for (int j = ytid; j <  N/NPixel; j+=yoffset){
//                 for(int t=0;t<8;t++){
//                               int x = (i + dir[t][0]);
//                               int y = (j + dir[t][1]);
//                 // for (int dx = -Radius; dx <= Radius; dx++)   //Search over nearby 8 blocks and the target block 
//                 // {
//                 //     for (int dy = -Radius; dy <= Radius; dy++)
//                 //     {
//                         if (x*NPixel >= 0 && x*NPixel < N && y*NPixel >= 0 && y*NPixel < N)
//                         {
//                             int temp = 0;
//                             // Calculating the hamming distance
//                             HammingDistance(tmd, tm2d, i*NPixel*N + j*NPixel, ((i+dx)*NPixel*N) + ((j+dy)*NPixel), &temp, N);
//                             // printf("temp for block %d where dx is %d and dy is %d is %d\n",i* N/NPixel+j,dx,dy,temp);
//                             // for (int k = 0; k < NPixel; k++)
//                             //     for (int l = 0; l < NPixel; l++)
//                             //         temp = temp + int (tmd[(i+k)*N+(j+l)] xor tm2d[(i+k+dx*Radius)*N+(j+l+dy*Radius)]);
//                             // printf("temp is %d\n",temp);
//                             if (temp < targetvalue[i* N/NPixel + j])
//                             {
//                                 targetvalue[i* N/NPixel + j] = temp;
//                                 // printf("targetvalue[%d] is %d\n",i* N/NPixel+j,targetvalue[i* N/NPixel + j]);
//                                 Targetind[i* N/NPixel + j] = x* N/NPixel + y;
//                             }
//                         }
//                 }
//             }
//         }
// }

// __global__ void OpticalFlowCalculation(int *tmd, int *tm2d, int *Targetind, int *velocity, int N){
//         int xtid = threadIdx.x + blockIdx.x * blockDim.x;
//         int ytid = threadIdx.y + blockIdx.y * blockDim.y;
//         int xoffset = gridDim.x*blockDim.x;
//         int yoffset = gridDim.y*blockDim.y;
//         for (int i = xtid; i < (N/NPixel); i+=xoffset){
//             for (int j = ytid; j < (N/NPixel); j+=yoffset){
//                 int tm2daddri= (Targetind[i*(N/NPixel)+j]*NPixel)/N;
//                 int tm2daddrj= Targetind[i*(N/NPixel) + j]-tm2daddri*(N/NPixel);
//                 int tm2daddr= tm2daddri*NPixel*N+tm2daddrj*NPixel;
//                 OpticalFlow(tmd, tm2d, i*NPixel*N + j*NPixel, tm2daddr, velocity, N);
//             }
//         }
// }



// int main(int argc, char **argv)
// {
// int N = read_int( argc, argv, "-n", 32 ); // the number of particles
// NBlocks = int (N/NPixel);
// // binNum = int(N / NPixel); // Should be around sqrt(N/2)

// int *tmd = (int*) malloc( N * N * sizeof(int) );
// int *tm2d = (int*) malloc( N * N * sizeof(int) );
// int *velocity = (int*) malloc( N * N * sizeof(int) );

// int *Targetind = (int*) malloc( NBlocks * NBlocks * sizeof(int) );
// int *targetvalue = (int*) malloc( NBlocks * NBlocks * sizeof(int) );

// int * d_tmd,*d_tm2d, *d_velocity, *d_Targetind, *d_targetvalue;
// cudaMalloc((void **) &d_tmd, N * N * sizeof(int));
// cudaMalloc((void **) &d_tm2d, N * N * sizeof(int));
// cudaMalloc((void **) &d_velocity, N * N * sizeof(int));
// cudaMalloc((void **) &d_Targetind, NBlocks * NBlocks * sizeof(int));
// cudaMalloc((void **) &d_targetvalue, NBlocks * NBlocks * sizeof(int));

// init_events(N, tmd,tm2d);
// // /////////////////////////////////////////////////
// // ///////////////// just for test /////////////////
// // /////////////////////////////////////////////////
// //         printf("tmd is a 4*4 matrix in 2*2 block:\n");
// //         printf("|%d  %d | %d  %d|:\n",tmd[0],tmd[1],tmd[2],tmd[3]);
// //         printf("|%d  %d | %d  %d|:\n",tmd[4],tmd[5],tmd[6],tmd[7]);
// //         printf("---------------\n");
// //         printf("|%d  %d | %d  %d|:\n",tmd[8],tmd[9],tmd[10],tmd[11]);
// //         printf("|%d  %d | %d  %d|:\n",tmd[12],tmd[13],tmd[14],tmd[15]);

// //         printf("tm2d is a 4*4 matrix in 2*2 block:\n");
// //         printf("|%d  %d | %d  %d|:\n",tm2d[0],tm2d[1],tm2d[2],tm2d[3]);
// //         printf("|%d  %d | %d  %d|:\n",tm2d[4],tm2d[5],tm2d[6],tm2d[7]);
// //         printf("---------------\n");
// //         printf("|%d  %d | %d  %d|:\n",tm2d[8],tm2d[9],tm2d[10],tm2d[11]);
// //         printf("|%d  %d | %d  %d|:\n",tm2d[12],tm2d[13],tm2d[14],tm2d[15]);
// // /////////////////////////////////////////////////
// // /////////////////////////////////////////////////
// // /////////////////////////////////////////////////

// cudaDeviceSynchronize();
// double copy_time = read_timer();
// cudaMemcpy(d_tmd, tmd, N * N * sizeof(int), cudaMemcpyHostToDevice);
// cudaMemcpy(d_tm2d, tm2d, N * N * sizeof(int), cudaMemcpyHostToDevice);
// cudaDeviceSynchronize();
// copy_time = read_timer() - copy_time;

// double simulation_time = read_timer( );
// ////////////////////////////////////////////////////////////////////////////////     
// ////////////////////////////  Target Defining      ////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////
//         int threadNum = NUM_THREADS;
//         int blks = min(1024,(N*N + NUM_THREADS - 1) / NUM_THREADS);
//         int blockNum = blks;//min(512,(n+threadNum-1)/threadNum);

//         cudaMemset(d_Targetind, 0, NBlocks * NBlocks * sizeof(int));
//         cudaMemset(d_targetvalue, 0, NBlocks * NBlocks * sizeof(int));
//         TargetDefining<<<blockNum,threadNum>>>(d_tmd, d_tm2d, d_Targetind, d_targetvalue, N );
//         cudaMemcpy(Targetind, d_Targetind, NBlocks * NBlocks * sizeof(int), cudaMemcpyDeviceToHost);
//         cudaMemcpy(targetvalue, d_targetvalue, NBlocks * NBlocks * sizeof(int), cudaMemcpyDeviceToHost);

// // /////////////////////////////////////////////////
// // ///////////////// just for test /////////////////
// // /////////////////////////////////////////////////
// //         printf("Target Indices for 2*2 blocks:\n");
// //         printf("|%d | %d|:\n",Targetind[0],Targetind[1]);
// //         printf("|%d | %d|:\n",Targetind[2],Targetind[3]);
// //         printf("Target Values for 2*2 blocks:\n");
// //         printf("|%d | %d|:\n",targetvalue[0],targetvalue[1]);
// //         printf("|%d | %d|:\n",targetvalue[2],targetvalue[3]);
// // /////////////////////////////////////////////////
// // /////////////////////////////////////////////////
// // /////////////////////////////////////////////////

//         // for (int i = 0; i < N; i+=NPixel)
//         // {
//         //     for (int j = 0; j < N; j+=NPixel)
//         //     {
//         //     count = ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel;
//         //     Targetind[count] = i*N + j; // Assign block in t-2d as target
//         //     // printf("count= %d, ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel= %d\n", count, ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel);
//         //     HammingDistance(NPixel, N, i*N + j, i*N + j, tmd, tm2d, &targetvalue[count]);
//         //     }
//         // }

// ////////////////////////////////////////////////////////////////////////////////     
// ////////////////////////////  Optical Matching      ////////////////////////////
// ////////////////////////////////////////////////////////////////////////////////
//         cudaMemcpy(d_Targetind, Targetind, NBlocks * NBlocks * sizeof(int), cudaMemcpyHostToDevice);
//         cudaMemcpy(d_targetvalue, targetvalue, NBlocks * NBlocks * sizeof(int), cudaMemcpyHostToDevice);
//         OpticalMatching<<<blockNum,threadNum>>>(d_tmd, d_tm2d, d_Targetind, d_targetvalue, N);
//         cudaMemcpy(Targetind, d_Targetind, NBlocks * NBlocks * sizeof(int), cudaMemcpyDeviceToHost);
//         // cudaMemcpy(targetvalue, d_targetvalue, NBlocks * NBlocks * sizeof(int), cudaMemcpyDeviceToHost);
// // /////////////////////////////////////////////////
// // ///////////////// just for test /////////////////
// // /////////////////////////////////////////////////
// //         printf("New Target Indices for 2*2 blocks:\n");
// //         printf("|%d | %d|:\n",Targetind[0],Targetind[1]);
// //         printf("|%d | %d|:\n",Targetind[2],Targetind[3]);
// //         printf("New Target Values for 2*2 blocks:\n");
// //         printf("|%d | %d|:\n",targetvalue[0],targetvalue[1]);
// //         printf("|%d | %d|:\n",targetvalue[2],targetvalue[3]);
// // /////////////////////////////////////////////////
// // /////////////////////////////////////////////////
// // /////////////////////////////////////////////////

//         //  for (int i = 0; i < N; i+=NPixel)
//         // {
//         //     for (int j = 0; j < N; j+=NPixel)
//         //     {
//         //         count = ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel;
//         //         // int *temp = (int*) calloc((2*Radius +1)*(2*Radius +1), sizeof(int));
//         //         // Targetind[i*NBlocks+j] = (i * NBlocks + j);
//         //         for (int dx = -Radius; dx <= Radius; dx++)   //Search over nearby 8 blocks and the target block 
//         //         {
//         //             for (int dy = -Radius; dy <= Radius; dy++)
//         //             {
//         //                 if (i + dx*NPixel >= 0 && i + dx*NPixel < N && j + dy*NPixel >= 0 && j + dy*NPixel < N)
//         //                 {
//         //                     temp = 0;
//         //                     // Calculating the hamming distance
//         //                     HammingDistance(NPixel, N, i*N + j, (i+dx*Radius)*N + (j+dy*Radius), tmd, tm2d, &temp);
//         //                     // for (int k = 0; k < NPixel; k++)
//         //                     //     for (int l = 0; l < NPixel; l++)
//         //                     //         temp = temp + int (tmd[(i+k)*N+(j+l)] xor tm2d[(i+k+dx*Radius)*N+(j+l+dy*Radius)]);
//         //                     // printf("temp is %d\n",temp);
//         //                     if (temp < targetvalue[count])
//         //                     {
//         //                      targetvalue[count] = temp;
//         //                         Targetind[count] = (i + dx*NPixel)*N + (j + dy*NPixel);
//         //                         // int Targetind[count] = 0;
//         //                     }
//         //                 }
//         //             }
//         //         }
//         //     }
//         // }
// ///////////////////////////////////////////////////////////////////////////////////     
// ////////////////////////////////////  Optical Flow Calculation ////////////////////
// ///////////////////////////////////////////////////////////////////////////////////
//         cudaMemcpy(d_Targetind, Targetind, NBlocks * NBlocks * sizeof(int), cudaMemcpyHostToDevice);
//         cudaMemset(d_velocity, 0, N * N * sizeof(int) );
//         OpticalFlowCalculation<<<blockNum,threadNum>>>(d_tmd, d_tm2d, d_Targetind, d_velocity, N);
//         // cudaMemcpy(velocity, d_velocity, N * N * sizeof(int), cudaMemcpyDeviceToHost);
// // /////////////////////////////////////////////////
// // ///////////////// just for test /////////////////
// // /////////////////////////////////////////////////
// //         printf("Diff for 4*4 matrix and 2*2 blocks:\n");
// //         printf("|%d  %d | %d  %d|:\n",velocity[0],velocity[1],velocity[2],velocity[3]);
// //         printf("|%d  %d | %d  %d|:\n",velocity[4],velocity[5],velocity[6],velocity[7]);
// //         printf("---------------\n");
// //         printf("|%d  %d | %d  %d|:\n",velocity[8],velocity[9],velocity[10],velocity[11]);
// //         printf("|%d  %d | %d  %d|:\n",velocity[12],velocity[13],velocity[14],velocity[15]);
// // /////////////////////////////////////////////////
// // /////////////////////////////////////////////////
// // /////////////////////////////////////////////////
// //         for (int i = 0; i < N; i+=NPixel)
// //         {
// //             for (int j = 0; j < N; j+=NPixel)
// //             {
// //                 count = ((i+NPixel-1)/NPixel)*NBlocks + (j+NPixel-1)/NPixel;
// // // Calculating the Optical Flow
// //                 OpticalFlow(NPixel, N, i*N + j, Targetind[count], tmd, tm2d, velocity);                
// //                 // for (int k = 0; k < NPixel; k++)
// //                 //     for (int l = 0; l < NPixel; l++)
// //                 //         OpticalFlow[(i+k)*N+(j+l)] = tmd[(i+k)*N+(j+l)] - tm2d[Targetind[count]+(k*N+l)];
// //             }
// //         }
//    cudaDeviceSynchronize();
//    simulation_time = read_timer( ) - simulation_time;
//    printf( "N: %d, NPixels: %d, NBlocks: %d, Radius: %d, GPGPU ST: %g s, CT: %g s\n", N, NPixel, N/NPixel, Radius, simulation_time,copy_time);


// // int i, j, count = 0; 
// // int *B[N];
// //     for (i=0; i<M; i++)
// //          B[i] = (int *)malloc(N * sizeof(int));

//     // for (int i = 0; i <  N; i++) 
//     //   for (int j = 0; j < N; j++)
//     //      printf("tmd in position %d*%d is:%d\n",i, j, tmd[i*N+j]);
// cudaFree(d_tmd);
// cudaFree(d_tm2d);
// cudaFree(d_velocity);

// free( tmd );
// free( tm2d );
// free( velocity );
// free( Targetind );
// free( targetvalue );
// return 0;
// }
// //module load cuda
// //salloc -N 1 -t 01:30:00 -p gpu

