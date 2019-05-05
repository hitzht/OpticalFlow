#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include <time.h>
#include <sys/time.h>

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

void buildBins(vector<vector<particle_t> > &bins, particle_t* particles, int n, double gridSize)
{
    // gridSize = sqrt(n*density);
    int gs = int(gridSize/cutoff) +1;
    binSize = cutoff * 2;  
    binNum = int(gridSize / binSize)+1; // Should be around sqrt(N/2)

    // printf("Grid Size: %.4lf\n",gridSize);
    // printf("Bin Size: %.2lf\n",binSize);
    // printf("number of grids: %d\n",gs);
    // printf("Number of Bins: %d*%d\n",binNum,binNum);
    // Increase\Decrease binNum to be something like 2^k?
    
    bins.resize(binNum * binNum);

    for (int i = 0; i < n; i++)
    {
        int x = int(particles[i].x / binSize);
        int y = int(particles[i].y / binSize);
        // printf("x is %d while particles[i].x is %f\n",x,particles[i].x);
        // printf("y is %d while particles[i].y is %f\n",y,particles[i].y);
        //this adds a new elemen`t to a vector
        bins[x*binNum + y].push_back(particles[i]);
        // printf("size of bin: %lx\n",bins.size());
    }
}

//
//  benchmarking program
//
int main( int argc, char **argv )
{    
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 

    //
    //  process command line parameters
    //
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
    
    int n = read_int( argc, argv, "-n", 5 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );




//     int particle_per_proc = (n + n_proc - 1) / n_proc; // this is written in this way to show at least
//     // one particle belongs to each processor
//     // printf("particle_per_proc is:%d and number of processors is %d\n",particle_per_proc,n_proc);
//     int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
// // when the number of particles is more than processors similar intervals in partition_sizes
//     for( int i = 0; i < n_proc+1; i++ ){
//         partition_offsets[i] = min( i * particle_per_proc, n );
//         // printf("partition_offsets[i] where i is %d is:%d\n", i, partition_offsets[i]);
// }

//     int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
//     for( int i = 0; i < n_proc; i++ ){
//         partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
//         printf("partition_sizes[i] where i is %d is:%d\n", i, partition_sizes[i]);
//     }

    //
    // //  allocate storage for local partition
    // //
    // vector<vector<particle_t> > local_bins;
    // int nlocal = partition_sizes[rank];
    // particle_t *local = (particle_t*) malloc( nlocal * sizeof(particle_t) );    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    double gridSize = set_size( n );
    if( rank == 0 ){
        init_particles( n, particles );
    }
    // MPI_Bcast(
    // void* data,
    // int count,
    // MPI_Datatype datatype,
    // int root,
    // MPI_Comm communicator)
    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);

    //
    //  set up the data partitioning across processors
    //
    vector<vector<particle_t> > particle_bins;
    buildBins(particle_bins, particles, n, gridSize);
    // delete[] particles;
    // particles = NULL;


    // // int MPI_Scatterv(const void *sendbuf, const int *sendcounts, const int *displs,
    // //              MPI_Datatype sendtype, void *recvbuf, int recvcount,
    // //              MPI_Datatype recvtype,
    // //              int root, MPI_Comm comm)
    // // this partition_sizes and partition_offset helps sending massages to different processors at
    // // the same size
    // MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );

    int bins_per_proc = binNum / n_proc;

// This section defines access to particles in bins
    int bins_start = bins_per_proc * rank;
    int bins_end = bins_per_proc * (rank + 1);

    if (rank == n_proc - 1)
        bins_end = binNum;

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;
        // // 
        // //  collect all global data locally (not good idea to do)
        // //
        // //         int MPIAPI MPI_Allgatherv(
        // //   _In_  void         *sendbuf,
        // //         int          sendcount,
        // //         MPI_Datatype sendtype,
        // //   _Out_ void         *recvbuf,
        // //   _In_  int          *recvcounts,
        // //   _In_  int          *displs,
        // //         MPI_Datatype recvtype,
        // //         MPI_Comm     comm
        // // );
        // MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );

        // //
        // //  save current step if necessary (slightly different semantics than in other codes)
        // //
        // if( find_option( argc, argv, "-no" ) == -1 )
        //   if( fsave && (step%SAVEFREQ) == 0 )
        //     save( fsave, n, particles );

        //
        //  compute all forces
        //

        for (int i = bins_start; i < bins_end; ++i) {
            for (int j = 0; j < binNum; ++j) {
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
                            vector<particle_t>& vec2 = particle_bins[(i+dx) * binNum + j + dy];
                            // Comparing two Adjacent Vectors
                            for (int k = 0; k < vec.size(); k++)
                                for (int l = 0; l < vec2.size(); l++)
                                    apply_force(vec[k], vec2[l], &dmin, &davg, &navg);
                        }
                    }
                }
            }
        }

        // for( int i = 0; i < nlocal; i++ )
        // {
        //     local[i].ax = local[i].ay = 0;
        //     for (int j = 0; j < n; j++ )
        //         apply_force( local[i], particles[j], &dmin, &davg, &navg );
        // }
       vector<particle_t> intra_move;
       vector<particle_t> inter_move;

        for (int i = bins_start; i < bins_end; ++i) 
        {
            for (int j = 0; j < binNum; ++j) 
            {
                vector<particle_t>& vec = particle_bins[i * binNum + j];
                int tail = vec.size(), k = 0;
                for (; k < tail; ) {
                    move(vec[k]);
                    int x = int(vec[k].x / binSize);
                    int y = int(vec[k].y / binSize);
                    if (bins_start <= x && x < bins_end) {
                        if (x == i && y == j)
                            ++k;
                        else {
                            intra_move.push_back(vec[k]);
                            vec[k] = vec[--tail];
                        }
                    } else {
                        //int who = x / x_bins_per_proc;
                        inter_move.push_back(vec[k]);
                        vec[k] = vec[--tail];
                    }
                }
                vec.resize(k);
            }
        }
                // printf("inter_move is %f\n", inter_move[1].x);


////////////////////////////////////////////////////////////////
////////////////////////// Intra moving ////////////////////////
////////////////////////////////////////////////////////////////
        for (int i = 0; i < intra_move.size(); ++i) {
            int x = int(intra_move[i].x / binSize);
            int y = int(intra_move[i].y / binSize);
            particle_bins[x*binNum+y].push_back(intra_move[i]);
        }

////////////////////////////////////////////////////////////////
////////////////////////// Inter moving ////////////////////////
////////////////////////////////////////////////////////////////        
        if (rank != 0) {
            for (int i = bins_start - 1, j = 0; j < binNum; ++j) {
                vector<particle_t>& bin = particle_bins[i * binNum + j];
                // printf("inter_move is %f\n", inter_move[1].x);

                bin.clear();
            }
            for (int i = bins_start, j = 0; j < binNum; ++j) {
                vector<particle_t>& bin = particle_bins[i * binNum + j];
                inter_move.insert(inter_move.end(), bin.begin(), bin.end());
                bin.clear();
            }
        }

        if (rank != n_proc - 1) {
            for (int i = bins_end, j = 0; j < binNum; ++j) {
                vector<particle_t>& bin = particle_bins[i * binNum + j];
                bin.clear();
            }
            for (int i = bins_end - 1, j = 0; j < binNum; ++j) {
                vector<particle_t>& bin = particle_bins[i * binNum + j];
                inter_move.insert(inter_move.end(), bin.begin(), bin.end());
                bin.clear();
            }
        }
        vector<particle_t> incoming_move;
        int send_count = inter_move.size();
        int recv_counts[n_proc];

        // printf("worker: %d. MPI_Gather.\n", rank);
        MPI_Gather(&send_count, 1, MPI_INT, recv_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // now root knows recv_counts

        int displs[n_proc];
        int total_num = 0;

        if (rank == 0) {
            displs[0] = 0;
            for (int i = 1; i < n_proc; ++i) {
                displs[i] = displs[i-1] + recv_counts[i-1];
            }
            total_num = recv_counts[n_proc-1] + displs[n_proc-1];
            // printf("worker: %d, 1. %d / %d.\n", rank, total_, total_num);
            // assert(total_ == total_num);
            incoming_move.resize(total_num);
        }

        // now root knows total_num.
        //printf("worker: %d. MPI_Gatherv.\n", rank);

        MPI_Gatherv(inter_move.data(), send_count, PARTICLE, 
            incoming_move.data(), recv_counts, displs, PARTICLE, 
            0, MPI_COMM_WORLD);

        //printf("worker: %d. Classify.\n", rank);

        vector<vector<particle_t> > scatter_particles;
        scatter_particles.resize(n_proc);

        if (rank == 0) {
            for (int i = 0; i < incoming_move.size(); ++i) {
                int x = int(incoming_move[i].x / binSize);

                assert(incoming_move[i].x >= 0 && incoming_move[i].y >= 0 &&
                    incoming_move[i].x <= gridSize && incoming_move[i].y <= gridSize);

                int who = min(x / bins_per_proc, n_proc-1);
                scatter_particles[who].push_back(incoming_move[i]);

                int row = x % bins_per_proc;
                if (row == 0 && who != 0)
                    scatter_particles[who - 1].push_back(incoming_move[i]);
                if (row == bins_per_proc-1 && who != n_proc-1)
                    scatter_particles[who + 1].push_back(incoming_move[i]);
            }
            for (int i = 0; i < n_proc; ++i) {
                recv_counts[i] = scatter_particles[i].size();
            }
            displs[0] = 0;
            for (int i = 1; i < n_proc; ++i) {
                displs[i] = displs[i-1] + recv_counts[i-1];
            }
            // printf("worker: %d, 2. %d / %d.\n", rank, total_, displs[n_proc-1] + recv_counts[n_proc-1]);
            // assert(total_ == displs[n_proc-1] + recv_counts[n_proc-1]);
        }

        // printf("worker: %d. MPI_Scatter.\n", rank);
        send_count = 0;
        MPI_Scatter(recv_counts, 1, MPI_INT, &send_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        vector<particle_t> outgoing_move;
        outgoing_move.resize(send_count);

        vector<particle_t> scatter_particles_flatten;
        for (int i = 0; i < scatter_particles.size(); ++i) {
            scatter_particles_flatten.insert(scatter_particles_flatten.end(),
                scatter_particles[i].begin(), scatter_particles[i].end());
        }

        // printf("worker: %d. MPI_Scatterv.\n", rank);
        MPI_Scatterv(scatter_particles_flatten.data(), recv_counts, displs, PARTICLE, 
            outgoing_move.data(), send_count, PARTICLE, 0, MPI_COMM_WORLD);

        // int total__ = 0;
        // MPI_Reduce(&send_count, &total__, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        // if (rank == 0) {
        //     assert(total_ == total__);
        // }

        // printf("worker: %d. Bin.\n", rank);
        for (int i = 0; i < send_count; ++i) {
            particle_t &p = outgoing_move[i];
            assert(p.x >= 0 && p.y >= 0 && p.x <= gridSize && p.y <= gridSize);
            int x = p.x / binSize;
            int y = p.y / binSize;
            //printf("bin %d. x %d. y %d", x*bin_count + y, x, y);
            //printf(", size %ld.\n", bins[x*bin_count + y].size());
            particle_bins[x*binNum + y].push_back(p);
        }

        // //
        // //  move particles
        // //
        // for( int i = 0; i < nlocal; i++ )
        //     move( local[i] );

        if( find_option( argc, argv, "-no" ) == -1 )
        {

          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
 
          if (rank == 0){
            //
            // Computing statistical data
            //
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) absmin = rdmin;
          }
        }


    }
    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
      printf( "n = %d, simulation time = %g seconds", n, simulation_time);

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
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
  
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
