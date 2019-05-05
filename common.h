#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include <vector>
#include <math.h>
using namespace std;

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

// struct bin_v
// {
//   double binSize;
//   double gridSize;
//   int binNum;
// };


//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (0.001*cutoff)
#define dt      0.0005

//
// particle data structure
//
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;

//
//  timing routines
//
double read_timer();
typedef std::vector<particle_t> bin_t;

//
//  simulation routines
//
double set_size( int n );
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
// void apply_force( particle_t &particle, particle_t &neighbor);

void move( particle_t &p );
//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );
//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

///////////////////////////////// Grid declaration /////////////////////////////
typedef struct linkedlist_t // typedef is used to give a symbolic name for the existing name
// in a c program. This is similar to defining alias for the command.
// 
{
	linkedlist_t * next;
	particle_t * value;
}linkedlist_t; 


typedef struct grid_t
{
	int size;
	// the reason for **grid is the fact that grid is a two dimensional array
	linkedlist_t ** grid;
}grid_t;


//
// grid routines
//

void grid_init(grid_t & grid1, int gridsize);
void grid_add(grid_t & grid1, particle_t * particle);
// void buildBins(vector<bin_t> &bins, particle_t* particles, int n, double gridSize);
//
// Calculate the grid coordinate from a real coordinate
//
inline static int grid_coord(double c)
{
    return (int)floor(c / cutoff);
}
inline static int grid_coord_flat(int size, double x, double y)
{
    return grid_coord(x) * size + grid_coord(y);
}

#endif