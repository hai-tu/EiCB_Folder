#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <string>

using namespace std;

/* ************************************************************************ */
/* helper function: freeMatrix: frees memory of the matrix                  */
/* ************************************************************************ */
static void freeMatrix(double **Matrix) {
  free(Matrix[0]);
  free(Matrix);
}

/* ************************************************************************ */
/* helper function: allocateMatrix: allocates memory for a matrix           */
/* ************************************************************************ */
static double **allocateMatrix(int x, int y) {
  int i;

  double *data_layer = (double *)malloc(x * y * sizeof(double));
  double **resultmatrix = (double **)malloc(x * sizeof(double *));

  if (data_layer == NULL || resultmatrix == NULL) {
    printf("ERROR ALLOCATING MEMORY\n");
    return NULL;
  }

  for (i = 0; i < x; i++) {
    resultmatrix[i] = data_layer + i * y;
  }

  return resultmatrix;
}

int main(int argc, char *argv[])

//
//  Purpose:
//
//    MAIN is the main program for HEATED_PLATE_OPENMP.
//
//  Discussion:
//
//    This code solves the steady state heat equation on a rectangular region.
//
//    The sequential version of this program needs approximately
//    18/epsilon iterations to complete.
//
//
//    The physical region, and the boundary conditions, are suggested
//    by this diagram;
//
//                   W = 0
//             +------------------+
//             |                  |
//    W = 100  |                  | W = 100
//             |                  |
//             +------------------+
//                   W = 100
//
//    The region is covered with a grid of M by N nodes, and an M by N
//    array W is used to record the temperature.  The correspondence between
//    array indices and locations in the region is suggested by giving the
//    indices of the four corners:
//
//                  I = 0
//          [0][0]-------------[0][N-1]
//             |                  |
//      J = 0  |                  |  J = N-1
//             |                  |
//        [M-1][0]-----------[M-1][N-1]
//                  I = M-1
//
//    The steady state solution to the discrete heat equation satisfies the
//    following condition at an interior grid point:
//
//      W[Central] = (1/4) * ( W[North] + W[South] + W[East] + W[West] )
//
//    where "Central" is the index of the grid point, "North" is the index
//    of its immediate neighbor to the "north", and so on.
//
//    Given an approximate solution of the steady state heat equation, a
//    "better" solution is given by replacing each interior point by the
//    average of its 4 neighbors - in other words, by using the condition
//    as an ASSIGNMENT statement:
//
//      W[Central]  <=  (1/4) * ( W[North] + W[South] + W[East] + W[West] )
//
//    If this process is repeated often enough, the difference between
//    successive estimates of the solution will go to zero.
//
//    This program carries out such an iteration, using a tolerance specified by
//    the user.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    24 October 2019
//
//  Author:
//
//    Original C version by Michael Quinn.
//    C++ version by John Burkardt, slightly adapted by Tim Jammer.
//
//  Reference:
//
//    Michael Quinn,
//    Parallel Programming in C with MPI and OpenMP,
//    McGraw-Hill, 2004,
//    ISBN13: 978-0071232654,
//    LC: QA76.73.C15.Q55.
//
//  Local parameters:
//
//    Local, double DIFF, the norm of the change in the solution from one
//    iteration to the next.
//
//    Local, double MEAN, the average of the boundary values, used to initialize
//    the values of the solution in the interior.
//
//    Local, double U[M][N], the solution at the previous iteration.
//
//    Local, double W[M][N], the solution computed at the latest iteration.
//
{
#define M 500
#define N 500

  omp_set_num_threads(6);
  double diff;
  double epsilon = 0.001;
  int i;
  int iterations;
  int iterations_print;
  int j;
  double mean;
  double **u = allocateMatrix(M, N);
  double **w = allocateMatrix(M, N);
  double wtime;

  cout << "\n";
  cout << "HEATED_PLATE_OPENMP\n";
  cout << "  C++/OpenMP version\n";
  cout
      << "  A program to solve for the steady state temperature distribution\n";
  cout << "  over a rectangular plate.\n";
  cout << "\n";
  cout << "  Spatial grid of " << M << " by " << N << " points.\n";
  cout << "  The iteration will be repeated until the change is <= " << epsilon
       << "\n";
  cout << "  Number of processors available = " << omp_get_num_procs() << "\n";
  cout << "  Number of threads =              " << omp_get_max_threads()
       << "\n";
  //
  //  Set the boundary values, which don't change.
  //
  
  mean = 0.0;

  #pragma omp parallel
  {
  #pragma omp for
  for (i = 1; i < M - 1; i++) {
    w[i][0] = 100.0;
  }
  #pragma omp for
  for (i = 1; i < M - 1; i++) {
    w[i][N - 1] = 100.0;
  }
  #pragma omp for
  for (j = 0; j < N; j++) {
    w[M - 1][j] = 100.0;
  }
  #pragma omp for
  for (j = 0; j < N; j++) {
    w[0][j] = 0.0;
  }
  //}
  //
  //  Average the boundary values, to come up with a reasonable
  //  initial value for the interior.
  //
  #pragma omp for reduction(+: mean)
  for (i = 1; i < M - 1; i++) {
    mean = mean + w[i][0] + w[i][N - 1];
  }
  #pragma omp for reduction(+: mean)
  for (j = 0; j < N; j++) {
    mean = mean + w[M - 1][j] + w[0][j];
  }
  }

  mean = mean / (double)(2 * M + 2 * N - 4);
  cout << "\n";
  cout << "  MEAN = " << mean << "\n";
  //
  //  Initialize the interior solution to the mean value.
  //

  #pragma omp for private(j)
  for (i = 1; i < M - 1; i++) {
    for (j = 1; j < N - 1; j++) {
      w[i][j] = mean;
    }
  }
  //
  //  iterate until the  new solution W differs from the old solution U
  //  by no more than EPSILON.
  //
  iterations = 0;
  iterations_print = 1;
  cout << "\n";
  cout << " Iteration  Change\n";
  cout << "\n";
  wtime = omp_get_wtime();
  diff = epsilon;

  while (epsilon <= diff) {
    //
    //  Save the old solution in U.
    //
    #pragma omp parallel for private(j)
    for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++) {
        u[i][j] = w[i][j];
      }
    }
    //
    //  Determine the new estimate of the solution at the interior points.
    //  The new solution W is the average of north, south, east and west
    //  neighbors.
    #pragma omp parallel for private(j)
    for (i = 1; i < M - 1; i++) {
      for (j = 1; j < N - 1; j++) {
        w[i][j] = (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1]) / 4.0;
      }
    }

    //  This initialization must be executed by one thread only
    //#pragma omp single
    diff = 0.0;

    #pragma omp parallel for private(j) reduction(max: diff)
    for (i = 1; i < M - 1; i++) {
      for (j = 1; j < N - 1; j++) {
        if (diff < fabs(w[i][j] - u[i][j])) {
          diff = fabs(w[i][j] - u[i][j]);
        }
      }
    }
    //}

    iterations++;
    if (iterations == iterations_print) {
      cout << "  " << setw(8) << iterations << "  " << diff << "\n";
      iterations_print = 2 * iterations_print;
    }
  }
  //}
  wtime = omp_get_wtime() - wtime;

  cout << "\n";
  cout << "  " << setw(8) << iterations << "  " << diff << "\n";
  cout << "\n";
  cout << "  Error tolerance achieved.\n";
  cout << "  Wallclock time = " << wtime << "\n";
  //
  //  Terminate.
  //
  freeMatrix(u);
  freeMatrix(w);

  cout << "\n";
  cout << "HEATED_PLATE_OPENMP:\n";
  cout << "  Normal end of execution.\n";

  return 0;

#undef M
#undef N
}
