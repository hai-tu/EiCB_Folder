#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <assert.h>

#include <omp.h>

#define STATESIZE 64

// returns random float in range -1 to 1
double rand_dbl(struct random_data *random_buf) {
  int rand_int;
  random_r(random_buf, &rand_int);
  return (double)rand_int / RAND_MAX * 2.0 - 1.0;
}

int main(int argc, char *argv[]) {
  // reading the number of points to use from the argument
  unsigned long num_points = 100000;
  if (argc > 1) {
    num_points = atol(argv[1]);
  }
  double x, y, circle_points = 0, square_points = 0;

  printf("Enter the number of cases: ");
  scanf("%ld", &num_points);
  //printf("%d\n", a);
  #pragma omp parallel
  {
  // be aware that each thread should have its own random state
  struct random_data *random_buf = calloc(1, sizeof(struct random_data));
  char *state_buf = calloc(STATESIZE, sizeof(char));
  // to be clear: random_data and state_buf are both part of the random state
  assert(random_buf != NULL && state_buf != NULL);
  // threads should use a different random seed
  unsigned int seed = time(NULL);
  initstate_r(seed, state_buf, STATESIZE, random_buf);
  

  // you may get a random point like this:
  #pragma omp for private(x, y) reduction(+:circle_points)//reduction(+:square_points) 
  for(int i=0;i<num_points;i++){
  	x = rand_dbl(random_buf);
  	y = rand_dbl(random_buf);

  	if (x*x + y*y <= 1)
  		circle_points += 1;

  	//square_points += 1;

  }
  free(random_buf);
  free(state_buf);

  }
  //double x = rand_dbl(random_buf);
  //double y = rand_dbl(random_buf);

  //double pi = 3141; // replace this with the actual computed value

  double pi = (4*circle_points)/num_points;
  //pi = (4*circle_points)/square_points;
  printf("Pi is approx. %f\n", pi);

  //free(random_buf);
  //free(state_buf);

  return 0;
}
