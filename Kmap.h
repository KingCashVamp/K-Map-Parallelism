/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

#include <omp.h>
#include <cstdint>
#include <vector>

#define MAX_DIM = 6;

struct Term {
  /* Define the data structure for min or max Term here. */ 
  int x, y, z, i , j, k;
};
