#include <iostream>
#include <random>
#include <cstdlib>
#include <thread>

#include <stdlib.h>
#include <time.h>

#include "matrix_strassen.h"

int main(int argc, char** argv) {
  /*Matrix x({{1, 2, 3}, {4, 5, 6}});
  Matrix y({{1, 4}, {2, 5}, {3, 6}});
  x.dump();
  y.dump();
  Matrix z(x*y);
  z.dump();
  return 0;*/
  srand (time(NULL));
  size_t m = 1024;
  size_t n = m;
  size_t k = m;
  size_t count = 10;
  Matrix a(m, n);
  Matrix b(n, k);
  for (size_t ctr = 0; ctr < count; ++ctr) {
    for (size_t i = 0; i < m; ++i) {
      for (size_t j = 0; j < m; ++j) {
          a.set(i, j, rand());
          b.set(i, j, rand());
        }
      }
    Matrix c = a*b;
  }
  return 0;
}
