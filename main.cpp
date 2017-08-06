#include <iostream>
#include <random>
#include <cstdlib>
#include <thread>
#include <chrono>
#include <vector>

#include <stdlib.h>
#include <time.h>

#include "matrix_strassen.h"

int main(int argc, char** argv) {
  /*Matrix x({{0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}});
  Matrix y({{1, 2, 1}, {0, 3, 2}, {1, 0, 3}, {2, 1, 0}});
  x.dump();
  std::cout << "*" << std::endl;
  y.dump();
  Matrix z(x*y);
  std::cout << "=" << std::endl;
  z.dump();
  return 0;*/
  srand (time(NULL));
  std::vector<size_t> sizes({32, 64, 100, 128, 256, 512, 1024});
  for (auto it = sizes.begin(); it != sizes.end(); ++it) {
    size_t m = *it;
    size_t n = m;
    size_t k = m;
    size_t count = 100;
    Matrix a(m, n);
    Matrix b(n, k);
    std::chrono::microseconds sum(0);
    for (size_t ctr = 0; ctr < count; ++ctr) {
      for (size_t i = 0; i < m; ++i) {
        for (size_t j = 0; j < m; ++j) {
            a.set(i, j, rand());
            b.set(i, j, rand());
          }
        }
      std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
      Matrix c = a*b;
      std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
      sum += std::chrono::duration_cast<std::chrono::microseconds> (end - start);
    }
    double multiplication_time_ms = static_cast<double>(sum.count()) / (count * 1000);
    std::cout << "M = " << m << ": time = " << multiplication_time_ms << std::endl;
  }
  return 0;
}
