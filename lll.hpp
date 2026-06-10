#pragma once

#include <vector>

using Matrix = std::vector<std::vector<double>>;

struct LLLResult {
  Matrix basis;    // LLL-редукований базис
  int swaps;       // кількість swap-ів (умова Ловаса)
  double time_ms;  // час роботи в мілісекундах
};

// Повертає LLL-редукований базис разом зі статистикою.
// B     - матриця, рядки якої - базисні вектори
// delta - параметр Ловаса, має бути в (0.25, 1)
LLLResult lll_reduce(Matrix B, double delta);
