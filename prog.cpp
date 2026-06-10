#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "lll.hpp"

// Коефіцієнт Адамара: відношення об'єму паралелепіпеда до добутку норм векторів
// Значення близьке до 1 означає, що базис близький до ортогонального
static double hadamard_ratio(const Matrix& B) {
  int m = (int)B.size();
  double vol2 = 1.0;

  Matrix gso = B;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < i; ++j) {
      double n2 = 0.0;
      for (double x : gso[j]) n2 += x * x;
      if (n2 < 1e-10) continue;
      double mu = 0.0;
      for (size_t k = 0; k < B[i].size(); ++k) mu += B[i][k] * gso[j][k];
      mu /= n2;
      for (size_t k = 0; k < gso[i].size(); ++k) gso[i][k] -= mu * gso[j][k];
    }
    double n2 = 0.0;
    for (double x : gso[i]) n2 += x * x;
    vol2 *= n2;
  }

  double prod_norms2 = 1.0;
  for (const auto& row : B) {
    double n2 = 0.0;
    for (double x : row) n2 += x * x;
    prod_norms2 *= n2;
  }

  if (prod_norms2 < 1e-30) return 0.0;
  // Hadamard = (det / prod||b_i||)^(1/m) = (vol2 / prod_norms2)^(1/(2m))
  return std::pow(vol2 / prod_norms2, 1.0 / (2.0 * m));
}

// Генерує випадкову цілочисельну повнорангову матрицю n x n.
// Елементи з діапазону [-range, range].
static Matrix random_full_rank(int n, int range, std::mt19937& rng) {
  std::uniform_int_distribution<int> dist(-range, range);
  Matrix B;
  while (true) {
    B.assign(n, std::vector<double>(n));
    for (auto& row : B)
      for (auto& x : row) x = dist(rng);

    // перевірка повного рангу: GSO не має нульових векторів
    Matrix gso = B;
    bool full_rank = true;
    for (int i = 0; i < n && full_rank; ++i) {
      for (int j = 0; j < i; ++j) {
        double n2 = 0.0;
        for (double x : gso[j]) n2 += x * x;
        if (n2 < 1e-10) continue;
        double mu = 0.0;
        for (int k = 0; k < n; ++k) mu += B[i][k] * gso[j][k];
        mu /= n2;
        for (int k = 0; k < n; ++k) gso[i][k] -= mu * gso[j][k];
      }
      double n2 = 0.0;
      for (double x : gso[i]) n2 += x * x;
      if (n2 < 1e-10) {
        full_rank = false;
      }
    }
    if (full_rank) return B;
    // матриця виявилась виродженою - генеруємо ще раз
  }
}

int main() {
  const int N_MATRICES = 50;
  const int N = 30;
  const int RANGE = 10;  // елементи з [-10, 10]
  const double DELTAS[] = {0.50, 0.75, 0.90, 0.95, 0.99};
  const int N_DELTAS = (int)(sizeof(DELTAS) / sizeof(DELTAS[0]));

  std::mt19937 rng(42);

  std::cout << std::fixed << std::setprecision(4);
  std::cout << std::setw(4) << "mat" << std::setw(6) << "delta" << std::setw(10) << "time(ms)" << std::setw(8) << "swaps" << std::setw(12)
            << "hadamard" << std::setw(14) << "||b1||"
            << "\n";
  std::cout << std::string(54, '-') << "\n";

  for (int mi = 0; mi < N_MATRICES; ++mi) {
    Matrix B = random_full_rank(N, RANGE, rng);

    for (int di = 0; di < N_DELTAS; ++di) {
      double delta = DELTAS[di];

      try {
        LLLResult r = lll_reduce(B, delta);

        double h = hadamard_ratio(r.basis);
        double nb1 = 0.0;
        for (double x : r.basis[0]) nb1 += x * x;
        nb1 = std::sqrt(nb1);

        std::cout << std::setw(4) << (mi + 1) << std::setw(6) << delta << std::setw(10) << r.time_ms << std::setw(8) << r.swaps << std::setw(12) << h
                  << std::setw(14) << nb1 << "\n";
      } catch (const std::exception& e) {
        std::cerr << "mat=" << (mi + 1) << " delta=" << delta << " error: " << e.what() << "\n";
      }
    }
  }

  return 0;
}
