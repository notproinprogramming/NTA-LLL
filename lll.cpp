#include "lll.hpp"

#include <chrono>
#include <cmath>
#include <stdexcept>

// Скалярний добуток двох векторів
static double dot(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0.0;
  for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// Норма вектора в квадраті
static double norm2(const std::vector<double>& v) { return dot(v, v); }

// Побудова ортогонального базису Грама-Шмідта.
// gso[i] - i-й ортогоналізований вектор,
// mu[i][j] - коефіцієнти Грама-Шмідта (i > j)
static void gram_schmidt(const Matrix& B, Matrix& gso, Matrix& mu) {
  int m = (int)B.size();
  gso = B;
  mu.assign(m, std::vector<double>(m, 0.0));

  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < i; ++j) {
      double n2 = norm2(gso[j]);
      if (n2 < 1e-10) continue;
      mu[i][j] = dot(B[i], gso[j]) / n2;
      for (size_t k = 0; k < B[i].size(); ++k) gso[i][k] -= mu[i][j] * gso[j][k];
    }
  }
}

// Перерахунок одного рядка GSO після зміни B[i].
static void update_gso_row(const Matrix& B, Matrix& gso, Matrix& mu, int i) {
  int m = (int)B.size();
  gso[i] = B[i];
  for (int j = 0; j < i; ++j) {
    double n2 = norm2(gso[j]);
    if (n2 < 1e-10) {
      mu[i][j] = 0.0;
      continue;
    }
    mu[i][j] = dot(B[i], gso[j]) / n2;
    for (size_t k = 0; k < B[i].size(); ++k) gso[i][k] -= mu[i][j] * gso[j][k];
  }
  // оновити mu[l][i] для рядків нижче
  for (int l = i + 1; l < m; ++l) {
    double n2 = norm2(gso[i]);
    if (n2 < 1e-10) {
      mu[l][i] = 0.0;
      continue;
    }
    mu[l][i] = dot(B[l], gso[i]) / n2;
  }
}

LLLResult lll_reduce(Matrix B, double delta) {
  if (delta <= 0.25 || delta >= 1.0) throw std::invalid_argument("delta must be in (0.25, 1)");

  int m = (int)B.size();

  LLLResult res;
  res.swaps = 0;

  if (m == 0) {
    res.basis = B;
    res.time_ms = 0.0;
    return res;
  }

  auto t0 = std::chrono::high_resolution_clock::now();

  Matrix gso, mu;
  gram_schmidt(B, gso, mu);  // крок 1 будуємо базис Грама-Шмідта

  int k = 1;
  while (k < m) {
    // крок 2.1 - умова мінімізації проекцій
    for (int j = k - 1; j >= 0; --j) {
      if (std::abs(mu[k][j]) > 0.5) {
        long long q = (long long)std::round(mu[k][j]);
        for (size_t i = 0; i < B[k].size(); ++i) B[k][i] -= q * B[j][i];
        update_gso_row(B, gso, mu, k);
      }
    }

    // крок 2.2 - умова Ловаса
    double lhs = norm2(gso[k]);
    double rhs = (delta - mu[k][k - 1] * mu[k][k - 1]) * norm2(gso[k - 1]);

    if (lhs < rhs) {
      std::swap(B[k - 1], B[k]);
      update_gso_row(B, gso, mu, k - 1);
      update_gso_row(B, gso, mu, k);
      ++res.swaps;
      k = std::max(k - 1, 1);
    } else {
      k++;
    }
  }

  auto t1 = std::chrono::high_resolution_clock::now();
  res.time_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
  res.basis = std::move(B);
  return res;
}
