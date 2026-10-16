#include <memory>
#include <toppra/constraint.hpp>
#include <toppra/geometric_path/piecewise_poly_path.hpp>
#include <toppra/solver/seidel.hpp>
#include <toppra/toppra.hpp>

#include "gtest/gtest.h"

namespace {
using namespace toppra;

/// Bounds either u or x by [low, high] at every grid point.
class BoxBound : public LinearConstraint {
 public:
  BoxBound(bool control, value_type low, value_type high)
      : LinearConstraint(0, 0, true, control, !control),
        m_control(control), m_bound(low, high) {}

 private:
  void computeParams_impl(const GeometricPath&, const Vector&, Vectors&, Vectors&,
                          Vectors&, Matrices&, Vectors&, Bounds& ubound,
                          Bounds& xbound) override {
    for (auto& bound : m_control ? ubound : xbound) bound = m_bound;
  }

  bool m_control;
  Bound m_bound;
};

std::shared_ptr<PiecewisePolyPath> straightPath() {
  Matrix coefficients(2, 1);
  coefficients << 1, 0;
  return std::make_shared<PiecewisePolyPath>(Matrices{coefficients},
                                             std::vector<value_type>{0, 1});
}

/// Maximizes u at stage 0, where x = 4 and the next x is at most 10. Without
/// other constraints, u reaches (10 - 4) / (2 * 0.25) = 12.
Vector maximizeControl(Solver& solver) {
  Vector g(2), solution(2);
  g << -1, 0;
  EXPECT_TRUE(solver.solveStagewiseOptim(0, Matrix(), g, Bound{4, 4}, Bound{0, 10},
                                         solution));
  return solution;
}
}  // namespace

TEST(SolverInitialization, HonorsControlBounds) {
  solver::Seidel solver;
  solver.initialize({std::make_shared<BoxBound>(true, -2, 2)}, straightPath(),
                    Vector::LinSpaced(5, 0, 1));
  const Vector solution = maximizeControl(solver);
  EXPECT_NEAR(solution[0], 2, 1e-12);
  EXPECT_EQ(solution[1], 4);
}
