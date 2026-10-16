#include <limits>
#include <memory>
#include <toppra/algorithm/toppra.hpp>
#include <toppra/constraint.hpp>
#include <toppra/constraint/linear_joint_acceleration.hpp>
#include <toppra/constraint/linear_joint_velocity.hpp>
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

/// Bounds u at every stage, by 0.75 or by a value that changes with the stage.
class StageControlLimit : public LinearConstraint {
 public:
  explicit StageControlLimit(bool constant)
      : LinearConstraint(1, 1, constant, false, false) {}

 private:
  void computeParams_impl(const GeometricPath&, const Vector& times, Vectors& a,
                          Vectors& b, Vectors& c, Matrices& F, Vectors& g, Bounds&,
                          Bounds&) override {
    for (std::size_t i = 0; i < a.size(); ++i) {
      a[i].setOnes();
      b[i].setZero();
      c[i].setZero();
    }
    for (std::size_t i = 0; i < F.size(); ++i) {
      F[i].setOnes();
      g[i].setConstant(constantF() ? 0.75 : (i + 1.0) / times.size());
    }
  }
};

/// Gives the tests access to the constraint records.
class InspectableSeidel : public solver::Seidel {
 public:
  std::size_t linearCount() const { return m_constraintsParams.lin.size(); }
  std::size_t boxCount() const { return m_constraintsParams.box.size(); }
};

std::shared_ptr<PiecewisePolyPath> straightPath() {
  Matrix coefficients(2, 1);
  coefficients << 1, 0;
  return std::make_shared<PiecewisePolyPath>(Matrices{coefficients},
                                             std::vector<value_type>{0, 1});
}

/// A path with as many joints as given, which changes with the trial.
std::shared_ptr<PiecewisePolyPath> curvedPath(int joints, int trial) {
  Matrix coefficients(3, joints);
  for (int j = 0; j < joints; ++j)
    coefficients.col(j) << 0.1 * (j + 1), 1 + 0.2 * j, 0.1 * trial;
  return std::make_shared<PiecewisePolyPath>(Matrices{coefficients},
                                             std::vector<value_type>{0, 1});
}

/// A straight path of length 1 with acceleration limit 1.
std::unique_ptr<algorithm::TOPPRA> straightProblem(value_type velocity = 2) {
  LinearConstraintPtrs constraints{
      std::make_shared<constraint::LinearJointVelocity>(Vector::Constant(1, -velocity),
                                                        Vector::Constant(1, velocity)),
      std::make_shared<constraint::LinearJointAcceleration>(Vector::Constant(1, -1),
                                                            Vector::Constant(1, 1))};
  auto problem = std::make_unique<algorithm::TOPPRA>(constraints, straightPath());
  problem->solver(std::make_shared<solver::Seidel>());
  problem->setGridpoints(Vector::LinSpaced(51, 0, 1));
  return problem;
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

// Every initialization appended the constraint records to those of the
// previous one. A re-initialized solver solved with each constraint several
// times, and read beyond the records of a different grid.
TEST(SolverInitialization, ReplacesChangedPathsGridsDimensionsAndConstraints) {
  auto reused = std::make_shared<InspectableSeidel>();
  for (int trial = 0; trial < 12; ++trial) {
    SCOPED_TRACE(trial);
    const int joints = 1 + trial % 3;
    const int accelerationConstraints = trial % 3;
    const auto discretization = trial % 2 ? DiscretizationType::Interpolation
                                          : DiscretizationType::Collocation;
    const Vector velocity = Vector::Constant(joints, 1 + 0.1 * trial);
    LinearConstraintPtrs constraints{
        std::make_shared<constraint::LinearJointVelocity>(-velocity, velocity)};
    for (int i = 0; i < accelerationConstraints; ++i) {
      const Vector acceleration = Vector::Constant(joints, 0.5 + 0.2 * trial + i);
      constraints.push_back(std::make_shared<constraint::LinearJointAcceleration>(
          -acceleration, acceleration));
    }
    for (auto& constraint : constraints) constraint->discretizationType(discretization);
    const auto path = curvedPath(joints, trial);
    const Vector gridpoints = Vector::LinSpaced(trial % 2 ? 9 : 33, 0, 1);

    algorithm::TOPPRA reusing(constraints, path);
    reusing.solver(reused);
    reusing.setGridpoints(gridpoints);
    ASSERT_EQ(reusing.computePathParametrization(0, 0), ReturnCode::OK)
        << reusing.getErrorMessage();
    EXPECT_EQ(reused->linearCount(), static_cast<std::size_t>(accelerationConstraints));
    EXPECT_EQ(reused->boxCount(), 1u);

    algorithm::TOPPRA fresh(constraints, path);
    fresh.solver(std::make_shared<solver::Seidel>());
    fresh.setGridpoints(gridpoints);
    ASSERT_EQ(fresh.computePathParametrization(0, 0), ReturnCode::OK);
    EXPECT_EQ(reusing.getParameterizationData().parametrization,
              fresh.getParameterizationData().parametrization);
  }
}

TEST(SolverInitialization, ReplacesBoxBoundsOfAnotherKind) {
  InspectableSeidel solver;
  const auto path = straightPath();
  const Vector times = Vector::LinSpaced(5, 0, 1);
  solver.initialize({std::make_shared<BoxBound>(false, 0, 1)}, path, times);
  solver.initialize({std::make_shared<BoxBound>(true, -2, 2)}, path, times);
  EXPECT_NEAR(maximizeControl(solver)[0], 2, 1e-12);
  // The x bounds [2, 3] exclude x = 4, and no control bound is left.
  solver.initialize({std::make_shared<BoxBound>(false, 2, 3)}, path, times);
  Vector g(2), solution(2);
  g << -1, 0;
  EXPECT_FALSE(solver.solveStagewiseOptim(0, Matrix(), g, Bound{4, 4}, Bound{0, 10},
                                          solution));
  ASSERT_TRUE(solver.solveStagewiseOptim(0, Matrix(), g, Bound{2.5, 2.5}, Bound{0, 10},
                                         solution));
  EXPECT_NEAR(solution[0], 15, 1e-12);
  EXPECT_EQ(solver.boxCount(), 1u);
  EXPECT_EQ(solver.linearCount(), 0u);
}

TEST(SolverInitialization, ReplacesConstantAndStagewiseLinearRowsAfterGridChanges) {
  InspectableSeidel solver;
  const auto path = straightPath();
  for (int trial = 0; trial < 8; ++trial) {
    SCOPED_TRACE(trial);
    const bool constant = trial % 2;
    const int count = constant ? 9 : 33;
    auto limit = std::make_shared<StageControlLimit>(constant);
    limit->discretizationType(DiscretizationType::Collocation);
    solver.initialize({limit}, path, Vector::LinSpaced(count, 0, 1));
    Vector g(2), solution(2);
    g << -1, 0;
    ASSERT_TRUE(solver.solveStagewiseOptim(count - 2, Matrix(), g, Bound{0, 0},
                                           Bound{0, 10}, solution));
    EXPECT_NEAR(solution[0], constant ? 0.75 : (count - 1.0) / count, 1e-12);
    EXPECT_EQ(solution[1], 0);
    EXPECT_EQ(solver.linearCount(), 1u);
  }
}

TEST(SolverInitialization, RejectsInvalidInputsBeforeReplacingTheProblem) {
  InspectableSeidel solver;
  const auto path = straightPath();
  const Vector times = Vector::LinSpaced(5, 0, 1);
  const LinearConstraintPtrs constraints{std::make_shared<BoxBound>(true, -2, 2)};
  solver.initialize(constraints, path, times);

  Vector nan = times, repeated = times;
  nan[2] = std::numeric_limits<value_type>::quiet_NaN();
  repeated[2] = repeated[1];
  for (const Vector& invalid : {Vector(), nan, repeated})
    EXPECT_THROW(solver.initialize(constraints, path, invalid), std::invalid_argument);
  EXPECT_THROW(solver.initialize(constraints, nullptr, times), std::invalid_argument);
  EXPECT_THROW(solver.initialize({nullptr}, path, times), std::invalid_argument);

  EXPECT_EQ(solver.nbStages(), 4u);
  EXPECT_NEAR(maximizeControl(solver)[0], 2, 1e-12);
}

// setGridpoints changed the number of stages before validating the grid, so
// that a rejected grid left a problem whose grid and stage count disagree.
TEST(AlgorithmReinitialization, RejectedGridsLeaveTheProblemUnchanged) {
  auto problem = straightProblem();
  ASSERT_EQ(problem->computePathParametrization(0, 0), ReturnCode::OK);
  const Vector expected = problem->getParameterizationData().parametrization;

  Vector nan = Vector::LinSpaced(9, 0, 1), infinite = nan, repeated = nan, decreasing = nan;
  nan[4] = std::numeric_limits<value_type>::quiet_NaN();
  infinite[4] = std::numeric_limits<value_type>::infinity();
  repeated[4] = repeated[3];
  decreasing[4] = -1;
  const Vector otherInterval = Vector::LinSpaced(9, 0.1, 1);
  for (const Vector& rejected :
       {Vector(), Vector(Vector::Zero(1)), nan, infinite, repeated, decreasing,
        otherInterval}) {
    EXPECT_THROW(problem->setGridpoints(rejected), std::invalid_argument);
    ASSERT_EQ(problem->computePathParametrization(0, 0), ReturnCode::OK)
        << problem->getErrorMessage();
    EXPECT_EQ(problem->getParameterizationData().parametrization, expected);
  }
}

// The solver setter kept the algorithm initialized, so that a solver set
// after a computation was used without ever being initialized.
TEST(AlgorithmReinitialization, InitializesAReplacedSolver) {
  auto problem = straightProblem();
  ASSERT_EQ(problem->computePathParametrization(0, 0), ReturnCode::OK);
  const Vector expected = problem->getParameterizationData().parametrization;
  problem->solver(std::make_shared<solver::Seidel>());
  ASSERT_EQ(problem->computePathParametrization(0, 0), ReturnCode::OK)
      << problem->getErrorMessage();
  EXPECT_EQ(problem->getParameterizationData().parametrization, expected);
}

// An algorithm initialized its solver only once. When another algorithm
// initialized the same solver, the first one kept solving on the other's
// tables.
TEST(AlgorithmReinitialization, InitializesASharedSolverAgainAfterAnotherAlgorithm) {
  auto shared = std::make_shared<solver::Seidel>();
  auto fast = straightProblem(2), slow = straightProblem(0.5);
  fast->solver(shared);
  slow->solver(shared);
  ASSERT_EQ(fast->computePathParametrization(0, 0), ReturnCode::OK);
  const Vector expected = fast->getParameterizationData().parametrization;
  ASSERT_EQ(slow->computePathParametrization(0, 0), ReturnCode::OK);
  ASSERT_NE(slow->getParameterizationData().parametrization, expected);
  ASSERT_EQ(fast->computePathParametrization(0, 0), ReturnCode::OK);
  EXPECT_EQ(fast->getParameterizationData().parametrization, expected);
}
