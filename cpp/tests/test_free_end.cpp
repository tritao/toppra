#include <cmath>
#include <limits>
#include <memory>
#include <toppra/algorithm/toppra.hpp>
#include <toppra/constraint/linear_joint_acceleration.hpp>
#include <toppra/constraint/linear_joint_velocity.hpp>
#include <toppra/geometric_path/piecewise_poly_path.hpp>
#include <toppra/solver/seidel.hpp>

#include "gtest/gtest.h"

namespace {
using namespace toppra;

constexpr DiscretizationType kSchemes[] = {DiscretizationType::Collocation,
                                           DiscretizationType::Interpolation};

/// Times a joint moving along a straight path s in [begin, end], with velocity
/// limit 2 (x <= 4) and acceleration limit 1, on `gridpoints` gridpoints.
std::unique_ptr<algorithm::TOPPRA> line(value_type begin, value_type end,
                                        DiscretizationType scheme,
                                        int gridpoints = 51) {
  Matrix coefficients(2, 1);
  coefficients << 1, begin;
  auto path = std::make_shared<PiecewisePolyPath>(Matrices{coefficients},
                                                  std::vector<value_type>{begin, end});
  LinearConstraintPtrs constraints{
      std::make_shared<constraint::LinearJointVelocity>(Vector::Constant(1, -2),
                                                        Vector::Constant(1, 2)),
      std::make_shared<constraint::LinearJointAcceleration>(Vector::Constant(1, -1),
                                                            Vector::Constant(1, 1))};
  for (auto& constraint : constraints) constraint->discretizationType(scheme);
  auto problem = std::make_unique<algorithm::TOPPRA>(constraints, path);
  problem->solver(std::make_shared<solver::Seidel>());
  problem->setInitialXBounds(Bound{0, 100});
  problem->setGridpoints(Vector::LinSpaced(gridpoints, begin, end));
  return problem;
}

/// Times q1 = s + 0.2 s^2, q2 = 0.5 s^2 - 0.1 s^3 for s in [begin, end], with
/// unit velocity and acceleration limits, on a grid of spacing 0.01.
std::unique_ptr<algorithm::TOPPRA> curve(value_type begin, value_type end,
                                         DiscretizationType scheme) {
  // Coefficients about the window's own origin, highest degree first.
  Matrix coefficients(4, 2);
  coefficients << 0, -0.1,
                  0.2, 0.5 - 0.3 * begin,
                  1 + 0.4 * begin, begin - 0.3 * begin * begin,
                  begin + 0.2 * begin * begin, 0.5 * begin * begin - 0.1 * std::pow(begin, 3);
  auto path = std::make_shared<PiecewisePolyPath>(Matrices{coefficients},
                                                  std::vector<value_type>{begin, end});
  LinearConstraintPtrs constraints{
      std::make_shared<constraint::LinearJointVelocity>(Vector::Constant(2, -1),
                                                        Vector::Constant(2, 1)),
      std::make_shared<constraint::LinearJointAcceleration>(Vector::Constant(2, -1),
                                                            Vector::Constant(2, 1))};
  for (auto& constraint : constraints) constraint->discretizationType(scheme);
  auto problem = std::make_unique<algorithm::TOPPRA>(constraints, path);
  problem->solver(std::make_shared<solver::Seidel>());
  problem->setInitialXBounds(Bound{0, 100});
  const int stages = static_cast<int>(std::lround((end - begin) * 100));
  problem->setGridpoints(Vector::LinSpaced(stages + 1, begin, end));
  return problem;
}

const Vector& profile(const algorithm::TOPPRA& problem) {
  return problem.getParameterizationData().parametrization;
}

value_type maxDifference(const Vector& a, const Vector& b) {
  EXPECT_EQ(a.size(), b.size());
  return a.size() == b.size() ? (a - b).cwiseAbs().maxCoeff()
                              : std::numeric_limits<value_type>::infinity();
}
}  // namespace

TEST(FreeEnd, ReachesAFinalSpeedInsideTheIntervalRatherThanItsUpperBound) {
  for (auto scheme : kSchemes) {
    auto problem = line(0, 0.1, scheme);
    ASSERT_EQ(problem->computePathParametrizationFromSquaredSpeed(0, Bound{0, 4}),
              ReturnCode::OK);
    // Accelerating at the limit over 0.1 reaches x = 2 a s = 0.2, below 4.
    const auto& data = problem->getParameterizationData();
    EXPECT_NEAR(profile(*problem)[50], 0.2, 1e-12);
    for (int i = 0; i < 50; ++i) {
      const auto u = (data.parametrization[i + 1] - data.parametrization[i]) /
                     (2 * (data.gridpoints[i + 1] - data.gridpoints[i]));
      EXPECT_NEAR(u, 1, 1e-10) << "stage " << i;
    }
  }
}

TEST(FreeEnd, MatchesTheFixedEndParametrizationForAOnePointInterval) {
  for (auto scheme : kSchemes) {
    auto legacy = line(0, 1, scheme), squared = line(0, 1, scheme);
    ASSERT_EQ(legacy->computePathParametrization(0.25, 0.5), ReturnCode::OK);
    ASSERT_EQ(squared->computePathParametrizationFromSquaredSpeed(0.25, Bound{0.25, 0.25}),
              ReturnCode::OK);
    EXPECT_LE(maxDifference(profile(*legacy), profile(*squared)), 1e-12);
    EXPECT_EQ(profile(*squared)[50], 0.25);
  }
}

TEST(FreeEnd, WindowsMatchTheWholePathOnTheSameGrid) {
  for (auto scheme : kSchemes) {
    auto whole = line(0, 0.2, scheme, 101), left = line(0, 0.1, scheme),
         right = line(0.1, 0.2, scheme);
    const Bound stop{0, 0};
    ASSERT_EQ(right->computeControllableSetsFromSquaredSpeed(stop), ReturnCode::OK);
    const Bound boundary = right->getParameterizationData().controllable_sets.row(0);
    ASSERT_EQ(left->computePathParametrizationFromSquaredSpeed(0, boundary),
              ReturnCode::OK);
    const value_type reached = profile(*left)[50];
    ASSERT_EQ(right->computePathParametrizationFromSquaredSpeed(reached, stop),
              ReturnCode::OK);
    ASSERT_EQ(whole->computePathParametrizationFromSquaredSpeed(0, stop), ReturnCode::OK);
    EXPECT_LE(maxDifference(profile(*whole).head(51), profile(*left)), 1e-12);
    EXPECT_LE(maxDifference(profile(*whole).tail(51), profile(*right)), 1e-12);
  }
}

TEST(FreeEnd, CurvedWindowsMatchTheWholePathForAnIntervalOfFinalSpeeds) {
  for (auto scheme : kSchemes) {
    auto whole = curve(0, 2, scheme);
    const value_type start = 0.01;
    const Bound finalSpeeds{0.04, 0.09};
    ASSERT_EQ(whole->computePathParametrizationFromSquaredSpeed(start, finalSpeeds),
              ReturnCode::OK);

    const value_type edges[] = {0, 0.7, 1.3, 2};
    Bound interfaces[4];
    interfaces[3] = finalSpeeds;
    for (int w = 2; w >= 0; --w) {
      auto window = curve(edges[w], edges[w + 1], scheme);
      ASSERT_EQ(window->computeControllableSetsFromSquaredSpeed(interfaces[w + 1]),
                ReturnCode::OK);
      interfaces[w] = window->getParameterizationData().controllable_sets.row(0);
    }
    value_type reached = start;
    for (int w = 0; w < 3; ++w) {
      auto window = curve(edges[w], edges[w + 1], scheme);
      ASSERT_EQ(window->computePathParametrizationFromSquaredSpeed(reached,
                                                                   interfaces[w + 1]),
                ReturnCode::OK);
      const auto first = std::lround(edges[w] * 100);
      const Vector& expected = profile(*whole).segment(first, profile(*window).size());
      EXPECT_LE(maxDifference(expected, profile(*window)), 1e-8) << "window " << w;
      reached = profile(*window).tail<1>()[0];
    }
    EXPECT_GE(reached, finalSpeeds[0]);
    EXPECT_LE(reached, finalSpeeds[1]);
  }
}

TEST(FreeEnd, RestrictsTheIntervalToTheFinalGridpointsConstraints) {
  auto problem = line(0, 10, DiscretizationType::Interpolation);
  ASSERT_EQ(problem->computePathParametrizationFromSquaredSpeed(0, Bound{0, 100}),
            ReturnCode::OK);
  // The velocity limit 2 caps the final squared velocity at 4.
  EXPECT_NEAR(profile(*problem)[50], 4, 1e-12);
  EXPECT_EQ(problem->computeControllableSetsFromSquaredSpeed(Bound{5, 100}),
            ReturnCode::ERR_FAIL_CONTROLLABLE);
}

TEST(FreeEnd, KeepsTheLowerBoundAndRejectsUnreachableStartsAndInvalidInputs) {
  auto problem = line(0, 0.1, DiscretizationType::Interpolation);
  // Reaching x >= 1 from rest needs 0.5 of path; the start is not controllable.
  EXPECT_EQ(problem->computePathParametrizationFromSquaredSpeed(0, Bound{1, 2}),
            ReturnCode::ERR_FAIL_FORWARD_PASS);
  EXPECT_GT(problem->getParameterizationData().controllable_sets(0, 0), 0);
  const value_type nan = std::numeric_limits<value_type>::quiet_NaN();
  EXPECT_THROW(problem->computePathParametrizationFromSquaredSpeed(-1, Bound{0, 1}),
               std::invalid_argument);
  EXPECT_THROW(problem->computePathParametrizationFromSquaredSpeed(nan, Bound{0, 1}),
               std::invalid_argument);
  for (const Bound& invalid : {Bound{2, 1}, Bound{-1, 1}, Bound{0, nan}})
    EXPECT_THROW(problem->computeControllableSetsFromSquaredSpeed(invalid),
                 std::invalid_argument);
}

TEST(FreeEnd, SolvesASingleStage) {
  for (auto scheme : kSchemes) {
    auto problem = line(0, 0.1, scheme, 2);
    ASSERT_EQ(problem->computePathParametrizationFromSquaredSpeed(0, Bound{0, 4}),
              ReturnCode::OK);
    ASSERT_EQ(profile(*problem).size(), 2);
    EXPECT_NEAR(profile(*problem)[1], 0.2, 1e-12);
    EXPECT_EQ(problem->computePathParametrizationFromSquaredSpeed(0, Bound{1, 4}),
              ReturnCode::ERR_FAIL_FORWARD_PASS);
    EXPECT_EQ(problem->computeControllableSetsFromSquaredSpeed(Bound{5, 6}),
              ReturnCode::ERR_FAIL_CONTROLLABLE);
    ASSERT_EQ(problem->computePathParametrizationFromSquaredSpeed(0.1, Bound{0, 0}),
              ReturnCode::OK);
    EXPECT_EQ(profile(*problem)[1], 0);
  }
}
