#pragma once

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <memory>
#include <stdexcept>
#include <string>
#include <toppra/algorithm/toppra.hpp>
#include <toppra/constraint/linear_joint_acceleration.hpp>
#include <toppra/constraint/linear_joint_velocity.hpp>
#include <toppra/geometric_path/piecewise_poly_path.hpp>
#include <toppra/solver/seidel.hpp>

#include "corpus_manifest.hpp"
#include "gtest/gtest.h"

// Reads, writes and checks the replay files of corpus/ (see corpus/README.md).
namespace regression {
using namespace toppra;

struct RecordedCase {
  Matrices coefficients;
  std::vector<value_type> knots;
  Vector velocity, acceleration;
  value_type startSquared, endSquared;
  int stages;
};

inline value_type readFinite(std::istream& in) {
  value_type value;
  if (!(in >> value) || !std::isfinite(value))
    throw std::runtime_error("Invalid corpus number");
  return value;
}

inline RecordedCase readCase(std::istream& in) {
  in.imbue(std::locale::classic());
  std::string version;
  int spans = 0, joints = 0, stages = 0;
  if (!(in >> version >> spans >> joints >> stages) || version != "TOPPRA_REPLAY_V1" ||
      spans < 1 || spans > 128 || joints < 1 || joints > 64 || stages < 1 ||
      stages > 4096 || spans * stages > 4096)
    throw std::runtime_error("Invalid corpus header");
  RecordedCase c;
  c.stages = stages;
  c.startSquared = readFinite(in);
  c.endSquared = readFinite(in);
  if (c.startSquared < 0 || c.endSquared < 0)
    throw std::runtime_error("Negative corpus velocity");
  for (int i = 0; i <= spans; ++i) {
    const value_type knot = readFinite(in);
    if (i > 0 && knot <= c.knots.back())
      throw std::runtime_error("Corpus knots must increase");
    c.knots.push_back(knot);
  }
  c.velocity.resize(joints);
  c.acceleration.resize(joints);
  for (Vector* limits : {&c.velocity, &c.acceleration}) {
    for (int j = 0; j < joints; ++j) {
      (*limits)[j] = readFinite(in);
      if ((*limits)[j] <= 0) throw std::runtime_error("Corpus limits must be positive");
    }
  }
  for (int i = 0; i < spans; ++i) {
    int rows = 0;
    if (!(in >> rows) || rows < 1 || rows > 6)
      throw std::runtime_error("Invalid corpus polynomial degree");
    Matrix coefficients(rows, joints);
    for (int row = 0; row < rows; ++row)
      for (int j = 0; j < joints; ++j) coefficients(row, j) = readFinite(in);
    c.coefficients.push_back(coefficients);
  }
  std::string trailing;
  if (in >> trailing) throw std::runtime_error("Unexpected trailing corpus data");
  return c;
}

inline RecordedCase readFile(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("Cannot read corpus file: " + path);
  return readCase(in);
}

inline void writeRow(std::ostream& out, const Vector& row) {
  for (Eigen::Index i = 0; i < row.size(); ++i) out << (i ? " " : "") << row[i];
  out << '\n';
}

inline void writeCase(std::ostream& out, const RecordedCase& c) {
  out.imbue(std::locale::classic());
  out << std::setprecision(std::numeric_limits<value_type>::max_digits10);
  out << "TOPPRA_REPLAY_V1\n"
      << c.coefficients.size() << ' ' << c.velocity.size() << ' ' << c.stages << '\n'
      << c.startSquared << ' ' << c.endSquared << '\n';
  writeRow(out, Eigen::Map<const Vector>(c.knots.data(), c.knots.size()));
  writeRow(out, c.velocity);
  writeRow(out, c.acceleration);
  for (const auto& coefficients : c.coefficients) {
    out << coefficients.rows() << '\n';
    for (Eigen::Index row = 0; row < coefficients.rows(); ++row)
      writeRow(out, coefficients.row(row).transpose());
  }
}

/// A recorded case, set up as a TOPPRA problem on its grid.
struct Problem {
  explicit Problem(const RecordedCase& c) : recorded(c) {
    path = std::make_shared<PiecewisePolyPath>(c.coefficients, c.knots);
    constraints = {
        std::make_shared<constraint::LinearJointVelocity>(-c.velocity, c.velocity),
        std::make_shared<constraint::LinearJointAcceleration>(-c.acceleration,
                                                              c.acceleration)};
    for (auto& constraint : constraints)
      constraint->discretizationType(DiscretizationType::Collocation);
    const int spans = static_cast<int>(c.coefficients.size());
    gridpoints.resize(spans * c.stages + 1);
    for (int span = 0; span < spans; ++span)
      for (int k = 0; k < c.stages; ++k)
        gridpoints[span * c.stages + k] =
            c.knots[span] + (c.knots[span + 1] - c.knots[span]) * k / c.stages;
    gridpoints[spans * c.stages] = c.knots.back();
  }

  std::unique_ptr<algorithm::TOPPRA> algorithm() const {
    auto problem = std::make_unique<algorithm::TOPPRA>(constraints, path);
    problem->solver(std::make_shared<solver::Seidel>());
    problem->setGridpoints(gridpoints);
    problem->setInitialXBounds(Bound{0, 1e16});
    return problem;
  }

  /// Checks a parametrization against the limits of the recorded case, its
  /// velocities and its controllable sets, at every gridpoint.
  void checkProfile(const ParametrizationData& data) const {
    const auto N = gridpoints.size() - 1;
    ASSERT_EQ(data.parametrization.size(), N + 1);
    EXPECT_NEAR(data.parametrization[0], recorded.startSquared, 1e-12);
    EXPECT_NEAR(data.parametrization[N], recorded.endSquared, 1e-12);
    for (Eigen::Index i = 0; i <= N; ++i) {
      const value_type x = data.parametrization[i];
      ASSERT_TRUE(std::isfinite(x)) << "stage " << i;
      ASSERT_GE(x, 0) << "stage " << i;
      EXPECT_GE(x, data.controllable_sets(i, 0)) << "stage " << i;
      EXPECT_LE(x, data.controllable_sets(i, 1)) << "stage " << i;
      const Vector qs = path->eval_single(gridpoints[i], 1);
      for (Eigen::Index j = 0; j < qs.size(); ++j)
        EXPECT_LE(std::abs(qs[j]) * std::sqrt(x), recorded.velocity[j] * (1 + 1e-9))
            << "velocity at stage " << i << ", joint " << j;
      if (i == N) continue;
      const value_type u = (data.parametrization[i + 1] - x) /
                           (2 * (gridpoints[i + 1] - gridpoints[i]));
      const Vector qss = path->eval_single(gridpoints[i], 2);
      for (Eigen::Index j = 0; j < qs.size(); ++j)
        EXPECT_LE(std::abs(qs[j] * u + qss[j] * x), recorded.acceleration[j] * (1 + 1e-9))
            << "acceleration at stage " << i << ", joint " << j;
    }
  }

  const RecordedCase recorded;
  std::shared_ptr<PiecewisePolyPath> path;
  LinearConstraintPtrs constraints;
  Vector gridpoints;
};

inline value_type relativeDifference(const Vector& actual, const Vector& expected) {
  return (actual - expected).cwiseAbs().maxCoeff() /
         std::max(value_type{1}, expected.cwiseAbs().maxCoeff());
}

/// Parametrizes a recorded case with computePathParametrization and with
/// computePathParametrizationFromSquaredSpeed, and checks both. The squared
/// profile is returned in `profile`.
inline void checkCase(const RecordedCase& c, Vector* profile = nullptr) {
  const Problem problem(c);
  auto legacy = problem.algorithm();
  ASSERT_EQ(legacy->computePathParametrization(c.startSquared, std::sqrt(c.endSquared)),
            ReturnCode::OK)
      << legacy->getErrorMessage();
  {
    SCOPED_TRACE("computePathParametrization");
    problem.checkProfile(legacy->getParameterizationData());
  }

  auto squared = problem.algorithm();
  ASSERT_EQ(squared->computePathParametrizationFromSquaredSpeed(
                c.startSquared, Bound{c.endSquared, c.endSquared}),
            ReturnCode::OK)
      << squared->getErrorMessage();
  {
    SCOPED_TRACE("computePathParametrizationFromSquaredSpeed");
    problem.checkProfile(squared->getParameterizationData());
  }
  // The final velocities differ by the rounding of sqrt(endSquared)^2.
  EXPECT_LE(relativeDifference(squared->getParameterizationData().parametrization,
                               legacy->getParameterizationData().parametrization),
            1e-10);
  const Vector expected = squared->getParameterizationData().parametrization;

  // The controllable sets kept after releasing the solver give exactly the
  // same profile, on a solver that another problem used meanwhile.
  auto workspace = std::make_shared<solver::Seidel>();
  squared->releaseSolverWorkspace();
  {
    auto coarser = problem.algorithm();
    const Eigen::Index N = problem.gridpoints.size() - 1;
    Vector gridpoints((N + 1) / 2 + 1);
    for (Eigen::Index i = 0; i + 1 < gridpoints.size(); ++i)
      gridpoints[i] = problem.gridpoints[2 * i];
    gridpoints[gridpoints.size() - 1] = problem.gridpoints[N];
    coarser->setGridpoints(gridpoints);
    coarser->solver(workspace);
    ASSERT_EQ(coarser->computeControllableSetsFromSquaredSpeed(Bound{0, 0}),
              ReturnCode::OK)
        << coarser->getErrorMessage();
  }
  ASSERT_EQ(squared->computeForwardFromControllableSets(c.startSquared, workspace),
            ReturnCode::OK)
      << squared->getErrorMessage();
  EXPECT_EQ(squared->getParameterizationData().parametrization, expected);
  if (profile) *profile = expected;
}

inline void checkNamed(const char* name) {
  checkCase(readFile(std::string(kCorpusDirectory) + '/' + name + ".case"));
}
}  // namespace regression
