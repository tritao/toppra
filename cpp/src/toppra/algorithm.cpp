#include <toppra/algorithm.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include "toppra/toppra.hpp"

namespace toppra {

PathParametrizationAlgorithm::PathParametrizationAlgorithm(
    LinearConstraintPtrs constraints, const GeometricPathPtr &path)
    : m_constraints(std::move(constraints)), m_path(path) {};

ReturnCode PathParametrizationAlgorithm::computePathParametrization(value_type vel_start,
                                                                    value_type vel_end) {
  m_controllableValid = false;
  initialize();
  m_solver->setupSolver();
  Bound vel_ends;
  vel_ends.setConstant(vel_end);
  m_data.ret_code = computeControllableSets(vel_ends);
  if ((int)m_data.ret_code > 0) {
    return m_data.ret_code;
  }
  m_data.ret_code = computeForwardPass(vel_start);
  return m_data.ret_code;
};

namespace {
void checkFinalSquaredSpeed(const Bound& end_squared) {
  if (!end_squared.allFinite() || end_squared[0] < 0 || end_squared[0] > end_squared[1])
    throw std::invalid_argument(
        "Final squared path velocities should be finite, nonnegative and ordered.");
}
}  // namespace

ReturnCode PathParametrizationAlgorithm::computeControllableSetsFromSquaredSpeed(
    const Bound& end_squared) {
  checkFinalSquaredSpeed(end_squared);
  initialize();
  m_solver->setupSolver();
  Bound end = end_squared;
  m_data.ret_code = constrainFinalSquaredSpeed(end);
  if (m_data.ret_code != ReturnCode::OK) return m_data.ret_code;
  // Only now are the controllable sets changed.
  m_controllableValid = false;
  m_data.ret_code = computeControllableSetsSquared(end);
  m_controllableValid = m_data.ret_code == ReturnCode::OK;
  return m_data.ret_code;
}

ReturnCode PathParametrizationAlgorithm::computeForwardFromControllableSets(
    value_type start_squared, SolverPtr workspace) {
  if (!std::isfinite(start_squared) || start_squared < 0)
    throw std::invalid_argument(
        "Initial squared path velocity should be finite and nonnegative.");
  if (!m_controllableValid)
    throw std::logic_error(
        "No controllable sets of computeControllableSetsFromSquaredSpeed for the "
        "current configuration.");
  adoptSolverWorkspace(std::move(workspace));
  initialize();
  m_solver->setupSolver();
  m_data.ret_code = computeForwardPass(start_squared);
  return m_data.ret_code;
}

ReturnCode PathParametrizationAlgorithm::updateControllableSetsFromSquaredSpeed(
    const Bound& end_squared, SolverPtr workspace) {
  checkFinalSquaredSpeed(end_squared);
  if (!m_controllableValid)
    throw std::logic_error(
        "No controllable sets of computeControllableSetsFromSquaredSpeed for the "
        "current configuration.");
  adoptSolverWorkspace(std::move(workspace));
  initialize();
  m_solver->setupSolver();
  Bound end = end_squared;
  m_data.ret_code = constrainFinalSquaredSpeed(end);
  if (m_data.ret_code != ReturnCode::OK) return m_data.ret_code;
  m_controllableValid = false;
  m_data.ret_code = computeControllableSetsSquared(end, true);
  m_controllableValid = m_data.ret_code == ReturnCode::OK;
  return m_data.ret_code;
}

void PathParametrizationAlgorithm::prepareSolverWorkspace(SolverPtr workspace) {
  adoptSolverWorkspace(std::move(workspace));
  initialize();
}

void PathParametrizationAlgorithm::releaseSolverWorkspace() {
  m_solver.reset();
  m_initialized = false;
}

void PathParametrizationAlgorithm::adoptSolverWorkspace(SolverPtr workspace) {
  if (!workspace || workspace == m_solver) return;
  if (m_solver)
    throw std::logic_error("Release the solver before giving the algorithm another one.");
  m_solver = std::move(workspace);
  m_initialized = false;
}

ReturnCode PathParametrizationAlgorithm::computePathParametrizationFromSquaredSpeed(
    value_type start_squared, const Bound& end_squared) {
  if (!std::isfinite(start_squared) || start_squared < 0)
    throw std::invalid_argument(
        "Initial squared path velocity should be finite and nonnegative.");
  if (computeControllableSetsFromSquaredSpeed(end_squared) != ReturnCode::OK)
    return m_data.ret_code;
  // The forward pass checks the start itself, with the problem of its first
  // stage. Compared with the controllable set, a start at zero would be
  // rejected where rounding left a lower bound of 1e-16, for instance.
  m_data.ret_code = computeForwardPass(start_squared);
  return m_data.ret_code;
}

ReturnCode PathParametrizationAlgorithm::constrainFinalSquaredSpeed(Bound& end_squared) {
  // computeControllableSets constrains each gridpoint only through the stage
  // that leaves it: the final gridpoint's own constraints apply here.
  const Bound requested = end_squared;
  Matrix H;
  Vector g(2), solution(2);
  const Bound unused{0, 0};
  for (int bound : {1, 0}) {
    // Maximize, then minimize, x within the interval.
    g << 0, bound == 1 ? -1 : 1;
    if (!m_solver->solveStagewiseOptim(m_N, H, g, end_squared, unused, solution)) {
      m_errorStream << "No final squared velocity in " << requested
                    << " is admissible at the final gridpoint.\n";
      return ReturnCode::ERR_FAIL_CONTROLLABLE;
    }
    end_squared[bound] = std::max(end_squared[0], std::min(end_squared[1], solution[1]));
  }
  return ReturnCode::OK;
}

ReturnCode PathParametrizationAlgorithm::computeControllableSets(
    const Bound &vel_ends) {
  return computeControllableSetsSquared(vel_ends.cwiseAbs2());
}

ReturnCode PathParametrizationAlgorithm::computeControllableSetsSquared(
    const Bound &end_squared, bool stopWhereUnchanged) {
  TOPPRA_LOG_DEBUG("computeControllableSets");
  bool solver_ret;
  Vector g_upper{2}, g_lower{2}, solution;
  g_upper << 1e-9, -1;
  g_lower << -1e-9, 1;
  if (stopWhereUnchanged && m_data.controllable_sets.row(m_N) == end_squared)
    return ReturnCode::OK;
  m_data.controllable_sets.row(m_N) = end_squared;

  Matrix H;
  Bound x (m_initXBound), x_next;
  for (std::size_t i = m_N - 1; i != (std::size_t)-1; i--) {
    TOPPRA_LOG_DEBUG(i << ", " << m_N);
    Bound previous{0, 0};
    if (stopWhereUnchanged) previous = m_data.controllable_sets.row(i);
    x_next = m_data.controllable_sets.row(i + 1);
    solver_ret = m_solver->solveStagewiseOptim(i, H, g_upper, x, x_next, solution);

    if (!solver_ret) {
      TOPPRA_LOG_DEBUG("Fail: controllable, upper problem, idx: " << i);
      m_errorStream <<
        "Failed to compute controllable set, upper problem, at idx " << i <<
        "\ng_upper = " << g_upper.transpose() <<
        "\nx = " << x <<
        "\nx_next = " << x_next <<
        '\n';
      return ReturnCode::ERR_FAIL_CONTROLLABLE;
    }

    m_data.controllable_sets(i, 1) = solution[1];
    m_data.controllable_upper_controls(i) = solution[0];

    solver_ret = m_solver->solveStagewiseOptim(i, H, g_lower, x, x_next, solution);

    TOPPRA_LOG_DEBUG("down: " << solution.transpose());
    if (!solver_ret) {
      TOPPRA_LOG_DEBUG("Fail: controllable, lower problem, idx: " << i);
      m_errorStream <<
        "Failed to compute controllable set, lower problem, at idx " << i <<
        "\ng_lower = " << g_lower.transpose() <<
        "\nx = " << x <<
        "\nx_next = " << x_next <<
        '\n';
      return ReturnCode::ERR_FAIL_CONTROLLABLE;
    }

    // For whatever reason, sometimes the solver return negative
    // solution despite having a set of positve bounds. This readjust
    // if the solution is negative.
    m_data.controllable_sets(i, 0) = std::max(0.0, solution[1]);
    // The problems of the earlier stages depend only on this set: if it is
    // unchanged, so are theirs. This stage's upper control, which depends on
    // the next set too, has been updated.
    if (stopWhereUnchanged && m_data.controllable_sets(i, 0) == previous[0]
        && m_data.controllable_sets(i, 1) == previous[1])
      break;
  }
  return ReturnCode::OK;
}

ReturnCode PathParametrizationAlgorithm::computeFeasibleSets() {
  initialize();
  bool solver_ret;
  Vector g_upper{2}, g_lower{2}, solution;
  g_upper << 1e-9, -1;
  g_lower << -1e-9, 1;

  Matrix H;
  Bound x, x_next;
  x << 0, 100;
  x_next << 0, 100;
  for (std::size_t i = 0; i < m_N + 1; i++) {
    solver_ret = m_solver->solveStagewiseOptim(i, H, g_upper, x, x_next, solution);

    if (!solver_ret) {
      TOPPRA_LOG_DEBUG("Fail: feasible, upper problem, idx: " << i);
      m_errorStream <<
        "Failed to compute feasible set, upper problem, at idx " << i <<
        "\ng_upper = " << g_upper.transpose() <<
        "\nx = " << x <<
        "\nx_next = " << x_next <<
        '\n';
      return ReturnCode::ERR_FAIL_FEASIBLE;
    }

    m_data.feasible_sets(i, 1) = solution[1];

    solver_ret = m_solver->solveStagewiseOptim(i, H, g_lower, x, x_next, solution);

    if (!solver_ret) {
      TOPPRA_LOG_DEBUG("Fail: feasible, lower problem, idx: " << i);
      m_errorStream <<
        "Failed to compute feasible set, lower problem, at idx " << i <<
        "\ng_lower = " << g_lower.transpose() <<
        "\nx = " << x <<
        "\nx_next = " << x_next <<
        '\n';
      return ReturnCode::ERR_FAIL_FEASIBLE;
    }

    m_data.feasible_sets(i, 0) = solution[1];
  }
  return ReturnCode::OK;
}

void PathParametrizationAlgorithm::setGridpoints(const Vector& gridpoints)
{
  Bound I (m_path->pathInterval());
  auto s = gridpoints.size();
  if (s < 2 || !gridpoints.allFinite())
    throw std::invalid_argument("gridpoints should contain at least two finite values.");
  if (   std::fabs(gridpoints[0] - I(0)) > TOPPRA_NEARLY_ZERO
      || std::fabs(gridpoints[gridpoints.size()-1] - I(1)) > TOPPRA_NEARLY_ZERO)
    throw std::invalid_argument("Path interval does not match first and last gridpoints.");
  if (((gridpoints.tail(s-1) - gridpoints.head(s-1)).array() <= 0).any())
    throw std::invalid_argument("gridpoints should be strictly increasing.");
  // Change the problem only once the grid is known to be valid.
  m_controllableValid = false;
  m_data.gridpoints = gridpoints;
  m_N = s-1;
  m_initialized = false;
}

void PathParametrizationAlgorithm::initialize() {
  m_errorStream = std::stringstream();
  // A solver shared with another algorithm must be initialized again once the
  // other algorithm has initialized it with its own problem.
  if (m_initialized && m_solver->generation() == m_solverGeneration) return;
  if (!m_solver)
    throw std::logic_error("You must set a solver first.");
  Bound I (m_path->pathInterval());

  // Intiialize the gridpoints based on N
  if (m_data.gridpoints.size() == 0) {

    // If the number of gridpoint is set to 0, try to find automatically a reasonable gridpoint
    if (m_N == 0){
      m_data.gridpoints = m_path->proposeGridpoints();
      m_N = m_data.gridpoints.size() - 1;
    }

    // Otherwise use linearly spaced gridpoints
    else
      m_data.gridpoints = Vector::LinSpaced(m_N + 1, I(0), I(1));

    TOPPRA_LOG_DEBUG("gridpoints used:" << m_data.gridpoints.transpose());
  }
  else if (m_data.gridpoints.size() != m_N + 1)
    throw std::invalid_argument("number of gridpoints does not match attribute N.");

  auto s = m_data.gridpoints.size();
  m_data.parametrization.resize(s);
  m_data.controllable_sets.resize(s, 2);
  m_data.controllable_upper_controls.resize(s - 1);
  m_data.feasible_sets.resize(s, 2);
  m_solver->initialize(m_constraints, m_path, m_data.gridpoints);
  m_solverGeneration = m_solver->generation();
  m_initialized = true;
}

}  // namespace toppra
