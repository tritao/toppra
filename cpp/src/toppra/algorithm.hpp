#ifndef TOPPRA_ALGORITHM_HPP
#define TOPPRA_ALGORITHM_HPP

#include <stdexcept>
#include <sstream>
#include <toppra/constraint.hpp>
#include <toppra/geometric_path.hpp>
#include <toppra/solver.hpp>
#include <toppra/toppra.hpp>

namespace toppra {

/// Return code for Path Parametrization algorithm.
enum class ReturnCode {

  /// Success
  OK = 0,

  /// Unknown error
  ERR_UNKNOWN = 1,

  /// Fail during computing controllable sets. Problem might be infeasible.
  ERR_FAIL_CONTROLLABLE = 2,

  /// Fail during forward pass. Numerical error occured.
  ERR_FAIL_FORWARD_PASS = 3,

  /// Problem is not initialized
  ERR_UNINITIALIZED = 4,

  /// Fail to ocmpute feasible sets.
  ERR_FAIL_FEASIBLE = 5,
};

struct ParametrizationData {
  /// \brief Grid-points used for solving the discretized problem.
  /// The number of points must equal m_N + 1.
  Vector gridpoints;

  ///  Output parametrization (squared path velocity)
  Vector parametrization;

  Matrix controllable_sets;
  Matrix feasible_sets;

  /// Path acceleration with which the upper bound of each controllable set
  /// reaches the next controllable set, found while computing the sets.
  Vector controllable_upper_controls;

  ///  Return code of the algorithm.
  ReturnCode ret_code = ReturnCode::ERR_UNINITIALIZED;
};

/** \brief Base class for time parametrization algorithms.
 *
 */
class PathParametrizationAlgorithm {
 public:
  /** Construct the problem instance.
   *
   *  \param  constraints  List of constraints.
   *  \param  path  The geometric path.
   *
   */
  PathParametrizationAlgorithm(LinearConstraintPtrs constraints,
                               const GeometricPathPtr &path);

  /** \brief Set the level of discretization used by the solver.
   *
   * If is zero, will attempt to detect automatically the most
   * suitable grid.  
   * 
   * TODO: Automatic gridpoint selection is not implemented in the CPP
   * API.
   */
  void setN(int N) { m_N = N; m_initialized = false; m_controllableValid = false; };

  /** \brief Set the gridpoints (points with the path intervals)
   *
   * If not set manually, then N equally distributed points is used.
   */
  void setGridpoints(const Vector& gridpoints);

  /** \brief Set the LP/QP solver
   *
   * Default to \ref solver::qpOASESWrapper
   *
   * A solver can be shared by several algorithms, one computation at a time:
   * an algorithm initializes it again when another one has initialized it since.
   */
  void solver(SolverPtr solver) {
    m_solver.swap(solver);
    m_initialized = false;
    m_controllableValid = false;
  };

  /** \brief Get output or result of algorithm.
   */
  const ParametrizationData& getParameterizationData() const { return m_data; };

  /** Compute the time parametrization of the given path.
   *
   * \param vel_start Initial squared path velocity.
   * \param vel_end Final path velocity, not squared. See
   *        computePathParametrizationFromSquaredSpeed for an interval of final
   *        speeds and uniformly squared units.
   * \return Return code. When not ReturnCode::OK,
   *         check PathParametrizationAlgorithm::getErrorMessage
   */
  virtual ReturnCode computePathParametrization(value_type vel_start = 0,
                                                value_type vel_end = 0);

  /** Compute the controllable sets for an interval of final squared path
   * velocities, end_squared = [lower, upper], which is first restricted to the
   * squared velocities that the constraints admit at the final gridpoint. Only
   * the backward pass is computed.
   *
   * \throw std::invalid_argument if end_squared is not finite, nonnegative and
   *        ordered.
   * \return ReturnCode::ERR_FAIL_CONTROLLABLE if no velocity of the interval
   *         is admissible or the controllable sets cannot be computed.
   */
  ReturnCode computeControllableSetsFromSquaredSpeed(const Bound& end_squared);

  /** Compute the fastest parametrization from the squared path velocity
   * start_squared to any final squared path velocity in end_squared, which
   * is the last entry of the parametrization.
   *
   * A path can be parametrized in consecutive windows that share their end
   * gridpoints. Computed from the last window backward, the controllable set
   * at the first gridpoint of each window is the final interval of the
   * previous one; computed forward, each window starts at the velocity the
   * previous one reached. The parametrization equals that of the whole path
   * on the same grid.
   *
   * \throw std::invalid_argument if start_squared is not finite and
   *        nonnegative, or end_squared not as in
   *        computeControllableSetsFromSquaredSpeed.
   * \return ReturnCode::ERR_FAIL_FORWARD_PASS if start_squared is not
   *         controllable.
   */
  ReturnCode computePathParametrizationFromSquaredSpeed(value_type start_squared,
                                                       const Bound& end_squared);

  /** Compute only the forward pass, from the controllable sets of the last
   * successful computeControllableSetsFromSquaredSpeed or
   * computePathParametrizationFromSquaredSpeed of this algorithm.
   *
   * The controllable sets stay valid until the grid, N, the initial x bounds
   * or the solver are set. Changes of the path or constraint objects
   * themselves are not detected: they must not change while the controllable
   * sets are reused. With the same type of solver, the result is that of
   * computePathParametrizationFromSquaredSpeed.
   *
   * \param workspace A solver for an algorithm that has released its own
   *        (releaseSolverWorkspace). It is initialized with this algorithm's
   *        problem.
   * \throw std::logic_error if there are no such controllable sets, or if a
   *        workspace is given to an algorithm that has another solver.
   */
  ReturnCode computeForwardFromControllableSets(value_type start_squared,
                                                SolverPtr workspace = SolverPtr());

  /** Update the controllable sets for another interval of final squared
   * velocities, as computeControllableSetsFromSquaredSpeed would compute
   * them, from those of its last successful computation (or update).
   *
   * The backward pass stops at the first gridpoint whose controllable set is
   * unchanged, since the earlier ones then are too. The conditions of
   * computeForwardFromControllableSets apply, and so does its workspace.
   * When no velocity of the interval is admissible, ERR_FAIL_CONTROLLABLE is
   * returned and the controllable sets are kept; when the backward pass
   * fails, they are discarded.
   *
   * \throw std::invalid_argument as computeControllableSetsFromSquaredSpeed.
   * \throw std::logic_error as computeForwardFromControllableSets.
   */
  ReturnCode updateControllableSetsFromSquaredSpeed(const Bound& end_squared,
                                                   SolverPtr workspace = SolverPtr());

  /** Release the solver, keeping the controllable sets for
   * computeForwardFromControllableSets. Many algorithms can thus keep their
   * controllable sets while one solver, given to each in turn, holds the
   * tables of a single problem at a time.
   */
  void releaseSolverWorkspace();

  /** Compute the sets of feasible squared velocities.
   */
  ReturnCode computeFeasibleSets();

  /** Set initial bounds on \f$ \dot{s}^2.
   * This is helpfull when the solver encounters numerical issues.
   */
  void setInitialXBounds (const Bound& xbound)
  {
    m_initXBound = xbound;
    m_controllableValid = false;
  }

  /** Get the error message when PathParametrizationAlgorithm::computePathParametrization
   * failed.
   */
  std::string getErrorMessage() const {
    return m_errorStream.str();
  }

  virtual ~PathParametrizationAlgorithm() {}

 protected:
  /** \brief Select solver and gridpoints to use.
   *
   * This method implements a simple way to select gridpoints.
   */
  virtual void initialize();

  /** \brief Compute the forward pass.
   *
   * Derived class should provide a suitable forward pass function,
   * depending on the desired objective.
   */
  virtual ReturnCode computeForwardPass(value_type vel_start) = 0;

  /** Compute the sets of controllable squared path velocities.
   */
  ReturnCode computeControllableSets(const Bound &vel_ends);

  /** Compute the sets of controllable squared path velocities, from the
   * interval end_squared of final squared path velocities.
   */
  ReturnCode computeControllableSetsSquared(const Bound &end_squared,
                                           bool stopWhereUnchanged = false);

  /** Restrict the interval of final squared velocities to those that the
   * constraints admit at the final gridpoint.
   */
  ReturnCode constrainFinalSquaredSpeed(Bound &end_squared);

  /** Use workspace as the solver of an algorithm that released its own. */
  void adoptSolverWorkspace(SolverPtr workspace);

  /** To be implemented in child method. */
  LinearConstraintPtrs m_constraints;
  GeometricPathPtr m_path;
  SolverPtr m_solver;

  std::stringstream m_errorStream;

  /// Struct containing algorithm output.
  ParametrizationData m_data;

  /// \brief Number of segments in the discretized problems.
  /// See m_gridpoints for more information.
  int m_N = 100;

  int m_initialized = false;
  /// Solver::generation after this algorithm initialized its solver.
  std::size_t m_solverGeneration = 0;
  /// Whether the controllable sets are those of a successful
  /// computeControllableSetsFromSquaredSpeed for the current configuration.
  bool m_controllableValid = false;

  /** Set initial bounds on \f$ \dot{s}^2.
   * \sa setInitialXBounds
   * \todo The hard-coded bound below avoids numerical issues in LP / QP solvers
   * when \f$ x \f$ becomes too big. This issue should be addressed in the
   * solver wrapper themselfves as numerical behaviors is proper to each
   * individual solver.
   * See https://github.com/hungpham2511/toppra/issues/156
   */
  Bound m_initXBound = {0, 100};
};

}  // namespace toppra

#endif
