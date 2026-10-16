#include <toppra/solver.hpp>

#include <cmath>

#include <toppra/constraint.hpp>
#include <toppra/geometric_path.hpp>

#ifdef BUILD_WITH_qpOASES
#include <toppra/solver/qpOASES-wrapper.hpp>
#endif
#ifdef BUILD_WITH_GLPK
#include <toppra/solver/glpk-wrapper.hpp>
#endif
#include <toppra/solver/seidel.hpp>

namespace toppra {

SolverPtr Solver::createDefault() {
#ifdef BUILD_WITH_qpOASES
    return std::make_shared<solver::qpOASESWrapper>();
#elif BUILD_WITH_GLPK
    return std::make_shared<solver::GLPKWrapper>();
#else
    return std::make_shared<solver::Seidel>();
#endif
}

void Solver::initialize (const LinearConstraintPtrs& constraints, const GeometricPathPtr& path,
    const Vector& times)
{
  // Validate the inputs before replacing the current problem.
  if (!path)
    throw std::invalid_argument("Invalid path.");
  if (times.size() == 0 || !times.allFinite())
    throw std::invalid_argument("Invalid times.");
  for (Eigen::Index i = 1; i < times.size(); ++i)
    if (!(times[i] > times[i-1]) || !std::isfinite(times[i] - times[i-1]))
      throw std::invalid_argument("Invalid times.");
  for (const auto& constraint : constraints)
    if (!constraint)
      throw std::invalid_argument("Invalid constraint.");

  ++m_generation;
  m_constraints = constraints;
  m_path = path;
  m_times = times;
  m_N = times.size()-1;
  m_nV = 2;
  m_deltas = times.tail(m_N) - times.head(m_N);
  /// \todo assert that the time range of the path equals [ times[0], times[m_N] ].

  // Compute the constraints parameters, replacing those of the previous problem.
  // The records are resized rather than cleared, so that a solver initialized
  // again keeps the storage of their vectors and matrices.
  std::size_t nbLin = 0, nbBox = 0;
  for (const auto& constraint : m_constraints) {
    if (constraint->hasLinearInequalities()) ++nbLin;
    if (constraint->hasUbounds() || constraint->hasXbounds()) ++nbBox;
  }
  m_constraintsParams.lin.resize(nbLin);
  m_constraintsParams.box.resize(nbBox);
  std::size_t iLin = 0, iBox = 0;
  LinearConstraintParams emptyLinParams;
  BoxConstraintParams emptyBoxParams;
  for (std::size_t c = 0; c < m_constraints.size(); ++c) {
    LinearConstraint* lin = m_constraints[c].get();
    LinearConstraintParams* lparam;
    if (lin->hasLinearInequalities()) {
      lparam = &m_constraintsParams.lin[iLin++];
      lparam->cid = c;
    } else
      lparam = &emptyLinParams;
    BoxConstraintParams* bparam;
    if (lin->hasUbounds() || lin->hasXbounds()) {
      bparam = &m_constraintsParams.box[iBox++];
      bparam->cid = c;
      // The record may hold bounds of the other kind, from another constraint.
      if (!lin->hasUbounds()) bparam->u.clear();
      if (!lin->hasXbounds()) bparam->x.clear();
    } else
      bparam = &emptyBoxParams;
    lin->computeParams(*path, times,
        lparam->a, lparam->b, lparam->c, lparam->F, lparam->g,
        bparam->u, bparam->x);
  }
}

} // namespace toppra
