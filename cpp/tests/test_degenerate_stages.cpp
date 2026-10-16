// Feasible paths on which TOPPRA with the Seidel solver failed. Their data are
// in corpus/, which describes them; these tests replay them by name.
#include "regression_corpus.hpp"

TEST(DegenerateStages, ForwardPassAlongNearlyParallelConstraints) {
  regression::checkNamed("ForwardPassAlongNearlyParallelConstraints");
}

TEST(DegenerateStages, ControllableSetsAlongNearlyParallelConstraints) {
  regression::checkNamed("ControllableSetsAlongNearlyParallelConstraints");
}

TEST(DegenerateStages, ForwardPassFromTheUpperBoundOfAControllableSet) {
  regression::checkNamed("ForwardPassFromTheUpperBoundOfAControllableSet");
}

TEST(DegenerateStages, ForwardPassStepEndingOneRoundingBelowTheBound) {
  regression::checkNamed("ForwardPassStepEndingOneRoundingBelowTheBound");
}
