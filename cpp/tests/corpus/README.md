# Solver replay corpus

Every `.case` file in this directory is replayed by `RegressionCorpus.ReplayAll` in
`all_tests`, and so by the C++ CI jobs. CMake picks up added files when it builds;
no C++ change is needed. The `DegenerateStages` tests replay the cases they were
written for by name.

The cases are feasible paths on which the Seidel solver failed: nearly parallel
constraints in the forward and backward passes, a forward pass starting at the
upper bound of a controllable set, and a forward step ending one rounding below
that bound. `FixedCoordinateCancellation` was found by the generator below.

The checker evaluates, at every gridpoint, the joint velocity and the
(collocation) joint acceleration limits, finite nonnegative squared velocities,
the initial and final velocities and membership of the controllable sets. It
checks the parametrization of `computePathParametrization` (whose final velocity
is not squared) and of `computePathParametrizationFromSquaredSpeed`, and that the
two agree. These checks concern the gridpoints only, not the continuous
acceleration between them.

`RegressionCorpus.SeededEquivalentInputs` permutes the joints of each case and
changes their signs and units by powers of two, scaling the limits with them.
These are the same problem, with the same path-velocity profile. A unit change
that would overflow, underflow to zero or lose bits falls back to a sign change.
The generator is an explicit 32-bit linear congruential generator, reproducible
on every standard library. Each case has its own stream, derived from its file
name and the seed, so that adding, removing or moving other cases does not change
it. By default, CI runs 64 variants of every case, with seed `0x5eed0302`.

From the repository root:

```sh
cmake -S cpp -B cpp/build -DPYTHON_BINDINGS=OFF -DBUILD_TESTING=ON
cmake --build cpp/build --target all_tests
cpp/build/tests/all_tests --gtest_filter='DegenerateStages.*:RegressionCorpus.*'
TOPPRA_FUZZ_SEED=0x5eed0302 TOPPRA_FUZZ_CASES=4096 cpp/build/tests/all_tests --gtest_filter=RegressionCorpus.SeededEquivalentInputs
```

A failing generated variant is saved under the build's `tests/corpus-failures/`
directory, named after the seed, the case and the iteration, which the test
trace also prints. The C++ CI job uploads these files when it fails. To replay
one without adding it to the corpus:

```sh
TOPPRA_REPLAY_CASE=/path/to/failure.case cpp/build/tests/all_tests --gtest_filter=RegressionCorpus.ReplayAll
```

Once a failure is understood, copy its `.case` file here under a meaningful name
and describe it below; CI then replays it. Neither the tests nor CI change the
corpus by themselves.

## Cases

- `ForwardPassAlongNearlyParallelConstraints`: the forward pass failed at stage
  22, on a nearly parallel row whose projected offset was only rounding.
- `ControllableSetsAlongNearlyParallelConstraints`: the backward pass failed to
  compute the controllable set of stage 18, for the same geometry.
- `ForwardPassFromTheUpperBoundOfAControllableSet`: the forward pass started on
  the upper bound of a controllable set, at a joint whose path-acceleration
  coefficient was rounding (-1.8e-16).
- `ForwardPassStepEndingOneRoundingBelowTheBound`: the step into stage 32
  ended one rounding below an equally degenerate upper bound.
- `FixedCoordinateCancellation`: a permuted and rescaled variant of
  `ForwardPassFromTheUpperBoundOfAControllableSet`, in which a residual of
  3.55e-15 left by eliminating the fixed squared velocity became a false bound
  on the control. The file keeps the variant, independently of the generator.

## File format

Whitespace-separated numbers in the classic locale:

1. `TOPPRA_REPLAY_V1`
2. The number of spans, the number of joints and the number of stages per span.
3. The initial and final **squared** path velocities. The checker takes the
   square root of the final one for `computePathParametrization`.
4. `span_count + 1` increasing knots.
5. One positive velocity limit per joint.
6. One positive acceleration limit per joint.
7. For each span, the number of coefficient rows, then the rows, one
   coefficient per joint, highest degree first and local to the span, as
   `PiecewisePolyPath` expects.

The reader rejects non-finite values, invalid sizes or limits, missing data and
trailing data. It accepts up to 128 spans, 64 joints, 4096 stages in all and
degree five. Failure files are written with `max_digits10` digits; a test checks
that writing and reading back is exact, signed zeros included.
