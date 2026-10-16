#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>

#include "regression_corpus.hpp"

namespace {
using namespace regression;

constexpr std::size_t kCorpusSize = sizeof(kCorpusFiles) / sizeof(*kCorpusFiles);

uint32_t option(const char* name, uint32_t fallback) {
  const char* text = std::getenv(name);
  if (!text) return fallback;
  char* end = nullptr;
  errno = 0;
  const auto value = std::strtoull(text, &end, 0);
  if (errno || end == text || *end || value > UINT32_MAX)
    throw std::runtime_error(std::string("Invalid corpus option: ") + name);
  return static_cast<uint32_t>(value);
}

std::string caseName(const std::string& path) {
  const auto slash = path.find_last_of("/\\");
  return slash == std::string::npos ? path : path.substr(slash + 1);
}

/// The stream of a case depends on its file name only, not on its directory or
/// on the other cases.
uint32_t caseSeed(uint32_t seed, const std::string& path) {
  uint32_t hash = 2166136261u;  // FNV-1a
  for (unsigned char ch : caseName(path)) {
    hash ^= ch;
    hash *= 16777619u;
  }
  return seed ^ hash;
}

/// An explicit generator gives the same variants on every standard library.
uint32_t next(uint32_t& state) {
  state = state * 1664525u + 1013904223u;
  return state;
}

/// The same problem with its joints permuted, and their signs and units
/// changed by powers of two where that is exact.
RecordedCase equivalent(const RecordedCase& original, uint32_t& state) {
  RecordedCase c = original;
  const int joints = static_cast<int>(c.velocity.size());
  std::vector<int> order(joints);
  for (int j = 0; j < joints; ++j) order[j] = j;
  for (int j = joints - 1; j > 0; --j) std::swap(order[j], order[next(state) % (j + 1)]);
  for (int j = 0; j < joints; ++j) {
    const int from = order[j];
    const auto bits = next(state);
    value_type factor =
        std::ldexp((bits & 1) ? -1.0 : 1.0, static_cast<int>((bits >> 1) % 5) - 2);
    const auto exact = [&factor](value_type value) {
      const value_type scaled = value * factor;
      return std::isfinite(scaled) &&
             (value == 0 || (scaled != 0 && scaled / factor == value));
    };
    bool representable = exact(original.velocity[from]) && exact(original.acceleration[from]);
    for (const auto& coefficients : original.coefficients)
      for (Eigen::Index row = 0; row < coefficients.rows(); ++row)
        representable = representable && exact(coefficients(row, from));
    if (!representable) factor = factor < 0 ? -1 : 1;
    c.velocity[j] = original.velocity[from] * std::abs(factor);
    c.acceleration[j] = original.acceleration[from] * std::abs(factor);
    for (std::size_t span = 0; span < c.coefficients.size(); ++span)
      c.coefficients[span].col(j) = original.coefficients[span].col(from) * factor;
  }
  return c;
}

int failures() {
  const auto* result = ::testing::UnitTest::GetInstance()->current_test_info()->result();
  int count = 0;
  for (int i = 0; i < result->total_part_count(); ++i)
    if (result->GetTestPartResult(i).failed()) ++count;
  return count;
}

void saveFailure(const RecordedCase& c, uint32_t seed, const std::string& file,
                 uint32_t iteration) {
  const auto path = std::string(kFailureDirectory) + "/seed-" + std::to_string(seed) +
                    "-case-" + caseName(file) + "-iteration-" + std::to_string(iteration) +
                    ".case";
  std::ofstream out(path);
  if (out) {
    writeCase(out, c);
    out.close();
  }
  std::cerr << "Replay corpus failure: " << path << '\n';
  if (!out) {
    std::cerr << "Could not save the replay file; its input follows:\n";
    writeCase(std::cerr, c);
  }
}
}  // namespace

TEST(RegressionCorpus, ReplayAll) {
  ASSERT_GE(kCorpusSize, 4u);
  for (const char* path : kCorpusFiles) {
    SCOPED_TRACE(path);
    checkCase(readFile(path));
  }
  if (const char* path = std::getenv("TOPPRA_REPLAY_CASE")) {
    SCOPED_TRACE(path);
    checkCase(readFile(path));
  }
}

TEST(RegressionCorpus, SerializationPreservesExactInputs) {
  for (const char* path : kCorpusFiles) {
    SCOPED_TRACE(path);
    const RecordedCase c = readFile(path);
    std::stringstream text;
    writeCase(text, c);
    const RecordedCase restored = readCase(text);
    std::stringstream canonical;
    writeCase(canonical, restored);
    // Includes the sign of zeros and every digit needed for an exact round trip.
    EXPECT_EQ(canonical.str(), text.str());
    EXPECT_EQ(restored.knots, c.knots);
    EXPECT_EQ(restored.stages, c.stages);
    EXPECT_EQ(restored.startSquared, c.startSquared);
    EXPECT_EQ(restored.endSquared, c.endSquared);
    EXPECT_EQ(restored.velocity, c.velocity);
    EXPECT_EQ(restored.acceleration, c.acceleration);
    ASSERT_EQ(restored.coefficients.size(), c.coefficients.size());
    for (std::size_t i = 0; i < c.coefficients.size(); ++i)
      EXPECT_EQ(restored.coefficients[i], c.coefficients[i]);
  }
}

TEST(RegressionCorpus, RejectsMalformedReplayFiles) {
  for (const char* text : {
           "",
           "TOPPRA_REPLAY_V2 1 1 32",
           "TOPPRA_REPLAY_V1 -1 1 32",
           "TOPPRA_REPLAY_V1 1 0 32",
           "TOPPRA_REPLAY_V1 1 1 0",
           "TOPPRA_REPLAY_V1 1 1 4097",
           "TOPPRA_REPLAY_V1 1 1 32 -1 0",
           "TOPPRA_REPLAY_V1 1 1 32 0 0 1 0",
           "TOPPRA_REPLAY_V1 1 1 32 0 0 0 1 -1",
           "TOPPRA_REPLAY_V1 1 1 32 0 0 0 1 1 1 7",
           "TOPPRA_REPLAY_V1 1 1 32 0 0 0 1 1 1 2 nan 0",
           "TOPPRA_REPLAY_V1 1 1 32 0 0 0 1 1 1 2 1 0 trailing"}) {
    std::stringstream input(text);
    EXPECT_THROW(readCase(input), std::runtime_error) << text;
  }
}

TEST(RegressionCorpus, EquivalentMutationsStayRepresentableAtExtremeUnits) {
  RecordedCase c;
  c.knots = {0, 1};
  c.stages = 32;
  c.startSquared = c.endSquared = 0;
  c.velocity.resize(2);
  c.velocity << std::numeric_limits<value_type>::max(),
      std::numeric_limits<value_type>::denorm_min();
  c.acceleration = c.velocity;
  Matrix coefficients(2, 2);
  coefficients.row(0) = c.velocity.transpose();
  coefficients.row(1).setZero();
  c.coefficients.push_back(coefficients);
  uint32_t state = 0;
  for (int i = 0; i < 256; ++i) {
    const auto mutated = equivalent(c, state);
    EXPECT_TRUE(mutated.velocity.allFinite());
    EXPECT_TRUE(mutated.acceleration.allFinite());
    EXPECT_TRUE(mutated.coefficients[0].allFinite());
    EXPECT_TRUE((mutated.velocity.array() > 0).all());
    EXPECT_EQ(Vector(mutated.coefficients[0].row(0).cwiseAbs().transpose()),
              mutated.velocity);
  }
}

TEST(RegressionCorpus, CaseStreamsAreIndependentOfCheckoutLocationAndCorpusOrder) {
  for (uint32_t seed : {0u, 0xffffffffu}) {
    EXPECT_EQ(caseSeed(seed, "/one/Recorded.case"), caseSeed(seed, "C:\\two\\Recorded.case"));
    EXPECT_EQ(caseSeed(seed, "Recorded.case"), caseSeed(seed, "/other/Recorded.case"));
    EXPECT_NE(caseSeed(seed, "Recorded.case"), caseSeed(seed, "Another.case"));
  }
}

TEST(RegressionCorpus, SeededEquivalentInputs) {
  const auto seed = option("TOPPRA_FUZZ_SEED", 0x5eed0302u);
  const auto count = option("TOPPRA_FUZZ_CASES", 64);
  ASSERT_GT(count, 0u);
  ASSERT_LE(count, 100000u);
  for (const char* file : kCorpusFiles) {
    uint32_t state = caseSeed(seed, file);
    const RecordedCase original = readFile(file);
    Vector expected;
    checkCase(original, &expected);
    if (failures()) return;
    for (uint32_t iteration = 0; iteration < count; ++iteration) {
      SCOPED_TRACE(::testing::Message() << "seed=" << seed << " file=" << file
                                        << " iteration=" << iteration);
      const RecordedCase c = equivalent(original, state);
      const int before = failures();
      Vector actual;
      checkCase(c, &actual);
      if (actual.size() == expected.size())
        EXPECT_LE(relativeDifference(actual, expected), 1e-10);
      else
        ADD_FAILURE() << "The equivalent input has no complete profile";
      if (failures() != before) {
        saveFailure(c, seed, file, iteration);
        return;
      }
    }
  }
}
