#include "../../src/MacroPlacementSolver.h"

#include <unordered_map>

#include "gtest/gtest.h"

namespace mpl2 {

TEST(MacroPlacementSolver, LegalizesOverlappingMacrosInCore)
{
  MacroPlacementSolver solver;
  const SolverCore core{0, 0, 100, 100};
  const std::vector<SolverRect> movable{{0, 0, 0, 40, 40}, {1, 10, 10, 50, 50}};

  const auto result = solver.legalizeMacros(core, movable, {});

  ASSERT_EQ(result.size(), 2);
  std::unordered_map<int, SolverPlacement> by_id;
  for (const SolverPlacement& placement : result) {
    by_id[placement.id] = placement;
  }
  EXPECT_GE(by_id.at(0).x, 0);
  EXPECT_GE(by_id.at(0).y, 0);
  EXPECT_LE(by_id.at(0).x + 40, 100);
  EXPECT_LE(by_id.at(0).y + 40, 100);
  EXPECT_GE(by_id.at(1).x, 0);
  EXPECT_GE(by_id.at(1).y, 0);
  EXPECT_LE(by_id.at(1).x + 40, 100);
  EXPECT_LE(by_id.at(1).y + 40, 100);

  const bool separated_x = by_id.at(0).x + 40 <= by_id.at(1).x
                           || by_id.at(1).x + 40 <= by_id.at(0).x;
  const bool separated_y = by_id.at(0).y + 40 <= by_id.at(1).y
                           || by_id.at(1).y + 40 <= by_id.at(0).y;
  EXPECT_TRUE(separated_x || separated_y);
}

TEST(MacroPlacementSolver, ReportsInfeasibleLegalizationWithEmptyResult)
{
  MacroPlacementSolver solver;
  const SolverCore core{0, 0, 50, 50};
  const std::vector<SolverRect> movable{{0, 0, 0, 40, 40}, {1, 0, 0, 40, 40}};

  const auto result = solver.legalizeMacros(core, movable, {});

  EXPECT_TRUE(result.empty());
}

TEST(MacroPlacementSolver, PlacesClusterMacrosWithoutOverlap)
{
  MacroPlacementSolver solver;
  const std::vector<SolverMacro> macros{{0, 30, 30}, {1, 30, 30}, {2, 20, 20}};
  const std::vector<SolverWeightedNet> nets{{0, 1, 2}, {1, 2, 1}};

  const auto result = solver.placeMacrosInCluster(100, 100, macros, nets, false);

  ASSERT_EQ(result.size(), 3);
  std::unordered_map<int, SolverPlacement> by_id;
  for (const SolverPlacement& placement : result) {
    by_id[placement.id] = placement;
    EXPECT_GE(placement.x, 0);
    EXPECT_GE(placement.y, 0);
  }
  for (size_t i = 0; i < macros.size(); ++i) {
    for (size_t j = i + 1; j < macros.size(); ++j) {
      const auto& a = macros[i];
      const auto& b = macros[j];
      const auto& pa = by_id.at(a.id);
      const auto& pb = by_id.at(b.id);
      const bool separated_x = pa.x + a.width <= pb.x || pb.x + b.width <= pa.x;
      const bool separated_y = pa.y + a.height <= pb.y || pb.y + b.height <= pa.y;
      EXPECT_TRUE(separated_x || separated_y);
    }
  }
}


TEST(MacroPlacementSolver, OptimizesConnectedMacrosByCenterDistance)
{
  MacroPlacementSolver solver;
  const std::vector<SolverMacro> macros{{0, 80, 10}, {1, 10, 10}};
  const std::vector<SolverWeightedNet> nets{{0, 1, 100}};

  const auto result = solver.placeMacrosInCluster(100, 20, macros, nets, true);

  ASSERT_EQ(result.size(), 2);
  std::unordered_map<int, SolverPlacement> by_id;
  for (const SolverPlacement& placement : result) {
    by_id[placement.id] = placement;
  }
  const double wide_center_x = by_id.at(0).x + macros[0].width / 2.0;
  const double small_center_x = by_id.at(1).x + macros[1].width / 2.0;
  EXPECT_LE(std::abs(wide_center_x - small_center_x), 10.0);
}

TEST(MacroPlacementSolver, KeepsAllClusterMacrosWhenFineSolveIsInfeasible)
{
  MacroPlacementSolver solver;
  const std::vector<SolverMacro> macros{{0, 80, 80}, {1, 80, 80}};
  const std::vector<SolverWeightedNet> nets{{0, 1, 1}};

  const auto result = solver.placeMacrosInCluster(100, 100, macros, nets, false);

  EXPECT_EQ(result.size(), macros.size());
}

}  // namespace mpl2
