// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <cstdint>
#include <optional>
#include <unordered_map>
#include <vector>

namespace mpl2 {

struct SolverRect
{
  int id = 0;
  double lx = 0.0;
  double ly = 0.0;
  double ux = 0.0;
  double uy = 0.0;
};

struct SolverMacro
{
  int id = 0;
  double width = 0.0;
  double height = 0.0;
};

struct SolverWeightedNet
{
  int source = 0;
  int target = 0;
  double weight = 0.0;
};

struct SolverPlacement
{
  int id = 0;
  double x = 0.0;
  double y = 0.0;
};

struct SolverCore
{
  double lx = 0.0;
  double ly = 0.0;
  double ux = 0.0;
  double uy = 0.0;
};

class MacroPlacementSolver
{
 public:
  std::vector<SolverPlacement> legalizeMacros(
      const SolverCore& core,
      const std::vector<SolverRect>& movable_macros,
      const std::vector<SolverRect>& fixed_macros,
      double max_time_in_seconds = 10.0) const;

  std::vector<SolverPlacement> placeMacrosInCluster(
      double layout_width,
      double layout_height,
      const std::vector<SolverMacro>& macros,
      const std::vector<SolverWeightedNet>& nets,
      bool use_long_time_limit) const;

  std::vector<SolverPlacement> moveMacroGroupsToBoundary(
      const SolverCore& core,
      const std::vector<SolverMacro>& macro_groups,
      const std::unordered_map<int, SolverPlacement>& locations,
      const char* block_name) const;

 private:
  std::vector<SolverPlacement> coarseClusterSolve(
      double layout_width,
      double layout_height,
      const std::vector<SolverMacro>& macros,
      const std::vector<SolverWeightedNet>& nets,
      const std::unordered_map<int, SolverPlacement>& initial,
      double max_time_in_seconds) const;

  std::vector<SolverPlacement> fineClusterSolve(
      double layout_width,
      double layout_height,
      const std::vector<SolverMacro>& macros,
      const std::vector<SolverWeightedNet>& nets,
      const std::vector<SolverPlacement>& coarse_result,
      double max_time_in_seconds) const;
};

}  // namespace mpl2
