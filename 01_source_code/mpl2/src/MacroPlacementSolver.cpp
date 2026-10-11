// SPDX-License-Identifier: BSD-3-Clause

#include "MacroPlacementSolver.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/sat/model.h"
#include "ortools/util/sorted_interval_list.h"

namespace mpl2 {
namespace {

using operations_research::Domain;
using operations_research::sat::BoolVar;
using operations_research::sat::CpModelBuilder;
using operations_research::sat::CpSolverResponse;
using operations_research::sat::CpSolverStatus;
using operations_research::sat::IntVar;
using operations_research::sat::IntervalVar;
using operations_research::sat::LinearExpr;
using operations_research::sat::Model;
using operations_research::sat::NewSatParameters;
using operations_research::sat::SatParameters;
using operations_research::sat::SolutionIntegerValue;
using operations_research::sat::SolveCpModel;

constexpr int64_t kGridSize = 10;

int64_t roundToInt(const double value)
{
  return static_cast<int64_t>(std::llround(value));
}

void setTimeLimit(Model& model, const double max_time_in_seconds)
{
  SatParameters parameters;
  parameters.set_max_time_in_seconds(max_time_in_seconds);
  model.Add(NewSatParameters(parameters));
}

bool isFeasible(const CpSolverResponse& response)
{
  return response.status() == CpSolverStatus::OPTIMAL
         || response.status() == CpSolverStatus::FEASIBLE;
}

LinearExpr sumExprs(const std::vector<LinearExpr>& exprs)
{
  LinearExpr sum;
  for (const LinearExpr& expr : exprs) {
    sum += expr;
  }
  return sum;
}

std::unordered_map<int, SolverPlacement> toPlacementMap(
    const std::vector<SolverPlacement>& placements)
{
  std::unordered_map<int, SolverPlacement> map;
  for (const SolverPlacement& placement : placements) {
    map[placement.id] = placement;
  }
  return map;
}

std::unordered_map<int, SolverPlacement> greedyInitialPlacement(
    const double layout_width,
    const double layout_height,
    std::vector<SolverMacro> macros,
    const std::vector<SolverWeightedNet>& nets)
{
  std::unordered_map<int, double> incident_weight;
  for (const SolverWeightedNet& net : nets) {
    incident_weight[net.source] += net.weight;
    incident_weight[net.target] += net.weight;
  }
  std::sort(macros.begin(),
            macros.end(),
            [&incident_weight](const SolverMacro& lhs, const SolverMacro& rhs) {
              const double lhs_weight = incident_weight[lhs.id];
              const double rhs_weight = incident_weight[rhs.id];
              if (lhs_weight != rhs_weight) {
                return lhs_weight > rhs_weight;
              }
              const double lhs_area = lhs.width * lhs.height;
              const double rhs_area = rhs.width * rhs.height;
              if (lhs_area != rhs_area) {
                return lhs_area > rhs_area;
              }
              return lhs.id < rhs.id;
            });

  double x_pos = 0.0;
  double y_pos = layout_height;
  std::unordered_map<int, SolverPlacement> positions;
  for (const SolverMacro& macro : macros) {
    if (x_pos + macro.width > layout_width) {
      x_pos = 0.0;
      y_pos -= macro.height;
    }
    if (y_pos - macro.height < 0.0) {
      break;
    }
    positions[macro.id] = {macro.id, std::round(x_pos), std::round(y_pos - macro.height)};
    x_pos += macro.width;
  }
  return positions;
}


std::vector<SolverPlacement> completePlacementFallback(
    const std::vector<SolverMacro>& macros,
    const std::unordered_map<int, SolverPlacement>& initial)
{
  std::vector<SolverPlacement> fallback;
  fallback.reserve(macros.size());
  for (const SolverMacro& macro : macros) {
    const auto it = initial.find(macro.id);
    if (it != initial.end()) {
      fallback.push_back(it->second);
    } else {
      fallback.push_back({macro.id, 0.0, 0.0});
    }
  }
  return fallback;
}

int64_t gridSpan(const double value)
{
  return std::max<int64_t>(1, static_cast<int64_t>(std::ceil(value / kGridSize)));
}

bool overlapsWithFixed(const SolverPlacement& location,
                       const SolverMacro& shape,
                       const std::unordered_map<int, SolverMacro>& fixed_shapes,
                       const std::unordered_map<int, SolverPlacement>& fixed_locations,
                       const int64_t halo)
{
  const double lx = location.x - halo;
  const double ly = location.y - halo;
  const double ux = location.x + shape.width + halo;
  const double uy = location.y + shape.height + halo;
  for (const auto& [id, fixed_loc] : fixed_locations) {
    const SolverMacro& fixed_shape = fixed_shapes.at(id);
    const double flx = fixed_loc.x;
    const double fly = fixed_loc.y;
    const double fux = fixed_loc.x + fixed_shape.width;
    const double fuy = fixed_loc.y + fixed_shape.height;
    if (!(ux <= flx || lx >= fux || uy <= fly || ly >= fuy)) {
      return true;
    }
  }
  return false;
}

IntVar newAbsVar(CpModelBuilder& cp_model,
                 const LinearExpr& expr,
                 const int64_t upper_bound,
                 const std::string& name)
{
  IntVar value = cp_model.NewIntVar(Domain(0, std::max<int64_t>(0, upper_bound)))
                     .WithName(name);
  cp_model.AddAbsEquality(value, expr);
  return value;
}

IntVar addOverlapArea(CpModelBuilder& cp_model,
                      const IntVar& x,
                      const IntVar& y,
                      const SolverMacro& macro,
                      const SolverPlacement& fixed_loc,
                      const SolverMacro& fixed_shape,
                      const SolverCore& core,
                      const int64_t halo,
                      const int index)
{
  const int64_t core_r = roundToInt(core.ux);
  const int64_t core_t = roundToInt(core.uy);
  const int64_t width = roundToInt(macro.width);
  const int64_t height = roundToInt(macro.height);
  const int64_t fx = roundToInt(fixed_loc.x) - halo;
  const int64_t fy = roundToInt(fixed_loc.y) - halo;
  const int64_t fw = roundToInt(fixed_shape.width) + 2 * halo;
  const int64_t fh = roundToInt(fixed_shape.height) + 2 * halo;

  const IntVar end_x = cp_model.NewIntVar(Domain(0, core_r)).WithName("end_x_" + std::to_string(index));
  cp_model.AddEquality(end_x, x + width);
  const IntVar min_end_x = cp_model.NewIntVar(Domain(0, core_r)).WithName("min_end_x_" + std::to_string(index));
  cp_model.AddMinEquality(min_end_x, {end_x, fx + fw});
  const IntVar max_start_x = cp_model.NewIntVar(Domain(0, core_r)).WithName("max_start_x_" + std::to_string(index));
  cp_model.AddMaxEquality(max_start_x, {x, fx});
  const IntVar raw_w = cp_model.NewIntVar(Domain(-core_r, core_r)).WithName("raw_w_" + std::to_string(index));
  cp_model.AddEquality(raw_w, min_end_x - max_start_x);
  const IntVar overlap_w = cp_model.NewIntVar(Domain(0, core_r)).WithName("overlap_w_" + std::to_string(index));
  cp_model.AddMaxEquality(overlap_w, {raw_w, 0});

  const IntVar end_y = cp_model.NewIntVar(Domain(0, core_t)).WithName("end_y_" + std::to_string(index));
  cp_model.AddEquality(end_y, y + height);
  const IntVar min_end_y = cp_model.NewIntVar(Domain(0, core_t)).WithName("min_end_y_" + std::to_string(index));
  cp_model.AddMinEquality(min_end_y, {end_y, fy + fh});
  const IntVar max_start_y = cp_model.NewIntVar(Domain(0, core_t)).WithName("max_start_y_" + std::to_string(index));
  cp_model.AddMaxEquality(max_start_y, {y, fy});
  const IntVar raw_h = cp_model.NewIntVar(Domain(-core_t, core_t)).WithName("raw_h_" + std::to_string(index));
  cp_model.AddEquality(raw_h, min_end_y - max_start_y);
  const IntVar overlap_h = cp_model.NewIntVar(Domain(0, core_t)).WithName("overlap_h_" + std::to_string(index));
  cp_model.AddMaxEquality(overlap_h, {raw_h, 0});

  const IntVar area = cp_model.NewIntVar(Domain(0, std::max<int64_t>(0, width * height)))
                          .WithName("overlap_area_" + std::to_string(index));
  cp_model.AddMultiplicationEquality(area, overlap_w, overlap_h);
  return area;
}

SolverPlacement solveBoundaryStage(const SolverCore& core,
                                   const SolverMacro& macro,
                                   SolverPlacement location,
                                   const std::unordered_map<int, SolverMacro>& fixed_shapes,
                                   const std::unordered_map<int, SolverPlacement>& fixed_locations,
                                   const int64_t halo)
{
  const int64_t left_boundary = roundToInt(core.lx);
  const int64_t bottom_boundary = roundToInt(core.ly);
  const int64_t right_boundary = roundToInt(core.ux - macro.width);
  const int64_t top_boundary = roundToInt(core.uy - macro.height);

  CpModelBuilder first_model;
  const IntVar x1 = first_model.NewIntVar(Domain(left_boundary, right_boundary));
  const IntVar y1 = first_model.NewIntVar(Domain(bottom_boundary, top_boundary));
  const double left = location.x - core.lx;
  const double right = core.ux - location.x - macro.width;
  const double bottom = location.y - core.ly;
  const double top = core.uy - location.y - macro.height;
  char direction = 'L';
  double best = left;
  if (right < best) {
    best = right;
    direction = 'R';
  }
  if (bottom < best) {
    best = bottom;
    direction = 'B';
  }
  if (top < best) {
    direction = 'T';
  }

  const IntVar dx1 = first_model.NewIntVar(Domain(0, roundToInt(core.ux - core.lx)));
  const IntVar dy1 = first_model.NewIntVar(Domain(0, roundToInt(core.uy - core.ly)));
  if (direction == 'L') {
    first_model.AddEquality(dx1, x1 - roundToInt(core.lx));
    first_model.AddEquality(dy1, y1 - roundToInt(location.y));
    first_model.Minimize(100 * dx1 + dy1);
  } else if (direction == 'R') {
    first_model.AddEquality(dx1, roundToInt(core.ux) - x1 - roundToInt(macro.width));
    first_model.AddEquality(dy1, y1 - roundToInt(location.y));
    first_model.Minimize(100 * dx1 + dy1);
  } else if (direction == 'B') {
    first_model.AddEquality(dx1, x1 - roundToInt(location.x));
    first_model.AddEquality(dy1, y1 - roundToInt(core.ly));
    first_model.Minimize(dx1 + 100 * dy1);
  } else {
    first_model.AddEquality(dx1, x1 - roundToInt(location.x));
    first_model.AddEquality(dy1, roundToInt(core.uy) - y1 - roundToInt(macro.height));
    first_model.Minimize(dx1 + 100 * dy1);
  }

  Model first_sat;
  const CpSolverResponse first_response = SolveCpModel(first_model.Build(), &first_sat);
  if (isFeasible(first_response)) {
    location.x = SolutionIntegerValue(first_response, x1);
    location.y = SolutionIntegerValue(first_response, y1);
  }

  int64_t low_x = left_boundary;
  int64_t high_x = right_boundary;
  int64_t low_y = bottom_boundary;
  int64_t high_y = top_boundary;
  if (direction == 'L') {
    low_x = roundToInt(location.x);
    high_x = std::min<int64_t>(roundToInt(location.x + macro.width), right_boundary);
  } else if (direction == 'R') {
    low_x = std::max<int64_t>(left_boundary, roundToInt(location.x - 2 * macro.width));
    high_x = roundToInt(location.x);
  } else if (direction == 'B') {
    low_y = roundToInt(location.y);
    high_y = std::min<int64_t>(roundToInt(location.y + macro.height), top_boundary);
  } else {
    low_y = std::max<int64_t>(bottom_boundary, roundToInt(location.y - 2 * macro.height));
    high_y = roundToInt(location.y);
  }

  const double core_center_x = (core.lx + core.ux) / 2.0;
  const double core_center_y = (core.ly + core.uy) / 2.0;
  const double macro_center_x = location.x + macro.width / 2.0;
  const double macro_center_y = location.y + macro.height / 2.0;
  if (fixed_locations.empty()) {
    if (direction == 'L' || direction == 'R') {
      location.y = macro_center_y > core_center_y ? top_boundary : bottom_boundary;
    } else {
      location.x = macro_center_x > core_center_x ? right_boundary : left_boundary;
    }
    return location;
  }

  CpModelBuilder model;
  const IntVar x = model.NewIntVar(Domain(low_x, std::max(low_x, high_x)));
  const IntVar y = model.NewIntVar(Domain(low_y, std::max(low_y, high_y)));
  std::vector<LinearExpr> overlap_terms;
  int index = 0;
  for (const auto& [id, fixed_location] : fixed_locations) {
    overlap_terms.push_back(
        addOverlapArea(model, x, y, macro, fixed_location, fixed_shapes.at(id), core, halo, index++));
  }
  const int64_t macro_area = roundToInt(macro.width * macro.height);
  const IntVar total_overlap = model.NewIntVar(Domain(0, std::max<int64_t>(0, macro_area * fixed_locations.size())));
  model.AddEquality(total_overlap, sumExprs(overlap_terms));
  const IntVar dx_origin = newAbsVar(model, x - roundToInt(location.x), roundToInt(core.ux), "dx_origin");
  const IntVar dy_origin = newAbsVar(model, y - roundToInt(location.y), roundToInt(core.uy), "dy_origin");
  const IntVar origin_distance = model.NewIntVar(Domain(0, roundToInt(core.ux + core.uy)));
  model.AddEquality(origin_distance, dx_origin + dy_origin);

  if (overlapsWithFixed(location, macro, fixed_shapes, fixed_locations, halo)) {
    model.Minimize(1000 * total_overlap + origin_distance);
  } else {
    const IntVar corner_distance = model.NewIntVar(Domain(0, roundToInt(core.ux + core.uy)));
    if (direction == 'L' || direction == 'R') {
      if (macro_center_y >= core_center_y) {
        model.AddEquality(corner_distance, roundToInt(core.uy) - (y + roundToInt(macro.height)));
      } else {
        model.AddEquality(corner_distance, y - roundToInt(core.ly));
      }
    } else {
      if (macro_center_x >= core_center_x) {
        model.AddEquality(corner_distance, roundToInt(core.ux) - (x + roundToInt(macro.width)));
      } else {
        model.AddEquality(corner_distance, x - roundToInt(core.lx));
      }
    }
    model.Minimize(100 * corner_distance + origin_distance);
  }

  Model sat_model;
  const CpSolverResponse response = SolveCpModel(model.Build(), &sat_model);
  if (isFeasible(response)) {
    return {macro.id,
            static_cast<double>(SolutionIntegerValue(response, x)),
            static_cast<double>(SolutionIntegerValue(response, y))};
  }
  return location;
}

}  // namespace

std::vector<SolverPlacement> MacroPlacementSolver::legalizeMacros(
    const SolverCore& core,
    const std::vector<SolverRect>& movable_macros,
    const std::vector<SolverRect>& fixed_macros,
    const double max_time_in_seconds) const
{
  CpModelBuilder cp_model;
  std::unordered_map<int, IntVar> x_vars;
  std::unordered_map<int, IntVar> y_vars;
  auto no_overlap = cp_model.AddNoOverlap2D();

  const int64_t core_l = roundToInt(core.lx);
  const int64_t core_r = roundToInt(core.ux);
  const int64_t core_b = roundToInt(core.ly);
  const int64_t core_t = roundToInt(core.uy);

  for (const SolverRect& macro : movable_macros) {
    const int64_t width = roundToInt(macro.ux - macro.lx);
    const int64_t height = roundToInt(macro.uy - macro.ly);
    const IntVar x = cp_model.NewIntVar(Domain(core_l, core_r - width))
                         .WithName("x_" + std::to_string(macro.id));
    const IntVar y = cp_model.NewIntVar(Domain(core_b, core_t - height))
                         .WithName("y_" + std::to_string(macro.id));
    x_vars.emplace(macro.id, x);
    y_vars.emplace(macro.id, y);
    no_overlap.AddRectangle(cp_model.NewIntervalVar(x, width, x + width),
                            cp_model.NewIntervalVar(y, height, y + height));
  }

  for (const SolverRect& macro : movable_macros) {
    const int64_t width = roundToInt(macro.ux - macro.lx);
    const int64_t height = roundToInt(macro.uy - macro.ly);
    const IntVar& x = x_vars.at(macro.id);
    const IntVar& y = y_vars.at(macro.id);

    for (const SolverRect& fixed : fixed_macros) {
      const int64_t fx_l = roundToInt(fixed.lx);
      const int64_t fx_r = roundToInt(fixed.ux);
      const int64_t fy_b = roundToInt(fixed.ly);
      const int64_t fy_t = roundToInt(fixed.uy);

      const BoolVar left = cp_model.NewBoolVar();
      const BoolVar right = cp_model.NewBoolVar();
      const BoolVar down = cp_model.NewBoolVar();
      const BoolVar up = cp_model.NewBoolVar();

      cp_model.AddLessOrEqual(x + width, fx_l).OnlyEnforceIf(left);
      cp_model.AddGreaterThan(x + width, fx_l).OnlyEnforceIf(left.Not());
      cp_model.AddLessOrEqual(fx_r, x).OnlyEnforceIf(right);
      cp_model.AddGreaterThan(fx_r, x).OnlyEnforceIf(right.Not());
      cp_model.AddLessOrEqual(y + height, fy_b).OnlyEnforceIf(down);
      cp_model.AddGreaterThan(y + height, fy_b).OnlyEnforceIf(down.Not());
      cp_model.AddLessOrEqual(fy_t, y).OnlyEnforceIf(up);
      cp_model.AddGreaterThan(fy_t, y).OnlyEnforceIf(up.Not());
      cp_model.AddBoolOr({left, right, down, up});
    }
  }

  std::vector<LinearExpr> distances;
  for (const SolverRect& macro : movable_macros) {
    const int64_t width = roundToInt(macro.ux - macro.lx);
    const int64_t height = roundToInt(macro.uy - macro.ly);
    const IntVar dx = newAbsVar(cp_model,
                                x_vars.at(macro.id) - roundToInt(macro.lx),
                                core_r - width,
                                "dx_" + std::to_string(macro.id));
    const IntVar dy = newAbsVar(cp_model,
                                y_vars.at(macro.id) - roundToInt(macro.ly),
                                core_t - height,
                                "dy_" + std::to_string(macro.id));
    distances.push_back(dx + dy);
  }
  cp_model.Minimize(sumExprs(distances));

  Model model;
  setTimeLimit(model, max_time_in_seconds);
  const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);
  if (!isFeasible(response)) {
    return {};
  }

  std::vector<SolverPlacement> result;
  for (const SolverRect& macro : movable_macros) {
    result.push_back({macro.id,
                      static_cast<double>(SolutionIntegerValue(response, x_vars.at(macro.id))),
                      static_cast<double>(SolutionIntegerValue(response, y_vars.at(macro.id)))});
  }
  return result;
}

std::vector<SolverPlacement> MacroPlacementSolver::placeMacrosInCluster(
    const double layout_width,
    const double layout_height,
    const std::vector<SolverMacro>& macros,
    const std::vector<SolverWeightedNet>& nets,
    const bool use_long_time_limit) const
{
  const auto initial = greedyInitialPlacement(layout_width, layout_height, macros, nets);
  const double coarse_limit = use_long_time_limit ? 30.0 : (macros.size() > 10 ? 10.0 : 5.0);
  const double fine_limit = use_long_time_limit ? 30.0 : (macros.size() > 10 ? 10.0 : 5.0);
  const auto coarse = coarseClusterSolve(layout_width, layout_height, macros, nets, initial, coarse_limit);
  if (coarse.size() != macros.size()) {
    return completePlacementFallback(macros, initial);
  }
  return fineClusterSolve(layout_width, layout_height, macros, nets, coarse, fine_limit);
}

std::vector<SolverPlacement> MacroPlacementSolver::coarseClusterSolve(
    const double layout_width,
    const double layout_height,
    const std::vector<SolverMacro>& macros,
    const std::vector<SolverWeightedNet>& nets,
    const std::unordered_map<int, SolverPlacement>& initial,
    const double max_time_in_seconds) const
{
  CpModelBuilder cp_model;
  const int64_t coarse_width = roundToInt(layout_width / kGridSize);
  const int64_t coarse_height = roundToInt(layout_height / kGridSize);

  std::unordered_map<int, IntVar> x_vars;
  std::unordered_map<int, IntVar> y_vars;
  auto no_overlap = cp_model.AddNoOverlap2D();

  for (const SolverMacro& macro : macros) {
    const int64_t width = gridSpan(macro.width);
    const int64_t height = gridSpan(macro.height);
    const IntVar x = cp_model.NewIntVar(Domain(0, std::max<int64_t>(0, coarse_width - width)))
                         .WithName("x_" + std::to_string(macro.id));
    const IntVar y = cp_model.NewIntVar(Domain(0, std::max<int64_t>(0, coarse_height - height)))
                         .WithName("y_" + std::to_string(macro.id));
    x_vars.emplace(macro.id, x);
    y_vars.emplace(macro.id, y);
    no_overlap.AddRectangle(cp_model.NewIntervalVar(x, width, x + width),
                            cp_model.NewIntervalVar(y, height, y + height));

    const auto hint = initial.find(macro.id);
    if (hint != initial.end()) {
      cp_model.AddHint(x, roundToInt(hint->second.x / kGridSize));
      cp_model.AddHint(y, roundToInt(hint->second.y / kGridSize));
    }
  }

  std::vector<LinearExpr> distances;
  for (const SolverWeightedNet& net : nets) {
    if (x_vars.find(net.source) == x_vars.end() || x_vars.find(net.target) == x_vars.end()) {
      continue;
    }
    const auto source_macro = std::find_if(macros.begin(), macros.end(), [&net](const SolverMacro& macro) {
      return macro.id == net.source;
    });
    const auto target_macro = std::find_if(macros.begin(), macros.end(), [&net](const SolverMacro& macro) {
      return macro.id == net.target;
    });
    if (source_macro == macros.end() || target_macro == macros.end()) {
      continue;
    }
    const int64_t source_width = gridSpan(source_macro->width);
    const int64_t target_width = gridSpan(target_macro->width);
    const int64_t source_height = gridSpan(source_macro->height);
    const int64_t target_height = gridSpan(target_macro->height);
    const IntVar dx = newAbsVar(cp_model,
                                2 * x_vars.at(net.source) + source_width
                                    - (2 * x_vars.at(net.target) + target_width),
                                2 * coarse_width,
                                "dx_" + std::to_string(net.source) + "_" + std::to_string(net.target));
    const IntVar dy = newAbsVar(cp_model,
                                2 * y_vars.at(net.source) + source_height
                                    - (2 * y_vars.at(net.target) + target_height),
                                2 * coarse_height,
                                "dy_" + std::to_string(net.source) + "_" + std::to_string(net.target));
    distances.push_back(roundToInt(net.weight) * (dx + dy));
  }
  cp_model.Minimize(sumExprs(distances));

  Model model;
  setTimeLimit(model, max_time_in_seconds);
  const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);
  if (!isFeasible(response)) {
    return {};
  }

  std::vector<SolverPlacement> result;
  for (const SolverMacro& macro : macros) {
    result.push_back({macro.id,
                      static_cast<double>(SolutionIntegerValue(response, x_vars.at(macro.id))),
                      static_cast<double>(SolutionIntegerValue(response, y_vars.at(macro.id)))});
  }
  return result;
}

std::vector<SolverPlacement> MacroPlacementSolver::fineClusterSolve(
    const double layout_width,
    const double layout_height,
    const std::vector<SolverMacro>& macros,
    const std::vector<SolverWeightedNet>& nets,
    const std::vector<SolverPlacement>& coarse_result,
    const double max_time_in_seconds) const
{
  CpModelBuilder cp_model;
  const auto coarse_map = toPlacementMap(coarse_result);
  std::unordered_map<int, IntVar> x_vars;
  std::unordered_map<int, IntVar> y_vars;
  auto no_overlap = cp_model.AddNoOverlap2D();

  for (const SolverMacro& macro : macros) {
    const int64_t width = roundToInt(macro.width);
    const int64_t height = roundToInt(macro.height);
    int64_t cx = 0;
    int64_t cy = 0;
    const auto coarse_it = coarse_map.find(macro.id);
    if (coarse_it != coarse_map.end()) {
      cx = roundToInt(coarse_it->second.x * kGridSize);
      cy = roundToInt(coarse_it->second.y * kGridSize);
    }
    const int64_t max_x = std::max<int64_t>(0, roundToInt(layout_width) - width);
    const int64_t max_y = std::max<int64_t>(0, roundToInt(layout_height) - height);
    const int64_t low_x = std::max<int64_t>(cx - 2 * kGridSize, 0);
    const int64_t high_x = std::min<int64_t>(cx + 2 * kGridSize, max_x);
    const int64_t low_y = std::max<int64_t>(cy - 2 * kGridSize, 0);
    const int64_t high_y = std::min<int64_t>(cy + 2 * kGridSize, max_y);
    const IntVar x = cp_model.NewIntVar(Domain(low_x, std::max(low_x, high_x)))
                         .WithName("x_" + std::to_string(macro.id));
    const IntVar y = cp_model.NewIntVar(Domain(low_y, std::max(low_y, high_y)))
                         .WithName("y_" + std::to_string(macro.id));
    x_vars.emplace(macro.id, x);
    y_vars.emplace(macro.id, y);
    no_overlap.AddRectangle(cp_model.NewIntervalVar(x, width, x + width),
                            cp_model.NewIntervalVar(y, height, y + height));
  }

  std::vector<LinearExpr> distances;
  for (const SolverWeightedNet& net : nets) {
    if (x_vars.find(net.source) == x_vars.end() || x_vars.find(net.target) == x_vars.end()) {
      continue;
    }
    const auto source_macro = std::find_if(macros.begin(), macros.end(), [&net](const SolverMacro& macro) {
      return macro.id == net.source;
    });
    const auto target_macro = std::find_if(macros.begin(), macros.end(), [&net](const SolverMacro& macro) {
      return macro.id == net.target;
    });
    if (source_macro == macros.end() || target_macro == macros.end()) {
      continue;
    }
    const int64_t source_width = roundToInt(source_macro->width);
    const int64_t target_width = roundToInt(target_macro->width);
    const int64_t source_height = roundToInt(source_macro->height);
    const int64_t target_height = roundToInt(target_macro->height);
    const IntVar dx = newAbsVar(cp_model,
                                2 * x_vars.at(net.source) + source_width
                                    - (2 * x_vars.at(net.target) + target_width),
                                2 * roundToInt(layout_width),
                                "dx_" + std::to_string(net.source) + "_" + std::to_string(net.target));
    const IntVar dy = newAbsVar(cp_model,
                                2 * y_vars.at(net.source) + source_height
                                    - (2 * y_vars.at(net.target) + target_height),
                                2 * roundToInt(layout_height),
                                "dy_" + std::to_string(net.source) + "_" + std::to_string(net.target));
    distances.push_back(roundToInt(net.weight) * (dx + dy));
  }
  cp_model.Minimize(sumExprs(distances));

  Model model;
  setTimeLimit(model, max_time_in_seconds);
  const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);
  if (!isFeasible(response)) {
    if (coarse_result.size() != macros.size()) {
      return completePlacementFallback(macros, {});
    }
    std::vector<SolverPlacement> fallback;
    for (const SolverPlacement& placement : coarse_result) {
      fallback.push_back({placement.id, placement.x * kGridSize, placement.y * kGridSize});
    }
    return fallback;
  }

  std::vector<SolverPlacement> result;
  for (const SolverMacro& macro : macros) {
    result.push_back({macro.id,
                      static_cast<double>(SolutionIntegerValue(response, x_vars.at(macro.id))),
                      static_cast<double>(SolutionIntegerValue(response, y_vars.at(macro.id)))});
  }
  return result;
}

std::vector<SolverPlacement> MacroPlacementSolver::moveMacroGroupsToBoundary(
    const SolverCore& core,
    const std::vector<SolverMacro>& macro_groups,
    const std::unordered_map<int, SolverPlacement>& locations,
    const char* block_name) const
{
  int64_t halo = 10;
  const std::string name = block_name == nullptr ? "" : block_name;
  if (name == "mor1kx" || name == "eth_top" || name == "or1200_top") {
    halo = 2;
  }

  std::vector<SolverMacro> sorted_groups = macro_groups;
  std::sort(sorted_groups.begin(),
            sorted_groups.end(),
            [](const SolverMacro& lhs, const SolverMacro& rhs) {
              const double lhs_area = lhs.width * lhs.height;
              const double rhs_area = rhs.width * rhs.height;
              if (lhs_area == rhs_area) {
                return lhs.id < rhs.id;
              }
              return lhs_area > rhs_area;
            });

  std::unordered_map<int, SolverMacro> fixed_shapes;
  std::unordered_map<int, SolverPlacement> fixed_locations;
  for (const SolverMacro& group : sorted_groups) {
    const auto loc_it = locations.find(group.id);
    if (loc_it == locations.end()) {
      continue;
    }
    const SolverPlacement moved = solveBoundaryStage(
        core, group, loc_it->second, fixed_shapes, fixed_locations, halo);
    fixed_shapes[group.id] = group;
    fixed_locations[group.id] = moved;
  }

  if (sorted_groups.size() <= 1) {
    std::vector<SolverPlacement> result;
    for (const auto& [id, placement] : fixed_locations) {
      result.push_back(placement);
    }
    return result;
  }

  CpModelBuilder model;
  std::unordered_map<int, IntVar> x_vars;
  std::unordered_map<int, IntVar> y_vars;
  auto no_overlap = model.AddNoOverlap2D();
  bool first = true;
  int fixed_id = -1;
  for (const SolverMacro& group : sorted_groups) {
    const SolverPlacement& location = fixed_locations.at(group.id);
    if (first) {
      first = false;
      fixed_id = group.id;
      continue;
    }
    const int64_t width = roundToInt(group.width);
    const int64_t height = roundToInt(group.height);
    const IntVar x = model.NewIntVar(Domain(roundToInt(core.lx), roundToInt(core.ux) - width))
                         .WithName("x_" + std::to_string(group.id));
    const IntVar y = model.NewIntVar(Domain(roundToInt(core.ly), roundToInt(core.uy) - height))
                         .WithName("y_" + std::to_string(group.id));
    x_vars.emplace(group.id, x);
    y_vars.emplace(group.id, y);
    no_overlap.AddRectangle(model.NewIntervalVar(x - halo, width + 2 * halo, x + width + halo),
                            model.NewIntervalVar(y - halo, height + 2 * halo, y + height + halo));
    (void) location;
  }

  const SolverMacro& fixed_shape = fixed_shapes.at(fixed_id);
  const SolverPlacement& fixed_location = fixed_locations.at(fixed_id);
  const int64_t fx_l = roundToInt(fixed_location.x) - halo;
  const int64_t fx_r = roundToInt(fixed_location.x + fixed_shape.width) + halo;
  const int64_t fy_b = roundToInt(fixed_location.y) - halo;
  const int64_t fy_t = roundToInt(fixed_location.y + fixed_shape.height) + halo;
  for (const SolverMacro& group : sorted_groups) {
    if (group.id == fixed_id) {
      continue;
    }
    const int64_t width = roundToInt(group.width);
    const int64_t height = roundToInt(group.height);
    const IntVar& x = x_vars.at(group.id);
    const IntVar& y = y_vars.at(group.id);
    const BoolVar left = model.NewBoolVar();
    const BoolVar right = model.NewBoolVar();
    const BoolVar down = model.NewBoolVar();
    const BoolVar up = model.NewBoolVar();
    model.AddLessOrEqual(x + width, fx_l).OnlyEnforceIf(left);
    model.AddGreaterThan(x + width, fx_l).OnlyEnforceIf(left.Not());
    model.AddLessOrEqual(fx_r, x).OnlyEnforceIf(right);
    model.AddGreaterThan(fx_r, x).OnlyEnforceIf(right.Not());
    model.AddLessOrEqual(y + height, fy_b).OnlyEnforceIf(down);
    model.AddGreaterThan(y + height, fy_b).OnlyEnforceIf(down.Not());
    model.AddLessOrEqual(fy_t, y).OnlyEnforceIf(up);
    model.AddGreaterThan(fy_t, y).OnlyEnforceIf(up.Not());
    model.AddBoolOr({left, right, down, up});
  }

  const double core_width = core.ux - core.lx;
  const double core_height = core.uy - core.ly;
  std::vector<LinearExpr> distances;
  for (const SolverMacro& group : sorted_groups) {
    if (group.id == fixed_id) {
      continue;
    }
    const SolverPlacement& origin = fixed_locations.at(group.id);
    const int64_t width = roundToInt(group.width);
    const int64_t height = roundToInt(group.height);
    const IntVar dx = newAbsVar(model, x_vars.at(group.id) - roundToInt(origin.x), roundToInt(core.ux), "dx_" + std::to_string(group.id));
    const IntVar dy = newAbsVar(model, y_vars.at(group.id) - roundToInt(origin.y), roundToInt(core.uy), "dy_" + std::to_string(group.id));
    const int64_t weight = roundToInt((group.width / 20.0) * (group.height / 20.0));
    if (origin.x > core.lx + 0.2 * core_width && origin.x < core.ux - 0.2 * core_width
        && origin.y < core.ly + 0.2 * core_height) {
      distances.push_back((dx + 5 * dy) * weight);
    } else if (origin.x > core.lx + 0.2 * core_width && origin.x < core.ux - 0.2 * core_width
               && origin.y > core.uy - 0.2 * core_height) {
      distances.push_back((dx + 5 * dy) * weight);
    } else if (origin.y > core.ly + 0.2 * core_height && origin.y < core.uy - 0.2 * core_height
               && origin.x < core.lx + 0.2 * core_width) {
      distances.push_back((5 * dx + dy) * weight);
    } else if (origin.y > core.ly + 0.2 * core_height && origin.y < core.uy - 0.2 * core_height
               && origin.x > core.ux - 0.2 * core_width) {
      distances.push_back((5 * dx + dy) * weight);
    } else {
      distances.push_back((dx + dy) * weight);
    }
    (void) width;
    (void) height;
  }
  model.Minimize(sumExprs(distances));

  Model sat_model;
  setTimeLimit(sat_model, 10.0);
  const CpSolverResponse response = SolveCpModel(model.Build(), &sat_model);
  std::vector<SolverPlacement> result;
  result.push_back(fixed_location);
  for (const SolverMacro& group : sorted_groups) {
    if (group.id == fixed_id) {
      continue;
    }
    if (isFeasible(response)) {
      result.push_back({group.id,
                        static_cast<double>(SolutionIntegerValue(response, x_vars.at(group.id))),
                        static_cast<double>(SolutionIntegerValue(response, y_vars.at(group.id)))});
    } else {
      result.push_back(fixed_locations.at(group.id));
    }
  }
  return result;
}

}  // namespace mpl2
