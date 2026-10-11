// SPDX-License-Identifier: BSD-3-Clause

#include "hier_rtlmp.h"

#include <algorithm>
#include <array>
#include <deque>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

#include "db_sta/dbNetwork.hh"
#include "db_sta/dbSta.hh"
#include "object.h"
#include "odb/db.h"
#include "sta/Graph.hh"
#include "sta/Liberty.hh"
#include "sta/MinMax.hh"
#include "utl/Logger.h"

namespace mpl2 {

using utl::MPL;

namespace {

// Tunable knobs for the first ASP-DAC incremental refinement implementation.
// They are intentionally conservative because the flow starts from a legalized
// standard-cell placement; moving a macro too far can improve one local path but
// damage many already-placed cells around it.
constexpr double kArraySwapMinCostReductionRatio = 0.05;
constexpr double kFreeMoveMinCostReductionRatio = 0.10;
constexpr double kFreeMoveMaxDisplacementMicrons = 5.0;
constexpr int kFreeMoveMaxIterations = 18;
constexpr double kFreeMoveDisplacementPenalty = 2.0;
constexpr double kFreeMoveInitialStep = 0.10;
constexpr double kFreeMoveStepDecay = 0.10;
constexpr double kOrientMinHpwlReductionRatio = 0.05;

bool isHardMacro(odb::dbInst* inst)
{
  return inst != nullptr && inst->getMaster() != nullptr && inst->getMaster()->isBlock();
}

bool isFixedMacro(odb::dbInst* inst)
{
  const odb::dbPlacementStatus status = inst->getPlacementStatus();
  return status == odb::dbPlacementStatus::LOCKED
         || status == odb::dbPlacementStatus::FIRM;
}

bool isSignalNet(odb::dbNet* net)
{
  return net != nullptr && net->getSigType() == odb::dbSigType::SIGNAL;
}

bool isSequential(odb::dbInst* inst, sta::dbNetwork* network)
{
  if (inst == nullptr || network == nullptr) {
    return false;
  }
  const sta::LibertyCell* liberty_cell = network->libertyCell(inst);
  return liberty_cell != nullptr && liberty_cell->hasSequentials();
}

bool isSinkLike(odb::dbITerm* iterm)
{
  return iterm->getIoType() == odb::dbIoType::INPUT
         || iterm->getIoType() == odb::dbIoType::INOUT;
}

bool isDriverLike(odb::dbITerm* iterm)
{
  return iterm->getIoType() == odb::dbIoType::OUTPUT
         || iterm->getIoType() == odb::dbIoType::INOUT;
}

std::optional<float> pinSlack(odb::dbITerm* iterm,
                              sta::dbSta* sta,
                              sta::dbNetwork* network)
{
  if (iterm == nullptr || sta == nullptr || network == nullptr) {
    return std::nullopt;
  }

  // The caller is responsible for sourcing Liberty/SDC and estimating
  // parasitics after standard-cell placement.  Here we only query the current
  // OpenSTA graph, so the slack is exactly for the placement loaded in OpenDB.
  sta::Pin* pin = network->dbToSta(iterm);
  if (pin == nullptr || sta->graph() == nullptr) {
    return std::nullopt;
  }

  sta::Vertex* vertex = sta->graph()->pinLoadVertex(pin);
  if (vertex == nullptr) {
    vertex = sta->graph()->pinDrvrVertex(pin);
  }
  if (vertex == nullptr) {
    return std::nullopt;
  }
  return sta->vertexSlack(vertex, sta::MinMax::max());
}

std::optional<float> worstInstanceSlack(odb::dbInst* inst,
                                        sta::dbSta* sta,
                                        sta::dbNetwork* network)
{
  std::optional<float> worst;
  if (inst == nullptr) {
    return worst;
  }
  for (odb::dbITerm* iterm : inst->getITerms()) {
    const std::optional<float> slack = pinSlack(iterm, sta, network);
    if (!slack.has_value()) {
      continue;
    }
    if (!worst.has_value() || slack.value() < worst.value()) {
      worst = slack.value();
    }
  }
  return worst;
}

constexpr int kMaxEndpointTraceDepth = 30;
constexpr int kMaxEndpointTracePins = 20000;

void collectDownstreamRegisters(odb::dbInst* root_macro,
                                odb::dbITerm* driver_pin,
                                sta::dbNetwork* network,
                                std::set<odb::dbInst*>& found_regs,
                                std::set<odb::dbITerm*>& visited,
                                int)
{
  if (driver_pin == nullptr) {
    return;
  }

  // Use level-by-level BFS so each macro pin is associated with the nearest
  // sequential endpoints only.  This mirrors the paper's wording and avoids a
  // deep fanout cone overwhelming the refinement objective with distant regs.
  std::deque<std::pair<odb::dbITerm*, int>> queue;
  queue.push_back({driver_pin, 0});
  while (!queue.empty() && found_regs.empty() && visited.size() <= kMaxEndpointTracePins) {
    const auto [pin, depth] = queue.front();
    queue.pop_front();
    if (pin == nullptr || depth > kMaxEndpointTraceDepth || !visited.insert(pin).second) {
      continue;
    }
    odb::dbNet* net = pin->getNet();
    if (!isSignalNet(net)) {
      continue;
    }
    for (odb::dbITerm* sink_pin : net->getITerms()) {
      odb::dbInst* sink_inst = sink_pin->getInst();
      if (sink_inst == nullptr || sink_inst == root_macro || !isSinkLike(sink_pin)) {
        continue;
      }
      if (isSequential(sink_inst, network)) {
        found_regs.insert(sink_inst);
        continue;
      }
      if (isHardMacro(sink_inst)) {
        continue;
      }
      for (odb::dbITerm* next_pin : sink_inst->getITerms()) {
        if (isDriverLike(next_pin)) {
          queue.push_back({next_pin, depth + 1});
        }
      }
    }
  }
}

void collectUpstreamRegisters(odb::dbInst* root_macro,
                              odb::dbITerm* sink_pin,
                              sta::dbNetwork* network,
                              std::set<odb::dbInst*>& found_regs,
                              std::set<odb::dbITerm*>& visited,
                              int)
{
  if (sink_pin == nullptr) {
    return;
  }

  std::deque<std::pair<odb::dbITerm*, int>> queue;
  queue.push_back({sink_pin, 0});
  while (!queue.empty() && found_regs.empty() && visited.size() <= kMaxEndpointTracePins) {
    const auto [pin, depth] = queue.front();
    queue.pop_front();
    if (pin == nullptr || depth > kMaxEndpointTraceDepth || !visited.insert(pin).second) {
      continue;
    }
    odb::dbNet* net = pin->getNet();
    if (!isSignalNet(net)) {
      continue;
    }
    for (odb::dbITerm* driver_pin : net->getITerms()) {
      odb::dbInst* driver_inst = driver_pin->getInst();
      if (driver_inst == nullptr || driver_inst == root_macro || !isDriverLike(driver_pin)) {
        continue;
      }
      if (isSequential(driver_inst, network)) {
        found_regs.insert(driver_inst);
        continue;
      }
      if (isHardMacro(driver_inst)) {
        continue;
      }
      for (odb::dbITerm* prev_pin : driver_inst->getITerms()) {
        if (isSinkLike(prev_pin)) {
          queue.push_back({prev_pin, depth + 1});
        }
      }
    }
  }
}

std::vector<odb::dbInst*> collectMacros(odb::dbBlock* block)
{
  std::vector<odb::dbInst*> macros;
  if (block == nullptr) {
    return macros;
  }
  for (odb::dbInst* inst : block->getInsts()) {
    if (isHardMacro(inst)) {
      macros.push_back(inst);
    }
  }
  std::sort(macros.begin(), macros.end(), [](odb::dbInst* lhs, odb::dbInst* rhs) {
    return lhs->getName() < rhs->getName();
  });
  return macros;
}


int centerX(odb::dbInst* inst)
{
  return inst->getBBox()->getBox().xCenter();
}

int centerY(odb::dbInst* inst)
{
  return inst->getBBox()->getBox().yCenter();
}

bool nameLess(odb::dbInst* lhs, odb::dbInst* rhs)
{
  return lhs->getName() < rhs->getName();
}

bool hasArrayLikeAlignment(const std::vector<odb::dbInst*>& macros,
                           odb::dbBlock* block)
{
  if (macros.size() < 2 || block == nullptr) {
    return false;
  }

  // This is the lightweight ASP-DAC D-stage classifier. The thesis version can
  // later replace it with DBSCAN/outlier handling, but for now we only need a
  // stable way to separate memory-array macros from freely movable macros after
  // standard-cell placement. Same-master macros whose centers share a row or a
  // column are treated as array constrained.
  const int align_tolerance_dbu = block->micronsToDbu(2.0);
  for (size_t i = 0; i < macros.size(); i++) {
    for (size_t j = i + 1; j < macros.size(); j++) {
      const bool same_column = std::abs(centerX(macros[i]) - centerX(macros[j]))
                               <= align_tolerance_dbu;
      const bool same_row = std::abs(centerY(macros[i]) - centerY(macros[j]))
                            <= align_tolerance_dbu;
      if (same_column || same_row) {
        return true;
      }
    }
  }
  return false;
}


double macroCriticalityWeight(odb::dbInst* reg,
                              const std::map<odb::dbInst*, float>& reg_slack)
{
  const auto slack_it = reg_slack.find(reg);
  if (slack_it == reg_slack.end()) {
    return 1.0;
  }
  // Negative slack is already in nanoseconds.  Add a small floor so all traced
  // endpoints keep influence, while worse endpoints pull harder.
  return std::max(0.001, static_cast<double>(-slack_it->second));
}

double macroCostAtCenter(
    odb::dbInst* macro,
    double center_x,
    double center_y,
    const std::map<odb::dbInst*, std::set<odb::dbInst*>>& macro_to_regs,
    const std::map<odb::dbInst*, float>& reg_slack)
{
  const auto regs_it = macro_to_regs.find(macro);
  if (regs_it == macro_to_regs.end() || regs_it->second.empty()) {
    return 0.0;
  }

  double cost = 0.0;
  for (odb::dbInst* reg : regs_it->second) {
    const odb::Rect& reg_box = reg->getBBox()->getBox();
    const double dx = center_x - reg_box.xCenter();
    const double dy = center_y - reg_box.yCenter();
    cost += macroCriticalityWeight(reg, reg_slack) * (dx * dx + dy * dy);
  }
  return cost;
}

double macroCostAtOrigin(
    odb::dbInst* macro,
    int origin_x,
    int origin_y,
    const std::map<odb::dbInst*, std::set<odb::dbInst*>>& macro_to_regs,
    const std::map<odb::dbInst*, float>& reg_slack)
{
  const odb::Rect& box = macro->getBBox()->getBox();
  const double center_x = origin_x + box.dx() / 2.0;
  const double center_y = origin_y + box.dy() / 2.0;
  return macroCostAtCenter(macro, center_x, center_y, macro_to_regs, reg_slack);
}

float macroWorstSlack(odb::dbInst* macro,
                      const std::unordered_map<odb::dbInst*, std::vector<float>>& macro_slack)
{
  const auto slack_it = macro_slack.find(macro);
  if (slack_it == macro_slack.end() || slack_it->second.empty()) {
    return 0.0f;
  }
  return *std::min_element(slack_it->second.begin(), slack_it->second.end());
}

odb::Rect rectAtOrigin(odb::dbInst* macro, int origin_x, int origin_y)
{
  const odb::Rect& box = macro->getBBox()->getBox();
  return odb::Rect(origin_x, origin_y, origin_x + box.dx(), origin_y + box.dy());
}

bool hasMacroOverlapAt(odb::dbInst* macro,
                       const odb::Rect& candidate_box,
                       odb::dbBlock* block)
{
  for (odb::dbInst* inst : block->getInsts()) {
    if (inst == macro || !isHardMacro(inst) || !inst->isPlaced()) {
      continue;
    }
    if (candidate_box.overlaps(inst->getBBox()->getBox())) {
      return true;
    }
  }
  return false;
}

bool isLegalMacroBox(odb::dbInst* macro,
                     const odb::Rect& candidate_box,
                     odb::dbBlock* block)
{
  return block != nullptr && block->getCoreArea().contains(candidate_box)
         && !hasMacroOverlapAt(macro, candidate_box, block);
}

int clampOriginX(odb::dbInst* macro, int x, const odb::Rect& core)
{
  return std::clamp(x, core.xMin(), core.xMax() - static_cast<int>(macro->getBBox()->getDX()));
}

int clampOriginY(odb::dbInst* macro, int y, const odb::Rect& core)
{
  return std::clamp(y, core.yMin(), core.yMax() - static_cast<int>(macro->getBBox()->getDY()));
}

std::pair<double, double> weightedRegisterCenter(
    odb::dbInst* macro,
    const std::map<odb::dbInst*, std::set<odb::dbInst*>>& macro_to_regs,
    const std::map<odb::dbInst*, float>& reg_slack)
{
  const auto regs_it = macro_to_regs.find(macro);
  if (regs_it == macro_to_regs.end() || regs_it->second.empty()) {
    const odb::Rect& box = macro->getBBox()->getBox();
    return {static_cast<double>(box.xCenter()), static_cast<double>(box.yCenter())};
  }

  double sum_weight = 0.0;
  double target_x = 0.0;
  double target_y = 0.0;
  for (odb::dbInst* reg : regs_it->second) {
    const double weight = macroCriticalityWeight(reg, reg_slack);
    const odb::Rect& reg_box = reg->getBBox()->getBox();
    target_x += weight * reg_box.xCenter();
    target_y += weight * reg_box.yCenter();
    sum_weight += weight;
  }
  return {target_x / sum_weight, target_y / sum_weight};
}

bool hasMeaningfulCostReduction(double before, double after, double min_ratio)
{
  if (before <= 0.0) {
    return after < before;
  }
  return after <= before * (1.0 - min_ratio);
}

int clampNearInitial(int value, int initial, int max_delta)
{
  return std::clamp(value, initial - max_delta, initial + max_delta);
}

}  // namespace

bool HierRTLMP::runIncrementalMacroPlacement()
{
  return IncrementalMacroPlacement();
}

bool HierRTLMP::IncrementalMacroPlacement()
{
  logger_->info(MPL, 1026, "Starting incremental macro placement refinement.");
  IncreInitialize();
  CriticalMacroFilter();
  MacroClassify();

  logger_->info(MPL,
                1027,
                "Incremental macro placement diagnostics: {} critical macros, {} array groups, {} free macros.",
                critical_macro_to_regs_.size(),
                array_critical_macro_insts_.size(),
                free_critical_macro_insts_.size());

  if (critical_macro_to_regs_.empty()) {
    logger_->info(MPL, 1028, "No critical macros were selected; placement is unchanged.");
    return false;
  }

  const int accepted_swaps = TilingMacroExchange();
  const int moved_macros = FreeMacroOptimize();
  const int changed_orients = MacroOrientOptimizer();
  const bool changed = accepted_swaps > 0 || moved_macros > 0 || changed_orients > 0;
  logger_->info(MPL,
                1040,
                "Incremental macro placement changed macros: {} ({} swaps, {} moves, {} orientation changes).",
                changed ? "yes" : "no",
                accepted_swaps,
                moved_macros,
                changed_orients);
  return changed;
}

void HierRTLMP::IncreInitialize()
{
  critical_macro_to_regs_.clear();
  critical_regs_slack_.clear();
  macro_delay_fit.clear();
  tile_critical_macros_.clear();
  free_critical_macros_.clear();
  array_critical_macro_insts_.clear();
  free_critical_macro_insts_.clear();
  critical_macro_slack_.clear();
  macro_pin_regs_.clear();
  placed_macros_.clear();

  block_ = db_->getChip()->getBlock();
  const std::vector<odb::dbInst*> macros = collectMacros(block_);
  int fixed_count = 0;
  for (odb::dbInst* macro : macros) {
    const bool fixed = isFixedMacro(macro);
    placed_macros_[macro] = fixed;
    if (fixed) {
      fixed_count++;
    }
  }

  logger_->info(MPL,
                1029,
                "Initialized incremental macro placement on block {} with {} macros ({} fixed).",
                block_->getName(),
                macros.size(),
                fixed_count);
}

void HierRTLMP::CriticalMacroFilter()
{
  initializeSTA();
  const std::vector<odb::dbInst*> macros = collectMacros(block_);
  int macro_pin_count = 0;
  int connected_macro_pin_count = 0;
  int negative_macro_pin_count = 0;
  int candidate_endpoint_count = 0;
  int negative_endpoint_count = 0;
  int missing_slack_count = 0;

  for (odb::dbInst* macro : macros) {
    for (odb::dbITerm* macro_pin : macro->getITerms()) {
      macro_pin_count++;
      odb::dbNet* net = macro_pin->getNet();
      if (!isSignalNet(net)) {
        continue;
      }
      connected_macro_pin_count++;

      const std::optional<float> macro_pin_slack = pinSlack(macro_pin, sta_, network_);
      if (!macro_pin_slack.has_value()) {
        missing_slack_count++;
        continue;
      }
      if (macro_pin_slack.value() >= 0.0f) {
        continue;
      }
      negative_macro_pin_count++;
      critical_macro_slack_[macro].push_back(macro_pin_slack.value());

      std::set<odb::dbInst*> endpoints;
      std::set<odb::dbITerm*> visited_trace_pins;
      if (isDriverLike(macro_pin)) {
        findDownstreamRegisters(macro, macro_pin, endpoints, visited_trace_pins);
      }
      if (isSinkLike(macro_pin)) {
        findUpstreamRegisters(macro, macro_pin, endpoints, visited_trace_pins);
      }

      for (odb::dbInst* endpoint : endpoints) {
        candidate_endpoint_count++;
        const std::optional<float> endpoint_slack = worstInstanceSlack(endpoint, sta_, network_);
        if (!endpoint_slack.has_value()) {
          missing_slack_count++;
          continue;
        }
        if (endpoint_slack.value() < 0.0f) {
          negative_endpoint_count++;
          critical_macro_to_regs_[macro].insert(endpoint);
          macro_pin_regs_[macro][macro_pin].insert(endpoint);

          const auto slack_it = critical_regs_slack_.find(endpoint);
          if (slack_it == critical_regs_slack_.end()
              || endpoint_slack.value() < slack_it->second) {
            critical_regs_slack_[endpoint] = endpoint_slack.value();
          }
        }
      }
    }
  }

  // Keep the first implementation intentionally transparent: the counters tell
  // whether a run is using real post-placement STA data or merely seeing an
  // un-timed database.  If missing_slack_count is high, the Tcl flow probably
  // forgot to read Liberty/SDC and call estimate_parasitics after placement.
  logger_->info(MPL,
                1030,
                "Critical macro filter scanned {} macro pins ({} connected signal pins, {} negative macro pins, {} sequential endpoint candidates, {} negative endpoints, {} missing slack queries).",
                macro_pin_count,
                connected_macro_pin_count,
                negative_macro_pin_count,
                candidate_endpoint_count,
                negative_endpoint_count,
                missing_slack_count);
}

void HierRTLMP::MacroClassify()
{
  array_critical_macro_insts_.clear();
  free_critical_macro_insts_.clear();

  std::vector<odb::dbInst*> critical_macros;
  critical_macros.reserve(critical_macro_to_regs_.size());
  for (const auto& [macro, regs] : critical_macro_to_regs_) {
    if (macro != nullptr && !regs.empty()) {
      critical_macros.push_back(macro);
    }
  }
  std::sort(critical_macros.begin(), critical_macros.end(), nameLess);

  std::map<odb::dbMaster*, std::vector<odb::dbInst*>> macros_by_master;
  for (odb::dbInst* macro : critical_macros) {
    macros_by_master[macro->getMaster()].push_back(macro);
  }

  std::unordered_set<odb::dbInst*> array_macros;
  for (auto& [master, macros] : macros_by_master) {
    std::sort(macros.begin(), macros.end(), nameLess);
    if (master != nullptr && hasArrayLikeAlignment(macros, block_)) {
      array_critical_macro_insts_.push_back(macros);
      array_macros.insert(macros.begin(), macros.end());
    }
  }

  for (odb::dbInst* macro : critical_macros) {
    if (array_macros.find(macro) == array_macros.end()) {
      free_critical_macro_insts_.push_back(macro);
    }
  }

  // The following movement stages consume the dbInst-level groups. The old
  // HardMacro containers are left untouched for compatibility with the original
  // global RTLMP code above the Origin marker.
  logger_->info(MPL,
                1031,
                "Classified {} critical macros into {} array groups and {} free macros.",
                critical_macros.size(),
                array_critical_macro_insts_.size(),
                free_critical_macro_insts_.size());
}

int HierRTLMP::TilingMacroExchange()
{
  int accepted_swaps = 0;
  for (const std::vector<odb::dbInst*>& group : array_critical_macro_insts_) {
    if (group.size() < 2) {
      continue;
    }

    std::vector<odb::dbInst*> candidates = group;
    std::sort(candidates.begin(), candidates.end(), [&](odb::dbInst* lhs, odb::dbInst* rhs) {
      return macroWorstSlack(lhs, critical_macro_slack_)
             < macroWorstSlack(rhs, critical_macro_slack_);
    });

    // The paper limits array disturbance to at most half of the array.  Use the
    // worst-slack half as the active set, then accept only swaps that reduce the
    // timing-weighted distance objective for the two swapped macros.
    const size_t active_count = std::max<size_t>(1, group.size() / 2);
    candidates.resize(std::min(active_count, candidates.size()));

    std::set<odb::dbInst*> already_swapped;
    for (odb::dbInst* active_macro : candidates) {
      if (already_swapped.count(active_macro) != 0 || isFixedMacro(active_macro)) {
        continue;
      }

      odb::dbInst* best_partner = nullptr;
      double best_delta = 0.0;
      const odb::Point active_origin = active_macro->getOrigin();
      const double active_cost
          = macroCostAtOrigin(active_macro,
                              active_origin.x(),
                              active_origin.y(),
                              critical_macro_to_regs_,
                              critical_regs_slack_);

      for (odb::dbInst* partner : group) {
        if (partner == active_macro || isFixedMacro(partner)
            || already_swapped.count(partner) != 0) {
          continue;
        }
        const odb::Point partner_origin = partner->getOrigin();
        const double partner_cost
            = macroCostAtOrigin(partner,
                                partner_origin.x(),
                                partner_origin.y(),
                                critical_macro_to_regs_,
                                critical_regs_slack_);
        const double swapped_cost
            = macroCostAtOrigin(active_macro,
                                partner_origin.x(),
                                partner_origin.y(),
                                critical_macro_to_regs_,
                                critical_regs_slack_)
              + macroCostAtOrigin(partner,
                                  active_origin.x(),
                                  active_origin.y(),
                                  critical_macro_to_regs_,
                                  critical_regs_slack_);
        const double original_cost = active_cost + partner_cost;
        const double delta = swapped_cost - original_cost;
        if (hasMeaningfulCostReduction(original_cost, swapped_cost, kArraySwapMinCostReductionRatio)
            && delta < best_delta) {
          best_delta = delta;
          best_partner = partner;
        }
      }

      if (best_partner != nullptr) {
        const odb::Point active_origin = active_macro->getOrigin();
        const odb::Point partner_origin = best_partner->getOrigin();
        active_macro->setLocation(partner_origin.x(), partner_origin.y());
        best_partner->setLocation(active_origin.x(), active_origin.y());
        logger_->info(MPL,
                      1037,
                      "Array swap accepted: {} <-> {} (local cost delta {:.3e}).",
                      active_macro->getName(),
                      best_partner->getName(),
                      best_delta);
        already_swapped.insert(active_macro);
        already_swapped.insert(best_partner);
        accepted_swaps++;
      }
    }
  }

  logger_->info(MPL,
                1032,
                "Accepted {} array macro swaps across {} array groups.",
                accepted_swaps,
                array_critical_macro_insts_.size());
  return accepted_swaps;
}

int HierRTLMP::FreeMacroOptimize()
{
  int moved_macros = 0;
  const odb::Rect core = block_->getCoreArea();

  for (odb::dbInst* macro : free_critical_macro_insts_) {
    if (macro == nullptr || isFixedMacro(macro)) {
      continue;
    }

    const odb::Point initial_origin = macro->getOrigin();
    const odb::Rect initial_box = macro->getBBox()->getBox();
    const auto [target_x, target_y]
        = weightedRegisterCenter(macro, critical_macro_to_regs_, critical_regs_slack_);

    double current_center_x = initial_box.xCenter();
    double current_center_y = initial_box.yCenter();
    const double initial_cost = macroCostAtCenter(macro,
                                            current_center_x,
                                            current_center_y,
                                            critical_macro_to_regs_,
                                            critical_regs_slack_);
    double best_cost = initial_cost;
    int best_x = initial_origin.x();
    int best_y = initial_origin.y();
    const int max_move_dbu = block_->micronsToDbu(kFreeMoveMaxDisplacementMicrons);

    // Projected-gradient-style refinement: take damped steps toward the weighted
    // critical-register center, project into the core, and reject candidates that
    // overlap existing placed macros.  The displacement penalty is deliberately
    // strong in this first packaged version: post-placement macro moves disturb
    // standard-cell placement, so a move must be both small and clearly useful.
    for (int iter = 0; iter < kFreeMoveMaxIterations; iter++) {
      const double regularized_target_x
          = (target_x + kFreeMoveDisplacementPenalty * initial_box.xCenter())
            / (1.0 + kFreeMoveDisplacementPenalty);
      const double regularized_target_y
          = (target_y + kFreeMoveDisplacementPenalty * initial_box.yCenter())
            / (1.0 + kFreeMoveDisplacementPenalty);
      const double step = kFreeMoveInitialStep / (1.0 + kFreeMoveStepDecay * iter);
      current_center_x += step * (regularized_target_x - current_center_x);
      current_center_y += step * (regularized_target_y - current_center_y);

      int candidate_x = static_cast<int>(std::lround(current_center_x - initial_box.dx() / 2.0));
      int candidate_y = static_cast<int>(std::lround(current_center_y - initial_box.dy() / 2.0));
      candidate_x = clampNearInitial(candidate_x, initial_origin.x(), max_move_dbu);
      candidate_y = clampNearInitial(candidate_y, initial_origin.y(), max_move_dbu);
      candidate_x = clampOriginX(macro, candidate_x, core);
      candidate_y = clampOriginY(macro, candidate_y, core);

      const odb::Rect candidate_box = rectAtOrigin(macro, candidate_x, candidate_y);
      if (!isLegalMacroBox(macro, candidate_box, block_)) {
        continue;
      }

      const double candidate_cost
          = macroCostAtOrigin(macro,
                              candidate_x,
                              candidate_y,
                              critical_macro_to_regs_,
                              critical_regs_slack_);
      if (hasMeaningfulCostReduction(initial_cost, candidate_cost, kFreeMoveMinCostReductionRatio)
          && candidate_cost < best_cost) {
        best_cost = candidate_cost;
        best_x = candidate_x;
        best_y = candidate_y;
      }
    }

    if (best_x != initial_origin.x() || best_y != initial_origin.y()) {
      logger_->info(MPL,
                    1038,
                    "Free macro move accepted: {} ({:.3f}, {:.3f}) -> ({:.3f}, {:.3f}), local cost {:.3e} -> {:.3e}.",
                    macro->getName(),
                    block_->dbuToMicrons(initial_origin.x()),
                    block_->dbuToMicrons(initial_origin.y()),
                    block_->dbuToMicrons(best_x),
                    block_->dbuToMicrons(best_y),
                    initial_cost,
                    best_cost);
      macro->setLocation(best_x, best_y);
      moved_macros++;
    }
  }

  logger_->info(MPL,
                1033,
                "Moved {} free critical macros with projected-gradient refinement.",
                moved_macros);
  return moved_macros;
}

int HierRTLMP::MacroOrientOptimizer()
{
  int changed_orients = 0;
  const std::array<odb::dbOrientType, 4> candidate_orients = {
      odb::dbOrientType::R0,
      odb::dbOrientType::MX,
      odb::dbOrientType::MY,
      odb::dbOrientType::R180,
  };

  std::vector<odb::dbInst*> macros;
  for (const auto& group : array_critical_macro_insts_) {
    macros.insert(macros.end(), group.begin(), group.end());
  }
  macros.insert(macros.end(),
                free_critical_macro_insts_.begin(),
                free_critical_macro_insts_.end());
  std::sort(macros.begin(), macros.end(), nameLess);
  macros.erase(std::unique(macros.begin(), macros.end()), macros.end());

  for (odb::dbInst* macro : macros) {
    if (macro == nullptr || isFixedMacro(macro)) {
      continue;
    }

    const odb::dbOrientType original_orient = macro->getOrient();
    const odb::Point origin = macro->getOrigin();
    int64_t best_hpwl = ComputeMacroNetHpwl(macro);
    odb::dbOrientType best_orient = original_orient;

    for (const odb::dbOrientType& orient : candidate_orients) {
      macro->setOrient(orient);
      macro->setLocation(origin.x(), origin.y());
      if (!isLegalMacroBox(macro, macro->getBBox()->getBox(), block_)) {
        continue;
      }
      const int64_t hpwl = ComputeMacroNetHpwl(macro);
      if (hpwl <= static_cast<int64_t>((1.0 - kOrientMinHpwlReductionRatio) * best_hpwl)) {
        best_hpwl = hpwl;
        best_orient = orient;
      }
    }

    macro->setOrient(best_orient);
    macro->setLocation(origin.x(), origin.y());
    if (best_orient != original_orient) {
      logger_->info(MPL,
                    1039,
                    "Macro orientation accepted: {} {} -> {} (local HPWL {}).",
                    macro->getName(),
                    original_orient.getString(),
                    best_orient.getString(),
                    best_hpwl);
      changed_orients++;
    }
  }

  logger_->info(MPL,
                1036,
                "Changed {} macro orientations by local HPWL enumeration.",
                changed_orients);
  return changed_orients;
}

void HierRTLMP::initializeSTA()
{
  if (sta_ == nullptr) {
    logger_->warn(MPL, 1034, "OpenSTA handle is not available for incremental macro placement.");
    return;
  }

  // Recompute timing for the database as loaded by the Tcl flow.  This should
  // be called only after standard-cell placement and parasitic estimation so the
  // critical macro diagnosis reflects the placed design, not floorplan guesses.
  sta_->ensureGraph();
  sta_->updateTiming(false);
}

void HierRTLMP::analyzeMacroPins(const std::string& file_path,
                                 const std::string& output_file_path)
{
  logger_->info(MPL,
                1035,
                "Macro pin timing analysis is not implemented yet (input: {}, output: {}).",
                file_path,
                output_file_path);
}

void HierRTLMP::findDownstreamRegisters(odb::dbInst* macro,
                                        odb::dbITerm* iterm,
                                        std::set<odb::dbInst*>& found_regs,
                                        std::set<odb::dbITerm*>& visited)
{
  collectDownstreamRegisters(macro, iterm, network_, found_regs, visited, 0);
}

void HierRTLMP::findUpstreamRegisters(odb::dbInst* macro,
                                      odb::dbITerm* iterm,
                                      std::set<odb::dbInst*>& found_regs,
                                      std::set<odb::dbITerm*>& visited)
{
  collectUpstreamRegisters(macro, iterm, network_, found_regs, visited, 0);
}

void HierRTLMP::delayFitting(odb::dbInst* macro_inst)
{
  if (macro_inst != nullptr) {
    macro_delay_fit[macro_inst] = {0.0f, 0.0f, 0.0f};
  }
}

std::pair<float, float> HierRTLMP::freeMacroLP(odb::dbInst* inst, int)
{
  if (inst == nullptr) {
    return {0.0f, 0.0f};
  }
  const odb::Point origin = inst->getOrigin();
  return {static_cast<float>(block_->dbuToMicrons(origin.x())),
          static_cast<float>(block_->dbuToMicrons(origin.y()))};
}

void HierRTLMP::findBestOrient(odb::dbInst*,
                               std::unordered_map<odb::dbITerm*, std::set<odb::dbInst*>>)
{
}

void HierRTLMP::findBestOrient(odb::dbInst*)
{
}

int64_t HierRTLMP::ComputeMacroNetHpwl(odb::dbInst* macro)
{
  if (macro == nullptr) {
    return 0;
  }

  int64_t hpwl = 0;
  for (odb::dbITerm* iterm : macro->getITerms()) {
    odb::dbNet* net = iterm->getNet();
    if (!isSignalNet(net)) {
      continue;
    }
    bool has_pin = false;
    int min_x = std::numeric_limits<int>::max();
    int min_y = std::numeric_limits<int>::max();
    int max_x = std::numeric_limits<int>::min();
    int max_y = std::numeric_limits<int>::min();
    for (odb::dbITerm* net_iterm : net->getITerms()) {
      int x = 0;
      int y = 0;
      if (!net_iterm->getAvgXY(&x, &y)) {
        continue;
      }
      min_x = std::min(min_x, x);
      min_y = std::min(min_y, y);
      max_x = std::max(max_x, x);
      max_y = std::max(max_y, y);
      has_pin = true;
    }
    if (has_pin) {
      hpwl += static_cast<int64_t>(max_x - min_x) + static_cast<int64_t>(max_y - min_y);
    }
  }
  return hpwl;
}

}  // namespace mpl2
