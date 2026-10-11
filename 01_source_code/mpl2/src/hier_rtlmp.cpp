///////////////////////////////////////////////////////////////////////////
//
// BSD 3-Clause License
//
// Copyright (c) 2020, The Regents of the University of California
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its
//   contributors may be used to endorse or promote products derived from
//   this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
///////////////////////////////////////////////////////////////////////////////
#include "hier_rtlmp.h"
#include "MacroPlacementSolver.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <queue>
#include <thread>
#include "Mpl2Observer.h"
#include "SACoreHardMacro.h"
#include "SACoreSoftMacro.h"
#include "bus_synthesis.h"
#include "db_sta/dbNetwork.hh"
#include "object.h"
#include "odb/db.h"
#include "odb/util.h"
#include "par/PartitionMgr.h"
#include "sta/Liberty.hh"
#include "utl/Logger.h"

namespace mpl2 {

using std::string;

static string make_filename(const string& file_name)
{
  return std::regex_replace(file_name, std::regex("[\":<>\\|\\*\\?/]"), "--");
}

///////////////////////////////////////////////////////////
// Class HierRTLMP
using utl::MPL;

HierRTLMP::~HierRTLMP() = default;

// Constructors
HierRTLMP::HierRTLMP(sta::dbNetwork* network, odb::dbDatabase* db, sta::dbSta* sta, utl::Logger* logger, par::PartitionMgr* tritonpart)
    : network_(network),
      db_(db),
      logger_(logger),
      tritonpart_(tritonpart),
      sta_(sta),
      tree_(std::make_unique<PhysicalHierarchy>())
{
}

// Top Level Function
// The flow of our MacroPlacer is divided into 6 stages.
// 1) Multilevel Autoclustering:
//      Transform logical hierarchy into physical hierarchy.
// 2) Coarse Shaping -> Bottom - Up:
//      Determine the rough shape function for each cluster.
// 3) Fine Shaping -> Top - Down:
//      Refine the possible shapes of each cluster based on the fixed
//      outline and location of its parent cluster.
//      *This is executed within hierarchical placement method.
// 4) Hierarchical Macro Placement -> Top - Down
//      a) Placement of Clusters (one level at a time);
//      b) Placement of Macros (one macro cluster at a time).
// 5) Boundary Pushing
//      Push macro clusters to the boundaries of the design if they don't
//      overlap with either bundled IOs' blockages or other macros.
// 6) Orientation Improvement
//      Attempts macro flipping to improve WR.
void HierRTLMP::run()
{
  // 获取程序开始时间
  auto start = std::chrono::high_resolution_clock::now();

  // 统计固定宏模块的数量

  int total_macros = 0;
  int fixed_macros = 0;

  // block_->getInsts() 会返回 dbSet<dbInst>，直接用范围 for 循环即可
  for (odb::dbInst* inst : block_->getInsts()) {
      odb::dbMaster* master = inst->getMaster();
      
      // isBlock() 内部只是一个极快的枚举值对比
      if (master->isBlock()) { 
          total_macros++;
          odb::dbPlacementStatus status = inst->getPlacementStatus();
          
          // 同时检查 LOCKED 和 FIRM
          if (status == odb::dbPlacementStatus::LOCKED || 
              status == odb::dbPlacementStatus::FIRM) {
              fixed_macros++;
          }
      }
  }

  if (total_macros > 0 && fixed_macros > total_macros * 0.5) {
    // Incremental placement is intentionally disabled until the new algorithm
    // is reintroduced; use the global flow as the temporary fallback.
    GlobalMacroPlacement();
  } else {
      GlobalMacroPlacement();
  }

  // bool Incre_flag = false;
  // for (odb::dbInst* inst : block_->getInsts()) {
  //     odb::dbMaster* master = inst->getMaster();
  //     if (master->isBlock()) {  // 只检查宏
  //         if (inst->getPlacementStatus() == odb::dbPlacementStatus::LOCKED){
  //           Incre_flag = true;
  //           break;
  //         }
  //     }
  // }

  // if(Incre_flag){
  //   IncrementalMacroPlacement();
  // }else{
  //   GlobalMacroPlacement();
  // }

  // 计算运行时间
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;
  std::cout << "总运行时间: " << duration.count() << " 秒" << std::endl;
  clear();
}

void HierRTLMP::GlobalMacroPlacement(){
  runMultilevelAutoclustering();

  if (skip_macro_placement_) {
    logger_->info(MPL, 13, "Skipping macro placement.");
    return;
  }
  if (!tree_->has_std_cells) {
    resetSAParameters();
  }
  std::unique_ptr<Mpl2Observer> save_graphics;
  if (is_debug_only_final_result_) {
    save_graphics = std::move(graphics_);
  }
  
  runCoarseShaping();

  // 😀👆
  runNewMacroPlacement();
  // 🐂🐎

  // runHierarchicalMacroPlacement();
  // if (save_graphics) {
  //   graphics_ = std::move(save_graphics);
  //   graphics_->setMaxLevel(tree_->max_level);
  //   graphics_->drawResult();
  // }
  // Pusher pusher(logger_, tree_->root.get(), block_, boundary_to_io_blockage_);
  // pusher.pushMacrosToCoreBoundaries();
  // updateMacrosOnDb();
  // generateTemporaryStdCellsPlacement(tree_->root.get());
  // correctAllMacrosOrientation();
  // commitMacroPlacementToDb();
  // writeMacroPlacement(macro_placement_file_);

}

void HierRTLMP::init()
{
  block_ = db_->getChip()->getBlock();
}

////////////////////////////////////////////////////////////////////////
// Private functions
////////////////////////////////////////////////////////////////////////

// 添加的修改布局方法的函数
void HierRTLMP::runNewMacroPlacement(){

  // 更新连接和数据流
  for (auto& cluster : tree_->root->getChildren()) {
    clustering_engine_->updateInstancesAssociation(cluster.get());
  }
  clustering_engine_->updateConnections();
  clustering_engine_->updateDataFlow();

  // 对于顶层root,随机生成宽长比,返回id-宽长比对的map
  shapeGenerate();

  // 写聚类结果

  // 改成绝对路径

  
  stage_path_ = "/workspace/OpenROAD-flow-scripts/tools/OpenROAD/src/mpl2/temp/bloat";
  WriteClusterResult();

  // 调用replace进行扩散
  runReplaceTcl();

  // 将结果反馈回顶层簇的坐标
  rebackRePlaceClusterLocation();

  stage_path_ = "/workspace/OpenROAD-flow-scripts/tools/OpenROAD/src/mpl2/temp/first_lp";

  // 线性规划求解宏坐标
  InitializeMacroLocation();

  // 宏合法化
  stage_path_ = "/workspace/OpenROAD-flow-scripts/tools/OpenROAD/src/mpl2/temp/macro_legalize";
  macroLegalize(false);

  // 生成伪网络约束
  // clusterFakeNetStdGenerate();

  // 生成初始标准单元坐标，放在顶层簇的中心
  stdCellPlacement();
}

void HierRTLMP::macroLegalize(bool incre_flag)
{
  std::vector<Rect> space_placed;
  std::vector<Rect> fixed_macros;
  std::unordered_map<int, odb::dbInst*> order_inst_map;
  int order = 0;

  if (!incre_flag) {
    for (auto& [inst, hardmacro] : tree_->maps.inst_to_hard) {
      float x = hardmacro->getX();
      float y = hardmacro->getY();
      float width = hardmacro->getWidth();
      float height = hardmacro->getHeight();
      space_placed.emplace_back(x, y, x + width, y + height, false);
      order_inst_map[order] = inst;
      order++;
    }
  } else {
    for (auto hardmacro : free_critical_macros_) {
      odb::dbInst* inst = hardmacro->getInst();
      float x = hardmacro->getX();
      float y = hardmacro->getY();
      float width = hardmacro->getWidth();
      float height = hardmacro->getHeight();
      space_placed.emplace_back(x, y, x + width, y + height, false);
      order_inst_map[order] = inst;
      order++;
    }

    for (auto& [inst, fixed_macro] : tree_->maps.inst_to_hard) {
      auto it = std::find_if(free_critical_macros_.begin(),
                             free_critical_macros_.end(),
                             [&inst](HardMacro* ptr) {
                               return ptr->getInst() == inst;
                             });
      if (it == free_critical_macros_.end()) {
        float x = fixed_macro->getX();
        float y = fixed_macro->getY();
        float width = fixed_macro->getWidth();
        float height = fixed_macro->getHeight();
        fixed_macros.emplace_back(x, y, x + width, y + height, false);
      }
    }
  }

  const odb::Rect core_box = block_->getCoreArea();
  MacroPlacementSolver solver;
  SolverCore core{block_->dbuToMicrons(core_box.xMin()),
                  block_->dbuToMicrons(core_box.yMin()),
                  block_->dbuToMicrons(core_box.xMax()),
                  block_->dbuToMicrons(core_box.yMax())};
  std::vector<SolverRect> movable;
  for (int i = 0; i < static_cast<int>(space_placed.size()); ++i) {
    const Rect& rect = space_placed[i];
    movable.push_back({i, rect.xMin(), rect.yMin(), rect.xMax(), rect.yMax()});
  }
  std::vector<SolverRect> fixed;
  for (int i = 0; i < static_cast<int>(fixed_macros.size()); ++i) {
    const Rect& rect = fixed_macros[i];
    fixed.push_back({i, rect.xMin(), rect.yMin(), rect.xMax(), rect.yMax()});
  }
  const auto placements = solver.legalizeMacros(core, movable, fixed);
  if (placements.size() != movable.size()) {
    logger_->warn(MPL, 1020, "Macro legalization failed to find a feasible non-overlap solution; preserving current macro locations.");
    return;
  }
  for (const SolverPlacement& placement : placements) {
    odb::dbInst* inst = order_inst_map[placement.id];
    auto hardmacro = tree_->maps.inst_to_hard.find(inst)->second.get();
    float new_x = placement.x;
    float new_y = placement.y;
    const float core_lx = static_cast<float>(block_->dbuToMicrons(core_box.xMin()));
    const float core_ux = static_cast<float>(block_->dbuToMicrons(core_box.xMax()));
    const float core_ly = static_cast<float>(block_->dbuToMicrons(core_box.yMin()));
    const float core_uy = static_cast<float>(block_->dbuToMicrons(core_box.yMax()));
    new_x = std::max(new_x, core_lx);
    new_x = std::min(new_x, core_ux - hardmacro->getRealWidth());
    new_y = std::max(new_y, core_ly);
    new_y = std::min(new_y, core_uy - hardmacro->getRealHeight());
    inst->setPlacementStatus(odb::dbPlacementStatus::PLACED);
    hardmacro->setLocation(std::make_pair(new_x, new_y));
    inst->setLocation(block_->micronsToDbu(new_x), block_->micronsToDbu(new_y));
  }

  for (auto& [inst, hard_macro] : tree_->maps.inst_to_hard) {
    auto loc = inst->getLocation();
    float x = loc.getX();
    float y = loc.getY();
    float width = hard_macro->getRealWidthDBU();
    float height = hard_macro->getRealHeightDBU();
    float new_x = x;
    float new_y = y;
    if (x < core_box.xMin()) {
      new_x = core_box.xMin();
    } else if (x + width > core_box.xMax()) {
      new_x = core_box.xMax() - width;
    }
    if (y < core_box.yMin()) {
      new_y = core_box.yMin();
    } else if (y + height > core_box.yMax()) {
      new_y = core_box.yMax() - height;
    }
    inst->setLocation(new_x, new_y);
  }
}

void HierRTLMP::InitializeMacroLocation(){
  // 求解一个簇之间宏的连接的线网
  auto macro_dataflow_nets = createMacroDataflowNet(tree_->root.get());
  stage_path_ = "/workspace/OpenROAD-flow-scripts/tools/OpenROAD/src/mpl2/temp/first_lp";
  // 一级线性规划，求解簇内宏的相对位置关系
  for(auto& cluster : tree_->root->getChildren()){
    if(cluster->isIOCluster() || cluster->getClusterType() == StdCellCluster){
      continue;
    }
    int num_macro = cluster->getMetrics().getNumMacro(); 
    if(num_macro<2){
      for(auto& hardmacro : cluster->getHardMacros()){
        odb::dbInst* inst = hardmacro->getInst();
        float cluster_center_x = cluster->getX() + id_shape_map_[cluster->getId()].first/2;
        float cluster_center_y = cluster->getY() + id_shape_map_[cluster->getId()].second/2;

        float width = hardmacro->getWidth();
        float height = hardmacro->getHeight();
        float x = cluster_center_x - width/2;
        float y = cluster_center_y - height/2;
        inst->setOrient(hardmacro->getOrientation());
        inst->setLocation(block_->micronsToDbu(x), block_->micronsToDbu(y));
        inst->setPlacementStatus(odb::dbPlacementStatus::PLACED);

        hardmacro->setLocation(std::make_pair(x,y));
      }
      continue;
    }
    
    FirstLinearSolver(cluster.get(), macro_dataflow_nets);
  }

  // 二级线性规划(将簇内的宏推向边缘)适用于所有的宏簇，包括单宏簇
  stage_path_ = "/workspace/OpenROAD-flow-scripts/tools/OpenROAD/src/mpl2/temp/second_lp";

  SecondLinearSolver();


  // 计算每个簇的宏的面积，按照降序排序
  // 自定义比较函数，使得multiset是按照降序排列
  // auto cmp = [](const std::pair<float, int>& a, const std::pair<float, int>& b) {
  //     return a.first > b.first; // 降序
  // };
  // std::multiset<std::pair<float, int>, decltype(cmp)> macro_area_cmp(cmp);

  // // 已经放好的宏组，生成fence region
  // std::vector<std::vector<float>> fence_region;

  // for(auto& cluster : tree_->root->getChildren()){
  //   if(cluster->isIOCluster() || cluster->getClusterType() == StdCellCluster){
  //     continue;
  //   }
  //   float sum_area = 0.0f;
  //   for(auto hardmacro : cluster->getHardMacros()){
  //     sum_area += hardmacro->getWidth() * hardmacro->getHeight();
  //   }
  //   int id = cluster->getId();
  //   macro_area_cmp.insert({sum_area,id});
  // }

  // for (auto cmp : macro_area_cmp){
  //   auto current = tree_->maps.id_to_cluster.find(cmp.second);
  //   Cluster* cluster_ptr = current->second;
  //   SecondLinearSolver(cluster_ptr, fence_region);
  // }

  // for(auto& cluster : tree_->root->getChildren()){
  //   if(cluster->isIOCluster() || cluster->getClusterType() == StdCellCluster){
  //     continue;
  //   }
  //   SecondLinearSolver(cluster.get());
  // }

}

void HierRTLMP::SecondLinearSolver(){
  std::unordered_map<int, std::vector<float>> macro_big_rects;
  std::unordered_map<int, std::pair<float, float>> big_rects_shape;
  for (auto& cluster : tree_->root->getChildren()) {
    if (cluster->isIOCluster() || cluster->getClusterType() == StdCellCluster) {
      continue;
    }

    std::vector<float> current_big_rect;
    auto frontmacro = cluster->getHardMacros()[0];
    current_big_rect.push_back(frontmacro->getRealX());
    current_big_rect.push_back(frontmacro->getRealY());
    current_big_rect.push_back(frontmacro->getRealX() + frontmacro->getRealWidth());
    current_big_rect.push_back(frontmacro->getRealY() + frontmacro->getRealHeight());
    for (auto& hardmacro : cluster->getHardMacros()) {
      float x = hardmacro->getRealX();
      float y = hardmacro->getRealY();
      float width = hardmacro->getRealWidth();
      float height = hardmacro->getRealHeight();
      current_big_rect[0] = std::min(current_big_rect[0], x);
      current_big_rect[1] = std::min(current_big_rect[1], y);
      current_big_rect[2] = std::max(current_big_rect[2], x + width);
      current_big_rect[3] = std::max(current_big_rect[3], y + height);
    }
    macro_big_rects[cluster->getId()] = current_big_rect;
    big_rects_shape[cluster->getId()]
        = {current_big_rect[2] - current_big_rect[0], current_big_rect[3] - current_big_rect[1]};
  }

  const odb::Rect core_box = block_->getCoreArea();
  std::unordered_map<int, std::pair<float, float>> big_rects_loc;
  for (auto& cluster : tree_->root->getChildren()) {
    if (cluster->isIOCluster() || cluster->getClusterType() == StdCellCluster) {
      continue;
    }
    int id = cluster->getId();
    float cluster_center_x = cluster->getX() + id_shape_map_[id].first / 2;
    float cluster_center_y = cluster->getY() + id_shape_map_[id].second / 2;
    float rect_x = cluster_center_x - big_rects_shape[id].first / 2;
    float rect_y = cluster_center_y - big_rects_shape[id].second / 2;
    if (rect_x < block_->dbuToMicrons(core_box.xMin())) {
      rect_x = block_->dbuToMicrons(core_box.xMin());
    } else if (rect_x + big_rects_shape[id].first > block_->dbuToMicrons(core_box.xMax())) {
      rect_x = block_->dbuToMicrons(core_box.xMax()) - big_rects_shape[id].first;
    }
    if (rect_y < block_->dbuToMicrons(core_box.yMin())) {
      rect_y = block_->dbuToMicrons(core_box.yMin());
    } else if (rect_y + big_rects_shape[id].second > block_->dbuToMicrons(core_box.yMax())) {
      rect_y = block_->dbuToMicrons(core_box.yMax()) - big_rects_shape[id].second;
    }
    big_rects_loc[id] = {rect_x, rect_y};
  }

  std::vector<SolverMacro> solver_groups;
  std::unordered_map<int, SolverPlacement> solver_locations;
  for (auto& [id, shape] : big_rects_shape) {
    solver_groups.push_back({id, shape.first, shape.second});
    solver_locations[id] = {id, big_rects_loc[id].first, big_rects_loc[id].second};
  }
  SolverCore core{block_->dbuToMicrons(core_box.xMin()),
                  block_->dbuToMicrons(core_box.yMin()),
                  block_->dbuToMicrons(core_box.xMax()),
                  block_->dbuToMicrons(core_box.yMax())};
  MacroPlacementSolver solver;
  const auto placements = solver.moveMacroGroupsToBoundary(
      core, solver_groups, solver_locations, block_->getName().c_str());

  std::unordered_map<int, std::pair<float, float>> macro_move;
  for (const SolverPlacement& placement : placements) {
    macro_move[placement.id] = {placement.x - macro_big_rects[placement.id][0],
                                placement.y - macro_big_rects[placement.id][1]};
  }

  std::unordered_map<int, std::pair<float, float>> boundary_move;
  std::unordered_map<odb::dbInst*, std::pair<float, float>> macro_loc;
  for (auto& cluster : tree_->root->getChildren()) {
    if (cluster->isIOCluster() || cluster->getClusterType() == StdCellCluster) {
      continue;
    }
    int id = cluster->getId();
    boundary_move[id] = {0, 0};
    for (auto& hardmacro : cluster->getHardMacros()) {
      float x = hardmacro->getRealX() + macro_move[id].first;
      float y = hardmacro->getRealY() + macro_move[id].second;
      odb::dbInst* inst = hardmacro->getInst();
      macro_loc[inst] = {x, y};
      if (x < block_->dbuToMicrons(core_box.xMin())) {
        boundary_move[id].first++;
      } else if (x + hardmacro->getRealWidth() > block_->dbuToMicrons(core_box.xMax())) {
        boundary_move[id].first--;
      }
      if (y < block_->dbuToMicrons(core_box.yMin())) {
        boundary_move[id].second++;
      } else if (y + hardmacro->getRealHeight() > block_->dbuToMicrons(core_box.yMax())) {
        boundary_move[id].second--;
      }
    }

    for (auto& hardmacro : cluster->getHardMacros()) {
      odb::dbInst* inst = hardmacro->getInst();
      float x = macro_loc[inst].first + boundary_move[id].first;
      float y = macro_loc[inst].second + boundary_move[id].second;
      hardmacro->setRealLocation(std::make_pair(x, y));
      inst->setLocation(block_->micronsToDbu(x), block_->micronsToDbu(y));
      inst->setPlacementStatus(odb::dbPlacementStatus::LOCKED);
    }
  }
}

void HierRTLMP::SecondLinearSolver(Cluster* cluster, std::vector<std::vector<float>>& fence_region)
{
  // This overload was already disabled by the caller. Keep it as a no-op so
  // mpl2 no longer shells out to the old second_lp executable if it is called
  // accidentally before the model is reintroduced in C++.
  (void) cluster;
  (void) fence_region;
}

void HierRTLMP::FirstLinearSolver(Cluster* cluster, std::vector<MacroDataflowNet> macro_dataflow_nets){
  std::pair<float, float> final_tiling(0, 0);
  for (auto& tiling : cluster->getMacroTilings()) {
    final_tiling.first = std::max(final_tiling.first, tiling.first);
    final_tiling.second = std::max(final_tiling.second, tiling.second);
  }
  final_tiling.first = std::max(final_tiling.first, id_shape_map_[cluster->getId()].first);
  final_tiling.second = std::max(final_tiling.second, id_shape_map_[cluster->getId()].second);

  std::pair<float, float> min_width_height(0.0f, 0.0f);
  for (auto hardmacro : cluster->getHardMacros()) {
    min_width_height.first += hardmacro->getWidth();
    min_width_height.second += hardmacro->getHeight();
  }
  final_tiling.first = std::max(final_tiling.first, min_width_height.first);
  final_tiling.second = std::max(final_tiling.second, min_width_height.second);
  final_tiling.first *= 2;
  final_tiling.second *= 2;

  const odb::Rect core_box = block_->getCoreArea();
  final_tiling.first = std::min(final_tiling.first, static_cast<float>(block_->dbuToMicrons(core_box.xMax() - core_box.xMin())));
  final_tiling.second = std::min(final_tiling.second, static_cast<float>(block_->dbuToMicrons(core_box.yMax() - core_box.yMin())));

  auto hardmacros = cluster->getHardMacros();
  auto first_master = hardmacros[0]->getInst()->getMaster();
  bool same_type = true;
  for (auto hardmacro : hardmacros) {
    if (hardmacro->getInst()->getMaster() != first_master) {
      same_type = false;
      break;
    }
  }
  if (same_type) {
    auto num_tiling = SameMacroBestTiling(cluster);
    final_tiling.first = (num_tiling.first + 0.5) * hardmacros[0]->getWidth();
    final_tiling.second = (num_tiling.second + 0.5) * hardmacros[0]->getHeight();
  }

  std::unordered_map<int, odb::dbInst*> order_inst;
  std::vector<SolverMacro> solver_macros;
  int order = 0;
  for (auto hardmacro : cluster->getHardMacros()) {
    solver_macros.push_back({order, hardmacro->getWidth(), hardmacro->getHeight()});
    order_inst.emplace(order, hardmacro->getInst());
    order++;
  }

  std::vector<HardMacro*> children_macro = cluster->getHardMacros();
  std::vector<std::vector<float>> id_weight;
  for (auto& net : macro_dataflow_nets) {
    int driver_id = 0;
    int load_id = 0;
    odb::dbInst* driver_inst = net.driver_inst;
    odb::dbInst* load_inst = net.load_inst;
    auto driver_macro = tree_->maps.inst_to_hard.find(driver_inst)->second.get();
    if (std::find(children_macro.begin(), children_macro.end(), driver_macro) == children_macro.end()) {
      continue;
    }
    auto load_macro = tree_->maps.inst_to_hard.find(load_inst)->second.get();
    if (std::find(children_macro.begin(), children_macro.end(), load_macro) == children_macro.end()) {
      continue;
    }
    for (auto& [macro_order, inst] : order_inst) {
      if (inst == driver_inst) {
        driver_id = macro_order;
      }
      if (inst == load_inst) {
        load_id = macro_order;
      }
    }
    id_weight.push_back({static_cast<float>(driver_id), static_cast<float>(load_id), net.weight});
  }

  std::vector<std::vector<float>> rows_to_remove;
  for (const auto& current : id_weight) {
    if (current[0] > current[1]) {
      rows_to_remove.push_back(current);
    }
  }
  std::vector<std::vector<float>> id_increse_weight = id_weight;
  for (const auto& row : rows_to_remove) {
    id_increse_weight.erase(std::remove_if(id_increse_weight.begin(),
                                           id_increse_weight.end(),
                                           [&row](const std::vector<float>& r) {
                                             return r.size() >= 2 && r[0] == row[0] && r[1] == row[1];
                                           }),
                            id_increse_weight.end());
  }
  for (const auto& row : rows_to_remove) {
    id_increse_weight.erase(std::remove_if(id_increse_weight.begin(),
                                           id_increse_weight.end(),
                                           [&row](const std::vector<float>& r) {
                                             return r.size() >= 2 && r[0] == row[0] && r[1] == row[1];
                                           }),
                            id_increse_weight.end());
    for (auto& r : id_increse_weight) {
      if (r.size() >= 3 && r[0] == row[1] && r[1] == row[0]) {
        r[2] += row[2];
      }
    }
  }

  std::vector<SolverWeightedNet> solver_nets;
  for (auto current : id_increse_weight) {
    solver_nets.push_back({static_cast<int>(current[0]), static_cast<int>(current[1]), current[2]});
  }

  const bool use_long_time_limit = block_->getName() == "bp_be_top"
                                   || block_->getName() == "bp_fe_top"
                                   || block_->getName() == "bp_multi_top";
  MacroPlacementSolver solver;
  const auto placements = solver.placeMacrosInCluster(
      final_tiling.first, final_tiling.second, solver_macros, solver_nets, use_long_time_limit);
  if (placements.size() != solver_macros.size()) {
    logger_->warn(MPL, 1024, "Cluster macro placement for cluster {} returned {} placements for {} macros; preserving existing locations.",
                  cluster->getId(), placements.size(), solver_macros.size());
    return;
  }
  for (const SolverPlacement& placement : placements) {
    odb::dbInst* inst = order_inst[placement.id];
    auto hardmacro = tree_->maps.inst_to_hard.find(inst)->second.get();
    float true_x = placement.x + cluster->getX();
    float true_y = placement.y + cluster->getY();
    hardmacro->setLocation(std::make_pair(true_x, true_y));
  }
}

std::pair<int,int> HierRTLMP::SameMacroBestTiling(Cluster* cluster){
  // 判断是否为质数的 Lambda
  auto isPrime = [](int n) -> bool {
      if (n < 2) return false;
      for (int i = 2; i * i <= n; i++) {
          if (n % i == 0) return false;
      }
      return true;
  };

  // 初步
  // float best_ratio = 1.0;

  // const odb::Rect core_box = block_->getCoreArea();
  // float best_ratio = (core_box.xMax()-core_box.xMin())/(core_box.yMax()-core_box.yMin());

  auto hardmacros = cluster->getHardMacros();
  int num_macro = hardmacros.size();

  if(isPrime(num_macro)){
    num_macro++;
  }
  float width = hardmacros[0]->getWidth();
  float height = hardmacros[0]->getHeight();

  int best_n = 1, best_m = num_macro;
  float min_ratio_diff = std::max(((best_n * width) / (best_m * height)),((best_m * height) / (best_n * width)));

  for(int n=1;n<sqrt(num_macro)+1;n++){
    if(num_macro % n == 0){
      int m = num_macro/n;
      float temp_ratio_diff = std::max(((n * width) / (m * height)),((m * height) / (n * width)));
      if(temp_ratio_diff < min_ratio_diff){
        min_ratio_diff = temp_ratio_diff;
        best_n = n;
        best_m = m;
      }

      // 交换测试
      int temp_n = n;
      int temp_m = m;
      int temp = temp_m;
      temp_m = temp_n;
      temp_n = temp;

      float reverse_temp_ratio_diff = std::max(((temp_n * width) / (temp_m * height)),((temp_m * height) / (temp_n * width)));
      if(reverse_temp_ratio_diff < min_ratio_diff){
        min_ratio_diff = reverse_temp_ratio_diff;
        best_n = temp_n;
        best_m = temp_m;
      }
    }
  }
  return {best_n,best_m};
}

std::vector<MacroDataflowNet> HierRTLMP::createMacroDataflowNet(Cluster* parent){
  // std::unordered_set<HardMacro*> hard_macros(parent->getHardMacros().begin(),parent->getHardMacros().end());
  
  std::unordered_map<std::pair<odb::dbInst*,odb::dbInst*>,MacroDataflowNet, pairHash2MacroNet> macro_dataflow_net_map;

  // std::vector<MacroDataflowNet> macro_dataflow_nets;

  for(auto& [driver_pin, load_inst_vectors]: clustering_engine_->getDataFlow().macro_pins_and_macros){    
    odb::dbInst* driver_inst = driver_pin->getInst();
    auto driver_it = tree_->maps.inst_to_hard.find(driver_inst);
    HardMacro* driver_macro = driver_it->second.get();
    // 如果parent中没有driver_macro，跳过
    if(std::find(parent->getHardMacros().begin(), parent->getHardMacros().end(), driver_macro) == parent->getHardMacros().end()){
      continue;
    }

    // 如果parent中没有driver_macro，跳过
    // if(hard_macros.find(driver_macro) == hard_macros.end()){
    //   continue;
    // }

    int hop = 0;
    for(auto& load_insts : load_inst_vectors){
      for(auto& load_inst : load_insts){  
        // 如果parent中没有load_macro，跳过
        auto load_it = tree_->maps.inst_to_hard.find(load_inst);
        if(load_it == tree_->maps.inst_to_hard.end()){
          continue;
        }
        HardMacro* load_macro = load_it->second.get();

        if(std::find(parent->getHardMacros().begin(), parent->getHardMacros().end(), load_macro) == parent->getHardMacros().end()){
          continue;
        }

        // if(hard_macros.find(load_macro) == hard_macros.end()){
        //   continue;
        // }

        // 使用哈希表查找是否已经存在对应的 MacroDataflowNet
        auto key = std::make_pair(driver_inst, load_inst);
        auto it = macro_dataflow_net_map.find(key);

        if (it != macro_dataflow_net_map.end()) {
          // 修改其权重
          it->second.weight += std::pow(0.5, hop);
        } else {
          // 新建一个 net
          float weight = std::pow(0.5, hop);
          macro_dataflow_net_map[key] = MacroDataflowNet(weight, driver_inst, load_inst);
        }
        // // 判断macro_dataflow_nets中是否有了driver_inst和load_inst对应的net
        // auto it = std::find_if(macro_dataflow_nets.begin(), macro_dataflow_nets.end(),
        //                        [driver_inst, load_inst](const MacroDataflowNet& net) {
        //                          return net.driver_inst == driver_inst && net.load_inst == load_inst;
        //                        });
        // if(it != macro_dataflow_nets.end()){          
        //   // 得到这个MacroDataflowNet，修改其权重
        //   it->weight += std::pow(0.5,hop);
        // }else{
        //   // 如果macro_dataflow_nets中没有driver_inst和load_inst对应的net，新建一个
        //   float weight = std::pow(0.5,hop);
        //   macro_dataflow_nets.emplace_back(weight,driver_inst,load_inst);
        //   // MacroDataflowNet net(weight,driver_inst,load_inst);
        // }
      }
      hop++;
    }
  }

  // 将哈希表中的值转换为 vector
  std::vector<MacroDataflowNet> macro_dataflow_nets;
  for (const auto& [key, net] : macro_dataflow_net_map) {
    macro_dataflow_nets.push_back(net);
  }
  // 消除driver_inst和load_inst相同的net
  macro_dataflow_nets.erase(std::remove_if(macro_dataflow_nets.begin(), macro_dataflow_nets.end(),
                                           [](const MacroDataflowNet& net) {
                                             return net.driver_inst == net.load_inst;
                                           }),
                           macro_dataflow_nets.end());
  return macro_dataflow_nets;

}

void HierRTLMP::removeHardMacroOverlap(){
  for(auto& cluster : tree_->root->getChildren()){

    if(cluster->isIOCluster() || cluster->getClusterType() == StdCellCluster){
      continue;
    }

    auto macros = cluster->getHardMacros();
    int size = macros.size();
    for (int i = 0; i < size; ++i) {
      for (int j = i + 1; j < size; ++j) {

        HardMacro* macroi = macros[i];
        HardMacro* macroj = macros[j];
        
        // 判断二者相对位置关系
        bool i_left = true;
        bool i_bottom = true;

        if(macroi->getX() >= macroj->getX()){
          i_left = false;
        }
        if(macroi->getY() >= macroj->getY()){
          i_bottom = false;
        }

        if(i_left){
          if(macroi->getRealX() + macroi->getRealWidth() >= macroj->getRealX()){
            float new_x = macroi->getRealX() +1;
            macroi->setRealX(new_x);
          }
        }else{
          if(macroj->getRealX() + macroj->getRealWidth() >= macroi->getRealX()){
            float new_x = macroj->getRealX() +1;
            macroi->setRealX(new_x);
          }

        }

        if(i_bottom){

        }else{

        }

      }
    }
  }


  std::string path = stage_path_;
  // 计算每个cluster的bounding box
  std::unordered_map<int,Rect> id_box_map;
  for(auto entry : id_shape_map_){
    int id = entry.first;
    float width = entry.second.first;
    float height = entry.second.second;
    auto current = tree_->maps.id_to_cluster.find(id);
    Cluster* cluster_ptr = current->second;
    id_box_map.emplace(id,Rect(cluster_ptr->getX(),cluster_ptr->getX(),cluster_ptr->getX()+width,cluster_ptr->getX()+height));
  }

  // 快速判断是否重叠的情况
  auto is_overlap = [](Rect& a,Rect& b){
    if(a.xMin()>=b.xMax() || a.xMax()<=b.xMin() || a.yMin()>=b.yMax() || a.yMax()<=b.yMin()){
      return false;
    }
    return true;
  };

  // 计算重叠面积和区域
  auto calculate_overlap = [](Rect& a,Rect& b){
    float x_overlap = std::max(0.0f,std::min(a.xMax(),b.xMax())-std::max(a.xMin(),b.xMin()));
    float y_overlap = std::max(0.0f,std::min(a.yMax(),b.yMax())-std::max(a.yMin(),b.yMin()));
    float overlap_area = x_overlap*y_overlap;
    // 以下两种声明方式不知道哪个对
    // Rect overlap_coods = {std::max(a.xMin(),b.xMin()),std::max(a.yMin(),b.yMin()),std::min(a.xMax(),b.xMax()),std::min(a.yMax(),b.yMax())};
    Rect overlap_coods(std::max(a.xMin(),b.xMin()),std::max(a.yMin(),b.yMin()),std::min(a.xMax(),b.xMax()),std::min(a.yMax(),b.yMax()));
    
    return std::make_pair(overlap_area,overlap_coods);
  };

  // 找宏簇和混合簇的重叠区域
  // 遍历得到重叠的cluster以及重叠区域每加入一个cluster，就与之前的比对重叠区域
  std::vector<int> before_cluster_id;
  std::vector<std::pair<float,Rect>> overlap_area;
  std::vector<std::pair<std::pair<int,int>,std::vector<std::pair<float,Rect>>>> all_overlap;
  for(auto& current_id_box : id_box_map){
    if(before_cluster_id.empty()){
      before_cluster_id.push_back(current_id_box.first);
      continue;
    }
    for(auto& before_id : before_cluster_id){
      if(is_overlap(current_id_box.second,id_box_map[before_id])){
        overlap_area.push_back(calculate_overlap(current_id_box.second,id_box_map[before_id]));
        all_overlap.push_back(std::make_pair(std::make_pair(current_id_box.first,before_id),overlap_area));
      }
    }
    before_cluster_id.push_back(current_id_box.first);
  }

  // 重叠区域的处理，根据重叠区域的大小，调整cluster的shape
  // 标准单元簇的shape不变，宏簇和混合簇的shape调整，遍历overlap中id.first就行

  for(auto& overlap : all_overlap){
    int current_id = overlap.first.first;
    auto current = tree_->maps.id_to_cluster.find(current_id);
    Cluster* current_cluster_ptr = current->second;
    // 不调整标准单元簇
    if(current_cluster_ptr->getClusterType() == StdCellCluster){
      continue;
    }
    // 调整宏簇和混合簇
    int before_id = overlap.first.second;
    // auto before = tree_->maps.id_to_cluster.find(before_id);
    // Cluster* before_cluster_ptr = before->second;

    // 计算重叠区域面积
    float current_area = id_shape_map_[current_id].first * id_shape_map_[current_id].second;
    float before_area = id_shape_map_[before_id].first * id_shape_map_[before_id].second;
    float overlap_area = overlap.second[0].first;
    float current_overlap_ratio = overlap_area/current_area;
    float before_overlap_ratio = overlap_area/before_area;

    // 重叠区域占比大于0.5，调整shape
    if(current_overlap_ratio>0.5){
      float changed_width = id_shape_map_[current_id].first * sqrt(1/current_overlap_ratio);
      float changed_height = id_shape_map_[current_id].second * sqrt(current_overlap_ratio);
      id_shape_map_[current_id].first = changed_width;
      id_shape_map_[current_id].second = changed_height;
    }else if(before_overlap_ratio>0.5){
      auto before = tree_->maps.id_to_cluster.find(before_id);
      Cluster* before_cluster_ptr = before->second;

      if(before_cluster_ptr->getClusterType() == StdCellCluster){
        continue;
      }
      // 修改shape先修改tiling
      id_shape_map_[before_id].first = id_shape_map_[current_id].first * sqrt(1/current_overlap_ratio);
      id_shape_map_[before_id].second = id_shape_map_[current_id].second * sqrt(current_overlap_ratio);
    }
  }
  // quadtrees::QuadTree<Rect> quadtree;

}

////////////////////////////////////////////////////////////////////////
// Shape Generate
////////////////////////////////////////////////////////////////////////
void HierRTLMP::shapeGenerate(){
  // std::unordered_map<int,std::pair<float,float>> id_shape_map;  // id-宽长比的map 
  // std::random_device rd;
  // std::mt19937 gen(rd());

  const odb::Rect core_box = block_->getCoreArea();
  const float core_width = block_->dbuToMicrons(core_box.xMax()) - block_->dbuToMicrons(core_box.xMin());
  const float core_height = block_->dbuToMicrons(core_box.yMax()) - block_->dbuToMicrons(core_box.yMin());
   
  // for (auto& cluster : tree_->root.get()->getChildren()) {
  //   if(cluster->isIOCluster() == true){
  //     continue;
  //   }
  //   shape.first = core_width+1;
  //   shape.second = core_height+1;
  //   int id = cluster->getId();
  //   if(cluster->getClusterType() == MixedCluster){
  //     while(shape.first>core_width || shape.second>core_height){
  //       // 从macro_tilings中随机挑选
  //       std::pair<float, float> tiling_shape = setRandomShapes(*cluster);
  //       // 加入标准单元面积
  //       float std_area = cluster->getMetrics().getStdCellArea();
  //       float mixed_area = std_area + tiling_shape.first * tiling_shape.second;
  //       // 按文章来说， macro_tilings超过宽长比不会超过0.33-3的界限
  //       // 计算当前tiling情况下可能的的宽长比范围
  //       std::pair<float,float> width_range(tiling_shape.first,mixed_area/tiling_shape.first);
  //       // std::pair<float,float> height_range(tiling_shape.second,mixed_area/tiling_shape.second);
  //       std::pair<float,float> height_range(std::min(tiling_shape.second,mixed_area/tiling_shape.second),std::max(tiling_shape.second,mixed_area/tiling_shape.second));
  //       // 宽长比的界限,ratio = width/height
  //       float val1 = 0.33f, val2 = 3.0f;
  //       // 最小宽长比：max(0.33,最小宽/最大高),最大宽长比：min(3,最大宽/最小高)
  //       std::pair<float,float> ratio_range(std::max(val1,width_range.first/height_range.second),std::min(val2,width_range.second/height_range.first));
  //       // 随机生成宽长比
  //       std::uniform_real_distribution<float> dis(ratio_range.first, ratio_range.second);
  //       float random_ratio = dis(gen);
  //       // 宽/高=ratio
  //       shape.first = sqrt(mixed_area*random_ratio);
  //       shape.second = sqrt(mixed_area/random_ratio);
  //       // std::pair<float,float> shape(sqrt(mixed_area*random_ratio),sqrt(mixed_area/random_ratio));
  //     }
  //     // id_shape_map.emplace(id,shape); 
  //   } else if(cluster->getClusterType() == StdCellCluster) {
  //     while(shape.first>core_width || shape.second>core_height){
  //       // 宽长比随机生成0.33-3
  //       float std_area = cluster->getMetrics().getStdCellArea();
  //       std::uniform_real_distribution<float> dis(0.33f, 3.0f);
  //       float random_ratio = dis(gen);
  //       // 宽长比随机会导致replace有时候无法收敛，因此暂且全部标准单元簇设置为正方形，这样确保收敛
  //       random_ratio = 1.0f;
  //       // width = sqrt(area*ratio), height = sqrt(area/ratio)
  //       shape.first = sqrt(std_area*random_ratio);
  //       shape.second = sqrt(std_area/random_ratio);
  //       // std::pair<float,float> shape(sqrt(std_area*random_ratio), sqrt(std_area/random_ratio));        
  //     }
  //     // id_shape_map.emplace(id,shape);
  //   } else if(cluster->getClusterType() == HardMacroCluster) {
  //     while(shape.first>core_width || shape.second>core_height){
  //       std::pair<float, float> shape_temp = setRandomShapes(*cluster);
  //       shape.first = shape_temp.first;
  //       shape.second = shape_temp.second;
  //     }
  //     // id_shape_map.emplace(id,shape);
  //   }
  //   id_shape_map_.emplace(id,shape);
  // }
  
  std::pair<float,float> shape;
  for (auto& cluster : tree_->root->getChildren()){
    if(cluster->isIOCluster()){
      continue;
    }
    int id = cluster->getId();
    float area = 0.0f;
    if(cluster->getClusterType() == MixedCluster){
      area += cluster->getStdCellArea();
      for(auto hardmacro : cluster->getHardMacros()){
        area += hardmacro->getWidth() * hardmacro->getHeight();
      }
    } else if(cluster->getClusterType() == StdCellCluster){
      area += cluster->getStdCellArea();
    } else if(cluster->getClusterType() == HardMacroCluster) {
      for(auto hardmacro : cluster->getHardMacros()){
        area += hardmacro->getWidth() * hardmacro->getHeight();
      }
    }
    const float core_aspect_ratio = core_height > 0.0f ? core_width / core_height : 1.0f;
    shape.first = sqrt(area * core_aspect_ratio);
    shape.second = sqrt(area / core_aspect_ratio);
    id_shape_map_.emplace(id, shape);
  }
    
  // 膨胀
  float sum_area = 0.0;
  for(auto entry : id_shape_map_){
    sum_area += entry.second.first * entry.second.second;
  }

  const float core_area = core_width * core_height;
  expand_ratio_ = sum_area > 0.0f ? core_area / sum_area : 1.0f;

  for(auto& entry : id_shape_map_){
    entry.second.first *= sqrt(expand_ratio_);
    entry.second.second *= sqrt(expand_ratio_);
  }

  // 膨胀完后，若出现超级大的超出core_width或者core_height簇，需要调整
  for(auto& entry : id_shape_map_){
    if(entry.second.first>core_width || entry.second.second>core_height){
      float entry_area = entry.second.first * entry.second.second;
      if(entry.second.first>core_width){
        entry.second.first = core_width-1;
        entry.second.second = entry_area/core_width;
      }else if (entry.second.second>core_height){
        entry.second.second = core_height-1;
        entry.second.first = entry_area/core_height;
      }
    }
  }

}


std::pair<float, float> HierRTLMP::setRandomShapes(Cluster& cluster){
  // 仅针对混合簇和宏簇的macro_tilings生成一个随机的宽高
  const auto& macroTilings = cluster.getMacroTilings();
  if (macroTilings.empty()) {
        throw std::runtime_error("No macro tilings available for the cluster.");
  }

  std::random_device rd;
  std::mt19937 gen(rd()); 
  std::uniform_int_distribution<> dis(0, macroTilings.size()-1);
  int random_number = dis(gen);

  float tiling_width = cluster.getMacroTilings()[random_number].first;
  float tiling_height = cluster.getMacroTilings()[random_number].second;
  // std::pair<float, float> shape(tiling_width, tiling_height);
  return std::make_pair(tiling_width, tiling_height);
}

////////////////////////////////////////////////////////////////////////
// Write Cluster
////////////////////////////////////////////////////////////////////////
void HierRTLMP::WriteClusterResult(){
  // std::string path = stage_path_;
  bool current_io = true;

  ioid_clusterid_vector_.clear();

  for(auto& IOcluster:tree_->root->getChildren()){
    if(!IOcluster->isIOCluster()){
      continue;
    }
    for(auto& connect : IOcluster->getConnection()){
      ioid_clusterid_vector_.push_back(std::make_pair(std::make_pair(IOcluster->getId(),connect.first),current_io));
    }
    current_io = !current_io;
  }

  WriteVerilog();  
  WriteLib();
  // WriteDef();
  WriteLef();  

  if(diverge_flag_ == false){
    WriteDef();
  }

}

void HierRTLMP::WriteVerilog(){
  std::string path = stage_path_;
  std::string line;
  std::string verilog_name = "/root_cluster_"+block_->getName()+".v";
  std::string verilog_path = path + verilog_name;
  std::ofstream write_verilog(verilog_path);
  if(!write_verilog.is_open()){
    std::cerr<<"cannot create file: "<<verilog_path<<std::endl;
    logger_->error(MPL, 1006, "cannot create verilog file");
  }
  // 模块头
  write_verilog<<"module root_cluster (\n";
  int pin_index = 0;
  for(auto& pin_connect : ioid_clusterid_vector_){

    if(pin_connect.second){  //true的话，这是INPUT PIN
      write_verilog<<"inp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second);
    }else{
      write_verilog<<"outp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second);
    }

    pin_index++;
    if(pin_index != ioid_clusterid_vector_.size()){
      write_verilog<<",\n";
    }else{
      write_verilog<<"\n";
    }
  }

  // 👌
  write_verilog<<");\n";
  // cluster_pin_num_map用来存储每个簇有几个连接,从而确定每个簇的端口数
  std::map<int,int> cluster_pin_num_map;
  // all_nets是存储所有顶层簇之间的连接
  // A set to store all the nets between top-level clusters.
  // The PairComparator ensures that the pairs are stored in a specific order.
  std::set<std::pair<int, int>, PairComparator> all_nets;
  for(auto& cluster:tree_->root->getChildren()){
    if(cluster->isIOCluster()){
      continue;
    }
    int connected_num = 0;
    for (const auto& pair : cluster->getConnection()) {
      auto it = tree_->maps.id_to_cluster.find(pair.first);
      if (it == tree_->maps.id_to_cluster.end()) {
        continue;
        // logger_->error(MPL, 1009, "Cannot find cluster in id_to_cluster map {}",pair.first);
      }  // 找不到这个连接对应的cluster，这里就不报错，直接跳过
      
      Cluster* clusterPtr = it->second; // 获取 Cluster 指针
      if(clusterPtr->isIOCluster() == false){
        if(clusterPtr->getParent() == tree_->root.get()){
          connected_num++;
          all_nets.emplace(cluster->getId(),pair.first);
        }
      }
      // if(tree_->maps.id_to_cluster[pair.first]->isIOCluster() == false){  //这个连接不是与IO，且是顶层
      //   if(tree_->maps.id_to_cluster[pair.first]->getParent() == tree_->root.get()){
      //     connected_num++;
      //     all_nets.emplace(cluster->getId(),pair.first);
      //   }
      // }
    }
    cluster_pin_num_map.emplace(cluster->getId(),connected_num);
  }
  // 写线网
  for(int net_ordered = 1;net_ordered<all_nets.size()+1;net_ordered++){
    line = "wire   w"+std::to_string(net_ordered)+';';
    write_verilog<<line<<std::endl;
  }
  // 写pin相关的线网
  for (auto& pin_connect : ioid_clusterid_vector_){
    if(pin_connect.second){  //true的话，这是INPUT PIN
      line = "input inp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+";";
    }else{
      line = "output outp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+";";
    }
    write_verilog<<line<<std::endl;
  }


  // 写连接关系
  // 写一个map,ID-vector,这个vector里就是依次的线网
  // 遍历all_nets,找到每个2 pin-nets对应的两个map,在这俩映射的vector后push_back
  std::map<int,std::vector<int>> cluster_wire_map;
  int wire_order = 1;
  for(auto net : all_nets){
    // 判断当前cluster_wire_map中是否有net对应的map中的那个数组
    // 没有就开一个ID-vector
    if(cluster_wire_map.find(net.first) == cluster_wire_map.end()){
      std::vector<int> temp_vector1;
      cluster_wire_map[net.first] = temp_vector1;
    }
    if(cluster_wire_map.find(net.second) == cluster_wire_map.end()){
      std::vector<int> temp_vector2;
      cluster_wire_map[net.second] = temp_vector2;
    }

    cluster_wire_map[net.first].push_back(wire_order);
    cluster_wire_map[net.second].push_back(wire_order);
    wire_order++;
  }

  // 写cluster-wire
  int write_order = 1;
  std::string str = std::to_string(write_order);
  for(auto& cluster:tree_->root->getChildren()){
    if(cluster->isIOCluster()){
      continue;
    }
    // 类型是cluster的name_
    std::string type = "fake_macro" + std::to_string(write_order);
    // std::string type = cluster->getName();
    // 名称是U+顺序
    // std::string name = 'U'+std::to_string(write_order);
    // 名称修改成U+ID
    std::string name = 'U'+std::to_string(cluster->getId());

    // pin数量
    int pin_num = cluster_wire_map[cluster->getId()].size();
    // int pin_num = cluster_pin_num_map[cluster->getId()];
    line = "  "+type+' '+name+" ( ";

    // 先加与IO连接的
    for(auto& pin_connect : ioid_clusterid_vector_){
      if(pin_connect.first.second == cluster->getId()){
        if(pin_connect.second){  //pin是input,那这个单元连接的也是input，命名为I_pin+连接的簇id
          std::string wire = "inp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second);
          line = line+" .I_pin_"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+"("+wire+"),";
        }else{
          std::string wire = "outp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second);
          line = line+" .O_pin_"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+"("+wire+"),";
        }
      }
    }

    for(int i=0;i<pin_num;i++){
      // 找此时的线网
      int wire = cluster_wire_map[cluster->getId()][i];
      if((i+1)%2==1){
        if(i<pin_num-1){
          line = line + " .I"+std::to_string((i+2)/2)+"(w"+std::to_string(wire)+"),";
        }else{
          line = line + " .I"+std::to_string((i+2)/2)+"(w"+std::to_string(wire)+")";
        }
      // line = line + " .I"+std::to_string((i+2)/2)+"(w"+std::to_string(wire)+"),";
      }else{
        if(i<pin_num-1){
          line = line + " .O"+std::to_string((i+1)/2)+"(w"+std::to_string(wire)+"),";
        }else{
          line = line + " .O"+std::to_string((i+1)/2)+"(w"+std::to_string(wire)+")";
        }
      // line = line + " .O"+std::to_string((i+1)/2)+"(w"+std::to_string(wire)+"),";
      }
    }

    line = line + " ) ;";
    write_verilog<<line<<std::endl;
    write_order++;
  }

  line = "endmodule";
  write_verilog<<line<<std::endl;
  write_verilog.close();
}

void HierRTLMP::WriteDef(){
  std::string path = stage_path_;
  // 写def
  // 创建新def文件
  std::string case_name = block_->getName();

  std::string def_name = "/root_cluster_"+block_->getName()+".def";
  std::string def_path = path + def_name;
  std::string origin_def_path;
  
  // 依据不同例子效仿对应的def文件
  if(case_name == "bp_fe_top"){
    origin_def_path = "/workspace/OpenROAD-flow-scripts/flow/results/nangate45/bp_fe/base/2_2_floorplan_io.def";
  }else if(case_name == "black_parrot"){
    origin_def_path = "/workspace/OpenROAD-flow-scripts/flow/results/nangate45/bp/base/2_2_floorplan_io.def";
  }else if(case_name == "superblue18"){
    origin_def_path = "/workspace/iccad2015.ot/sp-result-mpl2/superblue18/pin-placed.def";
  }else if(case_name == "ariane"){
    origin_def_path = "/workspace/ariane133/ariane133_fp.def";
  }else if(case_name == "swerv_wrapper"){
    origin_def_path = "/workspace/OpenROAD-flow-scripts/flow/results/nangate45/swerv_wrapper43/base/2_3_floorplan_macro.def";
    // origin_def_path = "/workspace/OpenROAD-flow-scripts/flow/results/nangate45/swerv_wrapper/base/2_3_floorplan_macro.def";
  }else if(case_name == "bp_be_top"){
    origin_def_path = "/workspace/OpenROAD-flow-scripts/flow/results/nangate45/bp_be/base/2_3_floorplan_macro.def";
  }else if (case_name == "bp_multi_top"){
    origin_def_path = "/workspace/OpenROAD-flow-scripts/flow/results/nangate45/bp_multi/base/2_3_floorplan_macro.def";
  }else if (case_name == "npu_top"){
    origin_def_path = "/workspace/OpenROAD-flow-scripts/flow/results/nangate45/isa_npu/base/2_3_floorplan_macro.def"; 
  }else if (case_name == "dft_top"){
    origin_def_path = "/workspace/OpenROAD-flow-scripts/flow/results/nangate45/dft68/base/2_3_floorplan_macro.def"; 
  }else if (case_name == "superblue1"){
    origin_def_path = "/workspace/iccad2015.ot/sp-result-mpl2/superblue1/pin-placed.def";
  }else if (case_name == "superblue3"){
    origin_def_path = "/workspace/iccad2015.ot/sp-result-mpl2/superblue3/pin-placed.def";
  }else if (case_name == "superblue4"){
    origin_def_path = "/workspace/iccad2015.ot/sp-result-mpl2/superblue4/pin-placed.def";
  }else if (case_name == "superblue5"){
    origin_def_path = "/workspace/iccad2015.ot/sp-result-mpl2/superblue5/pin-placed.def";
  }else if (case_name == "superblue7"){
    origin_def_path = "/workspace/iccad2015.ot/sp-result-mpl2/superblue7/pin-placed.def";
  }else if (case_name == "superblue10"){
    origin_def_path = "/workspace/iccad2015.ot/sp-result-mpl2/superblue10/pin-placed.def";
  }else if (case_name == "superblue16"){
    origin_def_path = "/workspace/iccad2015.ot/sp-result-mpl2/superblue16/pin-placed.def";
  }else if (case_name == "superblue18"){
    origin_def_path = "/workspace/iccad2015.ot/sp-result-mpl2/superblue18/pin-placed.def";
  }else if (case_name == "mor1kx"){
    origin_def_path = "/workspace/OpenROAD-flow-scripts/flow/results/nangate45/mor1kx/base/2_3_floorplan_macro.def";
  }else if (case_name == "or1200_top"){
    origin_def_path = "/workspace/OpenROAD-flow-scripts/flow/results/nangate45/or1200/base/2_3_floorplan_macro.def";
  }else if (case_name == "eth_top"){
    origin_def_path = "/workspace/OpenROAD-flow-scripts/flow/results/nangate45/ethernet/base/2_3_floorplan_macro.def";
  }else{
    std::cerr<<"cannot find the origin def file"<<std::endl;
    logger_->error(MPL, 1002, "cannot find the origin def file");
  }
  // 打开原始def进行读取,这里先以bp_fe为例子，debug时找不到把例子名称放哪了
  // COMPONENTS前照着bp_fe/2_2_floorplan_io.def写
  // Archive/reproduction Tcl files can provide the exact benchmark floorplan DEF.
  // This keeps the top-level root placement setup tied to the active design
  // instead of relying on stale absolute development paths.
  if (const char* env_origin_def = std::getenv("MPL2_ORIGIN_DEF_PATH")) {
    if (std::strlen(env_origin_def) > 0) {
      origin_def_path = env_origin_def;
    }
  }

  std::ifstream origin_def(origin_def_path);
  if(!origin_def.is_open()){
    std::cerr<<" cannot open file: "<<origin_def_path<<std::endl;
  }

  std::ofstream write_def(def_path);
  if(!write_def.is_open()){
    std::cerr<<"cannot create file: "<<def_path<<std::endl;
    logger_->error(MPL, 1007, "cannot create file");
  }

  std::string line;
  while(std::getline(origin_def,line)){
    if(line.rfind("COMPONENTS",0) == 0){
      break;
    }
    write_def<<line<<std::endl;
  }
  origin_def.close();

  // COMONENTS部分,这个重载是bloat的，直接不写
  // line = "COMPONENTS "+std::to_string(id_shape_map.size())+" ;";
  // write_def<<line<<std::endl;
  // // pin之后的部分
  // int write_order = 1;
  // for(auto& cluster : tree_->root->getChildren()){
  //   if(cluster->isIOCluster() == true){
  //     continue;
  //   }
  //   // int id_cluster = cluster->getId();
  //   line = "    - "+std::string(1, 'U')+std::to_string(write_order)+" "+"fake_macro"+std::to_string(write_order)+" ;";
  //   write_order++;
  //   write_def<<line<<std::endl;
  // }
  // write_def<<"END COMPONENTS"<<std::endl; 

  // 写PINS
  write_def<<"PINS "<<std::to_string(ioid_clusterid_vector_.size())<<" ;"<<std::endl;

  for(auto& pin_connect:ioid_clusterid_vector_){
    // 找id和位置
    auto it = tree_->maps.id_to_cluster.find(pin_connect.first.first);
    if(it == tree_->maps.id_to_cluster.end()){
      continue;
    }
    Cluster* clusterPtr = it->second;
    float x = clusterPtr->getSoftMacro()->getX();
    float y = clusterPtr->getSoftMacro()->getY();
    float width = clusterPtr->getSoftMacro()->getWidth();
    float height = clusterPtr->getSoftMacro()->getHeight();
    int x_pin = 0;
    int y_pin = 0;
    const int pin_offset = block_->micronsToDbu(0.07);

    if(height==0){
      if(y>1){  // 👆pin
        x_pin = block_->micronsToDbu(x + width / 2);
        y_pin = block_->micronsToDbu(y + height) - pin_offset;
      }else{  // 👇pin
        x_pin = block_->micronsToDbu(x + width / 2);
        y_pin = block_->micronsToDbu(y) + pin_offset;        
      }
    }else if(width==0){
      if(x>1){  // 👉pin
        x_pin = block_->micronsToDbu(x + width) - pin_offset;
        y_pin = block_->micronsToDbu(y + height / 2);
      }else{  // 👈pin
        x_pin = block_->micronsToDbu(x) + pin_offset;
        y_pin = block_->micronsToDbu(y + height / 2);
      }
    }
    if(pin_connect.second){  //true的话，这是INPUT PIN
      line = "    - inp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+
                  " + NET inp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+" + DIRECTION INPUT + USE SIGNAL\n"
             "      + PORT\n"
             "        + LAYER metal6 ( -140 -140 ) ( 140 140 )\n"
             "        + PLACED ( "+std::to_string(x_pin)+" "+std::to_string(y_pin)+" ) N ;";
    }else{
      line = "    - outp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+
                  " + NET outp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+" + DIRECTION OUTPUT + USE SIGNAL\n"
             "      + PORT\n"
             "        + LAYER metal6 ( -140 -140 ) ( 140 140 )\n"
             "        + PLACED ( "+std::to_string(x_pin)+" "+std::to_string(y_pin)+" ) N ;";
    }
    write_def<<line<<std::endl;
  }

  write_def<<"END PINS"<<std::endl;


  // 写Net，比较困难，直接写verilog让openroad生成
  write_def<<"END DESIGN"<<std::endl;

}

void HierRTLMP::WriteLef(){
  std::string path = stage_path_;
  // 写lef，仿照macro_only.lef写
  std::string lef_name = "/root_cluster_"+block_->getName()+".lef";
  std::string lef_path = path + lef_name;

  std::string start = "VERSION 5.6 ;\n"
                      "BUSBITCHARS \"[]\" ;\n"
                      "DIVIDERCHAR \"/\" ;\n";
  
  // PIN,这里要修改,不是简单的每个簇一个input一个output
  // 连接关系不对,不能简单的认为只有一个输入输出,要数连接,根据IO_cluster判断这个连接是input还是output
  // 所有线网都是二端线网

  // 由于金属层的不清楚，所以I1-I4和O1-O4和HM_100x400_4x4的一样，多的就再从1开始
  // 用来存储每个簇有几个连接,从而确定每个簇的IO数,写PIN
  // 和WriteVerilog一样设计的pin数目
  std::set<std::pair<int, int>, PairComparator> all_nets;
  for(auto& cluster:tree_->root->getChildren()){
    if(cluster->isIOCluster()){
      continue;
    }
    for (const auto& pair : cluster->getConnection()) {
      auto it = tree_->maps.id_to_cluster.find(pair.first);
      if (it == tree_->maps.id_to_cluster.end()) {
        continue;
      } 
      Cluster* clusterPtr = it->second; // 获取 Cluster 指针
      if(clusterPtr->isIOCluster() == false){
        if(clusterPtr->getParent() == tree_->root.get()){
          all_nets.emplace(cluster->getId(),pair.first);
        }
      }
    }
  }

  std::map<int,std::vector<int>> cluster_wire_map;
  int wire_order = 1;
  for(auto net : all_nets){
    if(cluster_wire_map.find(net.first) == cluster_wire_map.end()){
      std::vector<int> temp_vector1;
      cluster_wire_map[net.first] = temp_vector1;
    }
    if(cluster_wire_map.find(net.second) == cluster_wire_map.end()){
      std::vector<int> temp_vector2;
      cluster_wire_map[net.second] = temp_vector2;
    }
    cluster_wire_map[net.first].push_back(wire_order);
    cluster_wire_map[net.second].push_back(wire_order);
    wire_order++;
  }

  std::ofstream write_lef(lef_path);
  if(!write_lef.is_open()){
    std::cerr<<"cannot create file: "<<lef_path<<std::endl;
    logger_->error(MPL, 1008, "cannot create lef file");
  }
  write_lef<<start<<std::endl;

  std::string line;
  int write_order = 1;
  for(auto& cluster:tree_->root->getChildren()){
    if(cluster->isIOCluster()){
      continue;
    }
    line = "MACRO fake_macro"+std::to_string(write_order);
    write_lef<<line<<std::endl;
    write_lef<<"  CLASS BLOCK ;"<<std::endl;
    line = "  FOREIGN fake_macro"+std::to_string(write_order)+" 0 0 ;";
    write_lef<<line<<std::endl;
    line = "  ORIGIN 0 0 ;\n"
           "  SYMMETRY X Y ;";
    write_lef<<line<<std::endl;
    std::pair<float,float> shape = id_shape_map_[cluster->getId()];
    line = "  SIZE "+std::to_string(std::floor(shape.first))+" BY "+std::to_string(std::floor(shape.second))+" ;";
    write_lef<<line<<std::endl;

    // 写pin
    int num_pin = cluster_wire_map[cluster->getId()].size();
    // 与IO连接的pin
    for(auto& pin_connect : ioid_clusterid_vector_){
      if(pin_connect.first.second == cluster->getId()){
        if(pin_connect.second){  //pin是input,那这个单元连接的也是input，命名为I_pin+连接的簇id
          line = "  PIN I_pin_"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+
          "\n    DIRECTION INPUT ;"
          "\n    PORT"
          "\n      LAYER metal3 ;"
          "\n        RECT 0 210 0.140 210.070  ;"
          "\n    END"
          "\n  END I_pin_"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second);
        }else{
          line = "  PIN O_pin_"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+
          "\n    DIRECTION OUTPUT ;"
          "\n    PORT"
          "\n      LAYER metal3 ;"
          "\n        RECT 99.860 210 100 210.070  ;"
          "\n    END"
          "\n  END O_pin_"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second);
        }
        write_lef<<line<<std::endl;
      }
    }



    for(int i=1;i<num_pin+1;i++){
      if(i%2 == 1){
        line = "  PIN I"+std::to_string((i+1)/2)+
               "\n    DIRECTION INPUT ;"
               "\n    PORT"
               "\n      LAYER metal3 ;"
               "\n        RECT 0 210 0.140 210.070  ;"
               "\n    END"
               "\n  END I"+std::to_string((i+1)/2);
      }else{
        line = "  PIN O"+std::to_string(i/2)+
               "\n    DIRECTION OUTPUT ;"
               "\n    PORT"
               "\n      LAYER metal3 ;"
               "\n        RECT 99.860 210 100 210.070  ;"
               "\n    END"
               "\n  END O"+std::to_string(i/2);
      }
      write_lef<<line<<std::endl;
    }
    // write_lef<<pin<<std::endl;
    line = "  OBS\n"
           "      LAYER metal3 ;\n"
           "        RECT 0 0 100 100 ;\n"
           "  END\n";
    write_lef<<line;
    line = "END fake_macro"+std::to_string(write_order)+"\n";
    write_order++;
    write_lef<<line<<std::endl;
  }

  write_lef.close();
  // 后面macro_only.lef中DoubleHeightSite及后续的部分暂时不写
}

void HierRTLMP::WriteLib(){
  std::string path = stage_path_;
  // 写lib,前面和macro_only一样，这是另外的lib和lef，之前tcl里要写读Nangate45.lef
  std::string lib_name = "/root_cluster_"+block_->getName()+".lib";
  std::string lib_path = path + lib_name;
  std::string origin_lib_path = "/workspace/OpenROAD-flow-scripts/tools/OpenROAD/src/mpl2/test/testcases/macro_only.lib";
  std::ifstream origin_lib(origin_lib_path);
  if(!origin_lib.is_open()){
    std::cerr<<" cannot open file: "<<origin_lib_path<<std::endl;
    logger_->error(MPL, 1004, "cannot open lib file");
  }
  std::ofstream write_lib(lib_path);
  if(!write_lib.is_open()){
    std::cerr<<" cannot open file: "<<lib_path<<std::endl;
    logger_->error(MPL, 1005, "cannot open lib file");
  }
  std::string line;

  // cell前和macro_only.lib一样
  while(std::getline(origin_lib,line)){ 
    if(line.compare(0, 8, "    cell") == 0){
      break;
    }
    write_lib<<line<<std::endl;
  }
  origin_lib.close();
  // cell部分
  // 和WriteVerilog一样设计的pin数目
  std::set<std::pair<int, int>, PairComparator> all_nets;
  for(auto& cluster:tree_->root->getChildren()){
    if(cluster->isIOCluster()){
      continue;
    }
    for (const auto& pair : cluster->getConnection()) {
      auto it = tree_->maps.id_to_cluster.find(pair.first);
      if (it == tree_->maps.id_to_cluster.end()) {
        continue;
      } 
      Cluster* clusterPtr = it->second; // 获取 Cluster 指针
      if(clusterPtr->isIOCluster() == false){
        if(clusterPtr->getParent() == tree_->root.get()){
          all_nets.emplace(cluster->getId(),pair.first);
        }
      }
    }
  }

  std::map<int,std::vector<int>> cluster_wire_map;
  int wire_order = 1;
  for(auto net : all_nets){
    if(cluster_wire_map.find(net.first) == cluster_wire_map.end()){
      std::vector<int> temp_vector1;
      cluster_wire_map[net.first] = temp_vector1;
    }
    if(cluster_wire_map.find(net.second) == cluster_wire_map.end()){
      std::vector<int> temp_vector2;
      cluster_wire_map[net.second] = temp_vector2;
    }
    cluster_wire_map[net.first].push_back(wire_order);
    cluster_wire_map[net.second].push_back(wire_order);
    wire_order++;
  }

  int write_order = 1;
  // 写cell
  for(auto& cluster:tree_->root->getChildren()){
    if(cluster->isIOCluster()){
      continue;
    }
    line = "    cell(fake_macro"+std::to_string(write_order)+") {";
    write_order++;
    write_lib<<line<<std::endl;
    int area = std::floor(std::floor(id_shape_map_[cluster->getId()].first) * std::floor(id_shape_map_[cluster->getId()].second));
    line = "        area : "+std::to_string(area)+";";
    write_lib<<line<<std::endl;
    write_lib<<"        interface_timing : true;\n"
               "        is_macro : true;\n\n";
    int pin_num = cluster_wire_map[cluster->getId()].size();

    // 与IO连接的pin
    for(auto& pin_connect : ioid_clusterid_vector_){
      if(pin_connect.first.second == cluster->getId()){
        if(pin_connect.second){  //pin是input,那这个单元连接的也是input，命名为I_pin+连接的簇id
          line = "    pin (I_pin_"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+") {";
          write_lib<<line<<std::endl;
          write_lib<<"        direction : input;\n"
                   "    }\n";
        }else{
          line = "    pin (O_pin_"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+") {";
          write_lib<<line<<std::endl;
          write_lib<<"        direction : output;\n"
                   "    }\n";
        }
      }
    }

    for(int i=1;i<pin_num +1;i++){
      if(i%2==1){
        line = "    pin (I"+std::to_string((i+1)/2)+") {\n"
               "        direction : input;\n"
               "    }";
      }else{
        line = "    pin (O"+std::to_string(i/2)+") {\n"
               "        direction : output;\n"
               "    }";
      }
      write_lib<<line<<std::endl;
    }
    write_lib<<"    }\n"<<std::endl;
  }
  write_lib<<"}\n";
  write_lib.close();
}

////////////////////////////////////////////////////////////////////////
// run Replace
////////////////////////////////////////////////////////////////////////
void HierRTLMP::runReplaceTcl(){
  std::string path = stage_path_;

  // 计算要扩散的宏簇数
  int num = 0;
  for(auto& cluster : tree_->root->getChildren()){
    if(cluster->isIOCluster()){
      continue;
    }
    num++;
  }

  // 写tcl文件
  std::string tcl_name = "root_"+block_->getName()+".tcl";
  std::string write_tcl_path = path + "/"+tcl_name;
  std::ofstream write_tcl(write_tcl_path);
  if(!write_tcl.is_open()){
    std::cerr<<"cannot create file: "<<write_tcl_path<<std::endl;
    logger_->error(MPL, 1001, "cannot create tcl file");
  }
  std::string line;
  if(block_->getName() == "superblue4"){
    write_tcl<<"read_lib /workspace/iccad2015.ot/superblue4/superblue4_Late.lib"<<std::endl;
    write_tcl<<"read_lef /workspace/iccad2015.ot/superblue4/superblue4.lef"<<std::endl;
  }else{
    write_tcl<<"read_lef /workspace/OpenROAD-flow-scripts/flow/platforms/nangate45/lef/NangateOpenCellLibrary.tech.lef"<<std::endl;
  }
  // write_tcl<<"read_lef /workspace/OpenROAD-flow-scripts/flow/platforms/nangate45/lef/NangateOpenCellLibrary.tech.lef"<<std::endl;
  write_tcl<<"read_lef "+stage_path_+"/root_cluster_"+block_->getName()+".lef"<<std::endl;
  write_tcl<<"read_lib "+stage_path_+"/root_cluster_"+block_->getName()+".lib"<<std::endl;
  write_tcl<<"read_verilog "+stage_path_+"/root_cluster_"+block_->getName()+".v"<<std::endl;
  write_tcl<<"link_design \"root_cluster\""<<std::endl;
  write_tcl<<"read_def "+stage_path_+"/root_cluster_"+block_->getName()+".def -floorplan_initialize"<<std::endl;
  
  // if(num<5){
  //   write_tcl << "catch {global_placement -density 0.9 -init_wirelength_coef 0.05 -bin_grid_count 64 -overflow 0.1} errMsg" << std::endl;
  // }else{
  //   write_tcl << "catch {global_placement -density 0.9 -init_wirelength_coef 0.05 -bin_grid_count 64 -overflow 0.2} errMsg" << std::endl;
  // }

  // 使用 catch 捕获错误
  write_tcl << "catch {global_placement -density 0.9 -init_wirelength_coef 0.05 -bin_grid_count 64 -overflow 0.2} errMsg" << std::endl;
  write_tcl << "if { [string length $errMsg] > 0 } { puts \"WARNING: root_cluster placement failed with error: $errMsg\" } else { puts \"INFO: root_cluster placement executed successfully\" }" << std::endl;

  // 确保 write_def 被执行
  write_tcl << "write_def " + block_->getName() + "_root_placed.def" << std::endl;
  write_tcl << "exit" << std::endl;

  write_tcl.close();

  // 改变工作目录
  if (chdir(path.c_str()) != 0) {
      perror("chdir failed");
      return;
  }
  // std::string tcl_path = "root_"+block_->getName()+".tcl";
  std::string command = "../openroad "+tcl_name;

  // 用popen读取RePlace的输出
  FILE* pipe = popen(command.c_str(),"r");
  if(!pipe){
    std::cerr<<"Failed to execute RePlace" << std::endl;
    return;
  }

  char buffer[256];
  while(fgets(buffer,sizeof(buffer),pipe)!=nullptr){
    std::cout<<buffer;
    if(strstr(buffer,"GPL-0305")!=nullptr){
      diverge_flag_ = true;  // 如果发现 GPL-0305 错误，设置 flag
    }else{
      diverge_flag_ = false;
    }
  }
  const int replace_status = pclose(pipe);
  if (replace_status != 0) {
    logger_->warn(MPL, 1021, "RePlace command '{}' exited with status {}.", command, replace_status);
  }

  // 先进行测试，无论收不收敛，都shrink
  // shrinkReplace();

  // 如果出现无法收敛，就shrink,再replace
  if(diverge_flag_ == true){
    shrinkReplace();
  }

  // std::cout<<"diverge_flag_ is "<<diverge_flag_<<std::endl;
  // 如果flag是true，修改id_shape_map_和stage_path_,shrink后再进行replace
}

void HierRTLMP::shrinkReplace(){
  // 读取bloat的def文件
  std::string path = stage_path_;
  std::string def_path = path+"/"+block_->getName()+"_root_placed.def";
  std::ifstream def_file(def_path);
  if(!def_file.is_open()){
    std::cerr<<"cannot open def file:"<<def_path<<std::endl;
  }

  // lambda表达式来拆分字符串
  auto splitString = [](const std::string& input) {
    std::istringstream stream(input);
    std::vector<std::string> tokens;
    std::string token;
    // 逐个提取
    while (stream >> token) {
      tokens.push_back(token); // 将每个非空格的单词添加到数组中
    }
    return tokens; // 返回包含所有单词的数组
  };

  std::unordered_map<int,std::pair<float,float>> shrink_loc_map;

  std::string line;
  bool COMPONENT_flag = false;
  while(std::getline(def_file,line)){
    if(line.substr(0, 10) == "COMPONENTS"){
      COMPONENT_flag = true;
      continue;
    }
    if(line.substr(0, 14) == "END COMPONENTS"){
      COMPONENT_flag = false;
      break;
    }
    // 使用 Lambda 表达式拆分字符串
    if(COMPONENT_flag == true){
      std::vector<std::string> tokens = splitString(line);
      std::string num_str = tokens[1].substr(1);
      int id = std::stoi(num_str);
      shrink_loc_map[id] = std::make_pair(block_->dbuToMicrons(std::stoi(tokens[6])),block_->dbuToMicrons(std::stoi(tokens[7])));
    }
  }

  // 计算shrink后的坐标，中心位置不变
  for(auto& entry : shrink_loc_map){
    float shrink_width = id_shape_map_[entry.first].first/expand_ratio_;
    float shrink_height = id_shape_map_[entry.first].second/expand_ratio_;

    entry.second.first = entry.second.first + id_shape_map_[entry.first].first/2 - shrink_width/2;
    entry.second.second = entry.second.second + id_shape_map_[entry.first].second/2 - shrink_height/2;
  }

  // 写shrink的root_cluster的一系列文件，改id_shape_map_是为了写lef和lib,reback时恢复
  for(auto& entry : id_shape_map_){
    entry.second.first = entry.second.first/expand_ratio_;
    entry.second.second = entry.second.second/expand_ratio_;
  }

  stage_path_ = "/workspace/OpenROAD-flow-scripts/tools/OpenROAD/src/mpl2/temp/shrink";
  WriteClusterResult();
  WriteDef(shrink_loc_map);
  // 调用Replace得到结果
  path = stage_path_;

  // 写tcl文件
  std::string tcl_name = "root_"+block_->getName()+".tcl";
  std::string write_tcl_path = path + "/"+tcl_name;
  std::ofstream write_tcl(write_tcl_path);
  if(!write_tcl.is_open()){
    std::cerr<<"cannot create file: "<<write_tcl_path<<std::endl;
    logger_->error(MPL, 1010, "cannot create tcl file");
  }

  write_tcl<<"read_lef /workspace/OpenROAD-flow-scripts/flow/platforms/nangate45/lef/NangateOpenCellLibrary.tech.lef"<<std::endl;
  write_tcl<<"read_lef "+stage_path_+"/root_cluster_"+block_->getName()+".lef"<<std::endl;
  write_tcl<<"read_lib "+stage_path_+"/root_cluster_"+block_->getName()+".lib"<<std::endl;
  write_tcl<<"read_verilog "+stage_path_+"/root_cluster_"+block_->getName()+".v"<<std::endl;
  write_tcl<<"link_design \"root_cluster\""<<std::endl;
  write_tcl<<"read_def "+stage_path_+"/root_cluster_"+block_->getName()+".def -floorplan_initialize"<<std::endl;

  // 使用 catch 捕获错误
  write_tcl << "catch {global_placement -incremental -density 0.6 -init_wirelength_coef 0.05 -overflow 0.2} errMsg" << std::endl;
  write_tcl << "if { [string length $errMsg] > 0 } { puts \"WARNING: shrink root_cluster_placed failed with error: $errMsg\" } else { puts \"INFO: shrink root_cluster_placed executed successfully\" }" << std::endl;
  write_tcl << "write_def " + block_->getName() + "_root_placed.def" << std::endl;
  write_tcl << "exit" << std::endl;

  write_tcl.close();

  // 改变工作目录
  if (chdir(path.c_str()) != 0) {
      perror("chdir failed");
      return;
  }
  std::string command = "../openroad "+tcl_name;

  // 用popen读取RePlace的输出
  FILE* pipe = popen(command.c_str(),"r");
  if(!pipe){
    std::cerr<<"Failed to execute RePlace" << std::endl;
    return;
  }

  char buffer[256];

  while(fgets(buffer,sizeof(buffer),pipe)!=nullptr){
    std::cout<<buffer;
  }
  const int replace_status = pclose(pipe);
  if (replace_status != 0) {
    logger_->warn(MPL, 1022, "Shrink RePlace command '{}' exited with status {}.", command, replace_status);
  }
  
  std::cout<<"shrink Replace finished"<<std::endl;

}

void HierRTLMP::WriteDef(std::unordered_map<int,std::pair<float,float>> shrink_loc_map){
  std::string path = stage_path_;
  std::string case_name = block_->getName();

  std::string def_name = "/root_cluster_"+block_->getName()+".def";
  std::string def_path = path + def_name;

  // origin_def看bloat的结果
  std::string target = "shrink";
  std::string replacement = "bloat";
  std::string origin_def_path = stage_path_;
  size_t pos = origin_def_path.find(target);
  origin_def_path.replace(pos,target.length(),replacement);
  origin_def_path = origin_def_path+"/"+block_->getName()+"_root_placed.def";
  
  // lambda表达式来拆分字符串
  auto splitString = [](const std::string& input) {
    std::istringstream stream(input);
    std::vector<std::string> tokens;
    std::string token;
    // 逐个提取
    while (stream >> token) {
      tokens.push_back(token); // 将每个非空格的单词添加到数组中
    }
    return tokens; // 返回包含所有单词的数组
  };

  // 按照bloat的def来写
  std::ifstream origin_def(origin_def_path);
  if(!origin_def.is_open()){
    std::cerr<<" cannot open file: "<<origin_def_path<<std::endl;
    // logger_->error(MPL, 1003, "cannot open def file");
  }

  std::ofstream write_def(def_path);
  if(!write_def.is_open()){
    std::cerr<<"cannot create file: "<<def_path<<std::endl;
    logger_->error(MPL, 1011, "cannot create file");
  }

  // COMONENTS部分
  std::string line;
  bool components_flag = false;
  while(std::getline(origin_def,line)){
    if(line.rfind("COMPONENTS",0) == 0){
      components_flag = true;
      write_def<<line<<std::endl;
      continue;
    }
    if(line.rfind("END COMPONENTS",0) == 0){
      write_def<<line<<std::endl;
      break;
    }

    if(components_flag == true){
      std::vector<std::string> tokens = splitString(line);
      std::string num_str = tokens[1].substr(1);
      int id = std::stoi(num_str);
      float new_x = shrink_loc_map[id].first;
      float new_y = shrink_loc_map[id].second;

      tokens[6] = std::to_string(block_->micronsToDbu(new_x));
      tokens[7] = std::to_string(block_->micronsToDbu(new_y));

      line = "    ";

      for(auto token : tokens){
        line = line + token + " ";
      }
    }

    write_def<<line<<std::endl;
  }
  origin_def.close();

  // 写PINS
  write_def<<"PINS "<<std::to_string(ioid_clusterid_vector_.size())<<" ;"<<std::endl;

  for(auto& pin_connect:ioid_clusterid_vector_){
    // 找id和位置
    auto it = tree_->maps.id_to_cluster.find(pin_connect.first.first);
    if(it == tree_->maps.id_to_cluster.end()){
      continue;
    }
    Cluster* clusterPtr = it->second;
    float x = clusterPtr->getSoftMacro()->getX();
    float y = clusterPtr->getSoftMacro()->getY();
    float width = clusterPtr->getSoftMacro()->getWidth();
    float height = clusterPtr->getSoftMacro()->getHeight();
    int x_pin = 0;
    int y_pin = 0;
    const int pin_offset = block_->micronsToDbu(0.07);

    if(height==0){
      if(y>1){  // 👆pin
        x_pin = block_->micronsToDbu(x + width / 2);
        y_pin = block_->micronsToDbu(y + height) - pin_offset;
      }else{  // 👇pin
        x_pin = block_->micronsToDbu(x + width / 2);
        y_pin = block_->micronsToDbu(y) + pin_offset;        
      }
    }else if(width==0){
      if(x>1){  // 👉pin
        x_pin = block_->micronsToDbu(x + width) - pin_offset;
        y_pin = block_->micronsToDbu(y + height / 2);
      }else{  // 👈pin
        x_pin = block_->micronsToDbu(x) + pin_offset;
        y_pin = block_->micronsToDbu(y + height / 2);
      }
    }
    if(pin_connect.second){  //true的话，这是INPUT PIN
      line = "    - inp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+
                  " + NET inp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+" + DIRECTION INPUT + USE SIGNAL\n"
             "      + PORT\n"
             "        + LAYER metal6 ( -140 -140 ) ( 140 140 )\n"
             "        + PLACED ( "+std::to_string(x_pin)+" "+std::to_string(y_pin)+" ) N ;";
    }else{
      line = "    - outp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+
                  " + NET outp"+std::to_string(pin_connect.first.first)+"_"+std::to_string(pin_connect.first.second)+" + DIRECTION OUTPUT + USE SIGNAL\n"
             "      + PORT\n"
             "        + LAYER metal6 ( -140 -140 ) ( 140 140 )\n"
             "        + PLACED ( "+std::to_string(x_pin)+" "+std::to_string(y_pin)+" ) N ;";
    }
    write_def<<line<<std::endl;
  }

  write_def<<"END PINS"<<std::endl;
  write_def<<"END DESIGN"<<std::endl;

  write_def.close();

}

////////////////////////////////////////////////////////////////////////
// Result Record
////////////////////////////////////////////////////////////////////////
void HierRTLMP::rebackRePlaceClusterLocation(){

  // 这里stage_path_就能确保是bloat的还是shrink的
  std::string path = stage_path_;
  // 因为规模小，直接写读流程bp_root_placed.def
  std::string def_name = "/"+block_->getName()+"_root_placed.def";
  // std::string def_name = "/bp_root_placed.def";
  std::string def_path = path+def_name;
  std::ifstream def_file(def_path);
  if (!def_file.is_open()) {
    logger_->error(MPL, 1023, "Cannot open RePlace result DEF: {}", def_path);
  }

  // lambda表达式来拆分字符串,这个有错误?
  auto splitString = [](const std::string& input) {
    std::istringstream stream(input);
    std::vector<std::string> tokens;
    std::string token;
    // 逐个提取
    while (stream >> token) {
      tokens.push_back(token); // 将每个非空格的单词添加到数组中
    }
    return tokens; // 返回包含所有单词的数组
  };

  // 为所有簇创建soft_macro
  for (auto& cluster : tree_->root->getChildren()) {
    if (cluster->isIOCluster()) {  // ignore all the io clusters
      continue;
    }
    auto soft_macro = std::make_unique<SoftMacro>(cluster.get());
    // no memory leakage, beacuse we set the soft macro, the old one
    // will be deleted
    cluster->setSoftMacro(std::move(soft_macro));
  }

  std::string line;
  bool COMPONENT_flag = false;
  while(std::getline(def_file,line)){
    if(line.substr(0, 10) == "COMPONENTS"){
      COMPONENT_flag = true;
      continue;
    }
    if(line.substr(0, 14) == "END COMPONENTS"){
      COMPONENT_flag = false;
      break;
    }
    // std::string input = "    - U14 fake_macro2 + PLACED ( 20140 1436242 ) N ;";
    // 使用 Lambda 表达式拆分字符串
    if(COMPONENT_flag == true){
      std::vector<std::string> tokens = splitString(line);
      std::string num_str = tokens[1].substr(1);
      int id = std::stoi(num_str);

      auto current = tree_->maps.id_to_cluster.find(id);
      if(current == tree_->maps.id_to_cluster.end()){
        continue;
      }
      Cluster* cluster_ptr = current->second;
      cluster_ptr->setX(block_->dbuToMicrons(std::stoi(tokens[6])));
      cluster_ptr->setY(block_->dbuToMicrons(std::stoi(tokens[7])));
      // cluster_ptr->getX();
      // tree_->maps.id_to_cluster[id]->setX(std::stoi(tokens[7]));
      // tree_->maps.id_to_cluster[id]->setY(std::stoi(tokens[8]));
    }
  }

  def_file.close();

  // id_shape_map_的回归bloat情况
  if(diverge_flag_ == true){
    for(auto& entry : id_shape_map_){
      entry.second.first = entry.second.first * expand_ratio_;
      entry.second.second = entry.second.second * expand_ratio_;
    }
  }

  // 所有二级子簇视为在父簇的中心
  // for(auto& cluster:tree_->root->getChildren()){
  //   if(cluster->isIOCluster()){
  //     continue;
  //   }
  //   PlaceSubCluster(cluster.get());
  // }

}

void HierRTLMP::stdCellPlacement(){
  // 修正顶层簇的坐标
  for(auto& entry : id_shape_map_){
    int id = entry.first;
    float bloat_width = entry.second.first;
    float bloat_height = entry.second.second;
    float new_width = bloat_width/sqrt(expand_ratio_);
    float new_height = bloat_height/sqrt(expand_ratio_);

    auto current = tree_->maps.id_to_cluster.find(id);
    Cluster* cluster_ptr = current->second;
    // const SoftMacro* soft_macro = cluster_ptr->getSoftMacro();

    float bloat_x = cluster_ptr->getX();
    float bloat_y = cluster_ptr->getY();

    float new_x = bloat_x + bloat_width/2 - new_width/2;
    float new_y = bloat_y + bloat_height/2 - new_height/2;

    cluster_ptr->setX(new_x);
    cluster_ptr->setY(new_y);
    cluster_ptr->getSoftMacro()->setWidth(new_width);
    cluster_ptr->getSoftMacro()->setHeight(new_height);
    entry.second.first = new_width;
    entry.second.second = new_height;
  }

  // 生成次级簇的坐标
  for(auto& cluster:tree_->root->getChildren()){
    if(cluster->isIOCluster()){
      continue;
    }
    PlaceSubCluster(cluster.get());
  }

  generateTemporaryStdCellsPlacement(tree_->root.get());

}

void HierRTLMP::generateTemporaryStdCellsPlacement(Cluster* cluster)
{
  if (cluster->isLeaf() && cluster->getNumStdCell() != 0) {  // 
    for (odb::dbModule* module : cluster->getDbModules()) {
      setModuleStdCellsLocation(cluster, module);
    }

    for (odb::dbInst* leaf_std_cell : cluster->getLeafStdCells()) {
      setTemporaryStdCellLocation(cluster, leaf_std_cell);
    }
  } else {
    for (const auto& child : cluster->getChildren()) {
      generateTemporaryStdCellsPlacement(child.get());
    }
  }
}

void HierRTLMP::PlaceSubCluster(Cluster* parent){
  
  int parent_id = parent->getId();
  auto it = id_shape_map_.find(parent_id);
  if(it == id_shape_map_.end()){
    return;
  }
  std::pair<float,float> shape = it->second;

  // 将所有二级子簇视为在父簇的中心，先得到这个二级子簇的面积，视为正方形
  for(auto& cluster : parent->getChildren()){
    // 为每个子簇创建soft_macro
    auto soft_macro = std::make_unique<SoftMacro>(cluster.get());
    cluster->setSoftMacro(std::move(soft_macro));
    float sub_cluster_area = cluster->getMetrics().getStdCellArea();
    if(cluster->getClusterType() != ClusterType::StdCellCluster){
      for(auto& hardmacro : cluster->getHardMacros()){
        sub_cluster_area += hardmacro->getWidth() * hardmacro->getHeight();
      }
    }

    float sub_cluster_side = std::sqrt(sub_cluster_area);
    cluster->setX(parent->getX() + shape.first/2 - sub_cluster_side/2);
    cluster->setY(parent->getY() + shape.second/2 - sub_cluster_side/2);
  }

}

////////////////////////////////////////////////////////////////////////
// Origin
////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
// Interfaces for setting up options
// Options related to macro placement

void HierRTLMP::setAreaWeight(float weight)
{
  area_weight_ = weight;
}

void HierRTLMP::setOutlineWeight(float weight)
{
  outline_weight_ = weight;
}

void HierRTLMP::setWirelengthWeight(float weight)
{
  wirelength_weight_ = weight;
}

void HierRTLMP::setGuidanceWeight(float weight)
{
  guidance_weight_ = weight;
}

void HierRTLMP::setFenceWeight(float weight)
{
  fence_weight_ = weight;
}

void HierRTLMP::setBoundaryWeight(float weight)
{
  boundary_weight_ = weight;
}

void HierRTLMP::setNotchWeight(float weight)
{
  notch_weight_ = weight;
}

void HierRTLMP::setMacroBlockageWeight(float weight)
{
  macro_blockage_weight_ = weight;
}

void HierRTLMP::setGlobalFence(float fence_lx,
                               float fence_ly,
                               float fence_ux,
                               float fence_uy)
{
  global_fence_lx_ = fence_lx;
  global_fence_ly_ = fence_ly;
  global_fence_ux_ = fence_ux;
  global_fence_uy_ = fence_uy;
}

void HierRTLMP::setHaloWidth(float halo_width)
{
  tree_->halo_width = halo_width;
}

void HierRTLMP::setHaloHeight(float halo_height)
{
  tree_->halo_height = halo_height;
}

void HierRTLMP::setGuidanceRegions(
    const std::map<odb::dbInst*, Rect>& guidance_regions)
{
  guides_ = guidance_regions;
}

// Options related to clustering
void HierRTLMP::setNumBundledIOsPerBoundary(int num_bundled_ios)
{
  tree_->bundled_ios_per_edge = num_bundled_ios;
}

void HierRTLMP::setClusterSize(int max_num_macro,
                               int min_num_macro,
                               int max_num_inst,
                               int min_num_inst)
{
  tree_->base_max_macro = max_num_macro;
  tree_->base_min_macro = min_num_macro;
  tree_->base_max_std_cell = max_num_inst;
  tree_->base_min_std_cell = min_num_inst;
}

void HierRTLMP::setClusterSizeTolerance(float tolerance)
{
  tree_->cluster_size_tolerance = tolerance;
}

void HierRTLMP::setMaxNumLevel(int max_num_level)
{
  tree_->max_level = max_num_level;
}

void HierRTLMP::setClusterSizeRatioPerLevel(float coarsening_ratio)
{
  tree_->cluster_size_ratio = coarsening_ratio;
}

void HierRTLMP::setLargeNetThreshold(int large_net_threshold)
{
  tree_->large_net_threshold = large_net_threshold;
}

void HierRTLMP::setSignatureNetThreshold(int signature_net_threshold)
{
  tree_->min_net_count_for_connection = signature_net_threshold;
}

void HierRTLMP::setPinAccessThreshold(float pin_access_th)
{
  pin_access_th_ = pin_access_th;
  pin_access_th_orig_ = pin_access_th;
}

void HierRTLMP::setTargetUtil(float target_util)
{
  target_util_ = target_util;
}

void HierRTLMP::setTargetDeadSpace(float target_dead_space)
{
  target_dead_space_ = target_dead_space;
}

void HierRTLMP::setMinAR(float min_ar)
{
  min_ar_ = min_ar;
}

void HierRTLMP::setSnapLayer(int snap_layer)
{
  snap_layer_ = snap_layer;
}

void HierRTLMP::setReportDirectory(const char* report_directory)
{
  report_directory_ = report_directory;
}

// Transform the logical hierarchy into a physical hierarchy.
void HierRTLMP::runMultilevelAutoclustering()
{
  clustering_engine_ = std::make_unique<ClusteringEngine>(
      block_, network_, logger_, tritonpart_);

  // Set target structures
  clustering_engine_->setDesignMetrics(metrics_);
  clustering_engine_->setTree(tree_.get());

  clustering_engine_->run();

  if (!tree_->has_unfixed_macros) {
    skip_macro_placement_ = true;
    return;
  }

  if (graphics_) {
    graphics_->finishedClustering(tree_->root.get());
  }
}

void HierRTLMP::runHierarchicalMacroPlacement()
{
  if (graphics_) {
    graphics_->startFine();
  }

  adjustMacroBlockageWeight();
  if (logger_->debugCheck(MPL, "hierarchical_macro_placement", 1)) {
    reportSAWeights();
  }

  if (bus_planning_on_) {
    adjustCongestionWeight();
    runHierarchicalMacroPlacement(tree_->root.get());
  } else {
    runHierarchicalMacroPlacementWithoutBusPlanning(tree_->root.get());
  }
}

void HierRTLMP::resetSAParameters()
{
  pos_swap_prob_ = 0.2;
  neg_swap_prob_ = 0.2;
  double_swap_prob_ = 0.2;
  exchange_swap_prob_ = 0.2;
  flip_prob_ = 0.2;
  resize_prob_ = 0.0;
  fence_weight_ = 0.0;
  boundary_weight_ = 0.0;
  notch_weight_ = 0.0;
  macro_blockage_weight_ = 0.0;
}

void HierRTLMP::runCoarseShaping()
{
  setRootShapes();

  if (tree_->has_only_macros) {
    logger_->warn(MPL, 27, "Design has only macros!");
    tree_->root->setClusterType(HardMacroCluster);
    return;
  }

  if (graphics_) {
    graphics_->startCoarse();
  }

  calculateChildrenTilings(tree_->root.get());

  setIOClustersBlockages();
  setPlacementBlockages();
}

void HierRTLMP::setRootShapes()
{
  auto root_soft_macro = std::make_unique<SoftMacro>(tree_->root.get());

  const float core_lx
      = static_cast<float>(block_->dbuToMicrons(block_->getCoreArea().xMin()));
  const float root_lx = std::max(core_lx, global_fence_lx_);

  const float core_ly
      = static_cast<float>(block_->dbuToMicrons(block_->getCoreArea().yMin()));
  const float root_ly = std::max(core_ly, global_fence_ly_);

  const float core_ux
      = static_cast<float>(block_->dbuToMicrons(block_->getCoreArea().xMax()));
  const float root_ux = std::min(core_ux, global_fence_ux_);

  const float core_uy
      = static_cast<float>(block_->dbuToMicrons(block_->getCoreArea().yMax()));
  const float root_uy = std::min(core_uy, global_fence_uy_);

  const float root_area = (root_ux - root_lx) * (root_uy - root_ly);
  const float root_width = root_ux - root_lx;
  const std::vector<std::pair<float, float>> root_width_list
  
      = {std::pair<float, float>(root_width, root_width)};

  root_soft_macro->setShapes(root_width_list, root_area);
  root_soft_macro->setWidth(root_width);  // This will set height automatically
  root_soft_macro->setX(root_lx);
  root_soft_macro->setY(root_ly);
  tree_->root->setSoftMacro(std::move(root_soft_macro));
}

// Compare two intervals according to the product
static bool comparePairProduct(const std::pair<float, float>& p1,
                               const std::pair<float, float>& p2)
{
  return p1.first * p1.second < p2.first * p2.second;
}

// Determine the macro tilings within each cluster in a bottom-up manner.
// (Post-Order DFS manner)
// Coarse shaping:  In this step, we only consider the size of macros
// Ignore all the standard-cell clusters.
// At this stage, we assume the standard cell placer will automatically
// place standard cells in the empty space between macros.
void HierRTLMP::calculateChildrenTilings(Cluster* parent)
{
  // base case, no macros in current cluster
  if (parent->getNumMacro() == 0) {
    return;
  }

  debugPrint(logger_,
             MPL,
             "coarse_shaping",
             1,
             "Determine shapes for {}",
             parent->getName());

  // Current cluster is a hard macro cluster
  if (parent->getClusterType() == HardMacroCluster) {
    debugPrint(logger_,
               MPL,
               "coarse_shaping",
               1,
               "{} is a Macro cluster",
               parent->getName());
    calculateMacroTilings(parent);
    return;
  }

  if (!parent->getChildren().empty()) {
    debugPrint(logger_,
               MPL,
               "coarse_shaping",
               1,
               "Started visiting children of {}",
               parent->getName());

    // Recursively visit the children of Mixed Cluster
    for (auto& cluster : parent->getChildren()) {
      if (cluster->getNumMacro() > 0) {
        calculateChildrenTilings(cluster.get());
      }
    }

    debugPrint(logger_,
               MPL,
               "coarse_shaping",
               1,
               "Done visiting children of {}",
               parent->getName());
  }

  if (graphics_) {
    graphics_->setCurrentCluster(parent);
  }

  // if the current cluster is the root cluster,
  // the shape is fixed, i.e., the fixed die.
  // Thus, we do not need to determine the shapes for it
  // calculate macro tiling for parent cluster based on
  // the macro tilings of its children
  std::vector<SoftMacro> macros;
  for (auto& cluster : parent->getChildren()) {
    if (cluster->getNumMacro() > 0) {
      SoftMacro macro = SoftMacro(cluster.get());
      macro.setShapes(cluster->getMacroTilings(), true);  // force_flag = true
      macros.push_back(macro);
    }
  }
  // if there is only one soft macro
  // the parent cluster has the shape of the its child macro cluster
  if (macros.size() == 1) {
    for (auto& cluster : parent->getChildren()) {
      if (cluster->getNumMacro() > 0) {
        parent->setMacroTilings(cluster->getMacroTilings());
        return;
      }
    }
  }

  debugPrint(
      logger_, MPL, "coarse_shaping", 1, "Running SA to calculate tiling...");

  // call simulated annealing to determine tilings
  std::set<std::pair<float, float>> macro_tilings;  // <width, height>
  // the probability of all actions should be summed to 1.0.
  const float action_sum = pos_swap_prob_ + neg_swap_prob_ + double_swap_prob_
                           + exchange_swap_prob_ + resize_prob_;

  const Rect outline(0, 0, tree_->root->getWidth(), tree_->root->getHeight());

  const int num_perturb_per_step = (macros.size() > num_perturb_per_step_ / 10)
                                       ? macros.size()
                                       : num_perturb_per_step_ / 10;
  // we vary the outline of parent cluster to generate different tilings
  // we first vary the outline width while keeping outline height fixed
  // Then we vary the outline height while keeping outline width fixed
  // Vary the outline width
  std::vector<float> vary_factor_list{1.0};
  float vary_step = 1.0 / num_runs_;  // change the outline by based on num_runs
  for (int i = 1; i < num_runs_; i++) {
    vary_factor_list.push_back(1.0 - i * vary_step);
  }
  int remaining_runs = num_runs_;
  int run_id = 0;
  while (remaining_runs > 0) {
    SoftSAVector sa_batch;
    const int run_thread
        = graphics_ ? 1 : std::min(remaining_runs, num_threads_);
    for (int i = 0; i < run_thread; i++) {
      const Rect new_outline(0,
                             0,
                             outline.getWidth() * vary_factor_list[run_id++],
                             outline.getHeight());
      if (graphics_) {
        graphics_->setOutline(micronsToDbu(new_outline));
      }
      std::unique_ptr<SACoreSoftMacro> sa
          = std::make_unique<SACoreSoftMacro>(tree_->root.get(),
                                              new_outline,
                                              macros,
                                              1.0,     // area weight
                                              1000.0,  // outline weight
                                              0.0,     // wirelength weight
                                              0.0,     // guidance weight
                                              0.0,     // fence weight
                                              0.0,     // boundary weight
                                              0.0,     // macro blockage
                                              0.0,     // notch weight
                                              0.0,     // no notch size
                                              0.0,     // no notch size
                                              pos_swap_prob_ / action_sum,
                                              neg_swap_prob_ / action_sum,
                                              double_swap_prob_ / action_sum,
                                              exchange_swap_prob_ / action_sum,
                                              resize_prob_ / action_sum,
                                              init_prob_,
                                              max_num_step_,
                                              num_perturb_per_step,
                                              random_seed_,
                                              graphics_.get(),
                                              logger_);
      sa_batch.push_back(std::move(sa));
    }
    if (sa_batch.size() == 1) {
      runSA<SACoreSoftMacro>(sa_batch[0].get());
    } else {
      // multi threads
      std::vector<std::thread> threads;
      threads.reserve(sa_batch.size());
      for (auto& sa : sa_batch) {
        threads.emplace_back(runSA<SACoreSoftMacro>, sa.get());
      }
      for (auto& th : threads) {
        th.join();
      }
    }
    // add macro tilings
    for (auto& sa : sa_batch) {
      if (sa->isValid(outline)) {
        macro_tilings.insert(
            std::pair<float, float>(sa->getWidth(), sa->getHeight()));
      }
    }
    remaining_runs -= run_thread;
  }
  // vary the outline height while keeping outline width fixed
  remaining_runs = num_runs_;
  run_id = 0;
  while (remaining_runs > 0) {
    SoftSAVector sa_batch;
    const int run_thread
        = graphics_ ? 1 : std::min(remaining_runs, num_threads_);
    for (int i = 0; i < run_thread; i++) {
      const Rect new_outline(0,
                             0,
                             outline.getWidth(),
                             outline.getHeight() * vary_factor_list[run_id++]);
      if (graphics_) {
        graphics_->setOutline(micronsToDbu(new_outline));
      }
      std::unique_ptr<SACoreSoftMacro> sa
          = std::make_unique<SACoreSoftMacro>(tree_->root.get(),
                                              new_outline,
                                              macros,
                                              1.0,     // area weight
                                              1000.0,  // outline weight
                                              0.0,     // wirelength weight
                                              0.0,     // guidance weight
                                              0.0,     // fence weight
                                              0.0,     // boundary weight
                                              0.0,     // macro blockage
                                              0.0,     // notch weight
                                              0.0,     // no notch size
                                              0.0,     // no notch size
                                              pos_swap_prob_ / action_sum,
                                              neg_swap_prob_ / action_sum,
                                              double_swap_prob_ / action_sum,
                                              exchange_swap_prob_ / action_sum,
                                              resize_prob_ / action_sum,
                                              init_prob_,
                                              max_num_step_,
                                              num_perturb_per_step,
                                              random_seed_,
                                              graphics_.get(),
                                              logger_);
      sa_batch.push_back(std::move(sa));
    }
    if (sa_batch.size() == 1) {
      runSA<SACoreSoftMacro>(sa_batch[0].get());
    } else {
      // multi threads
      std::vector<std::thread> threads;
      threads.reserve(sa_batch.size());
      for (auto& sa : sa_batch) {
        threads.emplace_back(runSA<SACoreSoftMacro>, sa.get());
      }
      for (auto& th : threads) {
        th.join();
      }
    }
    // add macro tilings
    for (auto& sa : sa_batch) {
      if (sa->isValid(outline)) {
        macro_tilings.insert(
            std::pair<float, float>(sa->getWidth(), sa->getHeight()));
      }
    }
    remaining_runs -= run_thread;
  }
  std::vector<std::pair<float, float>> tilings(macro_tilings.begin(),
                                               macro_tilings.end());
  std::sort(tilings.begin(), tilings.end(), comparePairProduct);
  for (auto& shape : tilings) {
    debugPrint(logger_,
               MPL,
               "coarse_shaping",
               2,
               "width: {}, height: {}, aspect_ratio: {}, min_ar: {}",
               shape.first,
               shape.second,
               shape.second / shape.first,
               min_ar_);
  }
  // we do not want very strange tilings if we have choices
  std::vector<std::pair<float, float>> new_tilings;
  for (auto& tiling : tilings) {
    if (tiling.second / tiling.first >= min_ar_
        && tiling.second / tiling.first <= 1.0 / min_ar_) {
      new_tilings.push_back(tiling);
    }
  }
  // if there are valid tilings
  if (!new_tilings.empty()) {
    tilings = std::move(new_tilings);
  }
  // update parent
  parent->setMacroTilings(tilings);
  if (tilings.empty()) {
    logger_->error(MPL,
                   3,
                   "There are no valid tilings for mixed cluster: {}",
                   parent->getName());
  } else {
    std::string line
        = "The macro tiling for mixed cluster " + parent->getName() + "  ";
    for (auto& shape : tilings) {
      line += " < " + std::to_string(shape.first) + " , ";
      line += std::to_string(shape.second) + " >  ";
    }
    line += "\n";
    debugPrint(logger_, MPL, "coarse_shaping", 2, "{}", line);
  }
}

void HierRTLMP::calculateMacroTilings(Cluster* cluster)
{
  if (cluster->getClusterType() != HardMacroCluster) {
    return;
  }

  std::vector<HardMacro*> hard_macros = cluster->getHardMacros();
  std::set<std::pair<float, float>> macro_tilings;  // <width, height>

  if (hard_macros.size() == 1) {
    float width = hard_macros[0]->getWidth();
    float height = hard_macros[0]->getHeight();

    std::vector<std::pair<float, float>> tilings;

    tilings.emplace_back(width, height);
    cluster->setMacroTilings(tilings);

    debugPrint(logger_,
               MPL,
               "coarse_shaping",
               1,
               "{} has only one macro, set tiling according to macro with halo",
               cluster->getName());
    return;
  }

  if (cluster->isArrayOfInterconnectedMacros()) {
    setTightPackingTilings(cluster);
    return;
  }

  if (graphics_) {
    graphics_->setCurrentCluster(cluster);
  }

  // otherwise call simulated annealing to determine tilings
  // set the action probabilities
  const float action_sum = pos_swap_prob_ + neg_swap_prob_ + double_swap_prob_
                           + exchange_swap_prob_;

  const Rect outline(0, 0, tree_->root->getWidth(), tree_->root->getHeight());

  // update macros
  std::vector<HardMacro> macros;
  macros.reserve(hard_macros.size());
  for (auto& macro : hard_macros) {
    macros.push_back(*macro);
  }
  int num_perturb_per_step = (macros.size() > num_perturb_per_step_ / 10)
                                 ? macros.size()
                                 : num_perturb_per_step_ / 10;
  if (cluster->getParent() == nullptr) {
    num_perturb_per_step = (macros.size() > num_perturb_per_step_ / 5)
                               ? macros.size()
                               : num_perturb_per_step_ / 5;
  }

  // To generate different macro tilings, we vary the outline constraints
  // we first vary the outline width while keeping outline_height fixed
  // Then we vary the outline height while keeping outline_width fixed
  // We vary the outline of cluster to generate different tilings
  std::vector<float> vary_factor_list{1.0};
  float vary_step = 1.0 / num_runs_;  // change the outline by at most halfly
  for (int i = 1; i < num_runs_; i++) {
    vary_factor_list.push_back(1.0 - i * vary_step);
  }
  int remaining_runs = num_runs_;
  int run_id = 0;
  while (remaining_runs > 0) {
    HardSAVector sa_batch;
    const int run_thread
        = graphics_ ? 1 : std::min(remaining_runs, num_threads_);
    for (int i = 0; i < run_thread; i++) {
      const Rect new_outline(0,
                             0,
                             outline.getWidth() * vary_factor_list[run_id++],
                             outline.getHeight());
      if (graphics_) {
        graphics_->setOutline(micronsToDbu(new_outline));
      }
      std::unique_ptr<SACoreHardMacro> sa
          = std::make_unique<SACoreHardMacro>(new_outline,
                                              macros,
                                              1.0,     // area_weight
                                              1000.0,  // outline weight
                                              0.0,     // wirelength weight
                                              0.0,     // guidance
                                              0.0,     // fence weight
                                              pos_swap_prob_ / action_sum,
                                              neg_swap_prob_ / action_sum,
                                              double_swap_prob_ / action_sum,
                                              exchange_swap_prob_ / action_sum,
                                              0.0,  // no flip
                                              init_prob_,
                                              max_num_step_,
                                              num_perturb_per_step,
                                              random_seed_ + run_id,
                                              graphics_.get(),
                                              logger_);
      sa_batch.push_back(std::move(sa));
    }
    if (sa_batch.size() == 1) {
      runSA<SACoreHardMacro>(sa_batch[0].get());
    } else {
      // multi threads
      std::vector<std::thread> threads;
      threads.reserve(sa_batch.size());
      for (auto& sa : sa_batch) {
        threads.emplace_back(runSA<SACoreHardMacro>, sa.get());
      }
      for (auto& th : threads) {
        th.join();
      }
    }
    // add macro tilings
    for (auto& sa : sa_batch) {
      if (sa->isValid(outline)) {
        macro_tilings.insert(
            std::pair<float, float>(sa->getWidth(), sa->getHeight()));
      }
    }
    remaining_runs -= run_thread;
  }
  // change the outline height while keeping outline width fixed
  remaining_runs = num_runs_;
  run_id = 0;
  while (remaining_runs > 0) {
    HardSAVector sa_batch;
    const int run_thread
        = graphics_ ? 1 : std::min(remaining_runs, num_threads_);
    for (int i = 0; i < run_thread; i++) {
      const Rect new_outline(0,
                             0,
                             outline.getWidth(),
                             outline.getHeight() * vary_factor_list[run_id++]);
      if (graphics_) {
        graphics_->setOutline(micronsToDbu(new_outline));
      }
      std::unique_ptr<SACoreHardMacro> sa
          = std::make_unique<SACoreHardMacro>(new_outline,
                                              macros,
                                              1.0,     // area_weight
                                              1000.0,  // outline weight
                                              0.0,     // wirelength weight
                                              0.0,     // guidance
                                              0.0,     // fence weight
                                              pos_swap_prob_ / action_sum,
                                              neg_swap_prob_ / action_sum,
                                              double_swap_prob_ / action_sum,
                                              exchange_swap_prob_ / action_sum,
                                              0.0,
                                              init_prob_,
                                              max_num_step_,
                                              num_perturb_per_step,
                                              random_seed_ + run_id,
                                              graphics_.get(),
                                              logger_);
      sa_batch.push_back(std::move(sa));
    }
    if (sa_batch.size() == 1) {
      runSA<SACoreHardMacro>(sa_batch[0].get());
    } else {
      // multi threads
      std::vector<std::thread> threads;
      threads.reserve(sa_batch.size());
      for (auto& sa : sa_batch) {
        threads.emplace_back(runSA<SACoreHardMacro>, sa.get());
      }
      for (auto& th : threads) {
        th.join();
      }
    }
    // add macro tilings
    for (auto& sa : sa_batch) {
      if (sa->isValid(outline)) {
        macro_tilings.insert(
            std::pair<float, float>(sa->getWidth(), sa->getHeight()));
      }
    }
    remaining_runs -= run_thread;
  }

  // sort the tilings based on area
  std::vector<std::pair<float, float>> tilings(macro_tilings.begin(),
                                               macro_tilings.end());
  std::sort(tilings.begin(), tilings.end(), comparePairProduct);
  for (auto& shape : tilings) {
    debugPrint(logger_,
               MPL,
               "coarse_shaping",
               2,
               "width: {}, height: {}",
               shape.first,
               shape.second);
  }
  // we only keep the minimum area tiling since all the macros has the same size
  // later this can be relaxed.  But this may cause problems because the
  // minimizing the wirelength may leave holes near the boundary
  std::vector<std::pair<float, float>> new_tilings;
  for (auto& tiling : tilings) {
    if (tiling.first * tiling.second <= tilings[0].first * tilings[0].second) {
      new_tilings.push_back(tiling);
    }
  }
  tilings = std::move(new_tilings);
  // update parent
  cluster->setMacroTilings(tilings);
  if (tilings.empty()) {
    logger_->error(MPL,
                   4,
                   "No valid tilings for hard macro cluster: {}",
                   cluster->getName());
  }

  std::string line = "Tiling for hard cluster " + cluster->getName() + "  ";
  for (auto& shape : tilings) {
    line += " < " + std::to_string(shape.first) + " , ";
    line += std::to_string(shape.second) + " >  ";
  }
  line += "\n";
  debugPrint(logger_, MPL, "coarse_shaping", 2, "{}", line);
}

// Used only for arrays of interconnected macros.
void HierRTLMP::setTightPackingTilings(Cluster* macro_array)
{
  std::vector<std::pair<float, float>> tight_packing_tilings;

  int divider = 1;
  int columns = 0, rows = 0;

  while (divider <= macro_array->getNumMacro()) {
    if (macro_array->getNumMacro() % divider == 0) {
      columns = macro_array->getNumMacro() / divider;
      rows = divider;

      // We don't consider tilings for right angle rotation orientations,
      // because they're not allowed in our macro placer.
      tight_packing_tilings.emplace_back(
          columns * macro_array->getHardMacros().front()->getWidth(),
          rows * macro_array->getHardMacros().front()->getHeight());
    }

    ++divider;
  }

  macro_array->setMacroTilings(tight_packing_tilings);
}

void HierRTLMP::setIOClustersBlockages()
{
  if (!tree_->maps.bterm_to_inst.empty()) {
    return;
  }

  IOSpans io_spans = computeIOSpans();
  const float depth = computeIOBlockagesDepth(io_spans);

  const Rect root(tree_->root->getX(),
                  tree_->root->getY(),
                  tree_->root->getX() + tree_->root->getWidth(),
                  tree_->root->getY() + tree_->root->getHeight());

  // Note that the range can be larger than the respective core dimension.
  // As SA only sees what is inside its current outline, this is not a problem.
  if (io_spans[L].second > io_spans[L].first) {
    const Rect left_io_blockage(root.xMin(),
                                io_spans[L].first,
                                root.xMin() + depth,
                                io_spans[L].second);

    boundary_to_io_blockage_[L] = left_io_blockage;
    macro_blockages_.push_back(left_io_blockage);
  }

  if (io_spans[T].second > io_spans[T].first) {
    const Rect top_io_blockage(io_spans[T].first,
                               root.yMax() - depth,
                               io_spans[T].second,
                               root.yMax());

    boundary_to_io_blockage_[T] = top_io_blockage;
    macro_blockages_.push_back(top_io_blockage);
  }

  if (io_spans[R].second > io_spans[R].first) {
    const Rect right_io_blockage(root.xMax() - depth,
                                 io_spans[R].first,
                                 root.xMax(),
                                 io_spans[R].second);

    boundary_to_io_blockage_[R] = right_io_blockage;
    macro_blockages_.push_back(right_io_blockage);
  }

  if (io_spans[B].second > io_spans[B].first) {
    const Rect bottom_io_blockage(io_spans[B].first,
                                  root.yMin(),
                                  io_spans[B].second,
                                  root.yMin() + depth);

    boundary_to_io_blockage_[B] = bottom_io_blockage;
    macro_blockages_.push_back(bottom_io_blockage);
  }
}

// Determine the range of IOs in each boundary of the die.
HierRTLMP::IOSpans HierRTLMP::computeIOSpans()
{
  IOSpans io_spans;

  odb::Rect die = block_->getDieArea();

  // Initialize spans based on the dimensions of the die area.
  io_spans[L]
      = {block_->dbuToMicrons(die.yMax()), block_->dbuToMicrons(die.yMin())};
  io_spans[T]
      = {block_->dbuToMicrons(die.xMax()), block_->dbuToMicrons(die.xMin())};
  io_spans[R] = io_spans[L];
  io_spans[B] = io_spans[T];

  for (auto term : block_->getBTerms()) {
    if (term->getSigType().isSupply()) {
      continue;
    }

    int lx = std::numeric_limits<int>::max();
    int ly = std::numeric_limits<int>::max();
    int ux = 0;
    int uy = 0;

    for (const auto pin : term->getBPins()) {
      for (const auto box : pin->getBoxes()) {
        lx = std::min(lx, box->xMin());
        ly = std::min(ly, box->yMin());
        ux = std::max(ux, box->xMax());
        uy = std::max(uy, box->yMax());
      }
    }

    // Modify ranges based on the position of the IO pins.
    if (lx <= die.xMin()) {
      io_spans[L].first = std::min(
          io_spans[L].first, static_cast<float>(block_->dbuToMicrons(ly)));
      io_spans[L].second = std::max(
          io_spans[L].second, static_cast<float>(block_->dbuToMicrons(uy)));
    } else if (uy >= die.yMax()) {
      io_spans[T].first = std::min(
          io_spans[T].first, static_cast<float>(block_->dbuToMicrons(lx)));
      io_spans[T].second = std::max(
          io_spans[T].second, static_cast<float>(block_->dbuToMicrons(ux)));
    } else if (ux >= die.xMax()) {
      io_spans[R].first = std::min(
          io_spans[R].first, static_cast<float>(block_->dbuToMicrons(ly)));
      io_spans[R].second = std::max(
          io_spans[R].second, static_cast<float>(block_->dbuToMicrons(uy)));
    } else {
      io_spans[B].first = std::min(
          io_spans[B].first, static_cast<float>(block_->dbuToMicrons(lx)));
      io_spans[B].second = std::max(
          io_spans[B].second, static_cast<float>(block_->dbuToMicrons(ux)));
    }
  }

  return io_spans;
}

// The depth of IO clusters' blockages is generated based on:
// 1) How many vertical or horizontal boundaries have signal IO pins.
// 2) The total length of the io spans in all used boundaries.
float HierRTLMP::computeIOBlockagesDepth(const IOSpans& io_spans)
{
  float sum_length = 0.0;

  for (auto& [pin_access, length] : io_spans) {
    if (length.second > length.first) {
      sum_length += std::abs(length.second - length.first);
    }
  }

  float std_cell_area = 0.0;
  for (auto& cluster : tree_->root->getChildren()) {
    if (cluster->getClusterType() == StdCellCluster) {
      std_cell_area += cluster->getArea();
    }
  }

  const float macro_dominance_factor
      = tree_->macro_with_halo_area
        / (tree_->root->getWidth() * tree_->root->getHeight());
  const float depth = (std_cell_area / sum_length)
                      * std::pow((1 - macro_dominance_factor), 2);

  debugPrint(logger_,
             MPL,
             "coarse_shaping",
             1,
             "Bundled IO clusters blokaged depth = {}",
             depth);

  return depth;
}

void HierRTLMP::setPlacementBlockages()
{
  for (odb::dbBlockage* blockage : block_->getBlockages()) {
    odb::Rect bbox = blockage->getBBox()->getBox();

    Rect bbox_micron(block_->dbuToMicrons(bbox.xMin()),
                     block_->dbuToMicrons(bbox.yMin()),
                     block_->dbuToMicrons(bbox.xMax()),
                     block_->dbuToMicrons(bbox.yMax()));

    placement_blockages_.push_back(bbox_micron);
  }
}

void HierRTLMP::adjustCongestionWeight()
{
  std::string snap_layer_name;

  for (auto& macro : tree_->maps.inst_to_hard) {
    odb::dbMaster* master = macro.first->getMaster();
    for (odb::dbMTerm* mterm : master->getMTerms()) {
      if (mterm->getSigType() == odb::dbSigType::SIGNAL) {
        for (odb::dbMPin* mpin : mterm->getMPins()) {
          for (odb::dbBox* box : mpin->getGeometry()) {
            odb::dbTechLayer* layer = box->getTechLayer();
            snap_layer_name = layer->getName();
          }
        }
      }
    }

    break;
  }

  // update weight
  if (dynamic_congestion_weight_flag_) {
    std::vector<std::string> layers;
    for (odb::dbTechLayer* layer : db_->getTech()->getLayers()) {
      if (layer->getType() == odb::dbTechLayerType::ROUTING) {
        layers.push_back(layer->getName());
      }
    }
    snap_layer_ = 0;
    for (int i = 0; i < layers.size(); i++) {
      if (layers[i] == snap_layer_name) {
        snap_layer_ = i + 1;
        break;
      }
    }
    if (snap_layer_ <= 0) {
      congestion_weight_ = 0.0;
    } else {
      congestion_weight_ = 1.0 * snap_layer_ / layers.size();
    }
    debugPrint(logger_,
               MPL,
               "bus_planning",
               1,
               "Adjusting congestion weight to {} - Snap layer is {}",
               congestion_weight_,
               snap_layer_);
  }
}

// Cluster Placement Engine Starts ...........................................
// The cluster placement is done in a top-down manner
// (Preorder DFS)
// The magic happens at how we determine the size of pin access
// If the size of pin access is too large, SA cannot generate the valid macro
// placement If the size of pin access is too small, the effect of bus synthesis
// can be ignored. Here our trick is to determine the size of pin access and
// standard-cell clusters together. We assume the region occupied by pin access
// will be filled by standard cells. More specifically, we first determine the
// width of pin access based on number of connections passing the pin access,
// then we calculate the height of the pin access based on the available area,
// standard-cell area and macro area. Note that, here we allow the utilization
// of standard-cell clusters larger than 1.0.  If there is no standard cells,
// the area of pin access is 0. Another important trick is that we call SA two
// times. The first time is to determine the location of pin access. The second
// time is to determine the location of children clusters. We assume the
// summation of pin access size is equal to the area of standard-cell clusters
void HierRTLMP::runHierarchicalMacroPlacement(Cluster* parent)
{
  if (parent->getClusterType() == HardMacroCluster) {
    placeMacros(parent);
    return;
  }

  if (parent->isLeaf()) {  // Cover IO Clusters && Leaf Std Cells
    return;
  }

  if (graphics_) {
    graphics_->setCurrentCluster(parent);
  }

  for (auto& cluster : parent->getChildren()) {
    clustering_engine_->updateInstancesAssociation(cluster.get());
  }
  // Place children clusters
  // map children cluster to soft macro
  for (auto& cluster : parent->getChildren()) {
    if (cluster->isIOCluster()) {  // ignore all the io clusters
      continue;
    }
    auto soft_macro = std::make_unique<SoftMacro>(cluster.get());
    // no memory leakage, beacuse we set the soft macro, the old one
    // will be deleted
    cluster->setSoftMacro(std::move(soft_macro));
  }

  // The simulated annealing outline is determined by the parent's shape
  const Rect outline(parent->getX(),
                     parent->getY(),
                     parent->getX() + parent->getWidth(),
                     parent->getY() + parent->getHeight());

  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Working on children of cluster: {}, Outline "
             "{}, {}  {}, {}",
             parent->getName(),
             outline.xMin(),
             outline.yMin(),
             outline.getWidth(),
             outline.getHeight());

  // Suppose the region, fence, guide has been mapped to cooresponding macros
  // This step is done when we enter the Hier-RTLMP program
  std::map<std::string, int> soft_macro_id_map;  // cluster_name, macro_id
  std::map<int, Rect> fences;
  std::map<int, Rect> guides;
  std::vector<SoftMacro> macros;
  std::vector<BundledNet> nets;

  std::vector<Rect> placement_blockages;
  std::vector<Rect> macro_blockages;

  findOverlappingBlockages(macro_blockages, placement_blockages, outline);

  // We store the bundled io clusters to push them into the macros' vector
  // only after it is already populated with the clusters we're trying to
  // place. This will facilitate how we deal with fixed terminals in SA moves.
  std::vector<Cluster*> io_clusters;

  // Each cluster is modeled as Soft Macro
  // The fences or guides for each cluster is created by merging
  // the fences and guides for hard macros in each cluster
  for (auto& cluster : parent->getChildren()) {
    if (cluster->isIOCluster()) {
      io_clusters.push_back(cluster.get());
      continue;
    }
    // for other clusters
    soft_macro_id_map[cluster->getName()] = macros.size();
    auto soft_macro = std::make_unique<SoftMacro>(cluster.get());
    clustering_engine_->updateInstancesAssociation(
        cluster.get());  // we need this step to calculate nets
    macros.push_back(*soft_macro);
    cluster->setSoftMacro(std::move(soft_macro));
    // merge fences and guides for hard macros within cluster
    if (cluster->getClusterType() == StdCellCluster) {
      continue;
    }
    Rect fence, guide;
    fence.mergeInit();
    guide.mergeInit();
    const std::vector<HardMacro*> hard_macros = cluster->getHardMacros();
    for (auto& hard_macro : hard_macros) {
      if (fences_.find(hard_macro->getName()) != fences_.end()) {
        fence.merge(fences_[hard_macro->getName()]);
      }

      auto itr = guides_.find(hard_macro->getInst());
      if (itr != guides_.end()) {
        guide.merge(itr->second);
      }
    }

    // Calculate overlap with outline
    fence.relocate(
        outline.xMin(), outline.yMin(), outline.xMax(), outline.yMax());
    guide.relocate(
        outline.xMin(), outline.yMin(), outline.xMax(), outline.yMax());
    if (fence.isValid()) {
      // current macro id is macros.size() - 1
      fences[macros.size() - 1] = fence;
    }
    if (guide.isValid()) {
      // current macro id is macros.size() - 1
      guides[macros.size() - 1] = guide;
    }
  }

  if (graphics_) {
    graphics_->setGuides(guides);
  }

  clustering_engine_->updateConnections();
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Finished calculating connection");

  int number_of_pin_access = 0;
  // Handle the pin access
  // Get the connections between pin accesses
  if (parent->getParent() != nullptr) {
    // the parent cluster is not the root cluster
    // First step model each pin access as the a dummy softmacro (width = 0.0,
    // height = 0.0) In our simulated annealing engine, the dummary softmacro
    // will no effect on SA We have four dummy SoftMacros based on our
    // definition
    std::vector<Boundary> pins = {L, T, R, B};
    for (auto& pin : pins) {
      soft_macro_id_map[toString(pin)] = macros.size();
      macros.emplace_back(0.0, 0.0, toString(pin));

      ++number_of_pin_access;
    }
    // add the connections between pin accesses, for example, L to R
    for (auto& [src_pin, pin_map] : parent->getBoundaryConnection()) {
      for (auto& [target_pin, weight] : pin_map) {
        nets.emplace_back(soft_macro_id_map[toString(src_pin)],
                          soft_macro_id_map[toString(target_pin)],
                          weight);
      }
    }
  }

  int num_of_macros_to_place = static_cast<int>(macros.size());

  // add the virtual connections (the weight related to IOs and macros belong to
  // the same cluster)
  for (const auto& [cluster1, cluster2] : parent->getVirtualConnections()) {
    BundledNet net(
        soft_macro_id_map[tree_->maps.id_to_cluster[cluster1]->getName()],
        soft_macro_id_map[tree_->maps.id_to_cluster[cluster2]->getName()],
        tree_->virtual_weight);
    net.src_cluster_id = cluster1;
    net.target_cluster_id = cluster2;
    nets.push_back(net);
  }

  // convert the connections between clusters to SoftMacros
  for (auto& cluster : parent->getChildren()) {
    const int src_id = cluster->getId();
    const std::string src_name = cluster->getName();
    for (auto& [cluster_id, weight] : cluster->getConnection()) {
      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 2,
                 " Cluster connection: {} {} {} ",
                 cluster->getName(),
                 tree_->maps.id_to_cluster[cluster_id]->getName(),
                 weight);
      const std::string name = tree_->maps.id_to_cluster[cluster_id]->getName();
      if (soft_macro_id_map.find(name) == soft_macro_id_map.end()) {
        float new_weight = weight;
        if (macros[soft_macro_id_map[src_name]].isStdCellCluster()) {
          new_weight *= tree_->virtual_weight;
        }
        // if the cluster_id is out of the parent cluster
        BundledNet net(
            soft_macro_id_map[src_name],
            soft_macro_id_map[toString(parent->getPinAccess(cluster_id).first)],
            new_weight);
        net.src_cluster_id = src_id;
        net.target_cluster_id = cluster_id;
        nets.push_back(net);
      } else if (src_id > cluster_id) {
        BundledNet net(
            soft_macro_id_map[src_name], soft_macro_id_map[name], weight);
        net.src_cluster_id = src_id;
        net.target_cluster_id = cluster_id;
        nets.push_back(net);
      }
    }
  }
  // merge nets to reduce runtime
  mergeNets(nets);

  if (parent->getParent() != nullptr) {
    // update the size of each pin access macro
    // each pin access macro with have their fences
    // check the net connection, i.e., how many nets crossing the boundaries
    std::map<int, float> net_map;
    net_map[soft_macro_id_map[toString(L)]] = 0.0;
    net_map[soft_macro_id_map[toString(T)]] = 0.0;
    net_map[soft_macro_id_map[toString(R)]] = 0.0;
    net_map[soft_macro_id_map[toString(B)]] = 0.0;
    for (auto& net : nets) {
      if (net_map.find(net.terminals.first) != net_map.end()) {
        net_map[net.terminals.first] += net.weight;
      }
      if (net_map.find(net.terminals.second) != net_map.end()) {
        net_map[net.terminals.second] += net.weight;
      }
    }

    debugPrint(logger_,
               MPL,
               "hierarchical_macro_placement",
               1,
               "Total tracks per micron:  {}",
               1.0 / pin_access_net_width_ratio_);

    // determine the width of each pin access
    float l_size
        = net_map[soft_macro_id_map[toString(L)]] * pin_access_net_width_ratio_;
    float r_size
        = net_map[soft_macro_id_map[toString(R)]] * pin_access_net_width_ratio_;
    float b_size
        = net_map[soft_macro_id_map[toString(B)]] * pin_access_net_width_ratio_;
    float t_size
        = net_map[soft_macro_id_map[toString(T)]] * pin_access_net_width_ratio_;
    const std::vector<std::pair<float, float>> tilings
        = parent->getMacroTilings();
    // When the program enter stage, the tilings cannot be empty
    // const float ar = outline_height / outline_width;
    // float max_height = std::sqrt(tilings[0].first * tilings[0].second * ar);
    // float max_width = max_height / ar;
    // max_height = std::max(max_height, tilings[0].first);
    // max_width = std::max(max_width, tilings[0].second);
    float max_height = 0.0;
    float max_width = 0.0;
    for (auto& tiling : tilings) {
      if (tiling.first <= outline.getWidth()
          && tiling.second <= outline.getHeight()) {
        max_width = std::max(max_width, tiling.first);
        max_height = std::max(max_height, tiling.second);
      }
    }
    max_width = std::min(max_width, outline.getWidth());
    max_height = std::min(max_height, outline.getHeight());
    debugPrint(logger_,
               MPL,
               "hierarchical_macro_placement",
               1,
               "Pin access calculation "
               "max_width: {}, max_height : {}",
               max_width,
               max_height);
    l_size = std::min(l_size, max_height);
    r_size = std::min(r_size, max_height);
    t_size = std::min(t_size, max_width);
    b_size = std::min(b_size, max_width);
    // determine the height of each pin access
    max_width = outline.getWidth() - max_width;
    max_height = outline.getHeight() - max_height;
    // the area of standard-cell clusters
    float std_cell_area = 0.0;
    for (auto& cluster : parent->getChildren()) {
      if (cluster->getClusterType() == StdCellCluster) {
        std_cell_area += cluster->getArea();
      }
    }
    // calculate the depth based on area
    float sum_length = 0.0;
    int num_hor_access = 0;
    int num_ver_access = 0;
    if (l_size > 0.0) {
      num_hor_access += 1;
      sum_length += l_size;
    }
    if (r_size > 0.0) {
      num_hor_access += 1;
      sum_length += r_size;
    }
    if (t_size > 0.0) {
      num_ver_access += 1;
      sum_length += t_size;
    }
    if (b_size > 0.0) {
      num_ver_access += 1;
      sum_length += b_size;
    }
    max_width = num_hor_access > 0 ? max_width / num_hor_access : max_width;
    max_height = num_ver_access > 0 ? max_height / num_ver_access : max_height;
    const float depth = std_cell_area / sum_length;
    // update the size of pin access macro
    if (l_size > 0) {
      const float temp_width = std::min(max_width, depth);
      macros[soft_macro_id_map[toString(L)]]
          = SoftMacro(temp_width, l_size, toString(L));
      fences[soft_macro_id_map[toString(L)]]
          = Rect(0.0, 0.0, temp_width, outline.getWidth());

      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 1,
                 "Left - width : {}, height : {}",
                 l_size,
                 temp_width);
    }
    if (r_size > 0) {
      const float temp_width = std::min(max_width, depth);
      macros[soft_macro_id_map[toString(R)]]
          = SoftMacro(temp_width, r_size, toString(R));
      fences[soft_macro_id_map[toString(R)]]
          = Rect(outline.getWidth() - temp_width,
                 0.0,
                 outline.getWidth(),
                 outline.getHeight());

      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 1,
                 "Right - width : {}, height : {}",
                 r_size,
                 temp_width);
    }
    if (t_size > 0) {
      const float temp_height = std::min(max_height, depth);
      macros[soft_macro_id_map[toString(T)]]
          = SoftMacro(t_size, temp_height, toString(T));
      fences[soft_macro_id_map[toString(T)]]
          = Rect(0.0,
                 outline.getHeight() - temp_height,
                 outline.getWidth(),
                 outline.getHeight());

      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 1,
                 "Top - width : {}, height : {}",
                 t_size,
                 temp_height);
    }
    if (b_size > 0) {
      const float temp_height = std::min(max_height, depth);
      macros[soft_macro_id_map[toString(B)]]
          = SoftMacro(b_size, temp_height, toString(B));
      fences[soft_macro_id_map[toString(B)]]
          = Rect(0.0, 0.0, outline.getWidth(), temp_height);

      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 1,
                 "Bottom - width : {}, height : {}",
                 b_size,
                 temp_height);
    }
  }

  for (Cluster* io_cluster : io_clusters) {
    soft_macro_id_map[io_cluster->getName()] = macros.size();

    macros.emplace_back(
        std::pair<float, float>(io_cluster->getX() - outline.xMin(),
                                io_cluster->getY() - outline.yMin()),
        io_cluster->getName(),
        io_cluster->getWidth(),
        io_cluster->getHeight(),
        io_cluster);
  }

  // Write the connections between macros
  std::ofstream file;
  std::string file_name = parent->getName();
  file_name = make_filename(file_name);

  file_name = report_directory_ + "/" + file_name;
  file.open(file_name + "net.txt");
  for (auto& net : nets) {
    file << macros[net.terminals.first].getName() << "   "
         << macros[net.terminals.second].getName() << "   " << net.weight
         << std::endl;
  }
  file.close();

  // Call Simulated Annealing Engine to place children
  // set the action probabilities
  // the summation of probabilities should be one.
  const float action_sum = pos_swap_prob_ + neg_swap_prob_ + double_swap_prob_
                           + exchange_swap_prob_ + resize_prob_;
  // In our implementation, target_util and target_dead_space are different
  // target_util is used to determine the utilization for MixedCluster
  // target_dead_space is used to determine the utilization for
  // StandardCellCluster We vary the target utilization to generate different
  // tilings
  std::vector<float> target_utils{target_util_};
  std::vector<float> target_dead_spaces{target_dead_space_};
  // In our implementation, the utilization can be larger than 1.
  for (int i = 1; i < num_target_util_; i++) {
    target_utils.push_back(target_util_ + i * target_util_step_);
  }
  // In our implementation, the target_dead_space should be less than 1.0.
  // The larger the target dead space, the higher the utilization.
  for (int i = 1; i < num_target_dead_space_; i++) {
    if (target_dead_space_ + i * target_dead_space_step_ < 1.0) {
      target_dead_spaces.push_back(target_dead_space_
                                   + i * target_dead_space_step_);
    }
  }
  // Since target_util and target_dead_space are independent variables
  // the combination should be (target_util, target_dead_space_list)
  // target util has higher priority than target_dead_space
  std::vector<float> target_util_list;
  std::vector<float> target_dead_space_list;
  for (auto& target_util : target_utils) {
    for (auto& target_dead_space : target_dead_spaces) {
      target_util_list.push_back(target_util);
      target_dead_space_list.push_back(target_dead_space);
    }
  }
  // The number of perturbations in each step should be larger than the
  // number of macros
  const int num_perturb_per_step = (macros.size() > num_perturb_per_step_)
                                       ? macros.size()
                                       : num_perturb_per_step_;
  int remaining_runs = target_util_list.size();
  int run_id = 0;
  SACoreSoftMacro* best_sa = nullptr;
  SoftSAVector sa_containers;  // The owner of SACores objects.
  float best_cost = std::numeric_limits<float>::max();
  // To give consistency across threads we check the solutions
  // at a fixed interval independent of how many threads we are using.
  const int check_interval = 10;
  int begin_check = 0;
  int end_check = std::min(check_interval, remaining_runs);
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Start Simulated Annealing Core");
  while (remaining_runs > 0) {
    SoftSAVector sa_batch;
    const int run_thread
        = graphics_ ? 1 : std::min(remaining_runs, num_threads_);
    for (int i = 0; i < run_thread; i++) {
      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 1,
                 "Start Simulated Annealing (run_id = {})",
                 run_id);

      std::vector<SoftMacro> shaped_macros = macros;  // copy for multithread

      const float target_util = target_util_list[run_id];
      const float target_dead_space = target_dead_space_list[run_id++];

      debugPrint(logger_,
                 MPL,
                 "fine_shaping",
                 1,
                 "Starting adjusting shapes for children of {}. target_util = "
                 "{}, target_dead_space = {}",
                 parent->getName(),
                 target_util,
                 target_dead_space);

      if (!runFineShaping(parent,
                          shaped_macros,
                          soft_macro_id_map,
                          target_util,
                          target_dead_space)) {
        debugPrint(logger_,
                   MPL,
                   "fine_shaping",
                   1,
                   "Cannot generate feasible shapes for children of {}, sa_id: "
                   "{}, target_util: {}, target_dead_space: {}",
                   parent->getName(),
                   run_id,
                   target_util,
                   target_dead_space);
        continue;
      }
      debugPrint(logger_,
                 MPL,
                 "fine_shaping",
                 1,
                 "Finished adjusting shapes for children of cluster {}",
                 parent->getName());
      // Note that all the probabilities are normalized to the summation of 1.0.
      // Note that the weight are not necessaries summarized to 1.0, i.e., not
      // normalized.
      std::unique_ptr<SACoreSoftMacro> sa
          = std::make_unique<SACoreSoftMacro>(tree_->root.get(),
                                              outline,
                                              shaped_macros,
                                              area_weight_,
                                              outline_weight_,
                                              wirelength_weight_,
                                              guidance_weight_,
                                              fence_weight_,
                                              boundary_weight_,
                                              macro_blockage_weight_,
                                              notch_weight_,
                                              notch_h_th_,
                                              notch_v_th_,
                                              pos_swap_prob_ / action_sum,
                                              neg_swap_prob_ / action_sum,
                                              double_swap_prob_ / action_sum,
                                              exchange_swap_prob_ / action_sum,
                                              resize_prob_ / action_sum,
                                              init_prob_,
                                              max_num_step_,
                                              num_perturb_per_step,
                                              random_seed_,
                                              graphics_.get(),
                                              logger_);
      sa->setNumberOfMacrosToPlace(num_of_macros_to_place);
      sa->setCentralizationAttemptOn(true);
      sa->setFences(fences);
      sa->setGuides(guides);
      sa->setNets(nets);
      sa->addBlockages(placement_blockages);
      sa->addBlockages(macro_blockages);
      sa_batch.push_back(std::move(sa));
    }
    if (sa_batch.size() == 1) {
      runSA<SACoreSoftMacro>(sa_batch[0].get());
    } else {
      // multi threads
      std::vector<std::thread> threads;
      threads.reserve(sa_batch.size());
      for (auto& sa : sa_batch) {
        threads.emplace_back(runSA<SACoreSoftMacro>, sa.get());
      }
      for (auto& th : threads) {
        th.join();
      }
    }
    remaining_runs -= run_thread;
    // add macro tilings
    for (auto& sa : sa_batch) {
      sa_containers.push_back(std::move(sa));
    }
    while (sa_containers.size() >= end_check) {
      while (begin_check < end_check) {
        auto& sa = sa_containers[begin_check];
        if (sa->isValid() && sa->getNormCost() < best_cost) {
          best_cost = sa->getNormCost();
          best_sa = sa.get();
        }
        ++begin_check;
      }
      // add early stop mechanism
      if (best_sa || remaining_runs == 0) {
        break;
      }
      end_check = begin_check + std::min(check_interval, remaining_runs);
    }
    if (best_sa) {
      break;
    }
  }
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Finished Simulated Annealing Core");

  if (best_sa == nullptr) {
    debugPrint(logger_,
               MPL,
               "hierarchical_macro_placement",
               1,
               "SA Summary for cluster {}",
               parent->getName());

    for (auto i = 0; i < sa_containers.size(); i++) {
      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 1,
                 "sa_id: {}, target_util: {}, target_dead_space: {}",
                 i,
                 target_util_list[i],
                 target_dead_space_list[i]);

      sa_containers[i]->printResults();
    }

    logger_->error(MPL, 5, "Failed on cluster {}", parent->getName());
  }

  if (best_sa->centralizationWasReverted()) {
    best_sa->alignMacroClusters();
  }
  best_sa->fillDeadSpace();

  // update the clusters and do bus planning
  std::vector<SoftMacro> shaped_macros;
  best_sa->getMacros(shaped_macros);
  file.open(file_name + ".fp.txt.temp");
  for (auto& macro : shaped_macros) {
    file << macro.getName() << "   " << macro.getX() << "   " << macro.getY()
         << "   " << macro.getWidth() << "   " << macro.getHeight()
         << std::endl;
  }
  file.close();

  if (parent->getParent() != nullptr) {
    // ***********************************************************************
    // Now convert the area occupied by pin access macros to hard macro
    // blockages. Note at this stage, the size of each pin access macro is 0.0.
    // But you cannot remove the pin access macros. You still need these pin
    // access macros to maintain the connections
    // ***********************************************************************
    if (shaped_macros[soft_macro_id_map[toString(L)]].getWidth() > 0.0) {
      const float l_ly = shaped_macros[soft_macro_id_map[toString(L)]].getY();
      const float l_width
          = shaped_macros[soft_macro_id_map[toString(L)]].getWidth();
      const float l_height
          = shaped_macros[soft_macro_id_map[toString(L)]].getHeight();
      macro_blockages.emplace_back(0.0, l_ly, l_width, l_ly + l_height);
      shaped_macros[soft_macro_id_map[toString(L)]]
          = SoftMacro(std::pair<float, float>(0.0, l_ly),
                      toString(L),
                      0.0,
                      l_height,
                      nullptr);
    }
    if (shaped_macros[soft_macro_id_map[toString(R)]].getWidth() > 0.0) {
      const float r_ly = shaped_macros[soft_macro_id_map[toString(R)]].getY();
      const float r_width
          = shaped_macros[soft_macro_id_map[toString(R)]].getWidth();
      const float r_height
          = shaped_macros[soft_macro_id_map[toString(R)]].getHeight();
      macro_blockages.emplace_back(outline.getWidth() - r_width,
                                   r_ly,
                                   outline.getWidth(),
                                   r_ly + r_height);
      shaped_macros[soft_macro_id_map[toString(R)]]
          = SoftMacro(std::pair<float, float>(outline.getWidth(), r_ly),
                      toString(R),
                      0.0,
                      r_height,
                      nullptr);
    }
    if (shaped_macros[soft_macro_id_map[toString(T)]].getWidth() > 0.0) {
      const float t_lx = shaped_macros[soft_macro_id_map[toString(T)]].getX();
      const float t_width
          = shaped_macros[soft_macro_id_map[toString(T)]].getWidth();
      const float t_height
          = shaped_macros[soft_macro_id_map[toString(T)]].getHeight();
      macro_blockages.emplace_back(t_lx,
                                   outline.getHeight() - t_height,
                                   t_lx + t_width,
                                   outline.getHeight());
      shaped_macros[soft_macro_id_map[toString(T)]]
          = SoftMacro(std::pair<float, float>(t_lx, outline.getHeight()),
                      toString(T),
                      t_width,
                      0.0,
                      nullptr);
    }
    if (shaped_macros[soft_macro_id_map[toString(B)]].getWidth() > 0.0) {
      const float b_lx = shaped_macros[soft_macro_id_map[toString(B)]].getX();
      const float b_width
          = shaped_macros[soft_macro_id_map[toString(B)]].getWidth();
      const float b_height
          = shaped_macros[soft_macro_id_map[toString(B)]].getHeight();
      macro_blockages.emplace_back(b_lx, 0.0, b_lx + b_width, b_height);
      shaped_macros[soft_macro_id_map[toString(B)]]
          = SoftMacro(std::pair<float, float>(b_lx, 0.0),
                      toString(B),
                      b_width,
                      0.0,
                      nullptr);
    }

    // Exclude the pin access macros from the sequence pair now that
    // they were converted to macro blockages.
    num_of_macros_to_place -= number_of_pin_access;

    macros = shaped_macros;
    remaining_runs = target_util_list.size();
    run_id = 0;
    best_sa = nullptr;
    sa_containers.clear();  // Destroy SACores used for pin access.
    best_cost = std::numeric_limits<float>::max();
    begin_check = 0;
    end_check = std::min(check_interval, remaining_runs);
    debugPrint(logger_,
               MPL,
               "hierarchical_macro_placement",
               1,
               "Start Simulated Annealing Core");
    while (remaining_runs > 0) {
      SoftSAVector sa_batch;
      const int run_thread
          = graphics_ ? 1 : std::min(remaining_runs, num_threads_);
      for (int i = 0; i < run_thread; i++) {
        debugPrint(logger_,
                   MPL,
                   "hierarchical_macro_placement",
                   1,
                   "Start Simulated Annealing (run_id = {})",
                   run_id);

        std::vector<SoftMacro> shaped_macros = macros;  // copy for multithread

        const float target_util = target_util_list[run_id];
        const float target_dead_space = target_dead_space_list[run_id++];

        debugPrint(logger_,
                   MPL,
                   "fine_shaping",
                   1,
                   "Starting adjusting shapes for children of {}. target_util "
                   "= {}, target_dead_space = {}",
                   parent->getName(),
                   target_util,
                   target_dead_space);

        if (!runFineShaping(parent,
                            shaped_macros,
                            soft_macro_id_map,
                            target_util,
                            target_dead_space)) {
          debugPrint(
              logger_,
              MPL,
              "fine_shaping",
              1,
              "Cannot generate feasible shapes for children of {}, sa_id: "
              "{}, target_util: {}, target_dead_space: {}",
              parent->getName(),
              run_id,
              target_util,
              target_dead_space);
          continue;
        }
        debugPrint(logger_,
                   MPL,
                   "fine_shaping",
                   1,
                   "Finished adjusting shapes for children of cluster {}",
                   parent->getName());
        // Note that all the probabilities are normalized to the summation
        // of 1.0. Note that the weight are not necessaries summarized to 1.0,
        // i.e., not normalized.
        std::unique_ptr<SACoreSoftMacro> sa = std::make_unique<SACoreSoftMacro>(
            tree_->root.get(),
            outline,
            shaped_macros,
            area_weight_,
            outline_weight_,
            wirelength_weight_,
            guidance_weight_,
            fence_weight_,
            boundary_weight_,
            macro_blockage_weight_,
            notch_weight_,
            notch_h_th_,
            notch_v_th_,
            pos_swap_prob_ / action_sum,
            neg_swap_prob_ / action_sum,
            double_swap_prob_ / action_sum,
            exchange_swap_prob_ / action_sum,
            resize_prob_ / action_sum,
            init_prob_,
            max_num_step_,
            num_perturb_per_step,
            random_seed_,
            graphics_.get(),
            logger_);
        sa->setNumberOfMacrosToPlace(num_of_macros_to_place);
        sa->setCentralizationAttemptOn(true);
        sa->setFences(fences);
        sa->setGuides(guides);
        sa->setNets(nets);
        sa->addBlockages(placement_blockages);
        sa->addBlockages(macro_blockages);
        sa_batch.push_back(std::move(sa));
      }
      if (sa_batch.size() == 1) {
        runSA<SACoreSoftMacro>(sa_batch[0].get());
      } else {
        // multi threads
        std::vector<std::thread> threads;
        threads.reserve(sa_batch.size());
        for (auto& sa : sa_batch) {
          threads.emplace_back(runSA<SACoreSoftMacro>, sa.get());
        }
        for (auto& th : threads) {
          th.join();
        }
      }
      remaining_runs -= run_thread;
      // add macro tilings
      for (auto& sa : sa_batch) {
        sa_containers.push_back(std::move(sa));
      }
      while (sa_containers.size() >= end_check) {
        while (begin_check < end_check) {
          auto& sa = sa_containers[begin_check];
          if (sa->isValid() && sa->getNormCost() < best_cost) {
            best_cost = sa->getNormCost();
            best_sa = sa.get();
          }
          ++begin_check;
        }
        // add early stop mechanism
        if (best_sa || remaining_runs == 0) {
          break;
        }
        end_check = begin_check + std::min(check_interval, remaining_runs);
      }
      if (best_sa) {
        break;
      }
    }
    debugPrint(logger_,
               MPL,
               "hierarchical_macro_placement",
               1,
               "Finished Simulated Annealing Core");
    if (best_sa == nullptr) {
      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 1,
                 "SA Summary for cluster {}",
                 parent->getName());

      for (auto i = 0; i < sa_containers.size(); i++) {
        debugPrint(logger_,
                   MPL,
                   "hierarchical_macro_placement",
                   1,
                   "sa_id: {}, target_util: {}, target_dead_space: {}",
                   i,
                   target_util_list[i],
                   target_dead_space_list[i]);

        sa_containers[i]->printResults();
      }

      logger_->error(MPL, 6, "Failed on cluster {}", parent->getName());
    }

    if (best_sa->centralizationWasReverted()) {
      best_sa->alignMacroClusters();
    }
    best_sa->fillDeadSpace();

    // update the clusters and do bus planning
    best_sa->getMacros(shaped_macros);
  }

  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Finished Simulated Annealing for cluster: {}\n",
             parent->getName());
  best_sa->printResults();
  // write the cost function. This can be used to tune the temperature schedule
  // and cost weight
  best_sa->writeCostFile(file_name + ".cost.txt");
  // write the floorplan information
  file.open(file_name + ".fp.txt");
  for (auto& macro : shaped_macros) {
    file << macro.getName() << "   " << macro.getX() << "   " << macro.getY()
         << "   " << macro.getWidth() << "   " << macro.getHeight()
         << std::endl;
  }
  file.close();

  updateChildrenShapesAndLocations(parent, shaped_macros, soft_macro_id_map);

  // check if the parent cluster still need bus planning
  for (auto& child : parent->getChildren()) {
    if (child->getClusterType() == MixedCluster) {
      debugPrint(logger_,
                 MPL,
                 "bus_planning",
                 1,
                 "Calling bus planning for cluster {}",
                 parent->getName());
      callBusPlanning(shaped_macros, nets);
      break;
    }
  }

  updateChildrenRealLocation(parent, outline.xMin(), outline.yMin());
  sa_containers.clear();

  // Continue cluster placement on children
  for (auto& cluster : parent->getChildren()) {
    runHierarchicalMacroPlacement(cluster.get());
  }

  clustering_engine_->updateInstancesAssociation(parent);
}

// Merge nets to reduce runtime
void HierRTLMP::mergeNets(std::vector<BundledNet>& nets)
{
  std::vector<int> net_class(nets.size(), -1);
  for (int i = 0; i < nets.size(); i++) {
    if (net_class[i] == -1) {
      for (int j = i + 1; j < nets.size(); j++) {
        if (nets[i] == nets[j]) {
          net_class[j] = i;
        }
      }
    }
  }

  for (int i = 0; i < net_class.size(); i++) {
    if (net_class[i] > -1) {
      nets[net_class[i]].weight += nets[i].weight;
    }
  }

  std::vector<BundledNet> merged_nets;
  for (int i = 0; i < net_class.size(); i++) {
    if (net_class[i] == -1) {
      merged_nets.push_back(nets[i]);
    }
  }
  nets.clear();
  nets = std::move(merged_nets);

  if (graphics_) {
    graphics_->setBundledNets(nets);
  }
}

// Recommendation from the original implementation:
// For single level, increase macro blockage weight to
// half of the outline weight.
void HierRTLMP::adjustMacroBlockageWeight()
{
  if (tree_->max_level == 1) {
    float new_macro_blockage_weight = outline_weight_ / 2.0;
    debugPrint(logger_,
               MPL,
               "hierarchical_macro_placement",
               1,
               "Tree max level is {}, Changing macro blockage weight from {} "
               "to {} (half of the outline weight)",
               tree_->max_level,
               macro_blockage_weight_,
               new_macro_blockage_weight);
    macro_blockage_weight_ = new_macro_blockage_weight;
  }
}

void HierRTLMP::reportSAWeights()
{
  logger_->report("\nSimmulated Annealing Weights:\n");
  logger_->report("Area = {}", area_weight_);
  logger_->report("Outline = {}", outline_weight_);
  logger_->report("WL = {}", wirelength_weight_);
  logger_->report("Guidance = {}", guidance_weight_);
  logger_->report("Fence = {}", fence_weight_);
  logger_->report("Boundary = {}", boundary_weight_);
  logger_->report("Notch = {}", notch_weight_);
  logger_->report("Macro Blockage = {}\n", macro_blockage_weight_);
}

void HierRTLMP::runHierarchicalMacroPlacementWithoutBusPlanning(Cluster* parent)
{
  if (parent->getClusterType() == HardMacroCluster) {
    placeMacros(parent);
    return;
  }

  if (parent->isLeaf()) {  // Cover IO Clusters && Leaf Std Cells
    return;
  }

  if (graphics_) {
    graphics_->setCurrentCluster(parent);
  }

  for (auto& cluster : parent->getChildren()) {
    clustering_engine_->updateInstancesAssociation(cluster.get());
  }
  // Place children clusters
  // map children cluster to soft macro
  for (auto& cluster : parent->getChildren()) {
    if (cluster->isIOCluster()) {  // ignore all the io clusters
      continue;
    }
    auto soft_macro = std::make_unique<SoftMacro>(cluster.get());
    // no memory leakage, beacuse we set the soft macro, the old one
    // will be deleted
    cluster->setSoftMacro(std::move(soft_macro));
  }

  // The simulated annealing outline is determined by the parent's shape
  const Rect outline(parent->getX(),
                     parent->getY(),
                     parent->getX() + parent->getWidth(),
                     parent->getY() + parent->getHeight());

  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Working on children of cluster: {}, Outline "
             "{}, {}  {}, {}",
             parent->getName(),
             outline.xMin(),
             outline.yMin(),
             outline.getWidth(),
             outline.getHeight());

  // Suppose the region, fence, guide has been mapped to cooresponding macros
  // This step is done when we enter the Hier-RTLMP program
  std::map<std::string, int> soft_macro_id_map;  // cluster_name, macro_id
  std::map<int, Rect> fences;
  std::map<int, Rect> guides;
  std::vector<SoftMacro> macros;
  std::vector<BundledNet> nets;

  std::vector<Rect> placement_blockages;
  std::vector<Rect> macro_blockages;

  findOverlappingBlockages(macro_blockages, placement_blockages, outline);

  // We store the bundled io clusters to push them into the macros' vector
  // only after it is already populated with the clusters we're trying to
  // place. This will facilitate how we deal with fixed terminals in SA moves.
  std::vector<Cluster*> io_clusters;

  // Each cluster is modeled as Soft Macro
  // The fences or guides for each cluster is created by merging
  // the fences and guides for hard macros in each cluster
  for (auto& cluster : parent->getChildren()) {
    if (cluster->isIOCluster()) {
      io_clusters.push_back(cluster.get());
      continue;
    }
    // for other clusters
    soft_macro_id_map[cluster->getName()] = macros.size();
    auto soft_macro = std::make_unique<SoftMacro>(cluster.get());
    clustering_engine_->updateInstancesAssociation(
        cluster.get());  // we need this step to calculate nets
    macros.push_back(*soft_macro);
    cluster->setSoftMacro(std::move(soft_macro));
    // merge fences and guides for hard macros within cluster
    if (cluster->getClusterType() == StdCellCluster) {
      continue;
    }
    Rect fence, guide;
    fence.mergeInit();
    guide.mergeInit();
    const std::vector<HardMacro*> hard_macros = cluster->getHardMacros();
    for (auto& hard_macro : hard_macros) {
      if (fences_.find(hard_macro->getName()) != fences_.end()) {
        fence.merge(fences_[hard_macro->getName()]);
      }
      auto itr = guides_.find(hard_macro->getInst());
      if (itr != guides_.end()) {
        guide.merge(itr->second);
      }
    }

    // Calculate overlap with outline
    fence.relocate(
        outline.xMin(), outline.yMin(), outline.xMax(), outline.yMax());
    guide.relocate(
        outline.xMin(), outline.yMin(), outline.xMax(), outline.yMax());

    if (fence.isValid()) {
      // current macro id is macros.size() - 1
      fences[macros.size() - 1] = fence;
    }
    if (guide.isValid()) {
      // current macro id is macros.size() - 1
      guides[macros.size() - 1] = guide;
    }
  }

  if (graphics_) {
    graphics_->setGuides(guides);
  }

  const int num_of_macros_to_place = static_cast<int>(macros.size());

  for (Cluster* io_cluster : io_clusters) {
    soft_macro_id_map[io_cluster->getName()] = macros.size();

    macros.emplace_back(
        std::pair<float, float>(io_cluster->getX() - outline.xMin(),
                                io_cluster->getY() - outline.yMin()),
        io_cluster->getName(),
        io_cluster->getWidth(),
        io_cluster->getHeight(),
        io_cluster);
  }

  // model other clusters as fixed terminals
  if (parent->getParent() != nullptr) {
    std::queue<Cluster*> parents;
    parents.push(parent);
    while (parents.empty() == false) {
      auto frontwave = parents.front();
      parents.pop();
      for (auto& cluster : frontwave->getParent()->getChildren()) {
        if (cluster->getId() != frontwave->getId()) {
          // model this as a fixed softmacro
          soft_macro_id_map[cluster->getName()] = macros.size();
          macros.emplace_back(
              std::pair<float, float>(
                  cluster->getX() + cluster->getWidth() / 2.0 - outline.xMin(),
                  cluster->getY() + cluster->getHeight() / 2.0 - outline.yMin()),
              cluster->getName(),
              0.0,
              0.0,
              nullptr);
          debugPrint(
              logger_,
              MPL,
              "hierarchical_macro_placement",
              1,
              "fixed cluster : {}, lx = {}, ly = {}, width = {}, height = {}",
              cluster->getName(),
              cluster->getX(),
              cluster->getY(),
              cluster->getWidth(),
              cluster->getHeight());
        }
      }
      if (frontwave->getParent()->getParent() != nullptr) {
        parents.push(frontwave->getParent());
      }
    }
  }

  // update the connnection
  clustering_engine_->updateConnections();
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Finished calculating connection");
  clustering_engine_->updateDataFlow();
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Finished updating dataflow");

  // add the virtual connections (the weight related to IOs and macros belong to
  // the same cluster)
  for (const auto& [cluster1, cluster2] : parent->getVirtualConnections()) {
    BundledNet net(
        soft_macro_id_map[tree_->maps.id_to_cluster[cluster1]->getName()],
        soft_macro_id_map[tree_->maps.id_to_cluster[cluster2]->getName()],
        tree_->virtual_weight);
    net.src_cluster_id = cluster1;
    net.target_cluster_id = cluster2;
    nets.push_back(net);
  }

  // convert the connections between clusters to SoftMacros
  for (auto& cluster : parent->getChildren()) {
    const int src_id = cluster->getId();
    const std::string src_name = cluster->getName();
    for (auto& [cluster_id, weight] : cluster->getConnection()) {
      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 2,
                 " Cluster connection: {} {} {} ",
                 cluster->getName(),
                 tree_->maps.id_to_cluster[cluster_id]->getName(),
                 weight);
      const std::string name = tree_->maps.id_to_cluster[cluster_id]->getName();
      if (src_id > cluster_id) {
        BundledNet net(
            soft_macro_id_map[src_name], soft_macro_id_map[name], weight);
        net.src_cluster_id = src_id;
        net.target_cluster_id = cluster_id;
        nets.push_back(net);
      }
    }
  }
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Finished creating bundled connections");
  // merge nets to reduce runtime
  mergeNets(nets);

  // Write the connections between macros
  std::ofstream file;
  std::string file_name = parent->getName();
  file_name = make_filename(file_name);

  file_name = report_directory_ + "/" + file_name;
  file.open(file_name + ".net.txt");
  for (auto& net : nets) {
    file << macros[net.terminals.first].getName() << "   "
         << macros[net.terminals.second].getName() << "   " << net.weight
         << std::endl;
  }
  file.close();

  // Call Simulated Annealing Engine to place children
  // set the action probabilities
  // the summation of probabilities should be one.
  const float action_sum = pos_swap_prob_ + neg_swap_prob_ + double_swap_prob_
                           + exchange_swap_prob_ + resize_prob_;
  // In our implementation, target_util and target_dead_space are different
  // target_util is used to determine the utilization for MixedCluster
  // target_dead_space is used to determine the utilization for
  // StandardCellCluster We vary the target utilization to generate different
  // tilings
  std::vector<float> target_utils{target_util_};
  std::vector<float> target_dead_spaces{target_dead_space_};
  // In our implementation, the utilization can be larger than 1.
  for (int i = 1; i < num_target_util_; i++) {
    target_utils.push_back(target_util_ + i * target_util_step_);
  }
  // In our implementation, the target_dead_space should be less than 1.0.
  // The larger the target dead space, the higher the utilization.
  for (int i = 1; i < num_target_dead_space_; i++) {
    if (target_dead_space_ + i * target_dead_space_step_ < 1.0) {
      target_dead_spaces.push_back(target_dead_space_
                                   + i * target_dead_space_step_);
    }
  }
  // Since target_util and target_dead_space are independent variables
  // the combination should be (target_util, target_dead_space_list)
  // target util has higher priority than target_dead_space
  std::vector<float> target_util_list;
  std::vector<float> target_dead_space_list;
  for (auto& target_util : target_utils) {
    for (auto& target_dead_space : target_dead_spaces) {
      target_util_list.push_back(target_util);
      target_dead_space_list.push_back(target_dead_space);
    }
  }
  // The number of perturbations in each step should be larger than the
  // number of macros
  const int num_perturb_per_step = (macros.size() > num_perturb_per_step_)
                                       ? macros.size()
                                       : num_perturb_per_step_;
  int remaining_runs = target_util_list.size();
  int run_id = 0;
  SACoreSoftMacro* best_sa = nullptr;
  SoftSAVector sa_containers;  // The owner of SACore objects
  // To give consistency across threads we check the solutions
  // at a fixed interval independent of how many threads we are using.
  const int check_interval = 10;
  int begin_check = 0;
  int end_check = std::min(check_interval, remaining_runs);
  float best_cost = std::numeric_limits<float>::max();
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Start Simulated Annealing Core");
  while (remaining_runs > 0) {
    SoftSAVector sa_batch;
    const int run_thread
        = graphics_ ? 1 : std::min(remaining_runs, num_threads_);
    for (int i = 0; i < run_thread; i++) {
      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 1,
                 "Start Simulated Annealing (run_id = {})",
                 run_id);

      std::vector<SoftMacro> shaped_macros = macros;  // copy for multithread

      const float target_util = target_util_list[run_id];
      const float target_dead_space = target_dead_space_list[run_id++];

      debugPrint(logger_,
                 MPL,
                 "fine_shaping",
                 1,
                 "Starting adjusting shapes for children of {}. target_util = "
                 "{}, target_dead_space = {}",
                 parent->getName(),
                 target_util,
                 target_dead_space);

      if (!runFineShaping(parent,
                          shaped_macros,
                          soft_macro_id_map,
                          target_util,
                          target_dead_space)) {
        debugPrint(logger_,
                   MPL,
                   "fine_shaping",
                   1,
                   "Cannot generate feasible shapes for children of {}, sa_id: "
                   "{}, target_util: {}, target_dead_space: {}",
                   parent->getName(),
                   run_id,
                   target_util,
                   target_dead_space);
        continue;
      }
      debugPrint(logger_,
                 MPL,
                 "fine_shaping",
                 1,
                 "Finished generating shapes for children of cluster {}",
                 parent->getName());
      // Note that all the probabilities are normalized to the summation of 1.0.
      // Note that the weight are not necessaries summarized to 1.0, i.e., not
      // normalized.
      std::unique_ptr<SACoreSoftMacro> sa
          = std::make_unique<SACoreSoftMacro>(tree_->root.get(),
                                              outline,
                                              shaped_macros,
                                              area_weight_,
                                              outline_weight_,
                                              wirelength_weight_,
                                              guidance_weight_,
                                              fence_weight_,
                                              boundary_weight_,
                                              macro_blockage_weight_,
                                              notch_weight_,
                                              notch_h_th_,
                                              notch_v_th_,
                                              pos_swap_prob_ / action_sum,
                                              neg_swap_prob_ / action_sum,
                                              double_swap_prob_ / action_sum,
                                              exchange_swap_prob_ / action_sum,
                                              resize_prob_ / action_sum,
                                              init_prob_,
                                              max_num_step_,
                                              num_perturb_per_step,
                                              random_seed_,
                                              graphics_.get(),
                                              logger_);
      sa->setNumberOfMacrosToPlace(num_of_macros_to_place);
      sa->setCentralizationAttemptOn(true);
      sa->setFences(fences);
      sa->setGuides(guides);
      sa->setNets(nets);
      sa->addBlockages(placement_blockages);
      sa->addBlockages(macro_blockages);
      sa_batch.push_back(std::move(sa));
    }
    if (sa_batch.size() == 1) {
      runSA<SACoreSoftMacro>(sa_batch[0].get());
    } else {
      // multi threads
      std::vector<std::thread> threads;
      threads.reserve(sa_batch.size());
      for (auto& sa : sa_batch) {
        threads.emplace_back(runSA<SACoreSoftMacro>, sa.get());
      }
      for (auto& th : threads) {
        th.join();
      }
    }
    remaining_runs -= run_thread;
    // add macro tilings
    for (auto& sa : sa_batch) {
      sa_containers.push_back(std::move(sa));
    }
    while (sa_containers.size() >= end_check) {
      while (begin_check < end_check) {
        auto& sa = sa_containers[begin_check];
        if (sa->isValid() && sa->getNormCost() < best_cost) {
          best_cost = sa->getNormCost();
          best_sa = sa.get();
        }
        ++begin_check;
      }
      // add early stop mechanism
      if (best_sa || remaining_runs == 0) {
        break;
      }
      end_check = begin_check + std::min(check_interval, remaining_runs);
    }
    if (best_sa) {
      break;
    }
  }
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Finished Simulated Annealing Core");
  if (best_sa == nullptr) {
    debugPrint(logger_,
               MPL,
               "hierarchical_macro_placement",
               1,
               "SA Summary for cluster {}",
               parent->getName());

    for (auto i = 0; i < sa_containers.size(); i++) {
      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 1,
                 "sa_id: {}, target_util: {}, target_dead_space: {}",
                 i,
                 target_util_list[i],
                 target_dead_space_list[i]);

      sa_containers[i]->printResults();
    }

    runEnhancedHierarchicalMacroPlacement(parent);
  } else {
    if (best_sa->centralizationWasReverted()) {
      best_sa->alignMacroClusters();
    }
    best_sa->fillDeadSpace();

    std::vector<SoftMacro> shaped_macros;
    best_sa->getMacros(shaped_macros);
    file.open(file_name + ".fp.txt.temp");
    for (auto& macro : shaped_macros) {
      file << macro.getName() << "   " << macro.getX() << "   " << macro.getY()
           << "   " << macro.getWidth() << "   " << macro.getHeight()
           << std::endl;
    }
    file.close();
    debugPrint(logger_,
               MPL,
               "hierarchical_macro_placement",
               1,
               "Finished Simulated Annealing for cluster {}",
               parent->getName());
    best_sa->printResults();
    // write the cost function. This can be used to tune the temperature
    // schedule and cost weight
    best_sa->writeCostFile(file_name + ".cost.txt");
    // write the floorplan information
    file.open(file_name + ".fp.txt");
    for (auto& macro : shaped_macros) {
      file << macro.getName() << "   " << macro.getX() << "   " << macro.getY()
           << "   " << macro.getWidth() << "   " << macro.getHeight()
           << std::endl;
    }
    file.close();

    updateChildrenShapesAndLocations(parent, shaped_macros, soft_macro_id_map);

    updateChildrenRealLocation(parent, outline.xMin(), outline.yMin());
  }

  // Continue cluster placement on children
  for (auto& cluster : parent->getChildren()) {
    runHierarchicalMacroPlacementWithoutBusPlanning(cluster.get());
  }

  clustering_engine_->updateInstancesAssociation(parent);
}

// This function is used in cases with very high density, in which it may
// be very hard to generate a valid tiling for the clusters.
// Here, we may want to try setting the area of all standard-cell clusters to 0.
// This should be only be used in mixed clusters.
void HierRTLMP::runEnhancedHierarchicalMacroPlacement(Cluster* parent)
{
  if (parent->getClusterType() != MixedCluster) {
    return;
  }

  // We only run this enhanced macro placement version if there are no
  // further levels ahead in the current branch of the physical hierarchy.
  for (auto& cluster : parent->getChildren()) {
    if (cluster->getClusterType() == MixedCluster) {
      return;
    }
  }

  if (graphics_) {
    graphics_->setCurrentCluster(parent);
  }

  // Place children clusters
  // map children cluster to soft macro
  for (auto& cluster : parent->getChildren()) {
    if (cluster->isIOCluster()) {  // ignore all the io clusters
      continue;
    }
    auto soft_macro = std::make_unique<SoftMacro>(cluster.get());
    // no memory leakage, beacuse we set the soft macro, the old one
    // will be deleted
    cluster->setSoftMacro(std::move(soft_macro));
  }

  // The simulated annealing outline is determined by the parent's shape
  const Rect outline(parent->getX(),
                     parent->getY(),
                     parent->getX() + parent->getWidth(),
                     parent->getY() + parent->getHeight());

  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Working on children of cluster: {}, Outline "
             "{}, {}  {}, {}",
             parent->getName(),
             outline.xMin(),
             outline.yMin(),
             outline.getWidth(),
             outline.getHeight());

  // Suppose the region, fence, guide has been mapped to cooresponding macros
  // This step is done when we enter the Hier-RTLMP program
  std::map<std::string, int> soft_macro_id_map;  // cluster_name, macro_id
  std::map<int, Rect> fences;
  std::map<int, Rect> guides;
  std::vector<SoftMacro> macros;
  std::vector<BundledNet> nets;

  std::vector<Rect> placement_blockages;
  std::vector<Rect> macro_blockages;

  findOverlappingBlockages(macro_blockages, placement_blockages, outline);

  // We store the bundled io clusters to push them into the macros' vector
  // only after it is already populated with the clusters we're trying to
  // place. This will facilitate how we deal with fixed terminals in SA moves.
  std::vector<Cluster*> io_clusters;

  // Each cluster is modeled as Soft Macro
  // The fences or guides for each cluster is created by merging
  // the fences and guides for hard macros in each cluster
  for (auto& cluster : parent->getChildren()) {
    if (cluster->isIOCluster()) {
      io_clusters.push_back(cluster.get());
      continue;
    }
    // for other clusters
    soft_macro_id_map[cluster->getName()] = macros.size();
    auto soft_macro = std::make_unique<SoftMacro>(cluster.get());
    clustering_engine_->updateInstancesAssociation(
        cluster.get());  // we need this step to calculate nets
    macros.push_back(*soft_macro);
    cluster->setSoftMacro(std::move(soft_macro));
    // merge fences and guides for hard macros within cluster
    if (cluster->getClusterType() == StdCellCluster) {
      continue;
    }
    Rect fence, guide;
    fence.mergeInit();
    guide.mergeInit();
    const std::vector<HardMacro*> hard_macros = cluster->getHardMacros();
    for (auto& hard_macro : hard_macros) {
      if (fences_.find(hard_macro->getName()) != fences_.end()) {
        fence.merge(fences_[hard_macro->getName()]);
      }

      auto itr = guides_.find(hard_macro->getInst());
      if (itr != guides_.end()) {
        guide.merge(itr->second);
      }
    }

    // Calculate overlap with outline
    fence.relocate(
        outline.xMin(), outline.yMin(), outline.xMax(), outline.yMax());
    guide.relocate(
        outline.xMin(), outline.yMin(), outline.xMax(), outline.yMax());
    if (fence.isValid()) {
      // current macro id is macros.size() - 1
      fences[macros.size() - 1] = fence;
    }
    if (guide.isValid()) {
      // current macro id is macros.size() - 1
      guides[macros.size() - 1] = guide;
    }
  }

  if (graphics_) {
    graphics_->setGuides(guides);
  }

  const int macros_to_place = static_cast<int>(macros.size());

  for (Cluster* io_cluster : io_clusters) {
    soft_macro_id_map[io_cluster->getName()] = macros.size();

    macros.emplace_back(
        std::pair<float, float>(io_cluster->getX() - outline.xMin(),
                                io_cluster->getY() - outline.yMin()),
        io_cluster->getName(),
        io_cluster->getWidth(),
        io_cluster->getHeight(),
        io_cluster);
  }

  // model other clusters as fixed terminals
  if (parent->getParent() != nullptr) {
    std::queue<Cluster*> parents;
    parents.push(parent);
    while (parents.empty() == false) {
      auto frontwave = parents.front();
      parents.pop();
      for (auto& cluster : frontwave->getParent()->getChildren()) {
        if (cluster->getId() != frontwave->getId()) {
          // model this as a fixed softmacro
          soft_macro_id_map[cluster->getName()] = macros.size();
          macros.emplace_back(
              std::pair<float, float>(
                  cluster->getX() + cluster->getWidth() / 2.0 - outline.xMin(),
                  cluster->getY() + cluster->getHeight() / 2.0
                      - outline.yMin()),
              cluster->getName(),
              0.0,
              0.0,
              nullptr);
          debugPrint(
              logger_,
              MPL,
              "hierarchical_macro_placement",
              1,
              "fixed cluster : {}, lx = {}, ly = {}, width = {}, height = {}",
              cluster->getName(),
              cluster->getX(),
              cluster->getY(),
              cluster->getWidth(),
              cluster->getHeight());
        }
      }
      if (frontwave->getParent()->getParent() != nullptr) {
        parents.push(frontwave->getParent());
      }
    }
  }

  // update the connnection
  clustering_engine_->updateConnections();
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Finished calculating connection");
  clustering_engine_->updateDataFlow();
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Finished updating dataflow");

  // add the virtual connections (the weight related to IOs and macros belong to
  // the same cluster)
  for (const auto& [cluster1, cluster2] : parent->getVirtualConnections()) {
    BundledNet net(
        soft_macro_id_map[tree_->maps.id_to_cluster[cluster1]->getName()],
        soft_macro_id_map[tree_->maps.id_to_cluster[cluster2]->getName()],
        tree_->virtual_weight);
    net.src_cluster_id = cluster1;
    net.target_cluster_id = cluster2;
    nets.push_back(net);
  }

  // convert the connections between clusters to SoftMacros
  for (auto& cluster : parent->getChildren()) {
    const int src_id = cluster->getId();
    const std::string src_name = cluster->getName();
    for (auto& [cluster_id, weight] : cluster->getConnection()) {
      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 2,
                 " Cluster connection: {} {} {} ",
                 cluster->getName(),
                 tree_->maps.id_to_cluster[cluster_id]->getName(),
                 weight);
      const std::string name = tree_->maps.id_to_cluster[cluster_id]->getName();
      if (src_id > cluster_id) {
        BundledNet net(
            soft_macro_id_map[src_name], soft_macro_id_map[name], weight);
        net.src_cluster_id = src_id;
        net.target_cluster_id = cluster_id;
        nets.push_back(net);
      }
    }
  }
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Finished creating bundled connections");
  // merge nets to reduce runtime
  mergeNets(nets);

  // Write the connections between macros
  std::ofstream file;
  std::string file_name = parent->getName();
  file_name = make_filename(file_name);

  file_name = report_directory_ + "/" + file_name;
  file.open(file_name + ".net.txt");
  for (auto& net : nets) {
    file << macros[net.terminals.first].getName() << "   "
         << macros[net.terminals.second].getName() << "   " << net.weight
         << std::endl;
  }
  file.close();

  // Call Simulated Annealing Engine to place children
  // set the action probabilities
  // the summation of probabilities should be one.
  const float action_sum = pos_swap_prob_ + neg_swap_prob_ + double_swap_prob_
                           + exchange_swap_prob_ + resize_prob_;
  // In our implementation, target_util and target_dead_space are different
  // target_util is used to determine the utilization for MixedCluster
  // target_dead_space is used to determine the utilization for
  // StandardCellCluster We vary the target utilization to generate different
  // tilings
  std::vector<float> target_utils{1e6};
  std::vector<float> target_dead_spaces{0.99999};
  // Since target_util and target_dead_space are independent variables
  // the combination should be (target_util, target_dead_space_list)
  // target util has higher priority than target_dead_space
  std::vector<float> target_util_list;
  std::vector<float> target_dead_space_list;
  for (auto& target_util : target_utils) {
    for (auto& target_dead_space : target_dead_spaces) {
      target_util_list.push_back(target_util);
      target_dead_space_list.push_back(target_dead_space);
    }
  }
  // The number of perturbations in each step should be larger than the
  // number of macros
  const int num_perturb_per_step = (macros.size() > num_perturb_per_step_)
                                       ? macros.size()
                                       : num_perturb_per_step_;
  int remaining_runs = target_util_list.size();
  int run_id = 0;
  SACoreSoftMacro* best_sa = nullptr;
  SoftSAVector sa_containers;  // The owner of SACore objects
  float best_cost = std::numeric_limits<float>::max();
  // To give consistency across threads we check the solutions
  // at a fixed interval independent of how many threads we are using.
  const int check_interval = 10;
  int begin_check = 0;
  int end_check = std::min(check_interval, remaining_runs);
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Start Simulated Annealing Core");
  while (remaining_runs > 0) {
    SoftSAVector sa_batch;
    const int run_thread
        = graphics_ ? 1 : std::min(remaining_runs, num_threads_);
    for (int i = 0; i < run_thread; i++) {
      std::vector<SoftMacro> shaped_macros = macros;  // copy for multithread
      // determine the shape for each macro
      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 1,
                 "Start Simulated Annealing (run_id = {})",
                 run_id);
      const float target_util = target_util_list[run_id];
      const float target_dead_space = target_dead_space_list[run_id++];
      debugPrint(logger_,
                 MPL,
                 "fine_shaping",
                 1,
                 "Starting adjusting shapes for children of {}. target_util = "
                 "{}, target_dead_space = {}",
                 parent->getName(),
                 target_util,
                 target_dead_space);
      if (!runFineShaping(parent,
                          shaped_macros,
                          soft_macro_id_map,
                          target_util,
                          target_dead_space)) {
        debugPrint(logger_,
                   MPL,
                   "fine_shaping",
                   1,
                   "Cannot generate feasible shapes for children of {}, sa_id: "
                   "{}, target_util: {}, target_dead_space: {}",
                   parent->getName(),
                   run_id,
                   target_util,
                   target_dead_space);
        continue;
      }
      debugPrint(logger_,
                 MPL,
                 "fine_shaping",
                 1,
                 "Finished generating shapes for children of cluster {}",
                 parent->getName());
      // Note that all the probabilities are normalized to the summation of 1.0.
      // Note that the weight are not necessaries summarized to 1.0, i.e., not
      // normalized.
      std::unique_ptr<SACoreSoftMacro> sa
          = std::make_unique<SACoreSoftMacro>(tree_->root.get(),
                                              outline,
                                              shaped_macros,
                                              area_weight_,
                                              outline_weight_,
                                              wirelength_weight_,
                                              guidance_weight_,
                                              fence_weight_,
                                              boundary_weight_,
                                              macro_blockage_weight_,
                                              notch_weight_,
                                              notch_h_th_,
                                              notch_v_th_,
                                              pos_swap_prob_ / action_sum,
                                              neg_swap_prob_ / action_sum,
                                              double_swap_prob_ / action_sum,
                                              exchange_swap_prob_ / action_sum,
                                              resize_prob_ / action_sum,
                                              init_prob_,
                                              max_num_step_,
                                              num_perturb_per_step,
                                              random_seed_,
                                              graphics_.get(),
                                              logger_);
      sa->setNumberOfMacrosToPlace(macros_to_place);
      sa->setCentralizationAttemptOn(true);
      sa->setFences(fences);
      sa->setGuides(guides);
      sa->setNets(nets);
      sa->addBlockages(placement_blockages);
      sa->addBlockages(macro_blockages);
      sa_batch.push_back(std::move(sa));
    }
    if (sa_batch.size() == 1) {
      runSA<SACoreSoftMacro>(sa_batch[0].get());
    } else {
      // multi threads
      std::vector<std::thread> threads;
      threads.reserve(sa_batch.size());
      for (auto& sa : sa_batch) {
        threads.emplace_back(runSA<SACoreSoftMacro>, sa.get());
      }
      for (auto& th : threads) {
        th.join();
      }
    }
    remaining_runs -= run_thread;
    // add macro tilings
    for (auto& sa : sa_batch) {
      sa_containers.push_back(std::move(sa));
    }
    while (sa_containers.size() >= end_check) {
      while (begin_check < end_check) {
        auto& sa = sa_containers[begin_check];
        if (sa->isValid() && sa->getNormCost() < best_cost) {
          best_cost = sa->getNormCost();
          best_sa = sa.get();
        }
        ++begin_check;
      }
      // add early stop mechanism
      if (best_sa || remaining_runs == 0) {
        break;
      }
      end_check = begin_check + std::min(check_interval, remaining_runs);
    }
    if (best_sa) {
      break;
    }
  }
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Finished Simulated Annealing Core");
  if (best_sa == nullptr) {
    debugPrint(logger_,
               MPL,
               "hierarchical_macro_placement",
               1,
               "SA Summary for cluster {}",
               parent->getName());

    for (auto i = 0; i < sa_containers.size(); i++) {
      debugPrint(logger_,
                 MPL,
                 "hierarchical_macro_placement",
                 1,
                 "sa_id: {}, target_util: {}, target_dead_space: {}",
                 i,
                 target_util_list[i],
                 target_dead_space_list[i]);

      sa_containers[i]->printResults();
    }
    logger_->error(MPL, 40, "Failed on cluster {}", parent->getName());
  }

  if (best_sa->centralizationWasReverted()) {
    best_sa->alignMacroClusters();
  }
  best_sa->fillDeadSpace();

  // update the clusters and do bus planning
  std::vector<SoftMacro> shaped_macros;
  best_sa->getMacros(shaped_macros);
  file.open(file_name + ".fp.txt.temp");
  for (auto& macro : shaped_macros) {
    file << macro.getName() << "   " << macro.getX() << "   " << macro.getY()
         << "   " << macro.getWidth() << "   " << macro.getHeight()
         << std::endl;
  }
  file.close();

  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Finish Simulated Annealing for cluster {}",
             parent->getName());

  best_sa->printResults();
  // write the cost function. This can be used to tune the temperature schedule
  // and cost weight
  best_sa->writeCostFile(file_name + ".cost.txt");
  // write the floorplan information
  file.open(file_name + ".fp.txt");
  for (auto& macro : shaped_macros) {
    file << macro.getName() << "   " << macro.getX() << "   " << macro.getY()
         << "   " << macro.getWidth() << "   " << macro.getHeight()
         << std::endl;
  }
  file.close();

  updateChildrenShapesAndLocations(parent, shaped_macros, soft_macro_id_map);

  updateChildrenRealLocation(parent, outline.xMin(), outline.yMin());
}

// Verify the blockages' areas that have overlapped with current parent
// cluster. All the blockages will be converted to hard macros with fences.
void HierRTLMP::findOverlappingBlockages(std::vector<Rect>& macro_blockages,
                                         std::vector<Rect>& placement_blockages,
                                         const Rect& outline)
{
  for (auto& blockage : placement_blockages_) {
    computeBlockageOverlap(placement_blockages, blockage, outline);
  }

  for (auto& blockage : macro_blockages_) {
    computeBlockageOverlap(macro_blockages, blockage, outline);
  }

  if (graphics_) {
    graphics_->setOutline(micronsToDbu(outline));
    graphics_->setMacroBlockages(macro_blockages);
    graphics_->setPlacementBlockages(placement_blockages);
  }
}

void HierRTLMP::computeBlockageOverlap(std::vector<Rect>& overlapping_blockages,
                                       const Rect& blockage,
                                       const Rect& outline)
{
  const float b_lx = std::max(outline.xMin(), blockage.xMin());
  const float b_ly = std::max(outline.yMin(), blockage.yMin());
  const float b_ux = std::min(outline.xMax(), blockage.xMax());
  const float b_uy = std::min(outline.yMax(), blockage.yMax());

  if ((b_ux - b_lx > 0.0) && (b_uy - b_ly > 0.0)) {
    overlapping_blockages.emplace_back(b_lx - outline.xMin(),
                                       b_ly - outline.yMin(),
                                       b_ux - outline.xMin(),
                                       b_uy - outline.yMin());
  }
}

// Determine the shape of each cluster based on target utilization
// and target dead space.  In constrast to all previous works, we
// use two parameters: target utilization, target_dead_space.
// This is the trick part.  During our experiements, we found that keeping
// the same utilization of standard-cell clusters and mixed cluster will make
// SA very difficult to find a feasible solution.  With different utilization,
// SA can more easily find the solution. In our method, the target_utilization
// is used to determine the bloating ratio for mixed cluster, the
// target_dead_space is used to determine the bloating ratio for standard-cell
// cluster. The target utilization is based on tiling results we calculated
// before. The tiling results (which only consider the contribution of hard
// macros) will give us very close starting point.
bool HierRTLMP::runFineShaping(Cluster* parent,
                               std::vector<SoftMacro>& macros,
                               std::map<std::string, int>& soft_macro_id_map,
                               float target_util,
                               float target_dead_space)
{
  const float outline_width = parent->getWidth();
  const float outline_height = parent->getHeight();
  float pin_access_area = 0.0;
  float std_cell_cluster_area = 0.0;
  float std_cell_mixed_cluster_area = 0.0;
  float macro_cluster_area = 0.0;
  float macro_mixed_cluster_area = 0.0;
  // add the macro area for blockages, pin access and so on
  for (auto& macro : macros) {
    if (macro.getCluster() == nullptr) {
      pin_access_area += macro.getArea();  // get the physical-only area
    }
  }

  for (auto& cluster : parent->getChildren()) {
    if (cluster->isIOCluster()) {
      continue;  // IO clusters have no area
    }
    if (cluster->getClusterType() == StdCellCluster) {
      std_cell_cluster_area += cluster->getStdCellArea();
    } else if (cluster->getClusterType() == HardMacroCluster) {
      std::vector<std::pair<float, float>> shapes;
      for (auto& shape : cluster->getMacroTilings()) {
        if (shape.first < outline_width * (1 + conversion_tolerance_)
            && shape.second < outline_height * (1 + conversion_tolerance_)) {
          shapes.push_back(shape);
        }
      }
      if (shapes.empty()) {
        logger_->error(MPL,
                       7,
                       "Not enough space in cluster: {} for "
                       "child hard macro cluster: {}",
                       parent->getName(),
                       cluster->getName());
      }
      macro_cluster_area += shapes[0].first * shapes[0].second;
      cluster->setMacroTilings(shapes);
    } else {  // mixed cluster
      std_cell_mixed_cluster_area += cluster->getStdCellArea();
      std::vector<std::pair<float, float>> shapes;
      for (auto& shape : cluster->getMacroTilings()) {
        if (shape.first < outline_width * (1 + conversion_tolerance_)
            && shape.second < outline_height * (1 + conversion_tolerance_)) {
          shapes.push_back(shape);
        }
      }
      if (shapes.empty()) {
        logger_->error(MPL,
                       8,
                       "Not enough space in cluster: {} for "
                       "child mixed cluster: {}",
                       parent->getName(),
                       cluster->getName());
      }
      macro_mixed_cluster_area += shapes[0].first * shapes[0].second;
      cluster->setMacroTilings(shapes);
    }  // end for cluster type
  }

  // check how much available space to inflate for mixed cluster
  const float min_target_util
      = std_cell_mixed_cluster_area
        / (outline_width * outline_height - pin_access_area - macro_cluster_area
           - macro_mixed_cluster_area);
  if (target_util <= min_target_util) {
    target_util = min_target_util;  // target utilization for standard cells in
                                    // mixed cluster
  }
  // calculate the std_cell_util
  const float avail_space
      = outline_width * outline_height
        - (pin_access_area + macro_cluster_area + macro_mixed_cluster_area
           + std_cell_mixed_cluster_area / target_util);
  const float std_cell_util
      = std_cell_cluster_area / (avail_space * (1 - target_dead_space));
  // shape clusters
  if ((std_cell_cluster_area > 0.0 && avail_space < 0.0)
      || (std_cell_mixed_cluster_area > 0.0 && min_target_util <= 0.0)) {
    debugPrint(logger_,
               MPL,
               "fine_shaping",
               1,
               "No valid solution for children of {}"
               "avail_space = {}, min_target_util = {}",
               parent->getName(),
               avail_space,
               min_target_util);

    return false;
  }

  // set the shape for each macro
  for (auto& cluster : parent->getChildren()) {
    if (cluster->isIOCluster()) {
      continue;  // the area of IO cluster is 0.0
    }
    if (cluster->getClusterType() == StdCellCluster) {
      float area = cluster->getArea();
      float width = std::sqrt(area);
      float height = width;
      const float dust_threshold
          = 1.0 / macros.size();  // check if the cluster is the dust cluster
      const int dust_std_cell = 100;
      if ((width / outline_width <= dust_threshold
           && height / outline_height <= dust_threshold)
          || cluster->getNumStdCell() <= dust_std_cell) {
        width = 1e-3;
        height = 1e-3;
        area = width * height;
      } else {
        area = cluster->getArea() / std_cell_util;
        width = std::sqrt(area / min_ar_);
      }
      std::vector<std::pair<float, float>> width_list;
      width_list.emplace_back(width, area / width);
      macros[soft_macro_id_map[cluster->getName()]].setShapes(width_list, area);
    } else if (cluster->getClusterType() == HardMacroCluster) {
      macros[soft_macro_id_map[cluster->getName()]].setShapes(
          cluster->getMacroTilings());
      debugPrint(logger_,
                 MPL,
                 "fine_shaping",
                 2,
                 "hard_macro_cluster : {}",
                 cluster->getName());
      for (auto& shape : cluster->getMacroTilings()) {
        debugPrint(logger_,
                   MPL,
                   "fine_shaping",
                   2,
                   "    ( {} , {} ) ",
                   shape.first,
                   shape.second);
      }
    } else {  // Mixed cluster
      const std::vector<std::pair<float, float>> tilings
          = cluster->getMacroTilings();
      std::vector<std::pair<float, float>> width_list;
      // use the largest area
      float area = tilings[tilings.size() - 1].first
                   * tilings[tilings.size() - 1].second;
      area += cluster->getStdCellArea() / target_util;
      for (auto& shape : tilings) {
        if (shape.first * shape.second <= area) {
          width_list.emplace_back(shape.first, area / shape.second);
        }
      }

      debugPrint(logger_,
                 MPL,
                 "fine_shaping",
                 2,
                 "name:  {} area: {}",
                 cluster->getName(),
                 area);
      debugPrint(logger_, MPL, "fine_shaping", 2, "width_list :  ");
      for (auto& width : width_list) {
        debugPrint(logger_,
                   MPL,
                   "fine_shaping",
                   2,
                   " [  {} {}  ] ",
                   width.first,
                   width.second);
      }
      macros[soft_macro_id_map[cluster->getName()]].setShapes(width_list, area);
    }
  }

  return true;
}

// Call Path Synthesis to route buses
void HierRTLMP::callBusPlanning(std::vector<SoftMacro>& shaped_macros,
                                std::vector<BundledNet>& nets_old)
{
  std::vector<BundledNet> nets;
  for (auto& net : nets_old) {
    debugPrint(logger_, MPL, "bus_planning", 1, "net weight: {}", net.weight);
    if (net.weight > bus_net_threshold_) {
      nets.push_back(net);
    }
  }

  std::vector<int> soft_macro_vertex_id;
  std::vector<Edge> edge_list;
  std::vector<Vertex> vertex_list;
  if (!calNetPaths(shaped_macros,
                   soft_macro_vertex_id,
                   edge_list,
                   vertex_list,
                   nets,
                   congestion_weight_,
                   logger_)) {
    logger_->error(MPL, 9, "Bus planning has failed!");
  }
}

void HierRTLMP::placeMacros(Cluster* cluster)
{
  debugPrint(logger_,
             MPL,
             "hierarchical_macro_placement",
             1,
             "Place macros in cluster: {}",
             cluster->getName());

  if (graphics_) {
    graphics_->setCurrentCluster(cluster);
  }

  UniqueClusterVector macro_clusters;  // needed to calculate connections
  std::vector<HardMacro*> hard_macros = cluster->getHardMacros();
  std::vector<HardMacro> sa_macros;
  std::map<int, int> cluster_to_macro;
  std::set<odb::dbMaster*> masters;
  clustering_engine_->createTempMacroClusters(
      hard_macros, sa_macros, macro_clusters, cluster_to_macro, masters);

  const Rect outline(cluster->getX(),
                     cluster->getY(),
                     cluster->getX() + cluster->getWidth(),
                     cluster->getY() + cluster->getHeight());

  std::map<int, Rect> fences;
  std::map<int, Rect> guides;
  computeFencesAndGuides(hard_macros, outline, fences, guides);
  if (graphics_) {
    graphics_->setGuides(guides);
  }

  clustering_engine_->updateConnections();

  createFixedTerminals(outline, macro_clusters, cluster_to_macro, sa_macros);

  std::vector<BundledNet> nets
      = computeBundledNets(macro_clusters, cluster_to_macro);

  if (graphics_) {
    graphics_->setBundledNets(nets);
  }

  // Use exchange more often when there are more instances of a common
  // master.
  float exchange_swap_prob
      = exchange_swap_prob_ * 5
        * (1 - (masters.size() / (float) hard_macros.size()));

  // set the action probabilities (summation to 1.0)
  const float action_sum = pos_swap_prob_ * 10 + neg_swap_prob_ * 10
                           + double_swap_prob_ + exchange_swap_prob
                           + flip_prob_;

  float pos_swap_prob = pos_swap_prob_ * 10 / action_sum;
  float neg_swap_prob = neg_swap_prob_ * 10 / action_sum;
  float double_swap_prob = double_swap_prob_ / action_sum;
  exchange_swap_prob = exchange_swap_prob / action_sum;
  float flip_prob = flip_prob_ / action_sum;

  const int macros_to_place = static_cast<int>(hard_macros.size());

  int num_perturb_per_step = (macros_to_place > num_perturb_per_step_ / 10)
                                 ? macros_to_place
                                 : num_perturb_per_step_ / 10;

  SequencePair initial_seq_pair;
  if (cluster->isArrayOfInterconnectedMacros()) {
    setArrayTilingSequencePair(cluster, macros_to_place, initial_seq_pair);

    pos_swap_prob = 0.0f;
    neg_swap_prob = 0.0f;
    double_swap_prob = 0.0f;
    exchange_swap_prob = 0.95;
    flip_prob = 0.05;

    // Large arrays need more steps to properly converge.
    if (num_perturb_per_step > macros_to_place) {
      num_perturb_per_step *= 2;
    }
  }

  int run_id = 0;
  int remaining_runs = num_runs_;
  float best_cost = std::numeric_limits<float>::max();

  SACoreHardMacro* best_sa = nullptr;
  HardSAVector sa_containers;  // The owner of the SACore objects.

  while (remaining_runs > 0) {
    HardSAVector sa_batch;
    const int run_thread
        = graphics_ ? 1 : std::min(remaining_runs, num_threads_);

    for (int i = 0; i < run_thread; i++) {
      if (graphics_) {
        graphics_->setOutline(micronsToDbu(outline));
      }

      std::unique_ptr<SACoreHardMacro> sa = std::make_unique<SACoreHardMacro>(
          outline,
          sa_macros,
          area_weight_,
          outline_weight_ * (run_id + 1) * 10,
          wirelength_weight_ / (run_id + 1),
          guidance_weight_,
          fence_weight_,
          pos_swap_prob,
          neg_swap_prob,
          double_swap_prob,
          exchange_swap_prob,
          flip_prob,
          init_prob_,
          max_num_step_,
          num_perturb_per_step,
          random_seed_ + run_id,
          graphics_.get(),
          logger_);
      sa->setNumberOfMacrosToPlace(macros_to_place);
      sa->setNets(nets);
      sa->setFences(fences);
      sa->setGuides(guides);
      sa->setInitialSequencePair(initial_seq_pair);

      sa_batch.push_back(std::move(sa));

      run_id++;
    }
    if (sa_batch.size() == 1) {
      runSA<SACoreHardMacro>(sa_batch[0].get());
    } else {
      std::vector<std::thread> threads;
      threads.reserve(sa_batch.size());

      for (auto& sa : sa_batch) {
        threads.emplace_back(runSA<SACoreHardMacro>, sa.get());
      }

      for (auto& th : threads) {
        th.join();
      }
    }

    for (auto& sa : sa_batch) {
      SACoreWeights weights;
      weights.area = area_weight_;
      weights.outline = outline_weight_;
      weights.wirelength = wirelength_weight_;
      weights.guidance = guidance_weight_;
      weights.fence = fence_weight_;

      // Reset weights so we can compare the final costs.
      sa->setWeights(weights);
      if (sa->isValid(outline) && sa->getNormCost() < best_cost) {
        best_cost = sa->getNormCost();
        best_sa = sa.get();
      }

      sa_containers.push_back(std::move(sa));
    }

    remaining_runs -= run_thread;
  }

  if (best_sa == nullptr) {
    for (auto& sa : sa_containers) {
      sa->printResults();
    }
    logger_->error(MPL,
                   10,
                   "Macro placement failed for macro cluster: {}",
                   cluster->getName());
  } else {
    std::vector<HardMacro> best_macros;
    best_sa->getMacros(best_macros);
    best_sa->printResults();

    for (int i = 0; i < hard_macros.size(); i++) {
      *(hard_macros[i]) = best_macros[i];
    }
  }

  for (auto& hard_macro : hard_macros) {
    hard_macro->setX(hard_macro->getX() + outline.xMin());
    hard_macro->setY(hard_macro->getY() + outline.yMin());
  }

  clustering_engine_->clearTempMacroClusterMapping(macro_clusters);
  clustering_engine_->updateInstancesAssociation(cluster);
}

// Suppose we have a 2x2 array such as:
//      +-----++-----+
//      |     ||     |
//      |  0  ||  2  |
//      |     ||     |
//      +-----++-----+
//      +-----++-----+
//      |     ||     |
//      |  1  ||  3  |
//      |     ||     |
//      +-----++-----+
//  We can represent it through a sequence pair in which:
// 1) The positive sequence is made by going each column top
//    to bottom starting from the upper left corner: 0 1 2 3
// 2) The negative sequence is made by going each column bottom
//    to top starting from the lower left corner: 1 0 3 2
//
// Obs: Doing this will still keep the IO clusters at the end
// of the sequence pair, which is needed for the way SA handles it.
void HierRTLMP::setArrayTilingSequencePair(Cluster* cluster,
                                           const int macros_to_place,
                                           SequencePair& initial_seq_pair)
{
  // Set positive sequence
  for (int i = 0; i < macros_to_place; ++i) {
    initial_seq_pair.pos_sequence.push_back(i);
  }

  // Set negative sequence
  const int columns
      = cluster->getWidth() / cluster->getHardMacros().front()->getWidth();
  const int rows
      = cluster->getHeight() / cluster->getHardMacros().front()->getHeight();

  for (int i = 1; i <= columns; ++i) {
    for (int j = 1; j <= rows; j++) {
      initial_seq_pair.neg_sequence.push_back(rows * i - j);
    }
  }
}

void HierRTLMP::computeFencesAndGuides(
    const std::vector<HardMacro*>& hard_macros,
    const Rect& outline,
    std::map<int, Rect>& fences,
    std::map<int, Rect>& guides)
{
  for (int i = 0; i < hard_macros.size(); ++i) {
    if (fences_.find(hard_macros[i]->getName()) != fences_.end()) {
      fences[i] = fences_[hard_macros[i]->getName()];
      fences[i].relocate(
          outline.xMin(), outline.yMin(), outline.xMax(), outline.yMax());
    }
    auto itr = guides_.find(hard_macros[i]->getInst());
    if (itr != guides_.end()) {
      guides[i] = itr->second;
      guides[i].relocate(
          outline.xMin(), outline.yMin(), outline.xMax(), outline.yMax());
    }
  }
}

void HierRTLMP::createFixedTerminals(const Rect& outline,
                                     const UniqueClusterVector& macro_clusters,
                                     std::map<int, int>& cluster_to_macro,
                                     std::vector<HardMacro>& sa_macros)
{
  std::set<int> clusters_ids;

  for (auto& macro_cluster : macro_clusters) {
    for (auto [cluster_id, weight] : macro_cluster->getConnection()) {
      clusters_ids.insert(cluster_id);
    }
  }

  for (auto cluster_id : clusters_ids) {
    if (cluster_to_macro.find(cluster_id) != cluster_to_macro.end()) {
      continue;
    }
    auto& temp_cluster = tree_->maps.id_to_cluster[cluster_id];

    // model other cluster as a fixed macro with zero size
    cluster_to_macro[cluster_id] = sa_macros.size();
    sa_macros.emplace_back(
        std::pair<float, float>(
            temp_cluster->getX() + temp_cluster->getWidth() / 2.0
                - outline.xMin(),
            temp_cluster->getY() + temp_cluster->getHeight() / 2.0
                - outline.yMin()),
        temp_cluster->getName());
  }
}

std::vector<BundledNet> HierRTLMP::computeBundledNets(
    const UniqueClusterVector& macro_clusters,
    const std::map<int, int>& cluster_to_macro)
{
  std::vector<BundledNet> nets;

  for (auto& macro_cluster : macro_clusters) {
    const int src_id = macro_cluster->getId();

    for (auto [cluster_id, weight] : macro_cluster->getConnection()) {
      BundledNet net(
          cluster_to_macro.at(src_id), cluster_to_macro.at(cluster_id), weight);

      net.src_cluster_id = src_id;
      net.target_cluster_id = cluster_id;
      nets.push_back(net);
    }
  }

  return nets;
}

// Update the location based on the parent's perspective.
void HierRTLMP::updateChildrenShapesAndLocations(
    Cluster* parent,
    const std::vector<SoftMacro>& shaped_macros,
    const std::map<std::string, int>& soft_macro_id_map)
{
  for (auto& child : parent->getChildren()) {
    // IO clusters are fixed and have no area,
    // so their shapes and locations should not be updated
    if (!child->isIOCluster()) {
      *(child->getSoftMacro())
          = shaped_macros[soft_macro_id_map.at(child->getName())];
    }
  }
}

// Move children to their real location in the die area. The offset is the
// parent's origin.
void HierRTLMP::updateChildrenRealLocation(Cluster* parent,
                                           float offset_x,
                                           float offset_y)
{
  for (auto& child : parent->getChildren()) {
    child->setX(child->getX() + offset_x);
    child->setY(child->getY() + offset_y);
  }
}

// force-directed placement to generate guides for macros
// Attractive force and Repulsive force should be normalied separately
// Because their values can vary a lot.
void HierRTLMP::FDPlacement(std::vector<Rect>& blocks,
                            const std::vector<BundledNet>& nets,
                            float outline_width,
                            float outline_height,
                            const std::string& file_name)
{
  // following the ideas of Circuit Training
  logger_->report(
      "*****************************************************************");
  logger_->report("Start force-directed placement");
  const float ar = outline_height / outline_width;
  const std::vector<int> num_steps{1000, 1000, 1000, 1000};
  const std::vector<float> attract_factors{1.0, 100.0, 1.0, 1.0};
  const std::vector<float> repel_factors{1.0, 1.0, 100.0, 10.0};
  const float io_factor = 1000.0;
  const float max_size = std::max(outline_width, outline_height);
  std::mt19937 rand_gen(random_seed_);
  std::uniform_real_distribution<float> distribution(0.0, 1.0);

  auto calcAttractiveForce = [&](float attract_factor) {
    for (auto& net : nets) {
      const int& src = net.terminals.first;
      const int& sink = net.terminals.second;
      // float k = net.weight * attract_factor;
      float k = net.weight;
      if (blocks[src].fixed_flag == true || blocks[sink].fixed_flag == true
          || blocks[src].getWidth() < 1.0 || blocks[src].getHeight() < 1.0
          || blocks[sink].getWidth() < 1.0 || blocks[sink].getHeight() < 1.0) {
        k = k * io_factor;
      }
      const float x_dist
          = (blocks[src].getX() - blocks[sink].getX()) / max_size;
      const float y_dist
          = (blocks[src].getY() - blocks[sink].getY()) / max_size;
      const float dist = std::sqrt(x_dist * x_dist + y_dist * y_dist);
      const float f_x = k * x_dist * dist;
      const float f_y = k * y_dist * dist;
      blocks[src].addAttractiveForce(-1.0 * f_x, -1.0 * f_y);
      blocks[sink].addAttractiveForce(f_x, f_y);
    }
  };

  auto calcRepulsiveForce = [&](float repulsive_factor) {
    std::vector<int> macros(blocks.size());
    std::iota(macros.begin(), macros.end(), 0);
    std::sort(macros.begin(), macros.end(), [&](int src, int target) {
      return blocks[src].lx < blocks[target].lx;
    });
    // traverse all the macros
    auto iter = macros.begin();
    while (iter != macros.end()) {
      int src = *iter;
      for (auto iter_loop = ++iter; iter_loop != macros.end(); iter_loop++) {
        int target = *iter_loop;
        if (blocks[src].ux <= blocks[target].lx) {
          break;
        }
        if (blocks[src].ly >= blocks[target].uy
            || blocks[src].uy <= blocks[target].ly) {
          continue;
        }
        // ignore the overlap between clusters and IO ports
        if (blocks[src].getWidth() < 1.0 || blocks[src].getHeight() < 1.0
            || blocks[target].getWidth() < 1.0
            || blocks[target].getHeight() < 1.0) {
          continue;
        }
        if (blocks[src].fixed_flag == true
            || blocks[target].fixed_flag == true) {
          continue;
        }
        // apply the force from src to target
        if (src > target) {
          std::swap(src, target);
        }
        // check the overlap
        const float x_min_dist
            = (blocks[src].getWidth() + blocks[target].getWidth()) / 2.0;
        const float y_min_dist
            = (blocks[src].getHeight() + blocks[target].getHeight()) / 2.0;
        const float x_overlap
            = std::abs(blocks[src].getX() - blocks[target].getX()) - x_min_dist;
        const float y_overlap
            = std::abs(blocks[src].getY() - blocks[target].getY()) - y_min_dist;
        if (x_overlap <= 0.0 && y_overlap <= 0.0) {
          float x_dist
              = (blocks[src].getX() - blocks[target].getX()) / x_min_dist;
          float y_dist
              = (blocks[src].getY() - blocks[target].getY()) / y_min_dist;
          float dist = std::sqrt(x_dist * x_dist + y_dist * y_dist);
          const float min_dist = 0.01;
          if (dist <= min_dist) {
            x_dist = std::sqrt(min_dist);
            y_dist = std::sqrt(min_dist);
            dist = min_dist;
          }
          // const float f_x = repulsive_factor * x_dist / (dist * dist);
          // const float f_y = repulsive_factor * y_dist / (dist * dist);
          const float f_x = x_dist / (dist * dist);
          const float f_y = y_dist / (dist * dist);
          blocks[src].addRepulsiveForce(f_x, f_y);
          blocks[target].addRepulsiveForce(-1.0 * f_x, -1.0 * f_y);
        }
      }
    }
  };

  auto MoveBlock =
      [&](float attract_factor, float repulsive_factor, float max_move_dist) {
        for (auto& block : blocks) {
          block.resetForce();
        }
        if (attract_factor > 0) {
          calcAttractiveForce(attract_factor);
        }
        if (repulsive_factor > 0) {
          calcRepulsiveForce(repulsive_factor);
        }
        // normalization
        float max_f_a = 0.0;
        float max_f_r = 0.0;
        float max_f = 0.0;
        for (auto& block : blocks) {
          max_f_a = std::max(
              max_f_a,
              std::sqrt(block.f_x_a * block.f_x_a + block.f_y_a * block.f_y_a));
          max_f_r = std::max(
              max_f_r,
              std::sqrt(block.f_x_r * block.f_x_r + block.f_y_r * block.f_y_r));
        }
        max_f_a = std::max(max_f_a, 1.0f);
        max_f_r = std::max(max_f_r, 1.0f);
        // Move node
        // The move will be cancelled if the block will be pushed out of the
        // boundary
        for (auto& block : blocks) {
          const float f_x = attract_factor * block.f_x_a / max_f_a
                            + repulsive_factor * block.f_x_r / max_f_r;
          const float f_y = attract_factor * block.f_y_a / max_f_a
                            + repulsive_factor * block.f_y_r / max_f_r;
          block.setForce(f_x, f_y);
          max_f = std::max(
              max_f, std::sqrt(block.f_x * block.f_x + block.f_y * block.f_y));
        }
        max_f = std::max(max_f, 1.0f);
        for (auto& block : blocks) {
          const float x_dist
              = block.f_x / max_f * max_move_dist
                + (distribution(rand_gen) - 0.5) * 0.1 * max_move_dist;
          const float y_dist
              = block.f_y / max_f * max_move_dist
                + (distribution(rand_gen) - 0.5) * 0.1 * max_move_dist;
          block.move(x_dist, y_dist, 0.0f, 0.0f, outline_width, outline_height);
        }
      };

  // initialize all the macros
  for (auto& block : blocks) {
    block.makeSquare(ar);
    block.setLoc(outline_width * distribution(rand_gen),
                 outline_height * distribution(rand_gen),
                 0.0f,
                 0.0f,
                 outline_width,
                 outline_height);
  }

  // Iteratively place the blocks
  for (auto i = 0; i < num_steps.size(); i++) {
    // const float max_move_dist = max_size / num_steps[i];
    const float attract_factor = attract_factors[i];
    const float repulsive_factor = repel_factors[i];
    for (auto j = 0; j < num_steps[i]; j++) {
      MoveBlock(attract_factor,
                repulsive_factor,
                max_size / (1 + std::floor(j / 100)));
    }
  }
}

void HierRTLMP::setTemporaryStdCellLocation(Cluster* cluster,
                                            odb::dbInst* std_cell)
{
  const SoftMacro* soft_macro = cluster->getSoftMacro();

  if (!soft_macro) {
    return;
  }

  const int soft_macro_center_x_dbu
      = block_->micronsToDbu(soft_macro->getPinX());
  const int soft_macro_center_y_dbu
      = block_->micronsToDbu(soft_macro->getPinY());

  // Std cells are placed in the center of the cluster they belong to
  std_cell->setLocation(
      (soft_macro_center_x_dbu - std_cell->getBBox()->getDX() / 2),
      (soft_macro_center_y_dbu - std_cell->getBBox()->getDY() / 2));

  std_cell->setPlacementStatus(odb::dbPlacementStatus::PLACED);
}

void HierRTLMP::setModuleStdCellsLocation(Cluster* cluster,
                                          odb::dbModule* module)
{
  for (odb::dbInst* inst : module->getInsts()) {
    if (!inst->isCore()) {
      continue;
    }
    setTemporaryStdCellLocation(cluster, inst);
  }

  for (odb::dbModInst* mod_insts : module->getChildren()) {
    setModuleStdCellsLocation(cluster, mod_insts->getMaster());
  }
}

// Update the locations of std cells in odb using the locations that
// HierRTLMP estimates for the leaf standard clusters. This is needed
// for the orientation improvement step.
// void HierRTLMP::generateTemporaryStdCellsPlacement(Cluster* cluster)
// {
//   if (cluster->isLeaf() && cluster->getNumStdCell() != 0) {  // 
//     for (odb::dbModule* module : cluster->getDbModules()) {
//       setModuleStdCellsLocation(cluster, module);
//     }
//     for (odb::dbInst* leaf_std_cell : cluster->getLeafStdCells()) {
//       setTemporaryStdCellLocation(cluster, leaf_std_cell);
//     }
//   } else {
//     for (const auto& child : cluster->getChildren()) {
//       generateTemporaryStdCellsPlacement(child.get());
//     }
//   }
// }

float HierRTLMP::calculateRealMacroWirelength(odb::dbInst* macro)
{
  float wirelength = 0.0f;

  for (odb::dbITerm* iterm : macro->getITerms()) {
    if (iterm->getSigType() != odb::dbSigType::SIGNAL) {
      continue;
    }

    odb::dbNet* net = iterm->getNet();
    if (net != nullptr) {
      const odb::Rect bbox = net->getTermBBox();
      wirelength += block_->dbuToMicrons(bbox.dx() + bbox.dy());
    }
  }

  return wirelength;
}

void HierRTLMP::flipRealMacro(odb::dbInst* macro, const bool& is_vertical_flip)
{
  if (is_vertical_flip) {
    macro->setOrient(macro->getOrient().flipY());
  } else {
    macro->setOrient(macro->getOrient().flipX());
  }
}

void HierRTLMP::adjustRealMacroOrientation(const bool& is_vertical_flip)
{
  for (odb::dbInst* inst : block_->getInsts()) {
    if (!inst->isBlock()) {
      continue;
    }

    const float original_wirelength = calculateRealMacroWirelength(inst);
    odb::Point macro_location = inst->getLocation();

    // Flipping is done by mirroring the macro about the "Y" or "X" axis,
    // so, after flipping, we must manually set the location (lower-left corner)
    // again to move the macro back to the the position choosen by mpl2.
    flipRealMacro(inst, is_vertical_flip);
    inst->setLocation(macro_location.getX(), macro_location.getY());
    const float new_wirelength = calculateRealMacroWirelength(inst);

    debugPrint(logger_,
               MPL,
               "flipping",
               1,
               "Inst {} flip {} orig_WL {} new_WL {}",
               inst->getName(),
               is_vertical_flip ? "V" : "H",
               original_wirelength,
               new_wirelength);

    if (new_wirelength > original_wirelength) {
      flipRealMacro(inst, is_vertical_flip);
      inst->setLocation(macro_location.getX(), macro_location.getY());
    }
  }
}

void HierRTLMP::correctAllMacrosOrientation()
{
  // Apply vertical flip if necessary
  adjustRealMacroOrientation(true);

  // Apply horizontal flip if necessary
  adjustRealMacroOrientation(false);
}

void HierRTLMP::updateMacrosOnDb()
{
  for (const auto& [inst, hard_macro] : tree_->maps.inst_to_hard) {
    updateMacroOnDb(hard_macro.get());
  }
}

// We don't lock the macros here, because we'll attempt to improve
// orientation next.
void HierRTLMP::updateMacroOnDb(const HardMacro* hard_macro)
{
  odb::dbInst* inst = hard_macro->getInst();

  if (!inst) {
    return;
  }

  const int x = block_->micronsToDbu(hard_macro->getRealX());
  const int y = block_->micronsToDbu(hard_macro->getRealY());

  // Orientation must be set before location so we don't end up flipping
  // and misplacing the macro.
  inst->setOrient(hard_macro->getOrientation());
  inst->setLocation(x, y);
  inst->setPlacementStatus(odb::dbPlacementStatus::PLACED);
}

void HierRTLMP::commitMacroPlacementToDb()
{
  Snapper snapper(logger_);

  for (auto& [inst, hard_macro] : tree_->maps.inst_to_hard) {
    if (!inst) {
      continue;
    }

    snapper.setMacro(inst);
    snapper.snapMacro();

    inst->setPlacementStatus(odb::dbPlacementStatus::LOCKED);
  }
}

void HierRTLMP::setMacroPlacementFile(const std::string& file_name)
{
  macro_placement_file_ = file_name;
}

void HierRTLMP::writeMacroPlacement(const std::string& file_name)
{
  if (file_name.empty()) {
    return;
  }

  std::ofstream out(file_name);

  if (!out) {
    logger_->error(MPL, 11, "Cannot open file {}.", file_name);
  }

  out << odb::generateMacroPlacementString(block_);
}

void HierRTLMP::clear()
{
  tree_.reset();

  if (graphics_) {
    graphics_->eraseDrawing();
  }
}

void HierRTLMP::setBusPlanningOn(bool bus_planning_on)
{
  bus_planning_on_ = bus_planning_on;
}

void HierRTLMP::setDebug(std::unique_ptr<Mpl2Observer>& graphics)
{
  graphics_ = std::move(graphics);
}

void HierRTLMP::setDebugShowBundledNets(bool show_bundled_nets)
{
  graphics_->setShowBundledNets(show_bundled_nets);
}

void HierRTLMP::setDebugShowClustersIds(bool show_clusters_ids)
{
  graphics_->setShowClustersIds(show_clusters_ids);
}

void HierRTLMP::setDebugSkipSteps(bool skip_steps)
{
  graphics_->setSkipSteps(skip_steps);
}

void HierRTLMP::setDebugOnlyFinalResult(bool only_final_result)
{
  graphics_->setOnlyFinalResult(only_final_result);
  is_debug_only_final_result_ = only_final_result;
}

void HierRTLMP::setDebugTargetClusterId(const int target_cluster_id)
{
  graphics_->setTargetClusterId(target_cluster_id);
}

odb::Rect HierRTLMP::micronsToDbu(const Rect& micron_rect)
{
  return odb::Rect(block_->micronsToDbu(micron_rect.xMin()),
                   block_->micronsToDbu(micron_rect.yMin()),
                   block_->micronsToDbu(micron_rect.xMax()),
                   block_->micronsToDbu(micron_rect.yMax()));
}

//////// Pusher ////////

Pusher::Pusher(utl::Logger* logger,
               Cluster* root,
               odb::dbBlock* block,
               const std::map<Boundary, Rect>& boundary_to_io_blockage)
    : logger_(logger), root_(root), block_(block)
{
  core_ = block_->getCoreArea();
  setIOBlockages(boundary_to_io_blockage);
}

void Pusher::setIOBlockages(
    const std::map<Boundary, Rect>& boundary_to_io_blockage)
{
  for (const auto& [boundary, box] : boundary_to_io_blockage) {
    boundary_to_io_blockage_[boundary]
        = odb::Rect(block_->micronsToDbu(box.getX()),
                    block_->micronsToDbu(box.getY()),
                    block_->micronsToDbu(box.getX() + box.getWidth()),
                    block_->micronsToDbu(box.getY() + box.getHeight()));
  }
}

void Pusher::fetchMacroClusters(Cluster* parent,
                                std::vector<Cluster*>& macro_clusters)
{
  for (auto& child : parent->getChildren()) {
    if (child->getClusterType() == HardMacroCluster) {
      macro_clusters.push_back(child.get());

      for (HardMacro* hard_macro : child->getHardMacros()) {
        hard_macros_.push_back(hard_macro);
      }

    } else if (child->getClusterType() == MixedCluster) {
      fetchMacroClusters(child.get(), macro_clusters);
    }
  }
}

void Pusher::pushMacrosToCoreBoundaries()
{
  // Case in which the design has nothing but macros.
  if (root_->getClusterType() == HardMacroCluster) {
    return;
  }

  if (designHasSingleCentralizedMacroArray()) {
    return;
  }

  std::vector<Cluster*> macro_clusters;
  fetchMacroClusters(root_, macro_clusters);

  for (Cluster* macro_cluster : macro_clusters) {
    debugPrint(logger_,
               MPL,
               "boundary_push",
               1,
               "Macro Cluster {}",
               macro_cluster->getName());

    std::map<Boundary, int> boundaries_distance
        = getDistanceToCloseBoundaries(macro_cluster);

    if (logger_->debugCheck(MPL, "boundary_push", 1)) {
      logger_->report("Distance to Close Boundaries:");

      for (auto& [boundary, distance] : boundaries_distance) {
        logger_->report("{} {}", toString(boundary), distance);
      }
    }

    pushMacroClusterToCoreBoundaries(macro_cluster, boundaries_distance);
  }
}

bool Pusher::designHasSingleCentralizedMacroArray()
{
  int macro_cluster_count = 0;

  for (auto& child : root_->getChildren()) {
    switch (child->getClusterType()) {
      case MixedCluster:
        return false;
      case HardMacroCluster:
        ++macro_cluster_count;
        break;
      case StdCellCluster: {
        // Note: to check whether or not a std cell cluster is "tiny"
        // we use the area of its SoftMacro abstraction, because the
        // Cluster::getArea() will give us the actual std cell area
        // of the instances from that cluster.
        if (child->getSoftMacro()->getArea() != 0) {
          return false;
        }
      }
    }

    if (macro_cluster_count > 1) {
      return false;
    }
  }

  return true;
}

// We only group macros of the same size, so here we can use any HardMacro
// from the cluster to set the minimum distance from the respective
// boundary to trigger a push.
std::map<Boundary, int> Pusher::getDistanceToCloseBoundaries(
    Cluster* macro_cluster)
{
  std::map<Boundary, int> boundaries_distance;

  const odb::Rect cluster_box(
      block_->micronsToDbu(macro_cluster->getX()),
      block_->micronsToDbu(macro_cluster->getY()),
      block_->micronsToDbu(macro_cluster->getX() + macro_cluster->getWidth()),
      block_->micronsToDbu(macro_cluster->getY() + macro_cluster->getHeight()));

  HardMacro* hard_macro = macro_cluster->getHardMacros().front();

  Boundary hor_boundary_to_push;
  const int distance_to_left = std::abs(cluster_box.xMin() - core_.xMin());
  const int distance_to_right = std::abs(cluster_box.xMax() - core_.xMax());
  int smaller_hor_distance = 0;

  if (distance_to_left < distance_to_right) {
    hor_boundary_to_push = L;
    smaller_hor_distance = distance_to_left;
  } else {
    hor_boundary_to_push = R;
    smaller_hor_distance = distance_to_right;
  }

  const int hard_macro_width = hard_macro->getWidthDBU();
  if (smaller_hor_distance < hard_macro_width) {
    boundaries_distance[hor_boundary_to_push] = smaller_hor_distance;
  }

  Boundary ver_boundary_to_push;
  const int distance_to_top = std::abs(cluster_box.yMax() - core_.yMax());
  const int distance_to_bottom = std::abs(cluster_box.yMin() - core_.yMin());
  int smaller_ver_distance = 0;

  if (distance_to_bottom < distance_to_top) {
    ver_boundary_to_push = B;
    smaller_ver_distance = distance_to_bottom;
  } else {
    ver_boundary_to_push = T;
    smaller_ver_distance = distance_to_top;
  }

  const int hard_macro_height = hard_macro->getHeightDBU();
  if (smaller_ver_distance < hard_macro_height) {
    boundaries_distance[ver_boundary_to_push] = smaller_ver_distance;
  }

  return boundaries_distance;
}

void Pusher::pushMacroClusterToCoreBoundaries(
    Cluster* macro_cluster,
    const std::map<Boundary, int>& boundaries_distance)
{
  if (boundaries_distance.empty()) {
    return;
  }

  std::vector<HardMacro*> hard_macros = macro_cluster->getHardMacros();

  for (const auto& [boundary, distance] : boundaries_distance) {
    if (distance == 0) {
      continue;
    }

    for (HardMacro* hard_macro : hard_macros) {
      moveHardMacro(hard_macro, boundary, distance);
    }

    odb::Rect cluster_box(
        block_->micronsToDbu(macro_cluster->getX()),
        block_->micronsToDbu(macro_cluster->getY()),
        block_->micronsToDbu(macro_cluster->getX() + macro_cluster->getWidth()),
        block_->micronsToDbu(macro_cluster->getY()
                             + macro_cluster->getHeight()));

    moveMacroClusterBox(cluster_box, boundary, distance);

    // Check based on the shape of the macro cluster to avoid iterating each
    // of its HardMacros.
    if (overlapsWithHardMacro(cluster_box, hard_macros)
        || overlapsWithIOBlockage(cluster_box, boundary)) {
      debugPrint(logger_,
                 MPL,
                 "boundary_push",
                 1,
                 "Overlap found when moving {} to {}. Push reverted.",
                 macro_cluster->getName(),
                 toString(boundary));

      // Move back to original position.
      for (HardMacro* hard_macro : hard_macros) {
        moveHardMacro(hard_macro, boundary, (-distance));
      }
    }
  }
}

void Pusher::moveMacroClusterBox(odb::Rect& cluster_box,
                                 const Boundary boundary,
                                 const int distance)
{
  switch (boundary) {
    case NONE:
      return;
    case L: {
      cluster_box.moveDelta(-distance, 0);
      break;
    }
    case R: {
      cluster_box.moveDelta(distance, 0);
      break;
    }
    case T: {
      cluster_box.moveDelta(0, distance);
      break;
    }
    case B: {
      cluster_box.moveDelta(0, -distance);
      break;
    }
  }
}

void Pusher::moveHardMacro(HardMacro* hard_macro,
                           const Boundary boundary,
                           const int distance)
{
  switch (boundary) {
    case NONE:
      return;
    case L: {
      hard_macro->setXDBU(hard_macro->getXDBU() - distance);
      break;
    }
    case R: {
      hard_macro->setXDBU(hard_macro->getXDBU() + distance);
      break;
    }
    case T: {
      hard_macro->setYDBU(hard_macro->getYDBU() + distance);
      break;
    }
    case B: {
      hard_macro->setYDBU(hard_macro->getYDBU() - distance);
      break;
    }
  }
}

bool Pusher::overlapsWithHardMacro(
    const odb::Rect& cluster_box,
    const std::vector<HardMacro*>& cluster_hard_macros)
{
  for (const HardMacro* hard_macro : hard_macros_) {
    bool hard_macro_belongs_to_cluster = false;

    for (const HardMacro* cluster_hard_macro : cluster_hard_macros) {
      if (hard_macro == cluster_hard_macro) {
        hard_macro_belongs_to_cluster = true;
        break;
      }
    }

    if (hard_macro_belongs_to_cluster) {
      continue;
    }

    if (cluster_box.xMin() < hard_macro->getUXDBU()
        && cluster_box.yMin() < hard_macro->getUYDBU()
        && cluster_box.xMax() > hard_macro->getXDBU()
        && cluster_box.yMax() > hard_macro->getYDBU()) {
      return true;
    }
  }

  return false;
}

bool Pusher::overlapsWithIOBlockage(const odb::Rect& cluster_box,
                                    const Boundary boundary)
{
  if (boundary_to_io_blockage_.find(boundary)
      == boundary_to_io_blockage_.end()) {
    return false;
  }

  const odb::Rect box = boundary_to_io_blockage_.at(boundary);

  if (cluster_box.xMin() < box.xMax() && cluster_box.yMin() < box.yMax()
      && cluster_box.xMax() > box.xMin() && cluster_box.yMax() > box.yMin()) {
    return true;
  }

  return false;
}

//////// Snapper ////////

Snapper::Snapper(utl::Logger* logger) : logger_(logger), inst_(nullptr)
{
}

Snapper::Snapper(utl::Logger* logger, odb::dbInst* inst)
    : logger_(logger), inst_(inst)
{
}

void Snapper::snapMacro()
{
  snap(odb::dbTechLayerDir::VERTICAL);
  snap(odb::dbTechLayerDir::HORIZONTAL);
}

void Snapper::snap(const odb::dbTechLayerDir& target_direction)
{
  SameDirectionLayersData layers_data
      = computeSameDirectionLayersData(target_direction);

  int origin = target_direction == odb::dbTechLayerDir::VERTICAL
                   ? inst_->getOrigin().x()
                   : inst_->getOrigin().y();

  if (!layers_data.snap_layer) {
    // There are no pins to align with the track-grid.
    alignWithManufacturingGrid(origin);
    setOrigin(origin, target_direction);
    return;
  }

  odb::dbITerm* snap_pin = layers_data.layer_to_pin.at(layers_data.snap_layer);
  const LayerParameters& snap_layer_params
      = layers_data.layer_to_params.at(layers_data.snap_layer);

  if (!pinsAreAlignedWithTrackGrid(
          snap_pin, snap_layer_params, target_direction)) {
    // The idea here is to first align the origin of the macro with
    // the track-grid taking into account that the grid has a certain
    // offset with regards to (0,0) and, then, compensate the offset
    // of the pins themselves so that the lines of the grid cross
    // their center.
    origin
        = std::round(origin / snap_layer_params.pitch) * snap_layer_params.pitch
          + snap_layer_params.offset - snap_layer_params.pin_offset;

    alignWithManufacturingGrid(origin);
    setOrigin(origin, target_direction);
  }

  attemptSnapToExtraLayers(
      origin, layers_data, snap_layer_params, target_direction);
}

void Snapper::setOrigin(const int origin,
                        const odb::dbTechLayerDir& target_direction)
{
  if (target_direction == odb::dbTechLayerDir::VERTICAL) {
    inst_->setOrigin(origin, inst_->getOrigin().y());
  } else {
    inst_->setOrigin(inst_->getOrigin().x(), origin);
  }
}

SameDirectionLayersData Snapper::computeSameDirectionLayersData(
    const odb::dbTechLayerDir& target_direction)
{
  SameDirectionLayersData data;

  for (odb::dbITerm* iterm : inst_->getITerms()) {
    if (iterm->getSigType() != odb::dbSigType::SIGNAL) {
      continue;
    }

    for (odb::dbMPin* mpin : iterm->getMTerm()->getMPins()) {
      for (odb::dbBox* box : mpin->getGeometry()) {
        odb::dbTechLayer* layer = box->getTechLayer();
        if (layer->getDirection() == target_direction) {
          if (data.layer_to_pin.find(layer) != data.layer_to_pin.end()) {
            continue;
          }

          if (data.layer_to_pin.empty()) {
            data.snap_layer = layer;
          }

          data.layer_to_pin[layer] = iterm;
          data.layer_to_params[layer]
              = computeLayerParameters(layer, iterm, target_direction);
        }
      }
    }
  }

  return data;
}

LayerParameters Snapper::computeLayerParameters(
    odb::dbTechLayer* layer,
    odb::dbITerm* pin,
    const odb::dbTechLayerDir& target_direction)
{
  LayerParameters params;

  odb::dbBlock* block = inst_->getBlock();
  odb::dbTrackGrid* track_grid = block->findTrackGrid(layer);

  if (track_grid) {
    getTrackGrid(track_grid, params.offset, params.pitch, target_direction);
  } else {
    logger_->error(
        MPL, 39, "No track-grid found for layer {}", layer->getName());
  }

  params.pin_width = getPinWidth(pin, target_direction);
  params.lower_left_to_first_pin
      = getPinToLowerLeftDistance(pin, target_direction);

  // The distance between the pins and the lower-left corner
  // of the master of a macro instance may not be a multiple
  // of the track-grid, in these cases, we need to compensate
  // a small offset.
  int mterm_offset = 0;
  if (params.pitch != 0) {
    mterm_offset = params.lower_left_to_first_pin
                   - std::floor(params.lower_left_to_first_pin / params.pitch)
                         * params.pitch;
  }
  params.pin_offset = params.pin_width / 2 + mterm_offset;

  const odb::dbOrientType& orientation = inst_->getOrient();
  if (target_direction == odb::dbTechLayerDir::VERTICAL) {
    if (orientation == odb::dbOrientType::MY
        || orientation == odb::dbOrientType::R180) {
      params.pin_offset = -params.pin_offset;
    }
  } else if (orientation == odb::dbOrientType::MX
             || orientation == odb::dbOrientType::R180) {
    params.pin_offset = -params.pin_offset;
  }
  return params;
}

int Snapper::getPinWidth(odb::dbITerm* pin,
                         const odb::dbTechLayerDir& target_direction)
{
  int pin_width = 0;
  if (target_direction == odb::dbTechLayerDir::VERTICAL) {
    pin_width = pin->getBBox().dx();
  } else {
    pin_width = pin->getBBox().dy();
  }
  return pin_width;
}

void Snapper::getTrackGrid(odb::dbTrackGrid* track_grid,
                           int& origin,
                           int& step,
                           const odb::dbTechLayerDir& target_direction)
{
  // TODO: handle multiple patterns
  int count;
  if (target_direction == odb::dbTechLayerDir::VERTICAL) {
    track_grid->getGridPatternX(0, origin, count, step);
  } else {
    track_grid->getGridPatternY(0, origin, count, step);
  }
}

int Snapper::getPinToLowerLeftDistance(
    odb::dbITerm* pin,
    const odb::dbTechLayerDir& target_direction)
{
  int pin_to_origin = 0;

  odb::dbMTerm* mterm = pin->getMTerm();
  if (target_direction == odb::dbTechLayerDir::VERTICAL) {
    pin_to_origin = mterm->getBBox().xMin();
  } else {
    pin_to_origin = mterm->getBBox().yMin();
  }

  return pin_to_origin;
}

void Snapper::attemptSnapToExtraLayers(
    const int origin,
    const SameDirectionLayersData& layers_data,
    const LayerParameters& snap_layer_params,
    const odb::dbTechLayerDir& target_direction)
{
  const int total_number_of_layers
      = static_cast<int>(layers_data.layer_to_pin.size());
  const int total_attempts = 5;
  int remaining_attempts = total_attempts;

  int best_origin = origin;
  int best_number_of_snapped_layers = 1;
  int new_origin = origin;
  bool all_layers_snapped = false;

  while (remaining_attempts > 0) {
    const bool is_first_attempt = remaining_attempts == total_attempts;
    // The pins may already be aligned for all extra layers (in that case,
    // we don't need to move the macro at all), hence, we use the first
    // attempt to check if there are in fact any extra layers to snap.
    if (!is_first_attempt) {
      // Move one track ahead.
      new_origin += snap_layer_params.pitch;
      setOrigin(new_origin, target_direction);
    }

    int curr_number_of_snapped_layers = 1;
    for (const auto [layer, pin] : layers_data.layer_to_pin) {
      if (layer == layers_data.snap_layer) {
        continue;
      }

      if (pinsAreAlignedWithTrackGrid(
              pin, layers_data.layer_to_params.at(layer), target_direction)) {
        ++curr_number_of_snapped_layers;
      }

      if (curr_number_of_snapped_layers == total_number_of_layers) {
        all_layers_snapped = true;
        break;
      }
    }

    if (all_layers_snapped) {
      return;
    }

    if (curr_number_of_snapped_layers > best_number_of_snapped_layers) {
      best_number_of_snapped_layers = curr_number_of_snapped_layers;
      best_origin = new_origin;
    }

    --remaining_attempts;
  }

  setOrigin(best_origin, target_direction);
}

bool Snapper::pinsAreAlignedWithTrackGrid(
    odb::dbITerm* pin,
    const LayerParameters& layer_params,
    const odb::dbTechLayerDir& target_direction)
{
  int pin_center = target_direction == odb::dbTechLayerDir::VERTICAL
                       ? pin->getBBox().xCenter()
                       : pin->getBBox().yCenter();
  return (pin_center - layer_params.offset) % layer_params.pitch == 0;
}

void Snapper::alignWithManufacturingGrid(int& origin)
{
  const int manufacturing_grid
      = inst_->getDb()->getTech()->getManufacturingGrid();

  origin = std::round(origin / manufacturing_grid) * manufacturing_grid;
}

}  // namespace mpl2
