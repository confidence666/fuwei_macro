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

#pragma once

#include <limits>
#include <string>

#include "Mpl2Observer.h"
#include "clusterEngine.h"
#include "odb/dbTransform.h"
// 添加的库
#include <tcl.h>
#include <algorithm>
#include <unistd.h>
#include "absl/log/initialize.h"  // 日志初始化
#include <cmath>
// #include "ortools/linear_solver/linear_solver.h"

// 线性规划的库
#include <stdlib.h>
#include "ortools/sat/model.h"
#include "ortools/base/logging.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/cp_model.pb.h"
#include "ortools/sat/cp_model_solver.h"
#include "ortools/util/sorted_interval_list.h"



namespace odb {
class dbBlock;
class dbDatabase;
class dbInst;
class dbModule;
}  // namespace odb

namespace sta {
class dbNetwork;
class dbSta;
}  // namespace sta

namespace utl {
class Logger;
}

namespace par {
class PartitionMgr;
}

namespace mpl2 {
struct BundledNet;
class Cluster;
class HardMacro;
class Metrics;
struct Rect;
class SoftMacro;
class Snapper;
class SACoreSoftMacro;
class SACoreHardMacro;

// The parameters necessary to compute one coordinate of the new
// origin for aligning the macros' pins to the track-grid
struct LayerParameters
{
  int offset = 0;
  int pitch = 0;
  int pin_width = 0;
  int pin_offset = 0;
  int lower_left_to_first_pin = 0;
};

using LayersWithPinsMap = std::map<odb::dbTechLayer*, odb::dbITerm*>;
using LayerParametersMap = std::map<odb::dbTechLayer*, LayerParameters>;

// 添加的结构体
struct MacroDataflowNet{
  MacroDataflowNet() : weight(0.0), driver_inst(nullptr), load_inst(nullptr) {}
  // int id;
  float weight;
  odb::dbInst* driver_inst;
  odb::dbInst* load_inst;
  // std::set<odb::dbInst*> unorderd_items;
  MacroDataflowNet(float weight, odb::dbInst* driver_inst,odb::dbInst* load_inst) : 
    weight(weight), driver_inst(driver_inst), load_inst(load_inst) {}
};

  // KD树
  struct KDNode {
    int id;
    float x;
    float y;
    KDNode* left;
    KDNode* right;
    KDNode(int id, float x, float y) : id(id), x(x), y(y), left(nullptr), right(nullptr) {}
  };


// Hierarchical RTL-MP
// Support Multi-Level Clustering.
// Support designs with IO Pads.
// Support timing-driven macro placement (not enabled now, June, 2022)
// Identify RTL information based on logical hierarchy, connection signature and
// dataflow. Connection Signature is defined as the connection topology respect
// to other clusters. Dataflow is defined based on sequential graph.
// Timing-driven macro placement is based on uniform delay model.
class HierRTLMP
{
 public:
  HierRTLMP(sta::dbNetwork* network,
            odb::dbDatabase* db,
            sta::dbSta* sta,
            utl::Logger* logger,
            par::PartitionMgr* tritonpart);
  ~HierRTLMP();

  void init();
  void run();
  bool runIncrementalMacroPlacement();

  // Interfaces functions for setting options
  // Hierarchical Macro Placement Related Options
  void setGlobalFence(float fence_lx,
                      float fence_ly,
                      float fence_ux,
                      float fence_uy);
  void setHaloWidth(float halo_width);
  void setHaloHeight(float halo_height);
  void setGuidanceRegions(const std::map<odb::dbInst*, Rect>& guidance_regions);

  // Clustering Related Options
  void setNumBundledIOsPerBoundary(int num_bundled_ios);
  void setClusterSize(int max_num_macro,
                      int min_num_macro,
                      int max_num_inst,
                      int min_num_inst);
  void setClusterSizeTolerance(float tolerance);
  void setMaxNumLevel(int max_num_level);
  void setClusterSizeRatioPerLevel(float coarsening_ratio);
  void setLargeNetThreshold(int large_net_threshold);
  void setSignatureNetThreshold(int signature_net_threshold);
  void setAreaWeight(float area_weight);
  void setOutlineWeight(float outline_weight);
  void setWirelengthWeight(float wirelength_weight);
  void setGuidanceWeight(float guidance_weight);
  void setFenceWeight(float fence_weight);
  void setBoundaryWeight(float boundary_weight);
  void setNotchWeight(float notch_weight);
  void setMacroBlockageWeight(float macro_blockage_weight);
  void setPinAccessThreshold(float pin_access_th);
  void setTargetUtil(float target_util);
  void setTargetDeadSpace(float target_dead_space);
  void setMinAR(float min_ar);
  void setSnapLayer(int snap_layer);
  void setReportDirectory(const char* report_directory);
  void setDebug(std::unique_ptr<Mpl2Observer>& graphics);
  void setDebugShowBundledNets(bool show_bundled_nets);
  void setDebugShowClustersIds(bool show_clusters_ids);
  void setDebugSkipSteps(bool skip_steps);
  void setDebugOnlyFinalResult(bool only_final_result);
  void setDebugTargetClusterId(int target_cluster_id);
  void setBusPlanningOn(bool bus_planning_on);

  void setNumThreads(int threads) { num_threads_ = threads; }
  void setMacroPlacementFile(const std::string& file_name);
  void writeMacroPlacement(const std::string& file_name);

 private:
  using SoftSAVector = std::vector<std::unique_ptr<SACoreSoftMacro>>;
  using HardSAVector = std::vector<std::unique_ptr<SACoreHardMacro>>;
  using IOSpans = std::map<Boundary, std::pair<float, float>>;

////////////////////////////////////////////////////////////////////////
// TDMP的函数
////////////////////////////////////////////////////////////////////////
  void GlobalMacroPlacement();

  struct PairComparator {
    bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) const {
        return std::tie(a.first, a.second) < std::tie(b.first, b.second);
    }
  };  

  void runNewMacroPlacement();
  void rebackRePlaceClusterLocation();
  void removeHardMacroOverlap();

  void shapeGenerate();
  void WriteClusterResult();
  void WriteDef();
  void WriteDef(std::unordered_map<int,std::pair<float,float>> shrink_loc_map);
  void WriteLef();
  void WriteLib();
  void WriteVerilog();

  std::pair<float, float> setRandomShapes(Cluster& cluster);
  void stdCellPlacement();

  void InitializeMacroLocation();

  std::vector<std::pair<std::pair<int,int>,bool>> ioid_clusterid_vector_;//IOid,clusterid,是否是IO
  std::unordered_map<int,std::pair<float,float>> id_shape_map_;  //顶层簇的
  float expand_ratio_ = 1.0;
  std::string stage_path_;
  bool diverge_flag_ = false;

  // 运行顶层簇布局的函数
  void runReplaceTcl();
  void shrinkReplace();

  // 线性规划的函数
  void FirstLinearSolver(Cluster* cluster, std::vector<MacroDataflowNet> macro_dataflow_nets);
  std::pair<int,int> SameMacroBestTiling(Cluster* cluster);

  void SecondLinearSolver(Cluster* cluster,std::vector<std::vector<float>>& fence_region);
  void SecondLinearSolver();

  // 合法化
  void macroLegalize(bool incre_flag);

  // 伪网络约束
  void clusterFakeNetStdGenerate();
  // 二级簇的布局方法
  void PlaceSubCluster(Cluster* parent);
  std::vector<MacroDataflowNet> createMacroDataflowNet(Cluster* parent);
  // void SimpleSatProgram();
  // void insertKDtree();

  // 定义一个哈希函数用于 std::pair<odb::dbInst*, odb::dbInst*>
  struct pairHash2MacroNet {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
      auto hash1 = std::hash<T1>{}(p.first);
      auto hash2 = std::hash<T2>{}(p.second);
      return hash1 ^ hash2;
    }
  };



////////////////////////////////////////////////////////////////////////
// IncreMacro的函数
////////////////////////////////////////////////////////////////////////
bool IncrementalMacroPlacement(); 
void IncreInitialize();
// void CriticalMacroFilter();
void MacroClassify();
int TilingMacroExchange();
int FreeMacroOptimize();
// 新创建的函数
void CriticalMacroFilter();
struct DelaySample { 
  float dx;
  float dy;
  float delay;
  float weight;
}; // 用于delay函数拟合

std::map<odb::dbInst*,std::set<odb::dbInst*>> critical_macro_to_regs_;  // macro和连接的regs
std::map<odb::dbInst*,float> critical_regs_slack_; // reg的slack表
std::unordered_map<odb::dbInst*,std::vector<float>> macro_delay_fit; // 每个宏和寄存器拟合的delay参数
std::vector<std::vector<HardMacro*>> tile_critical_macros_; // tile的关键宏，一个vector中是能交换的
std::vector<HardMacro*> free_critical_macros_;  // free的关键宏

// Incremental placement is launched from a placed DEF/ODB after standard-cell
// placement, so it cannot depend on the global RTLMP HardMacro tree being
// rebuilt. Keep dbInst-level groups for the ASP-DAC refinement stages.
std::vector<std::vector<odb::dbInst*>> array_critical_macro_insts_;
std::vector<odb::dbInst*> free_critical_macro_insts_;

std::unordered_map<odb::dbInst*,std::vector<float>> critical_macro_slack_; // 统计关键宏的负slack的pin的slack

// 改成下面这个数据结构
std::unordered_map<odb::dbInst*, std::unordered_map<odb::dbITerm*, std::set<odb::dbInst*>>> macro_pin_regs_;
// macro-pin-reg


std::unordered_map<odb::dbInst*, bool> placed_macros_; // 不移动的已经locked的宏，在free中使用，避免重叠

void initializeSTA();
void analyzeMacroPins(const std::string& file_path,const std::string& output_file_path);
void findDownstreamRegisters(odb::dbInst* macro,odb::dbITerm* iterm,std::set<odb::dbInst*>& found_regs,std::set<odb::dbITerm*>& visited);
void findUpstreamRegisters(odb::dbInst* macro,odb::dbITerm* iterm,std::set<odb::dbInst*>& found_regs,std::set<odb::dbITerm*>& visited);
void delayFitting(odb::dbInst* macro_inst); // 拟合函数，每个宏一个delay的二次函数
std::pair<float,float> freeMacroLP(odb::dbInst* inst, int order);

int MacroOrientOptimizer();
void findBestOrient(odb::dbInst* macro,std::unordered_map<odb::dbITerm*,std::set<odb::dbInst*>> pin_regs);
void findBestOrient(odb::dbInst* macro);

int64_t ComputeMacroNetHpwl(odb::dbInst* macro);

// 需要更详细的pin和reg的连接

////////////////////////////////////////////////////////////////////////
// Origin
////////////////////////////////////////////////////////////////////////
  void runMultilevelAutoclustering();
  void runHierarchicalMacroPlacement();

  void resetSAParameters();

  void updateMacrosOnDb();
  void updateMacroOnDb(const HardMacro* hard_macro);
  void commitMacroPlacementToDb();
  void clear();
  void FDPlacement(std::vector<Rect>& blocks,
                   const std::vector<BundledNet>& nets,
                   float outline_width,
                   float outline_height,
                   const std::string& file_name);

  // Coarse Shaping
  void runCoarseShaping();
  void setRootShapes();
  void calculateChildrenTilings(Cluster* parent);
  void calculateMacroTilings(Cluster* cluster);
  void setTightPackingTilings(Cluster* macro_array);
  void setIOClustersBlockages();
  IOSpans computeIOSpans();
  float computeIOBlockagesDepth(const IOSpans& io_spans);
  void setPlacementBlockages();

  // Fine Shaping
  bool runFineShaping(Cluster* parent,
                      std::vector<SoftMacro>& macros,
                      std::map<std::string, int>& soft_macro_id_map,
                      float target_util,
                      float target_dead_space);

  // Hierarchical Macro Placement 1st stage: Cluster Placement
  void runHierarchicalMacroPlacement(Cluster* parent);
  void adjustMacroBlockageWeight();
  void reportSAWeights();
  void runHierarchicalMacroPlacementWithoutBusPlanning(Cluster* parent);
  void runEnhancedHierarchicalMacroPlacement(Cluster* parent);

  void findOverlappingBlockages(std::vector<Rect>& blockages,
                                std::vector<Rect>& placement_blockages,
                                const Rect& outline);
  void computeBlockageOverlap(std::vector<Rect>& overlapping_blockages,
                              const Rect& blockage,
                              const Rect& outline);
  void updateChildrenShapesAndLocations(
      Cluster* parent,
      const std::vector<SoftMacro>& shaped_macros,
      const std::map<std::string, int>& soft_macro_id_map);
  void updateChildrenRealLocation(Cluster* parent,
                                  float offset_x,
                                  float offset_y);
  void mergeNets(std::vector<BundledNet>& nets);

  // Hierarchical Macro Placement 2nd stage: Macro Placement
  void placeMacros(Cluster* cluster);
  void computeFencesAndGuides(const std::vector<HardMacro*>& hard_macros,
                              const Rect& outline,
                              std::map<int, Rect>& fences,
                              std::map<int, Rect>& guides);
  void createFixedTerminals(const Rect& outline,
                            const UniqueClusterVector& macro_clusters,
                            std::map<int, int>& cluster_to_macro,
                            std::vector<HardMacro>& sa_macros);
  std::vector<BundledNet> computeBundledNets(
      const UniqueClusterVector& macro_clusters,
      const std::map<int, int>& cluster_to_macro);
  void setArrayTilingSequencePair(Cluster* cluster,
                                  int macros_to_place,
                                  SequencePair& initial_seq_pair);

  // Orientation Improvement
  void generateTemporaryStdCellsPlacement(Cluster* cluster);
  void setModuleStdCellsLocation(Cluster* cluster, odb::dbModule* module);
  void setTemporaryStdCellLocation(Cluster* cluster, odb::dbInst* std_cell);

  void correctAllMacrosOrientation();
  float calculateRealMacroWirelength(odb::dbInst* macro);
  void adjustRealMacroOrientation(const bool& is_vertical_flip);
  void flipRealMacro(odb::dbInst* macro, const bool& is_vertical_flip);

  // Bus Planning
  void callBusPlanning(std::vector<SoftMacro>& shaped_macros,
                       std::vector<BundledNet>& nets_old);
  void adjustCongestionWeight();

  // Aux for conversion
  odb::Rect micronsToDbu(const Rect& micron_rect);

  sta::dbNetwork* network_ = nullptr;
  odb::dbDatabase* db_ = nullptr;
  odb::dbBlock* block_ = nullptr;
  utl::Logger* logger_ = nullptr;
  par::PartitionMgr* tritonpart_ = nullptr;
  sta::dbSta* sta_ = nullptr;
  std::unique_ptr<PhysicalHierarchy> tree_;

  std::unique_ptr<ClusteringEngine> clustering_engine_;

  // flag variables
  const bool dynamic_congestion_weight_flag_ = false;
  // Our experiments show that for most testcases, turn off bus planning
  // can generate better results.

  // We recommend that you turn on this flag for technology nodes with very
  // limited routing layers such as SkyWater130.  But for NanGate45,
  // ASASP7, you should turn off this option.
  bool bus_planning_on_ = false;

  // Parameters related to macro placement
  std::string report_directory_;
  std::string macro_placement_file_;

  // User can specify a global region for some designs
  float global_fence_lx_ = std::numeric_limits<float>::max();
  float global_fence_ly_ = std::numeric_limits<float>::max();
  float global_fence_ux_ = 0.0;
  float global_fence_uy_ = 0.0;

  const int num_runs_ = 10;    // number of runs for SA
  int num_threads_ = 10;       // number of threads
  const int random_seed_ = 0;  // random seed for deterministic

  float target_dead_space_ = 0.2;  // dead space for the cluster
  float target_util_ = 0.25;       // target utilization of the design
  const float target_dead_space_step_ = 0.05;  // step for dead space
  const float target_util_step_ = 0.1;         // step for utilization
  const float num_target_util_ = 10;
  const float num_target_dead_space_ = 20;

  float min_ar_ = 0.3;  // the aspect ratio range for StdCellCluster (min_ar_, 1
                        // / min_ar_)

  float pin_access_th_ = 0.1;  // each pin access is modeled as a SoftMacro
  float pin_access_th_orig_ = 0.1;
  float pin_access_net_width_ratio_
      = 0.1;  // define the ratio of number of connections
              // related to IOs to the range of these IO spans
  float notch_v_th_ = 10.0;
  float notch_h_th_ = 10.0;

  int snap_layer_ = 4;

  // SA related parameters
  // weight for different penalty
  float area_weight_ = 0.1;
  float outline_weight_ = 1.0;
  float wirelength_weight_ = 1.0;
  float guidance_weight_ = 10.0;
  float fence_weight_ = 10.0;
  float boundary_weight_ = 5.0;
  float notch_weight_ = 1.0;
  float macro_blockage_weight_ = 1.0;

  std::map<std::string, Rect> fences_;   // macro_name, fence
  std::map<odb::dbInst*, Rect> guides_;  // Macro -> Guidance Region
  std::vector<Rect> placement_blockages_;
  std::vector<Rect> macro_blockages_;
  std::map<Boundary, Rect> boundary_to_io_blockage_;

  // Fast SA hyperparameter
  float init_prob_ = 0.9;
  const int max_num_step_ = 2000;
  const int num_perturb_per_step_ = 500;

  // probability of each action
  float pos_swap_prob_ = 0.2;
  float neg_swap_prob_ = 0.2;
  float double_swap_prob_ = 0.2;
  float exchange_swap_prob_ = 0.2;
  float flip_prob_ = 0.4;
  float resize_prob_ = 0.4;

  // statistics of the design
  Metrics* metrics_ = nullptr;

  const int bus_net_threshold_ = 32;  // only for bus planning
  float congestion_weight_ = 0.5;     // for balance timing and congestion

  // since we convert from the database unit to the micrometer
  // during calculation, we may loss some accuracy.
  const float conversion_tolerance_ = 0.1;

  bool skip_macro_placement_ = false;

  std::unique_ptr<Mpl2Observer> graphics_;
  bool is_debug_only_final_result_{false};
};

class Pusher
{
 public:
  Pusher(utl::Logger* logger,
         Cluster* root,
         odb::dbBlock* block,
         const std::map<Boundary, Rect>& boundary_to_io_blockage);

  void pushMacrosToCoreBoundaries();

 private:
  void setIOBlockages(const std::map<Boundary, Rect>& boundary_to_io_blockage);
  bool designHasSingleCentralizedMacroArray();
  void pushMacroClusterToCoreBoundaries(
      Cluster* macro_cluster,
      const std::map<Boundary, int>& boundaries_distance);
  void fetchMacroClusters(Cluster* parent,
                          std::vector<Cluster*>& macro_clusters);
  std::map<Boundary, int> getDistanceToCloseBoundaries(Cluster* macro_cluster);
  void moveHardMacro(HardMacro* hard_macro, Boundary boundary, int distance);
  void moveMacroClusterBox(odb::Rect& cluster_box,
                           Boundary boundary,
                           int distance);
  bool overlapsWithHardMacro(
      const odb::Rect& cluster_box,
      const std::vector<HardMacro*>& cluster_hard_macros);
  bool overlapsWithIOBlockage(const odb::Rect& cluster_box, Boundary boundary);

  utl::Logger* logger_;

  Cluster* root_;
  odb::dbBlock* block_;
  odb::Rect core_;

  std::map<Boundary, odb::Rect> boundary_to_io_blockage_;
  std::vector<HardMacro*> hard_macros_;
};

struct SameDirectionLayersData
{
  LayersWithPinsMap layer_to_pin;
  LayerParametersMap layer_to_params;
  odb::dbTechLayer* snap_layer = nullptr;
};

class Snapper
{
 public:
  Snapper(utl::Logger* logger);
  Snapper(utl::Logger* logger, odb::dbInst* inst);

  void setMacro(odb::dbInst* inst) { inst_ = inst; }
  void snapMacro();

 private:
  void snap(const odb::dbTechLayerDir& target_direction);
  void alignWithManufacturingGrid(int& origin);
  void setOrigin(int origin, const odb::dbTechLayerDir& target_direction);
  bool pinsAreAlignedWithTrackGrid(odb::dbITerm* pin,
                                   const LayerParameters& layer_params,
                                   const odb::dbTechLayerDir& target_direction);

  SameDirectionLayersData computeSameDirectionLayersData(
      const odb::dbTechLayerDir& target_direction);
  LayerParameters computeLayerParameters(
      odb::dbTechLayer* layer,
      odb::dbITerm* pin,
      const odb::dbTechLayerDir& target_direction);
  void getTrackGrid(odb::dbTrackGrid* track_grid,
                    int& origin,
                    int& step,
                    const odb::dbTechLayerDir& target_direction);
  int getPinWidth(odb::dbITerm* pin,
                  const odb::dbTechLayerDir& target_direction);
  int getPinToLowerLeftDistance(odb::dbITerm* pin,
                                const odb::dbTechLayerDir& target_direction);
  void attemptSnapToExtraLayers(int origin,
                                const SameDirectionLayersData& layers_data,
                                const LayerParameters& snap_layer_params,
                                const odb::dbTechLayerDir& target_direction);

  utl::Logger* logger_;
  odb::dbInst* inst_;
};

}  // namespace mpl2
