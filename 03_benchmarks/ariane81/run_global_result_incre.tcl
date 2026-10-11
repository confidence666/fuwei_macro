# Auto-generated archive global flow for ariane81.
if {[info exists ::env(ASPDAC_MPL2_ARCHIVE_DIR)]} {
  set archive_dir $::env(ASPDAC_MPL2_ARCHIVE_DIR)
} else {
  set archive_dir "/workspace/archive/ASPDAC2026_Timing_Driven_MacroPlacer_mpl2_archive"
}
set result_dir "$archive_dir/04_results/ariane81/global_incre"
file mkdir $result_dir
set_thread_count 8
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/NangateOpenCellLibrary.tech.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/NangateOpenCellLibrary.macro.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/NangateOpenCellLibrary.macro.mod.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/NangateOpenCellLibrary.macro.rect.lef"
read_lef "$archive_dir/03_benchmarks/ariane81/inputs/flow/designs/nangate45/ariane81/lef/ariane_regfile.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_1024x32.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_128x116.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_128x256.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_128x32.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_2048x39.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_256x16.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_256x32.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_256x34.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_256x48.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_256x95.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_256x96.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_32x32.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_32x64.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_512x64.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_64x124.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_64x15.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_64x21.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_64x32.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_64x62.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_64x64.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_64x7.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/fakeram45_64x96.lef"
read_lef "$archive_dir/03_benchmarks/ariane81/inputs/flow/designs/nangate45/ariane81/lef/fifo_v3_1.lef"
read_lef "$archive_dir/03_benchmarks/ariane81/inputs/flow/designs/nangate45/ariane81/lef/fifo_v3_8.lef"
read_lib "$archive_dir/03_benchmarks/ariane81/inputs/flow/objects/nangate45/ariane81/base/lib/NangateOpenCellLibrary_typical.lib"
read_lib "$archive_dir/03_benchmarks/ariane81/inputs/flow/objects/nangate45/ariane81/base/lib/ariane_regfile.lib"
read_lib "$archive_dir/03_benchmarks/ariane81/inputs/flow/objects/nangate45/ariane81/base/lib/fakeram45_256x16.lib"
read_lib "$archive_dir/03_benchmarks/ariane81/inputs/flow/objects/nangate45/ariane81/base/lib/fakeram45_256x32.lib"
read_lib "$archive_dir/03_benchmarks/ariane81/inputs/flow/objects/nangate45/ariane81/base/lib/fakeram45_256x48.lib"
read_lib "$archive_dir/03_benchmarks/ariane81/inputs/flow/objects/nangate45/ariane81/base/lib/fifo_v3_1.lib"
read_lib "$archive_dir/03_benchmarks/ariane81/inputs/flow/objects/nangate45/ariane81/base/lib/fifo_v3_8.lib"
read_db "$archive_dir/04_results/ariane81/global/ariane81_global.odb"
read_sdc "$archive_dir/03_benchmarks/ariane81/inputs/flow/results/nangate45/ariane81/base/2_1_floorplan.sdc"

set_layer_rc -layer metal1 -resistance 5.4286e-03 -capacitance 7.41819E-02
set_layer_rc -layer metal2 -resistance 3.5714e-03 -capacitance 6.74606E-02
set_layer_rc -layer metal3 -resistance 3.5714e-03 -capacitance 8.88758E-02
set_layer_rc -layer metal4 -resistance 1.5000e-03 -capacitance 1.07121E-01
set_layer_rc -layer metal5 -resistance 1.5000e-03 -capacitance 1.08964E-01
set_layer_rc -layer metal6 -resistance 1.5000e-03 -capacitance 1.02044E-01
set_layer_rc -layer metal7 -resistance 1.8750e-04 -capacitance 1.10436E-01
set_layer_rc -layer metal8 -resistance 1.8750e-04 -capacitance 9.69714E-02
set_wire_rc -signal -layer metal3
set_wire_rc -clock -layer metal5

estimate_parasitics -placement
puts "METRIC_BEFORE_INCRE_WNS_BEGIN"
report_worst_slack -max
puts "METRIC_BEFORE_INCRE_WNS_END"
puts "METRIC_BEFORE_INCRE_TNS_BEGIN"
report_tns
puts "METRIC_BEFORE_INCRE_TNS_END"
set macros_changed [incremental_macro_placer]
if { !$macros_changed } {
  puts "METRIC_INCRE_NO_MACRO_CHANGE"
  write_def "$result_dir/ariane81_global_incre.def"
  write_db "$result_dir/ariane81_global_incre.odb"
  exit
}
write_def "$result_dir/ariane81_global_incre_macros.def"
write_db "$result_dir/ariane81_global_incre_macros.odb"
global_placement -timing_driven -skip_initial_place -incremental -density 0.7
detailed_placement
estimate_parasitics -placement
puts "METRIC_AFTER_INCRE_WNS_BEGIN"
report_worst_slack -max
puts "METRIC_AFTER_INCRE_WNS_END"
puts "METRIC_AFTER_INCRE_TNS_BEGIN"
report_tns
puts "METRIC_AFTER_INCRE_TNS_END"
write_def "$result_dir/ariane81_global_incre.def"
write_db "$result_dir/ariane81_global_incre.odb"
exit
