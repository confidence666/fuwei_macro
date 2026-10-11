# Auto-generated archive global flow for swerv_wrapper.
if {[info exists ::env(ASPDAC_MPL2_ARCHIVE_DIR)]} {
  set archive_dir $::env(ASPDAC_MPL2_ARCHIVE_DIR)
} else {
  set archive_dir "/workspace/archive/ASPDAC2026_Timing_Driven_MacroPlacer_mpl2_archive"
}
set result_dir "$archive_dir/04_results/swerv_wrapper/global"
file mkdir $result_dir
set_thread_count 8
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/NangateOpenCellLibrary.tech.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/NangateOpenCellLibrary.macro.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/NangateOpenCellLibrary.macro.mod.lef"
read_lef "$archive_dir/03_benchmarks/platforms/nangate45/lef/NangateOpenCellLibrary.macro.rect.lef"
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
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/NangateOpenCellLibrary_typical.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_1024x32.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_128x116.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_128x256.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_128x32.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_2048x39.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_256x16.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_256x32.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_256x34.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_256x48.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_256x95.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_256x96.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_32x32.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_32x64.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_512x64.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_64x124.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_64x15.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_64x21.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_64x32.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_64x62.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_64x64.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_64x7.lib"
read_lib "$archive_dir/03_benchmarks/platforms/nangate45/lib/fakeram45_64x96.lib"
read_db "$archive_dir/03_benchmarks/swerv_wrapper/inputs/flow/results/nangate45/swerv_wrapper/base/2_2_floorplan_io.odb"
read_sdc "$archive_dir/03_benchmarks/swerv_wrapper/inputs/flow/results/nangate45/swerv_wrapper/base/2_1_floorplan.sdc"

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

set ::env(MPL2_ORIGIN_DEF_PATH) "$archive_dir/03_benchmarks/swerv_wrapper/inputs/flow/results/nangate45/swerv_wrapper/base/2_3_floorplan_macro.def"
rtl_macro_placer -max_num_inst 30000 -min_num_inst 5000 -max_num_macro 12 -min_num_macro 4 -halo_width 10
write_def "$result_dir/swerv_wrapper_global_macro.def"
write_db "$result_dir/swerv_wrapper_global_macro.odb"
global_placement -timing_driven -density 0.7
detailed_placement
estimate_parasitics -placement
puts "METRIC_GLOBAL_WNS_BEGIN"
report_wns -digits 4
puts "METRIC_GLOBAL_WNS_END"
puts "METRIC_GLOBAL_TNS_BEGIN"
report_tns -digits 4
puts "METRIC_GLOBAL_TNS_END"
write_def "$result_dir/swerv_wrapper_global.def"
write_db "$result_dir/swerv_wrapper_global.odb"
exit
