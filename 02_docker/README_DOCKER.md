# Docker Reproduction Notes

Build the reproduction image from a Docker host:

```bash
cd ASPDAC2026_Timing_Driven_MacroPlacer_mpl2_archive
bash 02_docker/build_docker.sh
```

Run a benchmark with the rebuilt image.  The design argument is optional and defaults to `bp_be`:

```bash
bash 02_docker/run_global.sh bp_be
bash 02_docker/run_global_incre.sh bp_be
bash 02_docker/run_global.sh ariane81
bash 02_docker/run_global_incre.sh ariane81
```

The Dockerfile uses `openroad/ubuntu22.04-dev:latest`, checks out the recorded OpenROAD commit, applies `01_source_code/openroad_mpl2.patch`, and builds `/opt/OpenROAD/build/src/openroad`.
