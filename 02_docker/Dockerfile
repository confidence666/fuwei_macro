FROM openroad/ubuntu22.04-dev:latest

ARG NUM_THREADS=8
ARG OPENROAD_REPO=https://github.com/The-OpenROAD-Project/OpenROAD.git
ARG OPENROAD_COMMIT=7fb347f37d0068b4a9344d7e6ae1b668360e174e

WORKDIR /opt
RUN git clone --recursive ${OPENROAD_REPO} OpenROAD

WORKDIR /opt/OpenROAD
RUN git checkout ${OPENROAD_COMMIT} && git submodule update --init --recursive

COPY 01_source_code/openroad_mpl2.patch /opt/OpenROAD/openroad_mpl2.patch
RUN git apply /opt/OpenROAD/openroad_mpl2.patch

RUN ./etc/Build.sh -threads=${NUM_THREADS}

WORKDIR /workspace/archive
ENTRYPOINT ["/opt/OpenROAD/build/src/openroad"]
