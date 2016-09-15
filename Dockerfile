
FROM ubuntu:16.04

RUN apt-get -y update
RUN apt-get install -y \
  build-essential \
  gfortran \
  tcl \
  git \
  m4 \
  freeglut3 \
  doxygen \
  libblas-dev \
  liblapack-dev \
  libx11-dev \
  libnuma-dev \
  libcurl4-gnutls-dev \
  zlib1g-dev \
  libhwloc-dev

RUN apt-get install -y wget

ENV moose_env=moose-environment-5_Ubuntu-16.04_x86_64.deb
RUN wget http://mooseframework.org/static/media/uploads/files/$moose_env
RUN dpkg -i $moose_env

# Source MOOSE Environment

COPY . /moose/
WORKDIR /moose
ENV env_cmd="source /opt/moose/environments/moose_profile"
RUN /bin/bash -c "$env_cmd && ./scripts/update_and_rebuild_libmesh.sh"
WORKDIR /moose/test
RUN /bin/bash -c "$env_cmd && make -j7"
RUN /bin/bash -c "$env_cmd && ./run_tests -j7"
WORKDIR /

