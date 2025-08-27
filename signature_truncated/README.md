## Building Manually

This folder contains the full C implementation of the standard model signature from Section 5.1 using our truncated sampler.

Requirements:
- FLINT, which depends on GMP and MPFR
- AES-NI instructions
- cmake

FLINT is included as submodule to this repository, so please do a recursive pull.
To install FLINT, run the following commands (from the root directory of this README). It only needs to be installed once.
```shell
cd code/libs/flint
./bootstrap.sh 
./configure --enable-avx2 --disable-pthread CFLAGS="-O3 -Wall -march=native" CC="gcc-11"
make -j
make check
make install
```
To build our implementation, run the following commands (from the root directory of this README).
```shell
cd code
mkdir _build
cd _build
cmake ..
make -j
```

Then, running `./test` or `./bench` runs the tests or the benchmarks, respectively. Each change of the source code requires to run the last command (`make` or `make -j`).


## Building with Docker

We provide a Dockerfile to build a docker image corresponding to the implementation.
For information on how to install docker on your system visit the [docker docs](https://docs.docker.com/).
Depending on your setup, you may need to prefix the following commands with `sudo`.

### 1. Build the docker image
Navigate to the top-level directory which contains this README file. Then run (notice the trailing dot for the second command!)
```shell
docker build -t truncated-docker .
```
This takes some time as docker needs to pull the base-image for Ubuntu and then builds all dependencies
(i.e. GMP, MPFR and FLINT) as well as the test and benchmark executables from scratch.
In particular the `make check` instructions are *very* time consuming.

The libraries are built with standard options and fine-tuning is possible by adjusting the corresponding
`make` and `configure` commands in the Dockerfile. However, installation paths *must not* be changed.

Note that this step needs to be repeated for any changes in the code base or Dockerfile but not
for subsequent runs of the code.

### 2. Run the docker image
Issue the following command to start the docker container, where `YYY` is either `./test` or `./bench` to select
the tests or benchmarks respectively.
```shell
docker run -t truncated-docker YYY
```
You should now see the output of either the test or the benchmarking executables.

### 3. Stopping the running docker image
To get a list of running docker containers use `docker ps`. Then you can use `docker stop <container_id>` to stop
the execution. For more information see the [docker docs](https://docs.docker.com/engine/reference/builder/).

### 4. Cleaning up
To remove all stopped docker containers use `docker container prune`.
To remove docker images use `docker images` to get a list of all images and then `docker rm <image_id>` to remove
the image in question.