# Worst-Case Lattice Sampler with Truncated Gadgets and Applications
---

This repository contains the code by Corentin Jeudy for the paper **Worst-Case Lattice Sampler with Truncated Gadgets and Applications** by Corentin Jeudy and Olivier Sanders. The paper was published at Asiacrypt 2025 and is also available on the IACR ePrint Archive [here](https://eprint.iacr.org/2024/1952). Part of the code was taken from the implementation of [lattice anonymous credentials](https://github.com/Chair-for-Security-Engineering/lattice-anonymous-credentials) (by Sven Argo, Corentin Jeudy and Olivier Sanders) and [lattice blind signature](https://github.com/latticeblindsignature/lattice-blind-signature) (by Corentin Jeudy).  

The repository contains
- the parameter selection and security estimation scripts in python, and
- the implementation in C  

## Parameter Scripts

In the `/scripts` subdirectory, we provide a parameter estimation script which depends on the lattice-estimator.
The lattice-estimator is included as submodule, so please do a recursive pull.
```shell
git submodule update --init --recursive
```

You can then run the python scripts:  
- `parameters_standard_model_signature.py`: Parameter selection script for the standard model signature from [AGJLRS24] (optimized for standalone signature) and the one from Section 5.1 using our truncated sampler.  
- `parameters_group_signature.py`: Parameter selection script for the group signature from [LNPS21,LNP22] and the one from Section 5.2.1 using our truncated sampler (and/or the tricks from [AGJLRS24]).  
- `parameters_anonymous_credentials_AGJLRS24_bimodal.py`: Parameter selection script for the anonymous credentials from [AGJLRS24] (with or without the zero-knowledge optimizations of [LNP22,LN22]).  
- `parameters_anonymous_credentials.py`: Parameter selection script for the anonymous credentials with our truncated sampler (with or without the zero-knowledge optimizations of [LNP22,LN22]).  
- `parameters_blind_signature.py`: Parameter selection script for the blind signature from [JS24] and the one from Section 5.2.3 using our truncated sampler.  

**References**

[AGJLRS24] Argo, Güneysu, Jeudy, Land, Roux-Langlois, Sanders. _Practical Post-Quantum Signatures for Privacy_. CCS 2024  
[JS24] Jeudy, Sanders. _Improved Lattice Blind Signatures from Recycled Entropy_. IACR ePrint 2024/1289  
[LN22] Lyubashevsky, Nguyen. _BLOOM: Bimodal Lattice One-Out-of-Many Proofs and Applications_. ASIACRYPT 2022  
[LNP22] Lyubashevsky, Nguyen, Plançon. _Lattice-Based Zero-Knowledge Proofs and Applications: Shorter, Simpler, and More General_. CRYPTO 2022  
[LNPS21] Lyubashevsky, Nguyen, Plançon, Seiler. _Shorter Lattice-Based Group Signatures via ``Almost Free'' Encryption and Other Optimizations_. ASIACRYPT 2021

## C Implementations

The repository contains two C implementations:  
- `signature_full` folder contains the full implementation of the standard model signature from [AGJLRS24], with slightly optimized parameters for the standalone signature use-case.  
- `signature_truncated` folder contains the full implementation of the standard model signature from Section 5.1 using our truncated sampler.  

Both these implementations have their own README and Dockerfile. We nevertheless provide all the build instructions here for each signature. In the instructions, `X` must be replaced by either `full` or `truncated` depending on which implementation is built. 

Requirements:
- FLINT, which depends on GMP and MPFR
- AES-NI instructions
- cmake

FLINT is included as submodule to this repository, so please do a recursive pull.
To install FLINT, run the following commands (from the root directory of this repository). It only needs to be installed once.
```shell
cd signature_X/code/libs/flint
./bootstrap.sh 
./configure --enable-avx2 --disable-pthread CFLAGS="-O3 -Wall -march=native" CC="gcc-11"
make -j
make check
make install
```
To build our implementation, run the following commands (from the root directory of this repository).
```shell
cd signature_X/code
mkdir _build
cd _build
cmake ..
make -j
```

Then, running `./test` or `./bench` runs the tests or the benchmarks, respectively. Each change of the source code requires to run the last command (`make` or `make -j`).

## Building with Docker

We provide a Dockerfile to build a docker image corresponding to each implementation (either using our truncated sampler `/signature_truncated` or the full sampler `/signature_full`). We give the instructions to run our implementations with Docker, recalling that `X` should be replaced either by `full` or `truncated` depending on which implementation is built. 
For information on how to install docker on your system visit the [docker docs](https://docs.docker.com/).
Depending on your setup, you may need to prefix the following commands with `sudo`.

### 1. Build the docker image
Navigate to the top-level directory which contains this README file. Then run (notice the trailing dot for the second command!)
```shell
cd signature_X
docker build -t X-docker .
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
docker run -t X-docker YYY
```
You should now see the output of either the test or the benchmarking executables.

### 3. Stopping the running docker image
To get a list of running docker containers use `docker ps`. Then you can use `docker stop <container_id>` to stop
the execution. For more information see the [docker docs](https://docs.docker.com/engine/reference/builder/).

### 4. Cleaning up
To remove all stopped docker containers use `docker container prune`.
To remove docker images use `docker images` to get a list of all images and then `docker rm <image_id>` to remove
the image in question.