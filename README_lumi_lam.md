# Limited-Area Model spectral transforms on LUMI

## Usage (aka the most interesting part)

### (re)Compilation
source ~dadegrau2/deode/ectrans/lam/sources/ectrans/ENV_lumi
rm -rf ${BUILDDIR}/ectrans ${INSTALLDIR}/ectrans
mkdir -p ${BUILDDIR}/ectrans
cd ${BUILDDIR}/ectrans
ecbuild --prefix=${INSTALLDIR}/ectrans -Dfiat_ROOT=${INSTALLDIR}/fiat -DBUILD_SHARED_LIBS=OFF -DENABLE_FFTW=OFF -DENABLE_GPU=ON -DENABLE_OMPGPU=OFF -DENABLE_ACCGPU=ON -DENABLE_TESTS=OFF -DENABLE_GPU_AWARE_MPI=ON -DENABLE_CPU=OFF -DENABLE_ETRANS=ON ${SOURCEDIR}/ectrans
make -j16
make install

### test run

Allocate GPU resource with
    salloc --nodes=1 --ntasks-per-node=1 --cpus-per-task=1 --gpus-per-node=1 --account=project_462000140 --partition=standard-g --time=04:00:00 --mem=0
Recompile/run with 
    cd ${BASEDIR}/test/
    make -j16 -C ../build/ectrans/ install
    #srun ../install/ectrans/bin/ectrans-benchmark-gpu-sp-acc
    srun ../install/ectrans/bin/aatestprog-gpu-sp-acc

Note: recompiling like this may not be sufficient when modifying e.g. hip.cc files; do a full recompilation (above) in this case.

## Prerequisites: ecbuild and fiat

### ecbuild installation
    cd ${SOURCEDIR}
    git clone https://github.com/ecmwf/ecbuild.git
    cd ecbuild
    git checkout master
    sed -i -e "s/-Gfast//" cmake/compiler_flags/Cray_Fortran.cmake # remove obsolete switch -Gfast
    mkdir -p ${BUILDDIR}/ecbuild
    cd ${BUILDDIR}/ecbuild
    ${SOURCEDIR}/ecbuild/bin/ecbuild --prefix=${INSTALLDIR}/ecbuild ${SOURCEDIR}/ecbuild
    make
    make install

### fiat installation
With Cray compiler on lumi, one gets into trouble with OpenMP: for some reason, during linking the openmp library isn't found... This is solved by adding `${OpenMP_C_FLAGS}` in `programs/CMakeLists.txt`: `target_link_libraries( fiat-printbinding ${OpenMP_C_FLAGS} OpenMP::OpenMP_C )`
    cd ${SOURCEDIR}
    git clone https://github.com/ecmwf-ifs/fiat
    cd fiat
    git checkout main
	# ADD OpenMP::OpenMP_C here!
    rm -rf ${BUILDDIR}/fiat
    mkdir -p ${BUILDDIR}/fiat
    cd ${BUILDDIR}/fiat
    ecbuild -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=${INSTALLDIR}/fiat -DENABLE_MPI=ON -DBUILD_SHARED_LIBS=OFF -DENABLE_TESTS=OFF ${SOURCEDIR}/fiat
    make -j16
    make install

## History

### Different existing repositories contain elements for this:
* IAL/ectrans_withlam:withlam
integration of CPU lam sources in ectrans
* anmrde/ectrans:gpu_omp
OpenACC/(OpenMP)/hip branch of global with Cray compiler
* ddegrauwe/etrans:gpu_daand_lam
OpenACC/cufft branch for Nvidia compiler
* ddegrauwe/ectrans:gpu_daand_lam
modifs to global Nvidia-targeted code for lam:
    * LDGW switch to work on transposed arrays e.g. PGLAT in TRGTOL
    * generalizing FFTC plans to have stride
    * switch LALLOPERM2 (only used in etrans, not in ectrans)
* ddegrauwe/ectrans:gpu_lumi
OpenACC/rocfft branch for Cray compiler

### Merge plan:
1. Start from anmrde/ectrans:gpu_omp
2. integrate LAM (CPU) by *merging* in IAL/ectrans_withlam:withlam
3. move to GPU by *merging* ddegrauwe/ectrans:gpu_daand_lam
4. introduce optimizations by Lukas M.

### 1. Andreas' branch

When enabling GPU-aware MPI communications, Crya/Lumi complains about quite a lot of scalars not being present in OpenACC regions. Fixes for this were committed here.


### 2. Integrate LAM sources (CPU)

I created a new branch gpu-omp-daand-lam for this. `etrans` sources were taken from git@github.com:ACCORD-NWP/ectrans_withlam.git

Some incompatibilities due to version differences were solved as follows:
* NSTACK_MEMORY_TR isn't present in TPM_GEN (due to ectrans not being the latest IAL version). This was removed from eftdir_ctl, eftinv_ctl, eftdir_ctlad, eftinv_ctlad
* egath_spec was rewritten; not compatible with gath_spec_control that's in ectrans; I put back the version of egath_spec from cy43t2.
* same for edist_grid

### 3. LAM-GPU changes
Put the cpu-specific sources in `cpu`, and put OpenACC/cuFFT sources from Thomas B. in `gpu`.

Changes done for hipfft:
* remove all references to cuda, including tpm_fftc, fftc (cuda fft data type)
* removed LALLOPERM2 (assuming .FALSE.)
* removed LDGW, transposing ZGTF everywhere.