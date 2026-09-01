# FiPPS
Finite Element solver FiPPS²

# Manual

A manual (currently available in German only) describing the basic usage of FiPPS² is available under https://github.com/TUD-LFT/FiPPS2_manual

# Installation of required libraries

## FiPPS2

### 1. Setup `makeconfig.mk` & environment variables

Edit the `makeconfig.mk` in order to set VERSION and REVISION manually.

Assure that the proper MPI Fortran compiler is used (default: `mpif90`)

Set the environment variables to point to PETSc and SLEPc via

```
$ export PETSC_DIR=/mnt/appl/x86_64/petsc/3.22.5
$ export SLEPC_DIR=/mnt/appl/x86_64/slepc/3.22.2_mod
```

### 2. Build with `make`

In top level directory, build the program with `make` command.

The `FiPPS` executable contains all dependencies statically linked. It must be used together with `libuseroutput.so` on the `LD_LIBRARY_PATH`



## PETSc

### 1. Download current PETSc-Version

https://petsc.org/release/download/

### 2. Unpack .tar.gz

    tar xfvz [ARCHIVE].tar.gz

### 3. Set environment variables

Path, where PETSc can be found

    export PETSC_DIR=/btmpl/Software/PETSC/petsc-3.16.1

### 4. Configure

Configure - sequential MKL (recommended, as processes may get into conflict in OpenMP version)

```
module load compiler mkl mpich \
./configure \
  PETSC_ARCH=linux-gnu-intel \
  --with-mpi-dir=${MPI_ROOT} \
  --with-scalar-type=real \
  --with-debugging=0 \
  --with-fortran=1 \
  --with-blaslapack-dir=${MKLROOT} \
  --with-blacs=1 \
  --with-blacs-include=${MKLROOT}/include \
  --with-blacs-lib="-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl" \
  --with-scalapack=1 \
  --with-scalapack-lib="${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl" \
  --with-scalapack-include=${MKLROOT}/include/ \
  --download-metis=1 \
  --download-metis-use-doubleprecision=1 \
  --download-parmetis=1 \
  --download-mumps=1 \
  --prefix=/btmpl/Bibliotheken/petsc-3.22.5
```

Explanation of the options:

`PETSC_ARCH=linux-gnu-intel` - just a name (it is also the name of the folder in which the compiled files will be placed at first)

`--with-mpi-dir=${MPI_ROOT}` - Path in which MPI library can be found (MPI_ROOT has to be set in advance)

`--with-scalar-type=real` - no complex numbers

`--with-debugging=0` - Turn debugging off (therefore much faster)

`--with-fortran=1` - Also compile Fortran wrapper

`--with-blaslapack-dir="..."` - Defines BLAS/LAPACK pathes and libraries

`--with-blacs=1` - with BLACS

`--with-blacs-include=...` - Path to include files for BLACS

`--with-blacs-lib=...` - Path to libraries of BLACS

`--with-scalapack=1` - with SCALAPACK (for MUMPS)

`--with-scalapack-lib=...` - Path to libraries of SCALAPACK

`--with-scalapack-include=...` - Path to include files for SCALAPACK

`--download-metis=1`- Download and compile METIS

`--download-metis-use-doubleprecision=1`- Use METIS in the downloadable version with double precision

`--download-parmetis=1` - Download and compile PARMETIS (for MUMPS)

`--download-mumps=1` - Download and install MUMPS

`--prefix=...` - Path in which all necessary files are moved to when executing `make install`

### 5. Follow instructions

### 6. If required, install SLEPc

See section [SLEPc](#slepc)

## PETSc on Taurus

Load compiler, MPI, MKL and CMake:

```
module load iomkl/2018a
module load CMake
```

```
./configure \
--with-debugging=no \
--with-mpi-dir=${EBROOTOPENMPI} \
--with-fortran \
--with-scalar-type=real \
--with-blaslapack-dir=${MKLROOT} \
--with-blacs=1 \
--with-blacs-include=${MKLROOT}/include \
--with-blacs-lib=" -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl" \
--with-scalapack=1 \
--with-scalapack-include=${MKLROOT}/include \
--with-scalapack-lib=" ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl" \
--download-metis=1 \
--download-parmetis=1 \
--download-mumps=1 COPTFLAGS="-O3 -march=native" CXXOPTFLAGS="-O3 -march=native" FOPTFLAGS="-O3 -march=native" \
--prefix=/projects/geops2/PETSC/3.9.1_scs5
```

As sequential version

```
module load intel/2020a; \
./configure \
--with-debugging=no \
--with-mpi=0 \
--with-fortran \
--with-cc=icc \
--with-cxx=icpc \
--with-fc=ifort \
--with-scalar-type=real \
--with-blaslapack-dir=${MKLROOT} \
--with-blacs=1 \
--with-blacs-include=${MKLROOT}/include \
--with-blacs-lib=" -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl" \
--download-metis=1 \
--download-mumps=1 \
--with-mumps-serial=1 \
--CFLAGS="-O3 -xHost" \
--CXXFLAGS="-O3 -xHost" \
--FFLAGS="-O3 -xHost" \
--prefix=/projects/geops2/PETSC/3.16.1_seq
```



## SLEPc

### 1. Download current SLEPc-Version

https://slepc.upv.es/download/

### 2. Unpack .tar.gz

  tar xfvz [ARCHIVE].tar.gz

### 3. Customization

Customizations described in [Part 1](#part-1) to [Part 3](#part-3) ensure that the Cholesky decomposition of the matrix is not carried out again for each eigenvalue calculation. The Cholesky decomposition should be available after the static calculation.

#### Part 1

Change the method `STSetUp_Shift(ST st)` in the file `$SLEPC_DIR/src/sys/classes/st/impls/shift/shift.c` stored in folder `$SLEPC_DIR/src/st/impls/shift` from

```
PetscErrorCode STSetUp_Shift(ST st)
{
  PetscErrorCode ierr;
  PetscInt       k,nc,nmat=st->nmat;
  PetscScalar    *coeffs=NULL;

  PetscFunctionBegin;
  if (nmat>1) {
    ierr = STSetWorkVecs(st,1);CHKERRQ(ierr);
  }
  if (nmat>2) {  /* set-up matrices for polynomial eigenproblems */
    if (st->transform) {
      nc = (nmat*(nmat+1))/2;
      ierr = PetscMalloc1(nc,&coeffs);CHKERRQ(ierr);
      /* Compute coeffs */
      ierr = STCoeffs_Monomial(st,coeffs);CHKERRQ(ierr);
      /* T[n] = A_n */
      k = nmat-1;
      ierr = PetscObjectReference((PetscObject)st->A[k]);CHKERRQ(ierr);
      ierr = MatDestroy(&st->T[k]);CHKERRQ(ierr);
      st->T[k] = st->A[k];
      for (k=0;k<nmat-1;k++) {
        ierr = STMatMAXPY_Private(st,nmat>2?st->sigma:-st->sigma,0.0,k,coeffs?coeffs+((nmat-k)*(nmat-k-1))/2:NULL,PetscNot(st->state==ST_STATE_UPDATED),&st->T[k]);CHKERRQ(ierr);
      }
      ierr = PetscFree(coeffs);CHKERRQ(ierr);
      ierr = PetscObjectReference((PetscObject)st->T[nmat-1]);CHKERRQ(ierr);
      ierr = MatDestroy(&st->P);CHKERRQ(ierr);
      st->P = st->T[nmat-1];
      ierr = STKSPSetOperators(st,st->P,st->Pmat?st->Pmat:st->P);CHKERRQ(ierr);
    } else {
      for (k=0;k<nmat;k++) {
        ierr = PetscObjectReference((PetscObject)st->A[k]);CHKERRQ(ierr);
        ierr = MatDestroy(&st->T[k]);CHKERRQ(ierr);
        st->T[k] = st->A[k];
      }
    }
  }
  if (st->P) {
    ierr = KSPSetUp(st->ksp);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}
```

to 

```
PetscErrorCode STSetUp_Shift(ST st)
{
  PetscErrorCode ierr;
  PetscInt       k,nc,nmat=st->nmat;
  PetscScalar    *coeffs=NULL;

  PetscFunctionBegin;
  if (nmat>1) {
    ierr = STSetWorkVecs(st,1);CHKERRQ(ierr);
  }
  if (nmat>2) {  /* set-up matrices for polynomial eigenproblems */
    if (st->transform) {
      nc = (nmat*(nmat+1))/2;
      ierr = PetscMalloc1(nc,&coeffs);CHKERRQ(ierr);
      /* Compute coeffs */
      ierr = STCoeffs_Monomial(st,coeffs);CHKERRQ(ierr);
      /* T[n] = A_n */
      k = nmat-1;
      ierr = PetscObjectReference((PetscObject)st->A[k]);CHKERRQ(ierr);
      ierr = MatDestroy(&st->T[k]);CHKERRQ(ierr);
      st->T[k] = st->A[k];
      for (k=0;k<nmat-1;k++) {
        ierr = STMatMAXPY_Private(st,nmat>2?st->sigma:-st->sigma,0.0,k,coeffs?coeffs+((nmat-k)*(nmat-k-1))/2:NULL,PetscNot(st->state==ST_STATE_UPDATED),&st->T[k]);CHKERRQ(ierr);
      }
      ierr = PetscFree(coeffs);CHKERRQ(ierr);
      ierr = PetscObjectReference((PetscObject)st->T[nmat-1]);CHKERRQ(ierr);
      ierr = MatDestroy(&st->P);CHKERRQ(ierr);
      st->P = st->T[nmat-1];
      /*ierr = STKSPSetOperators(st,st->P,st->Pmat?st->Pmat:st->P);CHKERRQ(ierr);*/
    } else {
      for (k=0;k<nmat;k++) {
        ierr = PetscObjectReference((PetscObject)st->A[k]);CHKERRQ(ierr);
        ierr = MatDestroy(&st->T[k]);CHKERRQ(ierr);
        st->T[k] = st->A[k];
      }
    }
  }
  /*if (st->P) {
    ierr = KSPSetUp(st->ksp);CHKERRQ(ierr);
  }*/
  PetscFunctionReturn(0);
}
```

#### Part 2

Comment out `KSPReset` in the method `PetscErrorCode STReset(ST st)` in file `$SLEPC_DIR/src/sys/classes/st/interface/stfunc.c`. In the end, the method shall look as follows:

```
PetscErrorCode STReset(ST st)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (st) PetscValidHeaderSpecific(st,ST_CLASSID,1);
  if (!st) PetscFunctionReturn(0);
  STCheckNotSeized(st,1);
  if (st->ops->reset) { ierr = (*st->ops->reset)(st);CHKERRQ(ierr); }
  /*if (st->ksp) { ierr = KSPReset(st->ksp);CHKERRQ(ierr); }*/
  ierr = MatDestroyMatrices(PetscMax(2,st->nmat),&st->T);CHKERRQ(ierr);
  ierr = MatDestroyMatrices(PetscMax(2,st->nmat),&st->A);CHKERRQ(ierr);
  st->nmat = 0;
  ierr = PetscFree(st->Astate);CHKERRQ(ierr);
  ierr = MatDestroy(&st->Op);CHKERRQ(ierr);
  ierr = MatDestroy(&st->P);CHKERRQ(ierr);
  ierr = MatDestroy(&st->Pmat);CHKERRQ(ierr);
  ierr = VecDestroyVecs(st->nwork,&st->work);CHKERRQ(ierr);
  st->nwork = 0;
  ierr = VecDestroy(&st->wb);CHKERRQ(ierr);
  ierr = VecDestroy(&st->wht);CHKERRQ(ierr);
  ierr = VecDestroy(&st->D);CHKERRQ(ierr);
  st->state   = ST_STATE_INITIAL;
  st->opready = PETSC_FALSE;
  PetscFunctionReturn(0);
}
```

#### Part 3

Comment out `KSPDestroy` in the method `PetscErrorCode STDestroy(ST *st)` in file `$SLEPC_DIR/src/sys/classes/st/interface/stfunc.c`. In the end, the method shall look as follows:

```
PetscErrorCode STDestroy(ST *st)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*st) PetscFunctionReturn(0);
  PetscValidHeaderSpecific(*st,ST_CLASSID,1);
  if (--((PetscObject)(*st))->refct > 0) { *st = 0; PetscFunctionReturn(0); }
  ierr = STReset(*st);CHKERRQ(ierr);
  if ((*st)->ops->destroy) { ierr = (*(*st)->ops->destroy)(*st);CHKERRQ(ierr); }
  /*ierr = KSPDestroy(&(*st)->ksp);CHKERRQ(ierr);*/
  ierr = PetscHeaderDestroy(st);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
```

### 4. PETSc path

Path where PETSc is found. Settings are taken from there.

    export PETSC_DIR=/mnt/appl/x86_64/petsc/3.22.5

### 5. Path for SLEPc unpacking

Path to which SLEPc has been unpacked.

    export SLEPC_DIR=/btmpl/Software/PETSC/slepc-3.22.2

### 6. Configuration

    ./configure --prefix=/mnt/appl/x86_64/slepc/3.22.2_mod

Explanation of the option:
 
 `--prefix=/mnt/appl/x86_64/slepc/3.22.2_mod` - The SLEPC libraries are copied to this path after compiling with `make install`. Only the necessary files (without doc etc.) are transferred. (Saves disk space, especially with different versions).

### 7. Correction of includes

 Correct includes in the file `/mnt/appl/x86_64/slepc/3.22.2_mod/lib/slepc/conf/slepc_variables` from

    SLEPC_INCLUDE       = -I/mnt/appl/x86_64/slepc/3.22.2_mod/include -I/mnt/appl/x86_64/slepc/3.22.2_mod/include

to

    SLEPC_INCLUDE       = -I/mnt/appl/x86_64/slepc/3.22.2_mod/include -I${PETSC_DIR}/include
