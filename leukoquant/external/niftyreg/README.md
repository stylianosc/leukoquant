# NiftyReg (CUDA-enabled)

- `source/` — vendored, patched NiftyReg source. Pinned at commit `8538e7f~1`
  (`07e0558`, 2023-03-02) — the last commit before NiftyReg switched to the
  RNifti reimplementation of NIfTI I/O. Pinned specifically because it still
  uses the classic nifticlib-based `nifti1_io.h`/`nifti1.h`/`znzlib.h`, which
  is required for GIF to link against both NiftyReg and NiftySeg in the same
  binary without struct redefinition collisions (see `gif/README.md`).
  Patches applied on top of upstream are documented inline at each site
  (search for `[leukoquant patch]`).
- `gpu/` — **committed directly in the repo** (~40MB). This build statically
  links the CUDA math runtime (confirmed via `ldd` — no `libcudart`/
  `libcublas`/`libcusolver`/`libcusparse` needed), so unlike the `.sif`
  containers/GIF atlas database it never needed a Hugging Face download: the
  only unresolved CUDA dependency is a tiny `libcuda.so.1` driver stub. This
  is the actual runtime dependency for metrics/z-score/BaMoS's `--gpu` mode.
  `ensure_niftyreg_gpu()` (`leukoquant/utils/container_utils.py`) just
  sanity-checks it's present.
- `build_niftyreg_cuda.sh` — the build recipe that turns `source/` into what
  gets committed as `gpu/` (and separately, an input to `gif/build_gif.sh`).

## Building

Run on a GPU node (the CUDA sanity check in NiftyReg's own CMake needs a real
GPU to run against):

```
qsub -l gpu=true -pe gpu 1 -l tmem=8G,h_vmem=8G -l h_rt=1:0:0 \
  build_niftyreg_cuda.sh <install_prefix>
```

Produces `<install_prefix>/{bin,lib,include}`. This is what:
- gets copied into `gpu/` (`bin/`, `lib/`, `include/`, plus a `libcuda.so.1`
  stub copied from `$CUDA_HOME/lib64/stubs/libcuda.so`) and committed, and
- gets passed as `NIFTYREG_DIR` to `gif/build_gif.sh`.

Not something end users need to run themselves for normal `--gpu` use — the
committed `gpu/` is a prebuilt copy. Only needed when re-building NiftyReg
itself (e.g. after a source patch, or to support a new GPU architecture).
