# GIF

- `gif_source/` — vendored, patched GIF source (`KCL-BMEIS/gif`, `main`
  branch). Patches (documented inline, search for `[leukoquant patch]`):
  `-platf`/`-segPT` CLI flags added to select CPU/CUDA registration and the
  probability threshold, the `postprocessXML()` missing-return bug fixed
  (undefined behavior, the root cause of a real reproducible SIGSEGV), and
  Brain/TIV mask output restored (removed upstream in a private fork commit,
  `d6f47e2`, but still expected by downstream consumers -- see git log for
  full detail).
- `bin_gpu/` — the built `seg_GIF` binary, plus its bundled
  `libcudart.so.11.0` and a `libcuda.so.1` driver stub (fallback for nodes
  with no real GPU driver -- on a real GPU node, apptainer `--nv`'s real
  driver takes precedence automatically via normal dynamic-linker search
  order). Fully self-contained: `-platf 0` (CPU) and `-platf 1` (CUDA) both
  work from this one binary, no separate builds and no download needed.
- `GIF_200826.sh` — the runtime wrapper script `gif_workflow.smk` actually
  invokes.
- `build_gif.sh` — the build recipe for `bin_gpu/seg_GIF`.

## Building

Depends on two other builds having already been done first:

1. `leukoquant/external/niftyreg/build_niftyreg_cuda.sh` -> a NiftyReg
   install tree. Must be the CUDA-enabled build even if you only intend to
   run `-platf 0` -- there is only ever one `seg_GIF` binary.
2. `leukoquant/external/niftyseg/build_niftyseg.sh` -> a NiftySeg install
   tree.

Then, on a GPU node:

```
qsub -l gpu=true -pe gpu 1 -l tmem=8G,h_vmem=8G -l h_rt=1:0:0 \
  build_gif.sh <install_prefix> <niftyreg_install_dir> <niftyseg_install_dir> <eigen_include_dir>
```

Produces `<install_prefix>/bin/seg_GIF`. Copy that one binary into
`leukoquant/external/gif/bin_gpu/seg_GIF`, replacing what's there
(`libcudart.so.11.0`/`libcuda.so.1` alongside it don't need rebuilding --
they're copied CUDA runtime/driver-stub libraries, not part of this build).
