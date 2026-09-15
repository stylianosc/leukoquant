# NiftySeg (latest, vanilla upstream)

- `source_latest/` — vendored, unpatched latest NiftySeg source
  (`KCL-BMEIS/NiftySeg`, commit `16cf563`).
- `bin/` — the standard NiftySeg tools built from `source_latest/`
  (`seg_maths`, `seg_EM`, `seg_LabFusion`, `seg_LoAd`, `seg_stats`, etc).
  `BaMoS_WMH_*.sh` sources `seg_maths` from here (`PathNiftySeg`) rather than
  from `bamos/bin/` -- see `bamos/README.md` for why.
- `build_niftyseg.sh` — the build recipe.

## Building

```
build_niftyseg.sh <install_prefix> <eigen_include_dir>
```

`eigen_include_dir` must contain an `Eigen/` header subdirectory (only the
headers are needed, not `unsupported/`/`test/` -- NiftySeg's own
`ExternalProject_Add(Eigen)` download pulls thousands of files that can
exhaust a quota-limited filesystem's inode limit). `-DUSE_SYSTEM_EIGEN=ON`
bypasses that download entirely.

Produces `<install_prefix>/{bin,lib,include,priors}`. Only
`<install_prefix>/bin/*` is meant to be committed -- copy those into
`leukoquant/external/niftyseg/bin/`, replacing what's there.
