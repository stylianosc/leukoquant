# BaMoS

- `source/` — a private CMIC fork of NiftySeg 0.9.4, containing (among other
  things) two BaMoS-specific applications, `Seg_BiASM` and `Seg_Analysis`,
  that do not exist in mainline NiftySeg at all (confirmed via diff against
  `niftyseg/source_latest/seg-apps/`). `source/seg-apps/CMakeLists.txt` only
  builds those two -- everything else this fork's `seg-apps/` could build
  (`seg_maths`, `seg_EM`, `seg_LabFusion`, `seg_LoAd`, `seg_stats`) is a
  standard NiftySeg tool with no BaMoS-specific changes, so
  `BaMoS_WMH_*.sh` sources those from the separately built
  `leukoquant/external/niftyseg/bin/` instead (see `niftyseg/README.md`).
- `bin/` — the two BaMoS-specific binaries built from `source/`
  (`Seg_BiASM`, `Seg_Analysis`).
- `build_bamos.sh` — the build recipe.
- `scripts/` — `BaMoS_WMH_*.sh` and related processing scripts (not part of
  the build; these are what actually get run at pipeline time).
- `ICBM_Priors/` — prior images used by `BaMoS_WMH_*.sh`.

## Building

```
build_bamos.sh <install_prefix>
```

Produces `<install_prefix>/bin/{Seg_BiASM,Seg_Analysis}`. Copy those two
into `leukoquant/external/bamos/bin/`, replacing what's there.
