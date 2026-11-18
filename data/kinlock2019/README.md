# Kinlock (2019) plant interaction networks

Used by `read_kinlock` and the vignette *Plant communities* (`vignettes/05_plant_communities`).

## Provenance

| Item | Value |
|---|---|
| Title | Data from: A meta-analysis of plant interaction networks reveals competitive hierarchies as well as facilitation and intransitivity |
| Author | Nicole L. Kinlock |
| DOI | [10.5061/dryad.1sm06sp](https://doi.org/10.5061/dryad.1sm06sp) (Dryad version 3, 2020-10-28) |
| Paper | Kinlock N.L. (2019) *American Naturalist* 194:640–653, doi:10.1086/705293 |
| Licence | **CC0 1.0** (Dryad and Zenodo records) |
| Downloaded | 2026-11-13, from the Zenodo mirror of the Dryad record (zenodo.org record 4147030) |
| Integrity | The SHA-256 of the data file equals the digest Dryad publishes; both files are byte-identical to the Dryad v3 files |

## Files

SHA-256, size in bytes, file:

```
2c96dff7b8027b79aa7d5434f7e71db911e5483c1fc9cb5cf43537a02094e282  21203  Data_PlantInteractionNetworks.csv
a95c2e723a3133fd3414c378a46cc15b396d66ad2b8f9206b11f1877c1030284   3160  README_for_Data_PlantInteractionNetworks.csv
```

- `Data_PlantInteractionNetworks.csv`: 31 networks (one per study) × 27 columns; missing values
  are the string `NA`.
- `README_for_Data_PlantInteractionNetworks.csv`: the author's column dictionary (CR line
  endings).

## The interaction matrices

- `NetworkMonoCtrl` holds an n×n matrix of relative interaction indices (RII) against
  monoculture controls, flattened row by row. Entry (i, j) is the effect of species j on
  species i; the diagonal is 0 by construction. 25 of the 27 networks with this column are
  complete.
- `NetworkTrueCtrl` holds the same against plants grown alone (10 of 15 complete).
- RII = (P_with − P_without)/(P_with + P_without) lies in [−1, 1]; negative values are
  competition. RIIs are relative changes in performance, not per-capita Lotka–Volterra
  coefficients. `alpha_from_rii` reads them as such, an assumption it records in the matrix
  metadata. The files have no growth rates.
- The performance metric and its units differ between studies.
