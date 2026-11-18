# Adler et al. (2018) competition-coefficient database

Used by `read_adler` and the vignette *Plant communities* (`vignettes/05_plant_communities`).

## Provenance

| Item | Value |
|---|---|
| Title | Data from: Competition and coexistence in plant communities: intraspecific competition is stronger than interspecific competition |
| Authors | P.B. Adler, D. Smull, K.H. Beard, R.T. Choi, T. Furniss, A. Kulmatiski, J.M. Meiners, A.T. Tredennick, K.E. Veblen |
| DOI | [10.5061/dryad.q5mg97b](https://doi.org/10.5061/dryad.q5mg97b) (Dryad v1, 2019-05-21) |
| Paper | Adler et al. (2018) *Ecology Letters* 21:1319–1329, doi:10.1111/ele.13098 |
| Licence | **CC0 1.0** (Dryad and Zenodo records) |
| Downloaded | 2026-11-13: `datapluscode.zip` from the Zenodo mirror of the Dryad record (zenodo.org record 4999444), whose MD5 equals Dryad's |

Only the data table and its column descriptions are included here; the authors' R scripts are
in the original archive.

## Files

SHA-256, size in bytes, file (both extracted unchanged from `datapluscode.zip`, SHA-256
`4d44ecc136d136b48b3cecf610c1cc80e81252eabaea6944599a275da468adf1`):

```
f8969e508710321db098d434979139717512531539e0ec68aca1739183c21bab  306029  CompetitionRegressionData071517.csv
fa033575fb25d6613617e409efcd3e8b1d9939d7cd52724874657a85963dfbdd    2394  readme.txt
```

## The coefficients

- Long format: 1,184 rows, one coefficient per (study, focal species, competitor, treatment),
  with 41 studies.
- `competition coefficient` is the per-capita effect of the competitor on the focal species.
  Its sign convention is given per row in `Negative coefs mean competition` (`Yes`/`No`, or
  free text for a few studies).
- `read_adler` keeps (study, response, treatment) groups with every focal × competitor pair
  present and a `Yes`/`No` convention, puts signs in the Lotka–Volterra convention (negative =
  competition) and, by default, divides each row by its absolute intraspecific coefficient. It
  returns 39 matrices from 19 studies (2–6 species).
- Units are study-specific (survival, growth, biomass, population growth, …), so coefficients
  are comparable only within a study. The database has no growth rates.
