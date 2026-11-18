# Empirical data

Each reader records its conversions in the matrices' metadata. The data ship in `data/`, each
with a README giving its source, licence and checksums:

- `data/kinlock2019`: plant interaction networks of Kinlock (2019), CC0;
- `data/adler2018`: competition coefficients compiled by Adler et al. (2018), CC0;
- `data/stein2013`: the gut-microbiota model and data of Stein et al. (2013), CC BY 4.0
  (the article's supporting information) and MIT (the cLV repository's copy).

```@docs
read_kinlock
alpha_from_rii
read_adler
read_stein
```
