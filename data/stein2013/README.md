# Stein et al. (2013): gut microbiota of gnotobiotic mice under clindamycin

Two files are read by `read_stein(joinpath(pkgdir(EcologicalResilience), "data", "stein2013"))`, and
used by the vignette `06_gut_microbiota`:

- `derived/stein2013_glv.csv` holds the published generalised Lotka–Volterra fit.
- `clv_data_stein/raw_data.csv` holds the mouse time series.

## Provenance

| Item | Value |
|---|---|
| Paper | Stein R.R., Bucci V., Toussaint N.C., Buffie C.G., Rätsch G., Pamer E.G., Sander C., Xavier J.B. (2013) Ecological modeling from time-series inference: insight into dynamics and stability of intestinal microbiota. *PLoS Comput Biol* 9(12):e1003388. doi:10.1371/journal.pcbi.1003388. The underlying experiment is Buffie et al. (2012). |
| Source of the files | The cLV repository <https://github.com/tyjo/clv> (Joseph, Shenhav, Xavier, Halperin & Pe'er 2020, *PLoS Comput Biol*, doi:10.1371/journal.pcbi.1007917), commit `c184cd5f98a8f471f9e6131c88e7c9cf7717cf47` (2021-01-25), directory `data/stein/`. It redistributes the article's Dataset S1 (`Dataset_S1.xlsx`) and a CSV of its time series. |
| Downloaded | 2026-11-13 |
| Licence | The article and its Supporting Information are **CC BY 4.0** (PLoS); the data terms are those of the article's Dataset S1. The cLV repository's own files are **MIT**. Cite Stein et al. (2013) when using these data. |

## Files

SHA-256, size in bytes, path:

```
58ff99bc23c0524cdc38176086386e1cd75291d4f82682b1c9deacf0ac49ed20  1747  derived/stein2013_glv.csv
80093d52b91852615740b5d1c496c1056b2d06ca95652992f285fac680c079fd  1000  derived/extract_mmue.py
a141d18cb8578b7bfd01dfd92ddec8a02b5d83edc88cb49c3115dfeb2250f7a5  8424  clv_data_stein/raw_data.csv
```

### `derived/stein2013_glv.csv`: the published fit

- **Source:** extracted by `derived/extract_mmue.py` from sheet `MmuE` of Dataset S1 (SHA-256
  `7d0ad6623b1bed26f6fda4bcd26563682eb13e073410959bb90e178aebbc9077`, not bundled). Run it with
  `python3 -I extract_mmue.py Dataset_S1.xlsx stein2013_glv.csv`; it needs openpyxl.
- **Columns:**
  - `genus`;
  - the 11 interaction coefficients `M[genus, j]`, the effect of genus `j` on the row genus,
    with each column named by `j`;
  - `growth` (μ);
  - `susceptibility` (ε).
- **Row order:** as in the sheet, sorted by clindamycin susceptibility.
- **Model:** dxᵢ/dt = xᵢ (μᵢ + Σⱼ Mᵢⱼ xⱼ + εᵢ u(t)), where u(t) is the clindamycin signal.
- **Units:**
  - μ in d⁻¹;
  - M in (10¹¹ copies cm⁻³)⁻¹ d⁻¹;
  - ε in d⁻¹ per unit signal.
- **Regularisation:** Tikhonov, with λ* = (2.25, 9, 0.25) for (interactions, growth,
  susceptibility), from sheet `D_lambda`.

### `clv_data_stein/raw_data.csv`: the time series

- **Layout:** sheet `Y` of Dataset S1 without its two label rows, and with no trailing newline.
  It is wide and transposed: each row is a variable, and each of its 77 columns is a sample.
- **Rows:**
  - `Population`, `Replicate`, `ID` and `time (in d)`;
  - 11 taxa;
  - `Clindamycin signal`.
- **Units:** densities in 10¹¹ 16S rDNA copies per cm³. About 9% of the values are zero.
- **Design:** 9 mice in 3 populations.

| Population | Treatment | Mice (`ID`) | Days sampled |
|---|---|---|---|
| 1 | *C. difficile* spores only (no infection developed) | 2, 5, 8 | 0, 2, 6, 13 |
| 2 | one dose of clindamycin | 1, 4, 7 | 0, 2, 3, 4, 5, 6, 7, 9, 12, 16, 23 |
| 3 | clindamycin, then *C. difficile* spores the next day | 3, 6, 9 | as population 2; mouse 3 lacks day 9 |

The clindamycin signal is 1 at day 0 for populations 2 and 3, and 0 otherwise.

**Taxon order differs between the two files:** `read_stein` matches them by name.
