"""Extract the published gLV fit of Stein et al. (2013, PLoS Comput Biol 9:e1003388) from
sheet MmuE of their Dataset S1 into a tidy CSV: one row per genus (sorted by susceptibility, as
in the sheet), columns `genus`, the 11 interaction coefficients M[genus, j] (effect of genus j
on the row genus, named by j), `growth` (mu) and `susceptibility` (epsilon, clindamycin).
Usage: python3 -I extract_mmue.py Dataset_S1.xlsx stein2013_glv.csv"""
import csv, sys
import openpyxl

wb = openpyxl.load_workbook(sys.argv[1], read_only=True, data_only=True)
rows = list(wb["MmuE"].iter_rows(values_only=True))
header = rows[13]
genera = [g for g in header[1:12]]
body = rows[2:13]
assert [r[0] for r in body] == genera, "row and column genera differ"
with open(sys.argv[2], "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["genus"] + genera + ["growth", "susceptibility"])
    for r in body:
        w.writerow([r[0]] + [repr(float(x)) for x in r[1:12]] + [repr(float(r[12])), repr(float(r[13]))])
