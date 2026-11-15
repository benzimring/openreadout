# Provenance log — NMR processing, peak picking and integration

## 2026-09-23 — FID processing, phasing, baseline, peaks, integrals

**Corpus files used:** `nmrxiv-s275-1`, `-11`, `-13` … `-21` (TopSpin 4.3.0, CC0), `nmrxiv-s596-1`
(TOPSPIN 1.3, CC-BY-4.0), `nmrxiv-s846-50` (TopSpin 3.5, CC-BY-4.0), `nmrxiv-s501-7` and the
s501 VnmrJ directories (CC-BY-SA-4.0), `nmrxiv-s325` VnmrJ 4.0 (1H, 13C, DEPT), `nmrxiv-s837`
(TopSpin 3.7, FID-only), `nmrxiv-s200` JEOL Delta FIDs and their MestReNova JCAMP-DX exports
(CC0), `nmrxiv-s1243` JEOL FIDs, nmrglue test data `agilent_1d`, `nmrpy-test2`.

**Prior art consulted (read as documentation):**
- nmrglue (BSD-3-Clause, https://github.com/jjhelmus/nmrglue): `process/proc_base.py` (`em`, `gm`,
  `ps`, `fft` ordering, `fsh2`), `fileio/bruker.py` (`rm_dig_filter`, `bruker_dsp_table`),
  `fileio/jeol.py` (`X_OFFSET` as the carrier position in ppm), `fileio/varian.py`. Used to
  understand conventions; our implementation is our own.
- Chen, Weng, Goh, Garland, *An efficient algorithm for automatic phase correction of NMR spectra
  based on entropy minimization*, J. Magn. Reson. 158 (2002) 164–168 (ACME objective).
- Dietrich, Rüdel, Neumann, *Fast and precise automatic baseline correction of one- and
  two-dimensional NMR spectra*, J. Magn. Reson. 91 (1991) 1–11 (baseline recognition idea).

**What was inferred from what:**
- Bruker orientation and phase sign: processing the s275-1 FID and correlating with TopSpin's
  `pdata/1/1r` for the 16 combinations of spectrum reversal, phase sign, first-point handling and
  group-delay handling. Result: no conjugation, spectrum index k = FFT bin (N/2 − k) mod N (carrier
  at N/2; an earlier off-by-one ordering gave 0.4 correlation on 13C and was found by
  cross-correlating magnitudes: lag −1), correction e^{−i(PHC0 + PHC1·k/SI)}, first point × 0.5,
  group delay removed by rotation with the pre-response kept at negative time. Correlation
  1.000000 on s275-1/-11 and s846-50, scale 1.0006.
- `OFFSET = ((SFO1 − SF)·10⁶ + SW_h/2)/SF` reproduces the stored `OFFSET` (15.114770 vs 15.11477).
- Old consoles (s596-1: `AQ_mod` 1, DSPFVS 10, table group delay 59.083): the stored phases give
  an inverted spectrum and a residual first-order phase of 44°; not resolved (one file). Default
  phasing detects the inversion (negative fraction) and uses automatic phases instead.
- Varian: conjugating the FID puts TMS at 0.000 ppm and CHCl3 at 7.263 ppm in s325 1H with
  `reffrq` as the 0-ppm frequency; automatic phasing lands on −`rp` (−84.7 vs `rp` 85.0), and
  −`rp` gives ≤ 7 % negative signal on s325 1H/13C and s501 PROTON; `lp` could not be matched on
  the two files with non-zero `lp` (agilent_1d, nmrpy-test2).
- JEOL: the FID's first ~20 points are a digital-filter pre-response. The parameters `orders`
  (`"2 54 73"`) and `factors` (`"8  2"`) were read as stage count, filter orders and decimation
  factors; Σ (order−1)/2 / Π later factors = 19.656 points matches the observed onset, and the
  processed s200 qHNMR FID agrees with MestReNova's processed spectrum to 0.003 ppm (20/25 peaks).
  This reading of the two parameters is our inference from one instrument's files.
- DC offset: our unbaselined spectrum carries a constant the first point contributes that TopSpin's
  `1r` does not; the default first-order baseline removes it (peak comparisons use it).

## 2026-09-24 — rerun after the shared build directory was retired

All validation numbers were regenerated in a private target directory (the shared
`CARGO_BUILD_BUILD_DIR` could link another worktree's crates).

## 2026-09-26 — processing parity across 46 TopSpin experiments (Richard Zimring with Claude as assistant)

**Corpus files used:** 39 more Bruker experiments with TopSpin's own `pdata/1/1r` from 16 more
nmrXiv projects (P55, P72, P74, P80, P87, P88, P90, P98, P100, P101, P105, P107, P111, P122,
P146; CC-BY-4.0/CC0; `nmrxiv-parity/`, ids `nmrxiv-s<sample>-<experiment>`), with the 7 of
2026-09-23: TopSpin 3.5–4.3, 300–600 MHz, 1H, 13C (incl. DEPT and UDEFT), 31P, DSPFVS 10/20/21,
DIGMOD 1 and 3, FCOR 0/0.5/1, forward and backward linear prediction.

**Prior art consulted:** none new. nmrglue (BSD-3) is still only run as a comparison pipeline.

**What was measured and inferred (before changing the code):**
- With the stored phases our spectrum matched `1r` at correlation 1.000000 on every experiment
  acquired in DIGMOD 3 (integer `GRPDLY`), but only 0.9985–0.9993 on DIGMOD 1 experiments
  (`GRPDLY` 67.984–67.988). Fitting a residual phase against `1r` + i·`1i` gave, at the carrier,
  a zero-order offset of 1.1–2.8° with no first-order part, equal (within 0.4°) to
  180° × (⌈GRPDLY⌉ − GRPDLY) (2.13° for 67.9882, 2.48° for 67.9862, 2.84° for 67.9842). The same
  rule predicts the three "inverted" old-console experiments (DSPFVS 10, table delay 59.083 and
  similar): 180° × 0.917 = 165°, cos 165° = −0.966, the measured correlation (−0.963 … −0.966).
  Inference: the stored phases refer to a spectrum whose fractional group-delay correction
  carries a zero-order term π·(⌈G⌉ − G) relative to ours (ramp pivoted at the spectrum edge
  after a shift by ⌈G⌉ points instead of at the carrier after a shift by ⌊G⌋).
- FCOR (procs): the factor applied to the first FID point; FCOR 1 files kept the constant
  offset the halving removes, which showed as integral errors up to 48 % on the FID but not on
  `1r`. TDeff (procs): the number of FID points used (real + imaginary).

## 2026-09-26 — Magritek Spinsolve FIDs (Richard Zimring with Claude as assistant)

**Corpus files used:** `spinsolveproc-proton`, `spinsolveproc-t1`, the ten `spinsolve-zenodo20597567-*` series (the Spinsolve software's own processed spectra in `spectrum.pt1` / `*-Spectra.pt1`, one per FID row).

**Prior art consulted:** nmrglue 0.12 `fileio/spinsolve.py` (BSD-3-Clause): spectral width, observe frequency and carrier from `acqu.par`. Nothing else; no Magritek software or documentation.

**What was inferred from what:** fitting each FID row's transform to the plotted spectrum (a least-squares complex scale over the conjugated/non-conjugated, reversed/shifted variants): the plotted spectrum at ascending ppm point k equals the FFT of the conjugated FID (first point × 0.5, no apodization when `filter = "no"`) at frequency index k − n/2 times e^{+i·p0}, with p0 the first argument of `Phase(p0, p1)` in `processing.script` (= −`p0Phase` of `proc.par`), relative residual 2·10⁻⁷ on the 1D proton spectrum; `b1Freq` is the 0 ppm frequency. Hence `conjugate` = yes, stored φ0 = −p0, reference = `b1Freq`, no group delay. Unverified (no public example): a non-zero first-order phase, `filter = "yes"` apodization, `zf` > 1.

## 2026-10-11 — peak picking on water-suppressed metabolomics spectra

**Corpus files used:** `zenodo20308061-andrographis-mac50-3`, `zenodo20308061-andrographis-ult100-2` (Zenodo 20308061, CC-BY-4.0; all 27 experiments of the record were compared), `zenodo22712212-aiw-61` (Zenodo 22712212, CC-BY-4.0; 61 experiments compared), and the 1H experiments already in the corpus. **Prior art consulted:** none new.

**What was measured:** the highest point of TopSpin's `1r` in five windows against `analyze nmr-peaks --range-ppm` (`main_peak`), on the stored spectrum and on the FID processed here. The stored spectrum agrees everywhere except where two lines differ by less than 0.5 % in height (parabolic interpolation reorders them) and where the window cuts a line: a line just outside the window was reported because the index window was rounded outwards. Peaks are now kept only when their interpolated shift is inside the window. The FID processed here matches `1r` at correlation 0.99999 but sits on a baseline about 0.3 % of the tallest line lower than TopSpin's in these DIGMOD 3 data sets (the difference is in the first time-domain point, and no single first-point factor explains all files). This reorders weak lines of nearly equal height in 4 of 27 windows; it is listed as a known gap.

## 2026-10-11 — FIDs without stored referencing

**Corpus files used:** every Bruker experiment on disk with `pdata/1/procs` (nmrXiv and the two Zenodo bundles, about 290), comparing `procs` SF with `acqus` BF1. **Prior art consulted:** none.
**What was measured:** the referencing correction (SF − BF1)/SF is within 0.2 ppm of 0 in most experiments and 0.44 to 1.41 ppm (either sign) in 20. An FID processed without `procs` can only take BF1 as 0 ppm, so its shifts are off by that much. The FID processing now says so in a note when no `procs` gives SF.

## 2026-11-13 — cross peaks of stored 2-D spectra (`analyze nmr-cross-peaks`) (Ben Zimring with Claude as assistant)

**Corpus files used:** the 2D experiments with TopSpin's processed `pdata/1/2rr` in `nmrxiv-s73-bundle` (COSY, edited HSQC, HMBC), `nmrxiv-s572-bundle` (COSY, NOESY, edited HSQC, HMBC), `nmrxiv-s1407-bundle`, `nmrxiv-s1617-bundle` (COSY, HSQC, HMBC, TOCSY), `nmrxiv-s1517-bundle` (COSY, HMQC, HMBC, NOESY), `nmrxiv-s2331-bundle` (COSY, HSQC, HMBC, NOESY), `nmrxiv-s314-bundle` and `nmrxiv-s501-6` (ROESY, NOESY), all CC-BY-4.0 or CC-BY-SA-4.0; the JCAMP-DX 2D spectra `nmrxiv-s37-hsqc-jdx` (TopSpin export, CC-BY-4.0), `nmrxiv-s1-cosy-dx` (CC-BY-4.0) and `zenodo10800686-17-hsqc` (MestReNova export, CC0). Ground truth from the depositors, independent of our picker: the correlations picked in NMRium in `nmrxiv-s73-bundle/Nicotinic acid.nmrium` (4 HSQC, 3 COSY and 12 HMBC zones) and the assigned correlations in `nmrxiv-s572-bundle/Classics_Curcumin.nmredata.sdf` (NMReDATA tags `NMREDATA_2D_13C_1J_1H`, `NMREDATA_2D_13C_NJ_1H`, `NMREDATA_2D_1H_NJ_1H`, with the shifts of `NMREDATA_ASSIGNMENT`).

**Prior art consulted:** none. The detector is ours: points above a noise threshold, grouped into peaks by their topographic prominence (the drop from a maximum to the highest saddle that joins it to a higher maximum), which is the 2D form of the prominence rule `pick_peaks` already uses. NMRium and NMReDATA files were read as data (JSON and SDF text), not as code.

**What was inferred from what:**
- Orientation of a processed 2D spectrum: sweep (row) *j* of a Bruker `pdata/<n>` trace is the F1 point *j* of `extra.sweep_axis.axis` (`OFFSET` of `proc2s` at row 0, falling), and sample *k* is the F2 point *k* of `extra.axis`; checked by locating all 11 depositor HSQC correlations of nicotinic acid and curcumin.
- JCAMP-DX 2D spectra (`##NTUPLES= nD NMR SPECTRUM`): the axes are the header's `##FIRST=`/`##LAST=`/`##VAR_DIM=` per `##SYMBOL=` column, in Hz; the frequency that converts a dimension to ppm is its entry in the NTUPLES `##.OBSERVE FREQUENCY=` list when the list has one per dimension (MestReNova), else that dimension's `##$SF=` (TopSpin exports carry the acquisition parameters, then one `procs` block per dimension, F2 first, each with `##$SF=`, `##$AXNUC=` and `##$OFFSET=`).
- Experiment kind from the pulse program name stored in the file (`PULPROG`, `##.PULSE SEQUENCE=`), by the names seen in these files: `hsqc…` HSQC (`hsqced…` multiplicity-edited), `hmqc…` HMQC, `hmbc…` HMBC, `cosy…` COSY, `mlev…`/`dipsi…`/`tocsy…` TOCSY, `noesy…` NOESY, `roesy…` ROESY; otherwise from the two nuclei (heteronuclear or homonuclear). Tagged as inferred.
- Edited-HSQC signs: in the curcumin HSQC the depositor's own correlation list gives the enol CH (C10, 100.8 ppm) a negative intensity and the other CH and the OCH₃ positive ones, so the sign cannot be mapped to CH₂ without knowing how the spectrum was phased. Cross peaks report their sign only.

**What was measured (while choosing the detector's rules):**
- With one global noise SD, the depositor correlations were all found but COSY spectra listed 200 to 1100 cross peaks, most in the F2 columns of strong lines (t1 noise). A per-column noise lowered that; a column MAD ignores sparse non-uniform-sampling artifacts, so the column noise is the 80th percentile of |x − median| over 1.28. The 90th percentile lost two curcumin HMBC correlations in busy columns; the 80th kept them. Thresholding every point against its column noise split strong lines whose own F1 profile fills a column (synthetic 64-row spectra), and leaving all above-threshold points out of the column noise let t1 ridges through (curcumin COSY 240 cross peaks instead of 95); points are therefore grouped against the spectrum's SD, and each maximum is tested against its column's noise without its own rows (a sixteenth of F1 on each side, at least 8).
- Around strong HSQC lines, weaker cross peaks of 2 to 15 % of the line's height sit in its F2 column 1.1 to 2.1 F1 widths away (t1 truncation ripples; real HMBC correlations in the same column were 2.2 widths away and more) and in its F1 row in mirror pairs at ±52 Hz (`nmrxiv-s73`), ±80 Hz (`nmrxiv-s572`) and ±342 Hz (`nmrxiv-s37`) from it (decoupling sidebands). Hence `ripple_fraction` 0.2, `RIPPLE_WIDTHS` 2 and `SIDEBAND_HZ` 1000 with a mirror twin required.
- Phase-sensitive COSY cross peaks are antiphase, so multiplet components of both signs are merged into one cross peak.
- Homonuclear cross peaks: 11 of 15 depositor COSY correlations have a cross peak at the swapped position, 11 of 151 other off-diagonal cross peaks do; reported as `mirror_partner`, not used to filter.
- JCAMP-DX axes: in `nmrxiv-s37-hsqc-jdx` and `nmrxiv-s1-cosy-dx` (TopSpin exports) `##FIRST=`/`##LAST=` give F2 from 6403.99 and 2246.01 Hz down to 0, while the F2 `##$OFFSET=` is 12.71422 and 8.42538 ppm; with the offset the ethylbenzene CH₃ and aromatic CH land at 1.19 / 16.0 and 7.28 / 128.8 ppm and the COSY diagonal on F1 = F2, without it 3.29 and 0.94 ppm off. F1 agrees with its `##$OFFSET=` in both. In `zenodo10800686-17-hsqc` (MestReNova) the page labels step 94.4 Hz and the header 47.2 Hz over 512 F1 points; with the header the aromatic CH fall at 128 to 140 ppm and an OCH₃ at 54.5 ppm, with the labels at 26 to 50 and −121 ppm.
- Result on the depositors' correlations: 49 of 52 found (details in `docs/formats/nmr-processing.md`, "2D cross peaks").

## 2026-11-14 — cross peaks: files without `AXNUC`, long tables, wrong referencing (Ben Zimring with Claude as assistant)

**Why:** agents that ran `analyze nmr-cross-peaks` on 15 nmrXiv data sets outside the corpus met three problems: a TopSpin 1.3 data set without `##$AXNUC=`, tables that reached `max_cross_peaks`, and an HSQC whose `##$OFFSET=` values are near −200000 ppm.

**Corpus files used:** the stored 2-D spectra named in the 2026-11-13 entry, and copies of `nmrxiv-s73-bundle/5` that the tests edit in a temporary directory (`##$AXNUC=` removed, `##$OFFSET=` replaced). Ground truth as before: the depositors' correlations of nicotinic acid and curcumin.

**Also inspected, not added:** the HSQC spectra of the 15 nmrXiv samples of that benchmark (S41, S43, S47, S48, S273, S1258, S1266, S1272, S1275, S1393, S1491, S1492, S1564, S1569, S2308; CC-BY-4.0, CDLA-Permissive-1.0), for the number of cross peaks; of S48 also the `##$OFFSET=` values, and of S2308 the missing `##$AXNUC=`. The HSQC of S273 lists one methyl as eight cross peaks of alternating sign and was opened for that: its F1 carries truncation ripples further out than `RIPPLE_WIDTHS` reaches, which this change leaves as they are. No rule or constant comes from these spectra.

**Prior art consulted:** none.

**What was inferred from what:**
- A processed Bruker 2-D spectrum without `##$AXNUC=` in either dimension is still a chemical-shift spectrum when `proc2s` gives `OFFSET`, `SW_p` and `SF`. Its nuclei are `NUC1` of `acqus` (F2) and of `acqu2s` (F1).
