//! Cross peaks of stored 2D spectra (`analyze nmr-cross-peaks`) against correlations the
//! depositors picked themselves, which no picker of ours produced:
//!
//! - `nmrxiv-s73-bundle` (nicotinic acid, CC-BY-4.0): the zones of `Nicotinic acid.nmrium`
//!   (NMRium) for the edited HSQC (experiment 5), COSY (4) and HMBC (6);
//! - `nmrxiv-s572-bundle` (curcumin, CC-BY-4.0): the assigned correlations of
//!   `Classics_Curcumin.nmredata.sdf` (NMReDATA `NMREDATA_2D_13C_1J_1H`, `NMREDATA_2D_13C_NJ_1H`,
//!   `NMREDATA_2D_1H_NJ_1H`, at the shifts of `NMREDATA_ASSIGNMENT`) for the edited HSQC (5),
//!   HMBC (6) and COSY (3).
//!
//! The correlations are copied below as (F2, F1) shifts. A correlation is found when a cross
//! peak lies within 0.02 ppm in F2 and 0.6 ppm (13C) or 0.03 ppm (1H) in F1; COSY correlations
//! match either orientation. Recall is asserted, and the number of cross peaks listed is bounded
//! so that a change letting artifacts through shows up. Also checked: the JCAMP-DX 2D axes
//! (TopSpin's F2 from 0 Hz moved to its `$OFFSET`, MestReNova's F1 from the NTUPLES header),
//! that `analyze nmr-peaks` refuses a 2D spectrum, and, on copies of the nicotinic acid HSQC
//! made in a temporary directory, a data set without `##$AXNUC=` (as TopSpin 1.3 writes them)
//! and one whose `##$OFFSET=` values are far outside any spectrum.
//!
//! Run: `cargo test -p openreadout-corpus-tests --features corpus --test nmr_cross_peaks --
//! --nocapture` (`OPENREADOUT_CORPUS_DIR` overrides `corpus/files`).
#![cfg(feature = "corpus")]

use std::path::{Path, PathBuf};

use openreadout_signal::nmr::{
    self, CrossPeakReport, CrossPeakRequest, NmrRequest, analyze_cross_peaks,
};

#[path = "support/registry.rs"]
mod registry;

fn files_dir() -> PathBuf {
    std::env::var_os("OPENREADOUT_CORPUS_DIR").map_or_else(
        || Path::new(env!("CARGO_MANIFEST_DIR")).join("../../corpus/files"),
        PathBuf::from,
    )
}

/// The report for `rel` under the corpus directory, or `None` when it is not downloaded.
fn cross_peaks(rel: &str) -> Option<CrossPeakReport> {
    let path = files_dir().join(rel);
    if !path.exists() {
        eprintln!("skip {rel}: not downloaded (cargo xtask corpus fetch --tier standard)");
        return None;
    }
    let (_, mut ds) = registry::registry().open(&path).unwrap();
    let info = ds.info().unwrap();
    Some(analyze_cross_peaks(ds.as_mut(), &info, &CrossPeakRequest::default()).unwrap())
}

/// The report for the data set at `path`.
fn report(path: &Path) -> CrossPeakReport {
    let (_, mut ds) = registry::registry().open(path).unwrap();
    let info = ds.info().unwrap();
    analyze_cross_peaks(ds.as_mut(), &info, &CrossPeakRequest::default()).unwrap()
}

fn copy_dir(from: &Path, to: &Path) {
    std::fs::create_dir_all(to).unwrap();
    for entry in std::fs::read_dir(from).unwrap() {
        let entry = entry.unwrap();
        let target = to.join(entry.file_name());
        if entry.file_type().unwrap().is_dir() {
            copy_dir(&entry.path(), &target);
        } else {
            std::fs::copy(entry.path(), &target).unwrap();
        }
    }
}

/// A copy of the experiment `rel` in a temporary directory with each line of its processing
/// parameter files (`procs`, `proc2s`, `proc`, `proc2` of every `pdata/<n>`) passed through
/// `edit` (which gets the file name and the line; `None` drops the line). The files are edited
/// as bytes, since parameter files are not always UTF-8. `None` when the corpus file is not
/// downloaded.
fn edited_copy(
    rel: &str,
    edit: impl Fn(&str, &[u8]) -> Option<Vec<u8>>,
) -> Option<(tempfile::TempDir, PathBuf)> {
    let source = files_dir().join(rel);
    if !source.exists() {
        eprintln!("skip {rel}: not downloaded (cargo xtask corpus fetch --tier standard)");
        return None;
    }
    let dir = tempfile::tempdir().unwrap();
    let copy = dir.path().join("experiment");
    copy_dir(&source, &copy);
    for procno in std::fs::read_dir(copy.join("pdata")).unwrap() {
        let procno = procno.unwrap().path();
        for name in ["procs", "proc2s", "proc", "proc2"] {
            let file = procno.join(name);
            let Ok(bytes) = std::fs::read(&file) else {
                continue;
            };
            let mut out = Vec::with_capacity(bytes.len());
            for line in bytes.split(|b| *b == b'\n') {
                if let Some(l) = edit(name, line) {
                    out.extend_from_slice(&l);
                    out.push(b'\n');
                }
            }
            std::fs::write(&file, out).unwrap();
        }
    }
    Some((dir, copy))
}

/// Correlations of `truth` found among the cross peaks, and the ones missed.
fn recall(r: &CrossPeakReport, truth: &[(f64, f64)]) -> (usize, Vec<(f64, f64)>) {
    let homonuclear = r.f1.nucleus == r.f2.nucleus;
    let tol1 = if homonuclear { 0.03 } else { 0.6 };
    let near = |f2: f64, f1: f64| {
        r.cross_peaks
            .iter()
            .any(|p| (p.f2_ppm - f2).abs() <= 0.02 && (p.f1_ppm - f1).abs() <= tol1)
    };
    let mut missed = Vec::new();
    for &(f2, f1) in truth {
        if !(near(f2, f1) || homonuclear && near(f1, f2)) {
            missed.push((f2, f1));
        }
    }
    (truth.len() - missed.len(), missed)
}

struct Case {
    rel: &'static str,
    kind: &'static str,
    truth: &'static [(f64, f64)],
    /// Correlations that must be found.
    min_found: usize,
    /// Most cross peaks listed.
    max_listed: usize,
}

const NICOTINIC_HSQC: &[(f64, f64)] = &[
    (9.139, 149.869),
    (8.748, 152.313),
    (8.425, 137.79),
    (7.577, 123.702),
];
const NICOTINIC_COSY: &[(f64, f64)] = &[(8.751, 7.586), (8.43, 7.582), (9.141, 8.436)];
const NICOTINIC_HMBC: &[(f64, f64)] = &[
    (9.145, 166.239),
    (8.431, 166.252),
    (9.146, 152.13),
    (8.755, 149.661),
    (8.427, 149.828),
    (8.435, 152.059),
    (7.584, 152.157),
    (9.145, 137.664),
    (8.752, 137.767),
    (9.145, 127.141),
    (8.75, 123.75),
    (7.586, 127.168),
];
const CURCUMIN_HSQC: &[(f64, f64)] = &[
    (7.538, 140.716),
    (7.136, 123.126),
    (7.309, 111.226),
    (6.814, 115.672),
    (6.737, 121.034),
    (6.049, 100.763),
    (3.834, 55.706),
];
const CURCUMIN_HMBC: &[(f64, f64)] = &[
    (3.834, 147.969),
    (6.737, 100.763),
    (6.737, 126.361),
    (6.814, 126.361),
    (7.136, 111.226),
    (7.538, 111.226),
    (7.309, 123.126),
    (7.538, 123.126),
    (7.538, 121.034),
    (7.538, 126.361),
    (7.309, 126.361),
    (6.049, 121.034),
    (7.309, 140.716),
    (7.136, 140.716),
    (6.814, 147.969),
    (6.737, 183.137),
    (6.049, 183.137),
    (7.538, 183.137),
    (7.309, 149.146),
    (7.136, 149.146),
    (7.309, 147.969),
];
// the three correlations of the enol H (6.049 ppm) are long-range couplings at S/N 3 to 8
// against the t1 noise of their columns: not required
const CURCUMIN_COSY: &[(f64, f64)] = &[
    (7.538, 6.737),
    (7.136, 6.814),
    (6.737, 6.814),
    (6.049, 7.538),
    (6.049, 7.309),
    (6.049, 6.737),
];

const CASES: &[Case] = &[
    Case {
        rel: "nmrxiv-s73-bundle/5",
        kind: "hsqc",
        truth: NICOTINIC_HSQC,
        min_found: 4,
        max_listed: 12,
    },
    Case {
        rel: "nmrxiv-s73-bundle/4",
        kind: "cosy",
        truth: NICOTINIC_COSY,
        min_found: 3,
        max_listed: 60,
    },
    Case {
        rel: "nmrxiv-s73-bundle/6",
        kind: "hmbc",
        truth: NICOTINIC_HMBC,
        min_found: 12,
        max_listed: 30,
    },
    Case {
        rel: "nmrxiv-s572-bundle/5",
        kind: "hsqc",
        truth: CURCUMIN_HSQC,
        min_found: 7,
        max_listed: 35,
    },
    Case {
        rel: "nmrxiv-s572-bundle/6",
        kind: "hmbc",
        truth: CURCUMIN_HMBC,
        min_found: 19,
        max_listed: 80,
    },
    Case {
        rel: "nmrxiv-s572-bundle/3",
        kind: "cosy",
        truth: CURCUMIN_COSY,
        min_found: 3,
        max_listed: 200,
    },
];

#[test]
fn depositor_correlations_are_found() {
    for c in CASES {
        let Some(r) = cross_peaks(c.rel) else {
            continue;
        };
        let (found, missed) = recall(&r, c.truth);
        eprintln!(
            "{:22} {:5} found {found}/{} listed {:3} (noise SD {:.3e}) missed {missed:?}",
            c.rel,
            r.experiment.kind,
            c.truth.len(),
            r.cross_peaks.len(),
            r.noise_sd
        );
        assert_eq!(r.experiment.kind, c.kind, "{}", c.rel);
        assert_eq!(r.experiment.kind_source, "pulse_program", "{}", c.rel);
        assert!(found >= c.min_found, "{}: missed {missed:?}", c.rel);
        assert!(
            r.cross_peaks.len() <= c.max_listed,
            "{}: {} cross peaks listed",
            c.rel,
            r.cross_peaks.len()
        );
        if c.kind == "hsqc" {
            assert!(r.experiment.multiplicity_edited, "{}", c.rel);
            assert_eq!(r.f1.nucleus.as_deref(), Some("13C"), "{}", c.rel);
        }
        if c.kind == "cosy" {
            assert!(
                r.cross_peaks.iter().all(|p| p.on_diagonal.is_some()),
                "{}",
                c.rel
            );
        }
    }
}

#[test]
fn jcamp_dx_axes_are_placed() {
    // TopSpin export of an ethylbenzene HSQC: F2 is written from 0 Hz; the aromatic CH at 7.28 /
    // 128.8 ppm only lands there after the move to the F2 $OFFSET
    if let Some(r) = cross_peaks("nmrxiv-s37-hsqc-jdx.jdx") {
        assert!((r.f2.axis.first_ppm - 12.714_22).abs() < 1e-6);
        assert!(
            r.notes.iter().any(|n| n.contains("OFFSET")),
            "{:?}",
            r.notes
        );
        assert!(
            r.cross_peaks
                .iter()
                .any(|p| (p.f2_ppm - 7.283).abs() < 0.02 && (p.f1_ppm - 128.76).abs() < 0.6)
        );
        assert!(
            r.cross_peaks
                .iter()
                .any(|p| (p.f2_ppm - 1.189).abs() < 0.02 && (p.f1_ppm - 16.05).abs() < 0.6)
        );
    }
    // TopSpin export of a COSY: after the same move the diagonal lies on F1 = F2
    if let Some(r) = cross_peaks("nmrxiv-s1-cosy-dx.dx") {
        assert_eq!(r.experiment.kind, "cosy");
        let diagonal = r
            .cross_peaks
            .iter()
            .filter(|p| p.on_diagonal == Some(true))
            .count();
        assert!(diagonal >= 5, "{diagonal} diagonal peaks");
    }
    // MestReNova export: the F1 page labels are spaced twice the header's step; the header
    // gives the OCH3 at 3.69 / 54.5 ppm and aromatic CH near 128 ppm
    if let Some(r) = cross_peaks("zenodo10800686-17-hsqc.jdx") {
        assert!(
            (r.f1.axis.last_ppm() - (-9.545)).abs() < 0.01,
            "{:?}",
            r.f1.axis
        );
        assert!(
            r.cross_peaks
                .iter()
                .any(|p| (p.f2_ppm - 3.689).abs() < 0.02 && (p.f1_ppm - 54.5).abs() < 0.6)
        );
        assert!(
            r.cross_peaks
                .iter()
                .all(|p| p.f1_ppm > 100.0 || p.f2_ppm < 5.0)
        );
    }
}

#[test]
fn nmr_peaks_refuses_a_2d_spectrum() {
    let path = files_dir().join("nmrxiv-s73-bundle/5");
    if !path.exists() {
        return;
    }
    let (_, mut ds) = registry::registry().open(&path).unwrap();
    let info = ds.info().unwrap();
    let e = nmr::analyze(ds.as_mut(), &info, &NmrRequest::default()).unwrap_err();
    assert!(e.to_string().contains("2D spectrum"), "{e}");
    assert!(
        format!("{e:?}").contains("nmr-cross-peaks"),
        "the hint names the 2D command: {e:?}"
    );
}

#[test]
fn spectra_without_axnuc_take_the_acquisition_nuclei() {
    // TopSpin 1.3 writes no ##$AXNUC= in procs and proc2s
    let Some((_dir, copy)) = edited_copy("nmrxiv-s73-bundle/5", |_, line| {
        (!line.starts_with(b"##$AXNUC=")).then(|| line.to_vec())
    }) else {
        return;
    };
    let original = report(&files_dir().join("nmrxiv-s73-bundle/5"));
    let r = report(&copy);
    assert_eq!(r.f2.nucleus.as_deref(), Some("1H"));
    assert_eq!(r.f1.nucleus.as_deref(), Some("13C"));
    assert!(
        r.notes.iter().any(|n| n.contains("acquisition parameters")),
        "{:?}",
        r.notes
    );
    assert_eq!(r.experiment.kind, "hsqc");
    assert_eq!(r.cross_peaks.len(), original.cross_peaks.len());
    let (found, missed) = recall(&r, NICOTINIC_HSQC);
    assert_eq!(found, NICOTINIC_HSQC.len(), "missed {missed:?}");
    // and it is still a 2D spectrum to `analyze nmr-peaks`
    let (_, mut ds) = registry::registry().open(&copy).unwrap();
    let info = ds.info().unwrap();
    let e = nmr::analyze(ds.as_mut(), &info, &NmrRequest::default()).unwrap_err();
    assert!(e.to_string().contains("2D spectrum"), "{e}");
}

#[test]
fn offsets_outside_any_spectrum_get_a_note() {
    // a processed HSQC stored with ##$OFFSET= near -200000 ppm in both dimensions
    let Some((_dir, copy)) = edited_copy("nmrxiv-s73-bundle/5", |name, line| {
        Some(if line.starts_with(b"##$OFFSET=") {
            let offset = if name.starts_with("proc2") {
                -199_859.7
            } else {
                -200_016.8
            };
            format!("##$OFFSET= {offset}").into_bytes()
        } else {
            line.to_vec()
        })
    }) else {
        return;
    };
    let r = report(&copy);
    for dim in ["F2 axis", "F1 axis"] {
        assert!(
            r.notes
                .iter()
                .any(|n| n.contains(dim) && n.contains("referencing")),
            "{dim}: {:?}",
            r.notes
        );
    }
    let original = report(&files_dir().join("nmrxiv-s73-bundle/5"));
    assert!(
        !original.notes.iter().any(|n| n.contains("referencing")),
        "{:?}",
        original.notes
    );
}
