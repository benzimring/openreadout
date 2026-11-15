//! Cross peaks of stored 2D spectra: what `openreadout analyze nmr-cross-peaks` and the
//! `openreadout_nmr_cross_peaks` MCP tool run. `docs/formats/nmr-processing.md` ("2D cross
//! peaks") gives the rules in full, and the provenance log the measurements behind them.
//!
//! **Input.** A 2D spectrum the vendor software processed and stored: a trace with a
//! chemical-shift axis (F2) and one sweep per point of the indirect dimension F1 ([`axes`]).
//! Only the real part is read. FIDs (`ser`) are not processed here.
//!
//! **Noise.** [`noise_and_baseline`] takes the noise SD and the baseline from tiles of the
//! spectrum, so that tiles holding signal raise neither. t1 noise runs along F1, so each F2
//! point also has the noise of its own F1 column. A peak's S/N is measured against the noise
//! of its column without its own rows, and against the spectrum's SD when that is larger.
//!
//! **Peaks.** Positive and (unless `positive_only`) negative signal are searched separately.
//! Points above `min_snr` times the SD are visited from the tallest down and joined to their
//! visited 8-neighbors. Where two groups meet, the lower maximum stays a peak of its own when
//! its prominence is at least `min_prominence_snr` times its column's noise. This is the 2D
//! form of the prominence rule of [`pick_peaks`](super::peaks::pick_peaks). Positions come
//! from a parabola in each dimension, widths at half height from the F2 row and the F1 column.
//!
//! **Cross peaks.** Peaks of either sign closer than `merge_hz` in both dimensions are one
//! cross peak (the components of a multiplet). It is placed at its tallest peak, and `maxima`
//! counts its peaks.
//!
//! **Artifacts.** Truncation ripples of a taller cross peak are left out, and on heteronuclear
//! spectra its decoupling sidebands too (`drop_ripples`).
//!
//! **Homonuclear spectra.** Cross peaks on F1 = F2 are flagged `on_diagonal`. The others get
//! `mirror_partner`: whether a cross peak lies at the swapped position.
//!
//! **Experiment.** The kind comes from the pulse program name ([`experiment_kind`]), else from
//! the two nuclei, unless the request names it. A multiplicity-edited HSQC gives one
//! multiplicity group negative heights. Which group depends on how the spectrum was phased, so
//! only the sign is reported.

use std::collections::BTreeMap;

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use openreadout_core::model::{FileInfo, TraceInfo};
use openreadout_core::{Dataset, Error, Result};

use crate::nmr::MULTIPLET_GAP_HZ;
use crate::nmr::peaks::{half_width_of, median_in_place};
use crate::nmr::process::PpmAxis;
use crate::nmr::source::{self, is_fid, spectrum_axis};

/// Two-sided 80 % point of the standard normal distribution (`column_noise`).
const COLUMN_QUANTILE_Z: f64 = 1.281_551_6;
/// Tiles per dimension used by [`noise_and_baseline`].
const NOISE_TILES: usize = 32;
/// A width that could not be measured counts as this many points.
const UNMEASURED_WIDTH_POINTS: f64 = 4.0;
/// Largest spectrum read whole, in points (4096 × 4096).
pub const MAX_POINTS: usize = 1 << 24;
/// Default `merge_hz`: the multiplet gap of the 1-D integrals.
pub const DEFAULT_MERGE_HZ: f64 = MULTIPLET_GAP_HZ;
/// Default `ripple_fraction`.
pub const DEFAULT_RIPPLE_FRACTION: f64 = 0.2;
/// Truncation ripples are left out up to this many F1 widths from a cross peak.
pub const RIPPLE_WIDTHS: f64 = 2.0;
/// Decoupling sidebands are looked for up to this far from a cross peak in F2, Hz.
pub const SIDEBAND_HZ: f64 = 1000.0;
/// A table with more cross peaks than this gets a note with the counts at higher S/N.
pub const LONG_TABLE_ROWS: usize = 100;
/// The multiples of `min_snr` that note counts cross peaks at.
pub const SNR_LEVELS: [f64; 3] = [2.0, 5.0, 10.0];

/// What [`pick_cross_peaks`] found.
#[derive(Debug, Clone, Default)]
pub struct Picked {
    /// Cross peaks, high F2 shift first (at most `max_cross_peaks`).
    pub cross_peaks: Vec<CrossPeak>,
    /// Cross peaks found before `max_cross_peaks`.
    pub count: usize,
    /// Peaks left out as ripples or sidebands of a taller cross peak.
    pub ripples: usize,
    /// Cross peaks (before `max_cross_peaks`) with S/N of at least `min_snr` times each of
    /// [`SNR_LEVELS`].
    pub at_snr_levels: [usize; 3],
    /// Noise SD of the whole spectrum.
    pub noise_sd: f64,
    /// Baseline (median level).
    pub baseline: f64,
}

/// Thresholds and limits of cross-peak picking.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct CrossPeakOptions {
    /// Minimum S/N: height above the baseline over the noise of the peak's F1 column.
    pub min_snr: f64,
    /// Minimum prominence (height above the saddle to a higher peak), in noise SDs.
    pub min_prominence_snr: f64,
    /// Minimum height as a fraction of the tallest point of the searched region.
    pub min_height_fraction: f64,
    /// Search positive signal only (negative cross peaks are searched by default).
    pub positive_only: bool,
    /// Only pick between two F2 shifts, ppm.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub f2_range_ppm: Option<(f64, f64)>,
    /// Only pick between two F1 shifts, ppm.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub f1_range_ppm: Option<(f64, f64)>,
    /// Peaks closer than this in both dimensions are one cross peak, Hz.
    pub merge_hz: f64,
    /// Leave out peaks at most this fraction of the height of a taller cross peak in its F2
    /// column (truncation ripples) or F1 row (decoupling sidebands); 0 keeps them.
    pub ripple_fraction: f64,
    /// Keep at most this many cross peaks (the highest S/N).
    pub max_cross_peaks: usize,
}

impl Default for CrossPeakOptions {
    fn default() -> Self {
        CrossPeakOptions {
            min_snr: 10.0,
            min_prominence_snr: 5.0,
            min_height_fraction: 0.0,
            positive_only: false,
            f2_range_ppm: None,
            f1_range_ppm: None,
            merge_hz: DEFAULT_MERGE_HZ,
            ripple_fraction: DEFAULT_RIPPLE_FRACTION,
            max_cross_peaks: 500,
        }
    }
}

/// What [`analyze_cross_peaks`] should do.
#[derive(Debug, Clone, Default)]
pub struct CrossPeakRequest {
    /// Trace to use (default: the first stored 2D spectrum).
    pub trace: Option<u32>,
    /// Experiment kind to report instead of the one inferred from the file.
    pub experiment: Option<String>,
    /// Thresholds.
    pub options: CrossPeakOptions,
}

/// One dimension of a 2D spectrum.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct SpectrumDimension {
    /// Nucleus (`1H`, `13C`), when the file names it.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nucleus: Option<String>,
    /// Chemical-shift axis: point `k` at `first_ppm + k · step_ppm`.
    pub axis: PpmAxis,
}

/// What kind of 2D experiment the spectrum is.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct ExperimentInfo {
    /// `hsqc`, `hmqc`, `hmbc`, `cosy`, `tocsy`, `noesy`, `roesy`, or, when the pulse program
    /// does not say, `heteronuclear`, `homonuclear` or `unknown`.
    pub kind: String,
    /// A multiplicity-edited HSQC (pulse program name `hsqced…`): one multiplicity group has
    /// negative heights.
    pub multiplicity_edited: bool,
    /// The pulse program name stored in the file.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pulse_program: Option<String>,
    /// Where `kind` came from: `pulse_program`, `nuclei` or `request`.
    pub kind_source: String,
}

/// One cross peak.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct CrossPeak {
    /// Shift in the direct dimension, ppm.
    pub f2_ppm: f64,
    /// Shift in the indirect dimension, ppm.
    pub f1_ppm: f64,
    /// Stored value at the interpolated maximum (negative for a negative cross peak).
    pub height: f64,
    /// |height − baseline| over the noise of its F1 column (without its own rows; at least the
    /// spectrum's noise SD).
    pub snr: f64,
    /// Width at half height along F2, Hz (none when it cannot be measured).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub f2_width_hz: Option<f64>,
    /// Width at half height along F1, Hz.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub f1_width_hz: Option<f64>,
    /// F2 point of the maximum.
    pub f2_index: usize,
    /// F1 point (sweep) of the maximum.
    pub f1_index: usize,
    /// Peaks merged into this cross peak (multiplet components).
    pub maxima: u32,
    /// Homonuclear spectra: the cross peak lies on the diagonal F1 = F2.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub on_diagonal: Option<bool>,
    /// Homonuclear spectra, off the diagonal: another cross peak lies at the mirror position
    /// (F1 and F2 swapped). Real correlations usually have one; t1 noise rarely does.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mirror_partner: Option<bool>,
}

/// Output of `analyze nmr-cross-peaks`.
#[derive(Debug, Clone, Default, Serialize, Deserialize, JsonSchema)]
pub struct CrossPeakReport {
    /// The input path.
    pub path: String,
    /// Format id.
    pub format: String,
    /// Trace used (see `info` → `traces[]`).
    pub trace: u32,
    /// Its name.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub trace_name: Option<String>,
    /// The kind of experiment.
    pub experiment: ExperimentInfo,
    /// Direct dimension.
    pub f2: SpectrumDimension,
    /// Indirect dimension.
    pub f1: SpectrumDimension,
    /// Median level of the spectrum (heights are measured from it).
    pub baseline: f64,
    /// Robust noise standard deviation (units of `height`).
    pub noise_sd: f64,
    /// Thresholds used.
    pub options: CrossPeakOptions,
    /// Number of cross peaks found (before `max_cross_peaks`).
    pub cross_peak_count: usize,
    /// Cross peaks, high F2 shift first, then high F1 shift.
    pub cross_peaks: Vec<CrossPeak>,
    /// Remarks.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub notes: Vec<String>,
}

fn num(v: &Value) -> Option<f64> {
    source::num(Some(v))
}

fn clean_nucleus(s: &str) -> Option<String> {
    let s = s
        .trim()
        .trim_start_matches('^')
        .trim_matches(|c| c == '<' || c == '>');
    (!s.is_empty() && !s.eq_ignore_ascii_case("off")).then(|| s.to_string())
}

/// The dimensions of a JCAMP-DX 2D spectrum and the parameters that place its axes, as the
/// reader reports them in `extra.dimensions` and `extra.dimension_parameters`.
#[derive(Debug, Clone, Default, PartialEq)]
struct JcampHeader {
    /// `(first, last, points, unit)` of F2 from the NTUPLES header.
    f2: Option<(f64, f64, usize, String)>,
    /// The same for F1.
    f1: Option<(f64, f64, usize, String)>,
    /// F2 observe frequency, MHz, from a per-column `##.OBSERVE FREQUENCY=` list.
    f2_mhz: Option<f64>,
    /// F1 observe frequency, MHz, from the same list.
    f1_mhz: Option<f64>,
    /// F2 nucleus from a per-column `##.NUCLEUS=` list.
    f2_nucleus: Option<String>,
    /// F1 nucleus from the same list.
    f1_nucleus: Option<String>,
    /// `##$SF=` in file order (TopSpin exports: one `procs` block per dimension, F2 first).
    sf: Vec<f64>,
    /// `##$AXNUC=` in file order.
    axnuc: Vec<String>,
    /// `##$OFFSET=` in file order.
    offset: Vec<f64>,
}

/// [`JcampHeader`] of a JCAMP-DX trace: `extra.dimensions` holds the direct dimension (F2), then
/// the page variable (F1).
fn jcamp_header(t: &TraceInfo) -> JcampHeader {
    let dims = t.extra.get("dimensions").and_then(Value::as_array);
    let of = |i: usize, key: &str| dims?.get(i)?.get(key);
    let dim = |i: usize| {
        Some((
            of(i, "first").and_then(num)?,
            of(i, "last").and_then(num)?,
            usize::try_from(of(i, "size")?.as_u64()?).ok()?,
            of(i, "unit")
                .and_then(Value::as_str)
                .unwrap_or_default()
                .to_string(),
        ))
    };
    let nucleus = |i: usize| {
        of(i, "nucleus")
            .and_then(Value::as_str)
            .and_then(clean_nucleus)
    };
    let parameters = |key: &str| {
        t.extra
            .get("dimension_parameters")
            .and_then(|p| p.get(key))
            .and_then(Value::as_array)
            .map(Vec::as_slice)
            .unwrap_or_default()
    };
    JcampHeader {
        f2: dim(0),
        f1: dim(1),
        f2_mhz: of(0, "observe_frequency_mhz").and_then(num),
        f1_mhz: of(1, "observe_frequency_mhz").and_then(num),
        f2_nucleus: nucleus(0),
        f1_nucleus: nucleus(1),
        sf: parameters("spectrometer_frequency_mhz")
            .iter()
            .filter_map(num)
            .collect(),
        axnuc: parameters("nucleus")
            .iter()
            .filter_map(|v| v.as_str().map(str::to_string))
            .collect(),
        offset: parameters("offset_ppm").iter().filter_map(num).collect(),
    }
}

/// Whether a trace is a stored 2D chemical-shift spectrum: frequency domain, more than one
/// sweep and an F1 shift axis. A Bruker spectrum that names its F2 nucleus has to name the F1
/// nucleus too. Headers only.
pub fn is_2d_spectrum(t: &TraceInfo) -> bool {
    if is_fid(t) || t.sweep_count < 2 || spectrum_axis(t).is_none() {
        return false;
    }
    let Some(sa) = t.extra.get("sweep_axis") else {
        return false;
    };
    if let Some(ax) = sa.get("axis") {
        let ppm = ax
            .get("unit")
            .and_then(Value::as_str)
            .is_some_and(|u| u.eq_ignore_ascii_case("ppm"));
        let named = |v: Option<&Value>| v.and_then(Value::as_str).and_then(clean_nucleus).is_some();
        return ppm && (named(sa.get("nucleus")) || !named(t.extra.get("nucleus")));
    }
    sa.get("unit")
        .and_then(Value::as_str)
        .is_some_and(|u| u.eq_ignore_ascii_case("hz") || u.eq_ignore_ascii_case("ppm"))
        && sa
            .get("values")
            .and_then(Value::as_array)
            .is_some_and(|v| v.len() as u64 == u64::from(t.sweep_count))
}

/// The nuclei of the acquisition, F2 then F1: `nucleus` and the first of
/// `indirect_dimensions` of the data set's time-domain trace (Bruker `NUC1` of `acqus` and of
/// `acqu2s`).
fn acquisition_nuclei(traces: &[TraceInfo]) -> (Option<String>, Option<String>) {
    let Some(fid) = traces.iter().find(|t| is_fid(t)) else {
        return (None, None);
    };
    let name = |v: Option<&Value>| v.and_then(Value::as_str).and_then(clean_nucleus);
    let f1 = fid
        .extra
        .get("indirect_dimensions")
        .and_then(Value::as_array)
        .and_then(|d| d.first())
        .and_then(|d| name(d.get("nucleus")));
    (name(fid.extra.get("nucleus")), f1)
}

/// A generous range of chemical shifts for a nucleus, ppm: a spectrum of that nucleus overlaps
/// it unless its referencing is wrong. `None` for nuclei without an entry.
fn usual_shift_range(nucleus: &str) -> Option<(f64, f64)> {
    match nucleus.trim().to_ascii_uppercase().as_str() {
        "1H" | "2H" => Some((-50.0, 50.0)),
        "13C" => Some((-100.0, 400.0)),
        "15N" => Some((-600.0, 1200.0)),
        "19F" => Some((-600.0, 900.0)),
        "31P" => Some((-600.0, 800.0)),
        _ => None,
    }
}

/// A note for a dimension whose whole axis lies outside [`usual_shift_range`] of its nucleus.
fn implausible_axis(name: &str, d: &SpectrumDimension) -> Option<String> {
    let nucleus = d.nucleus.as_deref()?;
    let (lo, hi) = usual_shift_range(nucleus)?;
    let a = d.axis.first_ppm;
    let b = d.axis.ppm(d.axis.size.saturating_sub(1) as f64);
    let (min, max) = (a.min(b), a.max(b));
    (min.is_finite() && max.is_finite() && (max < lo || min > hi)).then(|| {
        format!(
            "the {name} axis runs from {max} to {min} ppm, outside the usual {nucleus} range ({lo} to {hi} ppm): the stored referencing is probably wrong, so check the {name} shifts"
        )
    })
}

fn not_2d(t: &TraceInfo) -> Error {
    Error::unsupported(
        "nmr-processing",
        format!("trace {} has no chemical-shift F1 axis", t.index),
        "Cross peaks need a stored 2D spectrum (Bruker pdata/<n> with 2rr, or a JCAMP-DX nD NMR SPECTRUM); `openreadout info` lists the traces and their extra.sweep_axis.",
    )
}

/// The F2 and F1 dimensions of a stored 2D spectrum, and notes on how they were placed.
/// The rules for JCAMP-DX axes are in `docs/formats/nmr-processing.md` ("2-D cross peaks").
pub fn axes(t: &TraceInfo) -> Result<(SpectrumDimension, SpectrumDimension, Vec<String>)> {
    let sa = t.extra.get("sweep_axis").ok_or_else(|| not_2d(t))?;
    let rows = t.sweep_count as usize;
    let f2_nucleus = t
        .extra
        .get("nucleus")
        .and_then(Value::as_str)
        .and_then(clean_nucleus);
    if let Some(ax) = sa.get("axis") {
        let f2 = spectrum_axis(t).ok_or_else(|| not_2d(t))?;
        let first = ax.get("first").and_then(num).ok_or_else(|| not_2d(t))?;
        let step = ax.get("step").and_then(num).ok_or_else(|| not_2d(t))?;
        let mhz = ax
            .get("spectrometer_frequency_mhz")
            .and_then(num)
            .or_else(|| sa.get("spectrometer_frequency_mhz").and_then(num))
            .unwrap_or(0.0);
        let f1_nucleus = sa
            .get("nucleus")
            .and_then(Value::as_str)
            .and_then(clean_nucleus);
        return Ok((
            SpectrumDimension {
                nucleus: f2_nucleus,
                axis: f2,
            },
            SpectrumDimension {
                nucleus: f1_nucleus,
                axis: PpmAxis {
                    first_ppm: first,
                    step_ppm: step,
                    size: rows,
                    reference_frequency_mhz: mhz,
                },
            },
            Vec::new(),
        ));
    }
    let h = jcamp_header(t);
    let (Some(d2), Some(d1)) = (h.f2.clone(), h.f1.clone()) else {
        return Err(not_2d(t));
    };
    let mut notes = Vec::new();
    let observe = t.extra.get("observe_frequency_mhz").and_then(num);
    let f2_mhz = h.f2_mhz.or_else(|| h.sf.first().copied()).or(observe);
    let f1_mhz = h.f1_mhz.or_else(|| h.sf.get(1).copied());
    let to_axis = |(first, last, n, unit): &(f64, f64, usize, String),
                   mhz: Option<f64>,
                   what: &str|
     -> Result<PpmAxis> {
        let n = *n;
        if n < 2 {
            return Err(not_2d(t));
        }
        let scale = if unit.eq_ignore_ascii_case("ppm") {
            1.0
        } else if unit.eq_ignore_ascii_case("hz") {
            let f = mhz.filter(|f| *f > 0.0).ok_or_else(|| {
                Error::unsupported(
                    "nmr-processing",
                    format!("this JCAMP-DX 2D spectrum names no {what} observe frequency"),
                    format!("The {what} values are in Hz and cannot be converted to ppm. Export the spectrum again with a ##.OBSERVE FREQUENCY= per dimension, or read the Bruker pdata directory."),
                )
            })?;
            1.0 / f
        } else {
            return Err(not_2d(t));
        };
        let (a, b) = (first * scale, last * scale);
        Ok(PpmAxis {
            first_ppm: a,
            step_ppm: (b - a) / (n - 1) as f64,
            size: n,
            reference_frequency_mhz: mhz.unwrap_or(0.0),
        })
    };
    let mut f2 = to_axis(&d2, f2_mhz, "F2")?;
    let mut f1 = to_axis(&d1, f1_mhz, "F1")?;
    if f1.size != rows {
        return Err(not_2d(t));
    }
    if h.offset.len() >= 2 {
        for (name, ax, off) in [("F2", &mut f2, h.offset[0]), ("F1", &mut f1, h.offset[1])] {
            if (ax.first_ppm - off).abs() > 2.0 * ax.step_ppm.abs() {
                notes.push(format!(
                    "{name} starts at {:.4} ppm by the JCAMP-DX ##FIRST= but at {off} ppm by the file's ##$OFFSET=; the axis was moved to the offset (TopSpin nD exports give F2 from 0 Hz; inferred)",
                    ax.first_ppm
                ));
                ax.first_ppm = off;
            }
        }
    }
    // the page labels should agree with the header's axis
    let pages: Vec<f64> = sa
        .get("values")
        .and_then(Value::as_array)
        .map(|v| v.iter().filter_map(num).collect())
        .unwrap_or_default();
    if let (Some(&p0), Some(&pn)) = (pages.first(), pages.last()) {
        // in the header's own unit
        let step = ((d1.1 - d1.0) / (rows - 1) as f64).abs();
        if (p0 - d1.0).abs() > 2.0 * step || (pn - d1.1).abs() > 2.0 * step {
            notes.push(format!(
                "the F1 page labels run from {p0} to {pn}, the NTUPLES header from {} to {}; the header is used",
                d1.0, d1.1
            ));
        }
    }
    let f1_nucleus = h
        .f1_nucleus
        .or_else(|| h.axnuc.get(1).and_then(|s| clean_nucleus(s)));
    let f2_nucleus = h
        .f2_nucleus
        .or_else(|| h.axnuc.first().and_then(|s| clean_nucleus(s)))
        .or(f2_nucleus);
    if h.f1_mhz.is_none() {
        notes.push(format!(
            "F1 shifts are Hz divided by the F1 frequency {} MHz from the file's second ##$SF= (inferred)",
            f1.reference_frequency_mhz
        ));
    }
    Ok((
        SpectrumDimension {
            nucleus: f2_nucleus,
            axis: f2,
        },
        SpectrumDimension {
            nucleus: f1_nucleus,
            axis: f1,
        },
        notes,
    ))
}

/// Experiment kind and multiplicity editing from a pulse program name (`hsqcedetgpsisp2.3` →
/// `("hsqc", true)`), by the names seen in corpus files; `None` when the name says nothing.
pub fn experiment_kind(pulse_program: &str) -> Option<(&'static str, bool)> {
    let p = pulse_program.trim().to_ascii_lowercase();
    // HSQC-TOCSY names contain `dipsi`/`mlev` too, so the heteronuclear kinds come first.
    if p.contains("hsqc") {
        return Some(("hsqc", p.contains("hsqced")));
    }
    for (marker, kind) in [
        ("hmqc", "hmqc"),
        ("hmbc", "hmbc"),
        ("tocsy", "tocsy"),
        ("mlev", "tocsy"),
        ("dipsi", "tocsy"),
        ("noesy", "noesy"),
        ("roesy", "roesy"),
        ("cosy", "cosy"),
    ] {
        if p.contains(marker) {
            return Some((kind, false));
        }
    }
    None
}

/// The pulse program name: the trace's own (`pulse_sequence` of JCAMP-DX) or any trace's
/// `pulse_program` (the Bruker time-domain trace).
fn pulse_program(info: &FileInfo, t: &TraceInfo) -> Option<String> {
    let get = |t: &TraceInfo| {
        ["pulse_program", "pulse_sequence"].iter().find_map(|k| {
            t.extra
                .get(*k)
                .and_then(Value::as_str)
                .map(|s| s.trim().trim_matches(|c| c == '<' || c == '>').to_string())
                .filter(|s| !s.is_empty())
        })
    };
    get(t).or_else(|| info.traces.iter().find_map(get))
}

const KINDS: [&str; 10] = [
    "hsqc",
    "hmqc",
    "hmbc",
    "cosy",
    "tocsy",
    "noesy",
    "roesy",
    "heteronuclear",
    "homonuclear",
    "unknown",
];

fn experiment(
    info: &FileInfo,
    t: &TraceInfo,
    f2: &SpectrumDimension,
    f1: &SpectrumDimension,
    requested: Option<&str>,
) -> Result<ExperimentInfo> {
    let pp = pulse_program(info, t);
    let from_pp = pp.as_deref().and_then(experiment_kind);
    if let Some(r) = requested {
        let r = r.trim().to_ascii_lowercase();
        let Some(k) = KINDS.iter().copied().find(|k| *k == r.as_str()) else {
            return Err(Error::Usage(format!(
                "experiment {r:?}: one of {}",
                KINDS.join(", ")
            )));
        };
        return Ok(ExperimentInfo {
            kind: k.to_string(),
            multiplicity_edited: k == "hsqc" && from_pp.is_some_and(|(_, e)| e),
            pulse_program: pp,
            kind_source: "request".into(),
        });
    }
    if let Some((kind, edited)) = from_pp {
        return Ok(ExperimentInfo {
            kind: kind.into(),
            multiplicity_edited: edited,
            pulse_program: pp,
            kind_source: "pulse_program".into(),
        });
    }
    let kind = match (&f2.nucleus, &f1.nucleus) {
        (Some(a), Some(b)) if a.eq_ignore_ascii_case(b) => "homonuclear",
        (Some(_), Some(_)) => "heteronuclear",
        _ => "unknown",
    };
    Ok(ExperimentInfo {
        kind: kind.into(),
        multiplicity_edited: false,
        pulse_program: pp,
        kind_source: "nuclei".into(),
    })
}

/// The trace to use: `requested` (checked) or the first stored 2D spectrum.
fn choose_trace(info: &FileInfo, requested: Option<u32>) -> Result<&TraceInfo> {
    if let Some(i) = requested {
        let t = info.traces.iter().find(|t| t.index == i).ok_or_else(|| {
            Error::Usage(format!(
                "trace {i} out of range (file has {} traces)",
                info.traces.len()
            ))
        })?;
        if !is_2d_spectrum(t) {
            return Err(Error::unsupported(
                "nmr-processing",
                format!("trace {i} is not a stored 2D spectrum"),
                "Pick a processed 2D spectrum (Bruker pdata/<n> with 2rr, a JCAMP-DX nD NMR SPECTRUM; see `openreadout info`). For 1-D spectra use `openreadout analyze nmr-peaks`.",
            ));
        }
        return Ok(t);
    }
    info.traces.iter().find(|t| is_2d_spectrum(t)).ok_or_else(|| {
        let has_ser = info
            .traces
            .iter()
            .any(|t| is_fid(t) && t.sweep_count > 1);
        let hint = if has_ser {
            "This data set has a 2D FID (ser) but no processed 2D spectrum; FIDs of 2D experiments are not processed here. Process it in the vendor software (or export JCAMP-DX) and read the result."
        } else {
            "Cross peaks need a stored 2D spectrum (Bruker pdata/<n> with 2rr, or a JCAMP-DX nD NMR SPECTRUM). For 1-D spectra use `openreadout analyze nmr-peaks`."
        };
        Error::unsupported(
            "nmr-processing",
            format!("this {} data set has no stored 2D spectrum", info.format.name),
            hint,
        )
    })
}

/// A stored 2D spectrum: `rows` F1 points × `cols` F2 points, row-major.
#[derive(Debug, Clone)]
pub struct Matrix {
    /// Values, `rows × cols`, F1 slowest.
    pub data: Vec<f64>,
    /// F1 points.
    pub rows: usize,
    /// F2 points.
    pub cols: usize,
}

impl Matrix {
    fn at(&self, r: usize, c: usize) -> f64 {
        self.data[r * self.cols + c]
    }
}

/// Read the real part of every sweep of `t`.
fn read_matrix(ds: &mut dyn Dataset, t: &TraceInfo, cols: usize) -> Result<Matrix> {
    let rows = t.sweep_count as usize;
    let total = rows.checked_mul(cols).filter(|n| *n <= MAX_POINTS);
    let Some(total) = total else {
        return Err(Error::unsupported(
            "nmr-processing",
            format!("a 2D spectrum of {rows} × {cols} points"),
            format!(
                "Cross peaks are picked on spectra of at most {MAX_POINTS} points; reprocess with a smaller SI."
            ),
        ));
    };
    let mut data = Vec::with_capacity(total);
    for r in 0..rows {
        let tr = ds.read_trace(t.index, r as u32, 0, cols as u64)?;
        let Some(ch) = tr.channels.into_iter().next() else {
            return Err(Error::corrupt(
                "nmr-processing",
                format!("trace {} sweep {r} has no channels", t.index),
            ));
        };
        if ch.len() < cols {
            return Err(Error::corrupt(
                "nmr-processing",
                format!(
                    "trace {} sweep {r} has {} points, the axis {cols}",
                    t.index,
                    ch.len()
                ),
            ));
        }
        data.extend_from_slice(&ch[..cols]);
    }
    Ok(Matrix { data, rows, cols })
}

/// Robust noise SD and baseline of `m` (module docs). Zero SD when the spectrum is too small or
/// constant.
pub fn noise_and_baseline(m: &Matrix) -> (f64, f64) {
    let tile = |n: usize| (n / NOISE_TILES).max(8).min(n.max(1));
    let (tr, tc) = (tile(m.rows), tile(m.cols));
    let mut sds = Vec::new();
    let mut medians = Vec::new();
    let mut buf = Vec::with_capacity(tr * tc);
    let mut r0 = 0;
    while r0 + tr <= m.rows {
        let mut c0 = 0;
        while c0 + tc <= m.cols {
            buf.clear();
            for r in r0..r0 + tr {
                let row = &m.data[r * m.cols + c0..r * m.cols + c0 + tc];
                buf.extend(row.iter().copied().filter(|v| v.is_finite()));
            }
            if buf.len() >= 16 {
                let med = median_in_place(&mut buf);
                for v in &mut buf {
                    *v = (*v - med).abs();
                }
                sds.push(1.4826 * median_in_place(&mut buf));
                medians.push(med);
            }
            c0 += tc;
        }
        r0 += tr;
    }
    if sds.is_empty() {
        return (0.0, 0.0);
    }
    sds.sort_by(f64::total_cmp);
    let sd = sds[(sds.len() - 1) / 4];
    let baseline = median_in_place(&mut medians);
    (sd, baseline)
}

/// Noise SD of `buf` from the 80th percentile of |x − median| (1.28 SD for Gaussian noise), at
/// least `floor`. Unlike the MAD it rises with artifacts that cover a fifth of a column.
fn quantile_noise(buf: &mut [f64], floor: f64) -> f64 {
    if buf.len() < 8 {
        return floor;
    }
    let med = median_in_place(buf);
    for v in buf.iter_mut() {
        *v = (*v - med).abs();
    }
    let k = (buf.len() * 80 / 100).min(buf.len() - 1);
    buf.select_nth_unstable_by(k, f64::total_cmp);
    (buf[k] / COLUMN_QUANTILE_Z).max(floor)
}

/// [`quantile_noise`] of every F1 column of `m` inside the windows. t1 noise and sampling
/// artifacts run along F1, so a column through a tall peak is noisier than the spectrum as a
/// whole.
fn column_noise(m: &Matrix, rows: (usize, usize), cols: (usize, usize), floor: f64) -> Vec<f64> {
    let mut buf = Vec::new();
    let mut col_sd = vec![floor; m.cols];
    for (c, sd) in col_sd.iter_mut().enumerate().take(cols.1 + 1).skip(cols.0) {
        buf.clear();
        buf.extend(
            (rows.0..=rows.1)
                .map(|r| m.at(r, c))
                .filter(|v| v.is_finite()),
        );
        *sd = quantile_noise(&mut buf, floor);
    }
    col_sd
}

/// [`quantile_noise`] of column `c` without the rows within `skip` of row `r`, so that a peak's
/// own F1 profile does not count against it.
fn noise_around(
    m: &Matrix,
    r: usize,
    c: usize,
    rows: (usize, usize),
    skip: usize,
    floor: f64,
) -> f64 {
    let buf: Vec<f64> = (rows.0..=rows.1)
        .filter(|&k| k.abs_diff(r) > skip)
        .map(|k| m.at(k, c))
        .filter(|v| v.is_finite())
        .collect();
    let mut buf = buf;
    quantile_noise(&mut buf, floor)
}

/// A maximum found by [`find_maxima`]: point, signed height above the baseline.
#[derive(Debug, Clone, Copy)]
struct Maximum {
    r: usize,
    c: usize,
    s: f64,
}

fn find(parent: &mut [u32], mut i: u32) -> u32 {
    while parent[i as usize] != i {
        let p = parent[i as usize];
        parent[i as usize] = parent[p as usize];
        i = p;
    }
    i
}

/// The thresholds of [`find_maxima`].
struct Limits<'a> {
    /// Noise of the F1 column of a point.
    noise: &'a dyn Fn(usize, usize) -> f64,
    /// Noise of that column without the point's own rows: a maximum is kept at `min_snr` times it.
    accept: &'a dyn Fn(usize, usize) -> f64,
    /// Noise SD of the spectrum: points at `min_snr` times it are grouped.
    sd: f64,
    min_snr: f64,
    /// A maximum stays separate at this many times `noise` above its saddle.
    min_prominence_snr: f64,
    /// Lowest height grouped.
    floor: f64,
}

/// Maxima of `sign · (m − baseline)` inside the row and column windows (module docs).
fn find_maxima(
    m: &Matrix,
    sign: f64,
    baseline: f64,
    lim: &Limits<'_>,
    rows: (usize, usize),
    cols: (usize, usize),
) -> Vec<Maximum> {
    let s = |r: usize, c: usize| sign * (m.at(r, c) - baseline);
    let mut cand: Vec<(usize, usize)> = Vec::new();
    for r in rows.0..=rows.1 {
        for c in cols.0..=cols.1 {
            let v = s(r, c);
            if v.is_finite() && v >= lim.floor && v >= lim.min_snr * lim.sd {
                cand.push((r, c));
            }
        }
    }
    cand.sort_by(|a, b| s(b.0, b.1).total_cmp(&s(a.0, a.1)));
    // slot of each visited point (index into `parent`)
    let mut slot: BTreeMap<(usize, usize), u32> = BTreeMap::new();
    let mut parent: Vec<u32> = Vec::with_capacity(cand.len());
    // the maximum of each group, kept at its root
    let mut top: Vec<Maximum> = Vec::with_capacity(cand.len());
    let mut kept: Vec<Maximum> = Vec::new();
    for &(r, c) in &cand {
        let here = Maximum { r, c, s: s(r, c) };
        let id = parent.len() as u32;
        parent.push(id);
        top.push(here);
        slot.insert((r, c), id);
        let mut roots: Vec<u32> = Vec::new();
        for dr in [-1i64, 0, 1] {
            for dc in [-1i64, 0, 1] {
                if dr == 0 && dc == 0 {
                    continue;
                }
                let (nr, nc) = (r as i64 + dr, c as i64 + dc);
                if nr < 0 || nc < 0 {
                    continue;
                }
                let Some(&n) = slot.get(&(nr as usize, nc as usize)) else {
                    continue;
                };
                let root = find(&mut parent, n);
                if !roots.contains(&root) {
                    roots.push(root);
                }
            }
        }
        if roots.is_empty() {
            continue;
        }
        // the group with the highest maximum absorbs the others and this point
        roots.sort_by(|a, b| top[*b as usize].s.total_cmp(&top[*a as usize].s));
        let main = roots[0];
        for &g in &roots[1..] {
            let lower = top[g as usize];
            if lower.s - here.s >= lim.min_prominence_snr * (lim.noise)(lower.r, lower.c) {
                kept.push(lower);
            }
            parent[g as usize] = main;
        }
        parent[id as usize] = main;
    }
    for (i, p) in parent.iter().enumerate() {
        if *p == i as u32 {
            kept.push(top[i]);
        }
    }
    kept.retain(|mx| mx.s >= lim.min_snr * (lim.accept)(mx.r, mx.c));
    kept
}

/// Offset in (−0.5, 0.5) of the parabola's vertex through (a, b, c) at (−1, 0, 1), and its
/// height; (0, b) when the three points do not curve downwards.
fn parabola(a: f64, b: f64, c: f64) -> (f64, f64) {
    let den = a - 2.0 * b + c;
    if den < 0.0 && a.is_finite() && c.is_finite() {
        let off = (0.5 * (a - c) / den).clamp(-0.5, 0.5);
        (off, b - 0.25 * (a - c) * off)
    } else {
        (0.0, b)
    }
}

/// Index window of a ppm range on `axis` (rounded outwards, clamped to the axis).
fn window(axis: &PpmAxis, range: Option<(f64, f64)>) -> (usize, usize) {
    let last = axis.size.saturating_sub(1);
    match range {
        Some((a, b)) => {
            let (i, j) = (axis.index_of(a), axis.index_of(b));
            let (i, j) = if i <= j { (i, j) } else { (j, i) };
            (
                i.floor().clamp(0.0, last as f64) as usize,
                j.ceil().clamp(0.0, last as f64) as usize,
            )
        }
        None => (0, last),
    }
}

fn inside(range: Option<(f64, f64)>, ppm: f64) -> bool {
    range.is_none_or(|(a, b)| (a.min(b)..=a.max(b)).contains(&ppm))
}

/// Cross peaks of the spectrum `m` on axes `f2` (columns) and `f1` (rows), with the noise SD and
/// baseline (module docs). `homonuclear` spectra get `on_diagonal` and `mirror_partner`, and no
/// sideband rule.
pub fn pick_cross_peaks(
    m: &Matrix,
    f2: &PpmAxis,
    f1: &PpmAxis,
    opts: &CrossPeakOptions,
    homonuclear: bool,
) -> Picked {
    let (sd, baseline) = noise_and_baseline(m);
    let empty = Picked {
        noise_sd: sd,
        baseline,
        ..Picked::default()
    };
    if m.rows < 2 || m.cols < 2 {
        return empty;
    }
    let rows = window(f1, opts.f1_range_ppm);
    let cols = window(f2, opts.f2_range_ppm);
    let mut tallest = 0.0_f64;
    for r in rows.0..=rows.1 {
        for c in cols.0..=cols.1 {
            let v = (m.at(r, c) - baseline).abs();
            if v.is_finite() {
                tallest = tallest.max(v);
            }
        }
    }
    if !(sd > 0.0) {
        return empty;
    }
    let col_sd = column_noise(m, rows, cols, sd);
    let noise = |_r: usize, c: usize| col_sd[c];
    // a peak's own rows: a sixteenth of F1 on each side, at least eight points
    let skip = (m.rows / 16).max(8);
    let accept = |r: usize, c: usize| noise_around(m, r, c, rows, skip, sd);
    let lim = Limits {
        noise: &noise,
        accept: &accept,
        sd,
        min_snr: opts.min_snr,
        min_prominence_snr: opts.min_prominence_snr,
        floor: opts.min_height_fraction * tallest,
    };
    let mut peaks: Vec<CrossPeak> = Vec::new();
    for sign in [1.0, -1.0] {
        if sign < 0.0 && opts.positive_only {
            break;
        }
        for mx in find_maxima(m, sign, baseline, &lim, rows, cols) {
            let s = |r: usize, c: usize| sign * (m.at(r, c) - baseline);
            let (r, c) = (mx.r, mx.c);
            // a maximum on the outermost row or column cannot be told from a rising edge
            if r == 0 || c == 0 || r + 1 == m.rows || c + 1 == m.cols {
                continue;
            }
            let row_n = |dc: i64| {
                let cc = c as i64 + dc;
                if cc < 0 || cc as usize >= m.cols {
                    f64::NAN
                } else {
                    s(r, cc as usize)
                }
            };
            let col_n = |dr: i64| {
                let rr = r as i64 + dr;
                if rr < 0 || rr as usize >= m.rows {
                    f64::NAN
                } else {
                    s(rr as usize, c)
                }
            };
            let (off2, h2) = parabola(row_n(-1), mx.s, row_n(1));
            let (off1, h1) = parabola(col_n(-1), mx.s, col_n(1));
            // both parabolas pass through the same maximum; add both rises
            let h = h2 + h1 - mx.s;
            let f2_ppm = f2.ppm(c as f64 + off2);
            let f1_ppm = f1.ppm(r as f64 + off1);
            if !(inside(opts.f2_range_ppm, f2_ppm) && inside(opts.f1_range_ppm, f1_ppm)) {
                continue;
            }
            let w2 = half_width_of(m.cols, c, h, |k| s(r, k)).map(|w| w * f2.hz_per_point());
            let w1 = half_width_of(m.rows, r, h, |k| s(k, c)).map(|w| w * f1.hz_per_point());
            peaks.push(CrossPeak {
                f2_ppm,
                f1_ppm,
                height: baseline + sign * h,
                snr: h / accept(r, c),
                // a parabola height above twice the sample would give a width below zero
                f2_width_hz: w2.filter(|w| w.is_finite() && *w > 0.0),
                f1_width_hz: w1.filter(|w| w.is_finite() && *w > 0.0),
                f2_index: c,
                f1_index: r,
                maxima: 1,
                on_diagonal: None,
                mirror_partner: None,
            });
        }
    }
    let cross = merge(peaks, f2, f1, opts.merge_hz);
    let (mut cross, ripples) =
        drop_ripples(cross, baseline, f2, f1, opts.ripple_fraction, !homonuclear);
    if homonuclear {
        let tol = |a: &PpmAxis| {
            let hz = if a.reference_frequency_mhz > 0.0 {
                opts.merge_hz / a.reference_frequency_mhz
            } else {
                0.0
            };
            (2.0 * a.step_ppm.abs()).max(hz)
        };
        let t = tol(f2).max(tol(f1));
        for p in &mut cross {
            p.on_diagonal = Some((p.f1_ppm - p.f2_ppm).abs() <= t);
        }
        let spots: Vec<(f64, f64, bool)> = cross
            .iter()
            .map(|p| (p.f2_ppm, p.f1_ppm, p.on_diagonal == Some(true)))
            .collect();
        for (i, p) in cross.iter_mut().enumerate() {
            if p.on_diagonal == Some(true) {
                continue;
            }
            p.mirror_partner = Some(spots.iter().enumerate().any(|(j, &(a, b, diag))| {
                j != i && !diag && (a - p.f1_ppm).abs() <= t && (b - p.f2_ppm).abs() <= t
            }));
        }
    }
    let count = cross.len();
    let at_snr_levels =
        SNR_LEVELS.map(|k| cross.iter().filter(|p| p.snr >= k * opts.min_snr).count());
    if cross.len() > opts.max_cross_peaks {
        cross.sort_by(|a, b| b.snr.total_cmp(&a.snr));
        cross.truncate(opts.max_cross_peaks);
    }
    cross.sort_by(|a, b| {
        b.f2_ppm
            .total_cmp(&a.f2_ppm)
            .then(b.f1_ppm.total_cmp(&a.f1_ppm))
    });
    Picked {
        cross_peaks: cross,
        count,
        ripples,
        at_snr_levels,
        noise_sd: sd,
        baseline,
    }
}

/// Leave out truncation ripples and decoupling sidebands of taller cross peaks (module docs).
/// Returns the kept cross peaks and the number left out.
fn drop_ripples(
    cross: Vec<CrossPeak>,
    baseline: f64,
    f2: &PpmAxis,
    f1: &PpmAxis,
    fraction: f64,
    heteronuclear: bool,
) -> (Vec<CrossPeak>, usize) {
    if !(fraction > 0.0) {
        return (cross, 0);
    }
    let amp = |p: &CrossPeak| (p.height - baseline).abs();
    let (hz2, hz1) = (f2.hz_per_point(), f1.hz_per_point());
    let (mhz2, mhz1) = (f2.reference_frequency_mhz, f1.reference_frequency_mhz);
    let mut order: Vec<CrossPeak> = cross;
    order.sort_by(|a, b| amp(b).total_cmp(&amp(a)));
    let mut kept: Vec<CrossPeak> = Vec::with_capacity(order.len());
    let mut dropped = 0;
    for (i, p) in order.iter().enumerate() {
        let ripple = kept.iter().any(|q| {
            if amp(p) > fraction * amp(q) {
                return false;
            }
            let w2 = q.f2_width_hz.unwrap_or(UNMEASURED_WIDTH_POINTS * hz2);
            let w1 = q.f1_width_hz.unwrap_or(UNMEASURED_WIDTH_POINTS * hz1);
            let in_row = |x: &CrossPeak| (x.f1_ppm - q.f1_ppm).abs() * mhz1 <= (0.5 * w1).max(hz1);
            let d2 = (p.f2_ppm - q.f2_ppm) * mhz2;
            let d1 = (p.f1_ppm - q.f1_ppm).abs() * mhz1;
            let column = d2.abs() <= w2.max(2.0 * hz2) && d1 <= RIPPLE_WIDTHS * w1.max(hz1);
            let sideband = heteronuclear
                && in_row(p)
                && d2.abs() > w2.max(2.0 * hz2)
                && d2.abs() <= SIDEBAND_HZ
                && order.iter().enumerate().any(|(j, x)| {
                    j != i
                        && amp(x) <= fraction * amp(q)
                        && in_row(x)
                        && ((x.f2_ppm - q.f2_ppm) * mhz2 + d2).abs() <= w2.max(3.0 * hz2)
                });
            column || sideband
        });
        if ripple {
            dropped += 1;
        } else {
            kept.push(p.clone());
        }
    }
    (kept, dropped)
}

/// Merge peaks closer than `merge_hz` in F2 and in F1 (at least one F1 point), of either sign
/// (COSY multiplets are antiphase), into cross peaks placed at their tallest member (single
/// linkage).
fn merge(peaks: Vec<CrossPeak>, f2: &PpmAxis, f1: &PpmAxis, merge_hz: f64) -> Vec<CrossPeak> {
    let n = peaks.len();
    let gap = |a: &PpmAxis, at_least_points: f64| {
        let hz = if a.reference_frequency_mhz > 0.0 && merge_hz > 0.0 {
            merge_hz / a.reference_frequency_mhz
        } else {
            0.0
        };
        hz.max(at_least_points * a.step_ppm.abs())
    };
    let g2 = gap(f2, 0.0);
    // adjacent F1 points always merge: a cross peak is wider than one F1 point
    let g1 = gap(f1, 1.0001);
    let mut parent: Vec<u32> = (0..n as u32).collect();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|a, b| peaks[*a].f2_ppm.total_cmp(&peaks[*b].f2_ppm));
    for (k, &i) in order.iter().enumerate() {
        for &j in &order[k + 1..] {
            let (a, b) = (&peaks[i], &peaks[j]);
            if b.f2_ppm - a.f2_ppm > g2 {
                break;
            }
            if (b.f1_ppm - a.f1_ppm).abs() <= g1 {
                let (ri, rj) = (find(&mut parent, i as u32), find(&mut parent, j as u32));
                if ri != rj {
                    parent[rj as usize] = ri;
                }
            }
        }
    }
    let mut groups: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for i in 0..n {
        let r = find(&mut parent, i as u32);
        groups.entry(r).or_default().push(i);
    }
    groups
        .into_values()
        .filter_map(|members| {
            let best = *members
                .iter()
                .max_by(|a, b| peaks[**a].snr.total_cmp(&peaks[**b].snr))?;
            let mut p = peaks[best].clone();
            p.maxima = members.len() as u32;
            Some(p)
        })
        .collect()
}

/// Pick the cross peaks of the stored 2D spectrum `req` asks for (module docs).
pub fn analyze_cross_peaks(
    ds: &mut dyn Dataset,
    info: &FileInfo,
    req: &CrossPeakRequest,
) -> Result<CrossPeakReport> {
    let opts = &req.options;
    let t = choose_trace(info, req.trace)?.clone();
    let (mut f2, mut f1, mut axis_notes) = axes(&t)?;
    let cols = (t.sample_count as usize).min(f2.axis.size);
    f2.axis.size = cols;
    if f2.nucleus.is_none() || f1.nucleus.is_none() {
        let (a2, a1) = acquisition_nuclei(&info.traces);
        let mut filled = Vec::new();
        for (name, dim, from) in [("F2", &mut f2, a2), ("F1", &mut f1, a1)] {
            if dim.nucleus.is_none()
                && let Some(n) = from
            {
                filled.push(format!("{name} {n}"));
                dim.nucleus = Some(n);
            }
        }
        if !filled.is_empty() {
            axis_notes.push(format!(
                "the processed spectrum names no nucleus for {}: taken from the acquisition parameters (inferred; TopSpin 1.3 writes no AXNUC)",
                filled.join(" and ")
            ));
        }
    }
    let exp = experiment(info, &t, &f2, &f1, req.experiment.as_deref())?;
    let homonuclear = match (&f2.nucleus, &f1.nucleus) {
        (Some(a), Some(b)) => a.eq_ignore_ascii_case(b),
        _ => matches!(
            exp.kind.as_str(),
            "cosy" | "tocsy" | "noesy" | "roesy" | "homonuclear"
        ),
    };
    let m = read_matrix(ds, &t, cols)?;
    let Picked {
        cross_peaks,
        count,
        ripples,
        at_snr_levels,
        noise_sd: sd,
        baseline,
    } = pick_cross_peaks(&m, &f2.axis, &f1.axis, opts, homonuclear);
    let mut notes: Vec<String> = axis_notes;
    notes.extend(implausible_axis("F2", &f2));
    notes.extend(implausible_axis("F1", &f1));
    if let Some(pct) = info.traces.iter().find_map(|t| {
        t.extra
            .get("non_uniform_sampling")
            .and_then(|n| n.get("amount_percent"))
            .and_then(num)
    }) {
        notes.push(format!(
            "non-uniform sampling ({pct:.0} % of the F1 increments): reconstruction artifacts can pass min_snr; raise it if weak cross peaks look doubtful"
        ));
    }
    if exp.multiplicity_edited {
        notes.push("multiplicity-edited HSQC: one multiplicity group has negative heights; which one (usually CH2) depends on how the spectrum was phased".into());
    }
    if ripples > 0 {
        notes.push(format!(
            "{ripples} peaks left out as truncation ripples or decoupling sidebands of a taller cross peak (ripple_fraction {}; 0 keeps them)",
            opts.ripple_fraction
        ));
    }
    if count > LONG_TABLE_ROWS {
        let [a, b, c] = at_snr_levels;
        let [ka, kb, kc] = SNR_LEVELS.map(|k| k * opts.min_snr);
        notes.push(format!(
            "{count} cross peaks pass min_snr {}: {a} have S/N of at least {ka}, {b} of {kb}, {c} of {kc}; raise min_snr or set min_height_fraction for a shorter table",
            opts.min_snr
        ));
    }
    if sd == 0.0 {
        notes.push("the noise SD is 0 (a constant or too small spectrum): no cross peaks".into());
    }
    if count > cross_peaks.len() {
        notes.push(format!(
            "{count} cross peaks found; the {} with the highest S/N are listed (max_cross_peaks)",
            cross_peaks.len()
        ));
    }
    Ok(CrossPeakReport {
        path: info.path.clone(),
        format: info.format.id.clone(),
        trace: t.index,
        trace_name: t.name.clone(),
        experiment: exp,
        f2,
        f1,
        baseline,
        noise_sd: sd,
        options: opts.clone(),
        cross_peak_count: count,
        cross_peaks,
        notes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A spectrum of `rows × cols` with Lorentzian-like cross peaks `(r, c, height, width)` on
    /// deterministic pseudo-random noise of SD about `noise`.
    fn synth(rows: usize, cols: usize, peaks: &[(f64, f64, f64, f64)], noise: f64) -> Matrix {
        let mut data = Vec::with_capacity(rows * cols);
        for r in 0..rows {
            for c in 0..cols {
                // sum of 12 uniforms in [-0.5, 0.5): approximately normal, SD 1
                let k = (r * cols + c) * 12;
                let g: f64 = (0..12)
                    .map(|j| crate::nmr::baseline::tests::noise(k + j))
                    .sum();
                let mut v = noise * g;
                for &(pr, pc, h, w) in peaks {
                    let (dr, dc) = ((r as f64 - pr) / w, (c as f64 - pc) / w);
                    v += h / ((1.0 + dr * dr) * (1.0 + dc * dc));
                }
                data.push(v);
            }
        }
        Matrix { data, rows, cols }
    }

    fn axis(first: f64, step: f64, size: usize, mhz: f64) -> PpmAxis {
        PpmAxis {
            first_ppm: first,
            step_ppm: step,
            size,
            reference_frequency_mhz: mhz,
        }
    }

    #[test]
    fn finds_positive_and_negative_cross_peaks() {
        // an "edited HSQC": 1H 10..0 ppm over 256 points (600 MHz), 13C 160..0 over 128
        let m = synth(
            128,
            256,
            &[
                (40.0, 60.0, 100.0, 1.5),
                (90.0, 200.0, -80.0, 1.5),
                (20.0, 30.0, 50.0, 1.5),
            ],
            1.0,
        );
        let f2 = axis(10.0, -10.0 / 255.0, 256, 600.0);
        let f1 = axis(160.0, -160.0 / 127.0, 128, 150.9);
        let r = pick_cross_peaks(&m, &f2, &f1, &CrossPeakOptions::default(), false);
        let (p, count, sd, base) = (r.cross_peaks, r.count, r.noise_sd, r.baseline);
        assert!((sd - 1.0).abs() < 0.2, "noise sd {sd}");
        assert!(base.abs() < 0.2, "baseline {base}");
        assert_eq!(count, 3, "{p:?}");
        let near = |r: f64, c: f64| {
            p.iter()
                .find(|x| {
                    (x.f1_ppm - f1.ppm(r)).abs() < 0.5 * f1.step_ppm.abs()
                        && (x.f2_ppm - f2.ppm(c)).abs() < 0.5 * f2.step_ppm.abs()
                })
                .unwrap_or_else(|| panic!("no cross peak at {r},{c}: {p:?}"))
        };
        assert!(near(40.0, 60.0).height > 90.0);
        assert!(near(90.0, 200.0).height < -70.0);
        // S/N is against the noise of the peak's own F1 column, which its tails raise
        assert!(near(20.0, 30.0).snr > 20.0);
        // high F2 shift first
        assert!(p.windows(2).all(|w| w[0].f2_ppm >= w[1].f2_ppm));
        let only_pos = CrossPeakOptions {
            positive_only: true,
            ..CrossPeakOptions::default()
        };
        let p = pick_cross_peaks(&m, &f2, &f1, &only_pos, false).cross_peaks;
        assert!(p.iter().all(|x| x.height > 0.0) && p.len() == 2);
    }

    #[test]
    fn multiplet_components_are_one_cross_peak() {
        // a doublet 7 Hz apart in F2 (600 MHz, 0.5 Hz per point → 14 points apart) is one cross
        // peak with two maxima; a peak 60 Hz away is another
        let f2 = axis(5.0, -0.5 / 600.0, 400, 600.0);
        let f1 = axis(80.0, -0.5, 64, 150.9);
        let m = synth(
            64,
            400,
            &[
                (30.0, 100.0, 100.0, 2.0),
                (30.0, 114.0, 100.0, 2.0),
                (30.0, 234.0, 100.0, 2.0),
            ],
            1.0,
        );
        let r = pick_cross_peaks(&m, &f2, &f1, &CrossPeakOptions::default(), false);
        let (p, count) = (r.cross_peaks, r.count);
        assert_eq!(count, 2, "{p:?}");
        assert_eq!(p.iter().map(|x| x.maxima).max(), Some(2));
        // merging switched off: three
        let none = CrossPeakOptions {
            merge_hz: 0.0,
            ..CrossPeakOptions::default()
        };
        let count = pick_cross_peaks(&m, &f2, &f1, &none, false).count;
        assert_eq!(count, 3);
    }

    #[test]
    fn shallow_dips_do_not_split_a_peak() {
        // two lines 3 points apart merge at a saddle above half height: one maximum (low noise,
        // so that noise on the long Lorentzian tails stays below the threshold)
        let f2 = axis(5.0, -0.01, 200, 600.0);
        let f1 = axis(80.0, -0.5, 64, 150.9);
        let m = synth(
            64,
            200,
            &[(30.0, 100.0, 100.0, 3.0), (30.0, 103.0, 95.0, 3.0)],
            0.05,
        );
        let opts = CrossPeakOptions {
            merge_hz: 0.0,
            ..CrossPeakOptions::default()
        };
        let r = pick_cross_peaks(&m, &f2, &f1, &opts, false);
        let (p, count) = (r.cross_peaks, r.count);
        assert_eq!(count, 1, "{p:?}");
    }

    #[test]
    fn diagonal_and_ranges() {
        // a COSY: 1H on both axes, 10..0 ppm, 200 × 200; a diagonal peak and a cross peak
        let a = axis(10.0, -10.0 / 199.0, 200, 400.0);
        let m = synth(
            200,
            200,
            &[(50.0, 50.0, 100.0, 1.5), (50.0, 120.0, 60.0, 1.5)],
            1.0,
        );
        let p = pick_cross_peaks(&m, &a, &a, &CrossPeakOptions::default(), true).cross_peaks;
        assert_eq!(p.len(), 2, "{p:?}");
        assert_eq!(p.iter().filter(|x| x.on_diagonal == Some(true)).count(), 1);
        let ranged = CrossPeakOptions {
            f2_range_ppm: Some((5.0, 0.0)),
            ..CrossPeakOptions::default()
        };
        let p = pick_cross_peaks(&m, &a, &a, &ranged, true).cross_peaks;
        assert_eq!(p.len(), 1);
        assert!(p[0].f2_ppm <= 5.0 && p[0].on_diagonal == Some(false));
    }

    #[test]
    fn noise_only_and_constant_spectra() {
        let f = axis(10.0, -0.1, 100, 400.0);
        let m = synth(100, 100, &[], 1.0);
        let p = pick_cross_peaks(&m, &f, &f, &CrossPeakOptions::default(), true).cross_peaks;
        assert!(p.is_empty(), "{p:?}");
        let flat = Matrix {
            data: vec![3.0; 100 * 100],
            rows: 100,
            cols: 100,
        };
        let r = pick_cross_peaks(&flat, &f, &f, &CrossPeakOptions::default(), true);
        let (p, sd, base) = (r.cross_peaks, r.noise_sd, r.baseline);
        assert!(p.is_empty() && sd == 0.0 && (base - 3.0).abs() < 1e-12);
        let tiny = Matrix {
            data: vec![1.0],
            rows: 1,
            cols: 1,
        };
        let p = pick_cross_peaks(&tiny, &f, &f, &CrossPeakOptions::default(), true).cross_peaks;
        assert!(p.is_empty());
    }

    #[test]
    fn ripples_and_sideband_pairs_are_left_out() {
        // an HSQC line with two 10 % F1 ripples in its column, a ±60 Hz pair of 5 % sidebands in
        // its row, and a real 15 % cross peak 40 F1 widths away in the same column
        let f2 = axis(8.0, -0.01, 300, 600.0);
        let f1 = axis(160.0, -0.5, 200, 150.9);
        let line = |r: f64, c: f64, h: f64| (r, c, h, 1.5);
        let m = synth(
            200,
            300,
            &[
                line(100.0, 150.0, 1000.0),
                line(96.0, 150.0, -100.0),
                line(104.0, 150.0, -100.0),
                line(100.0, 140.0, 50.0),
                line(100.0, 160.0, 50.0),
                line(20.0, 150.0, 150.0),
            ],
            1.0,
        );
        let r = pick_cross_peaks(&m, &f2, &f1, &CrossPeakOptions::default(), false);
        // the line and the real cross peak stay; the ripples and the sideband pair go
        assert_eq!(r.cross_peaks.len(), 2, "{:?}", r.cross_peaks);
        assert!(r.ripples >= 2, "{}", r.ripples);
        assert!(r.cross_peaks.iter().all(|x| x.height > 100.0));
        let keep_all = CrossPeakOptions {
            ripple_fraction: 0.0,
            ..CrossPeakOptions::default()
        };
        let r = pick_cross_peaks(&m, &f2, &f1, &keep_all, false);
        assert_eq!(r.ripples, 0);
        assert!(r.cross_peaks.len() >= 4, "{:?}", r.cross_peaks);
        // a single weak peak in the row has no mirror twin: kept
        let m = synth(
            200,
            300,
            &[line(100.0, 150.0, 1000.0), line(100.0, 140.0, 50.0)],
            1.0,
        );
        let r = pick_cross_peaks(&m, &f2, &f1, &CrossPeakOptions::default(), false);
        assert_eq!((r.ripples, r.cross_peaks.len()), (0, 2));
    }

    /// Like `synth`, for lines `(r, c, height, width, f1_phase_deg)` whose F1 profile is out of
    /// phase: absorption · cos + dispersion · sin.
    #[test]
    fn counts_at_higher_snr() {
        let m = synth(
            128,
            256,
            &[
                (40.0, 60.0, 400.0, 1.5),
                (90.0, 200.0, -80.0, 1.5),
                (20.0, 30.0, 30.0, 1.5),
            ],
            1.0,
        );
        let f2 = axis(10.0, -10.0 / 255.0, 256, 600.0);
        let f1 = axis(160.0, -160.0 / 127.0, 128, 150.9);
        let r = pick_cross_peaks(&m, &f2, &f1, &CrossPeakOptions::default(), false);
        assert_eq!(r.count, 3);
        // S/N of at least 20, 50 and 100 (min_snr 10 times SNR_LEVELS)
        assert_eq!(r.at_snr_levels, [3, 2, 1], "{:?}", r.cross_peaks);
    }

    fn trace(kind: &str, sweeps: u32, extra: serde_json::Value) -> TraceInfo {
        let mut e: BTreeMap<String, Value> = serde_json::from_value(extra).unwrap();
        e.insert("kind".into(), Value::from(kind));
        TraceInfo {
            index: 0,
            name: None,
            sample_rate_hz: 0.0,
            sample_count: 1024,
            sweep_count: sweeps,
            channels: Vec::new(),
            start_s: None,
            extra: e,
        }
    }

    #[test]
    fn spectra_without_a_stored_nucleus() {
        let ppm = |first: f64, mhz: f64| {
            serde_json::json!({"quantity": "chemical_shift", "unit": "ppm", "first": first,
                "step": -0.01, "size": 1024, "spectrometer_frequency_mhz": mhz})
        };
        // TopSpin 1.3: neither procs nor proc2s holds AXNUC
        let bare = trace(
            "processed_spectrum",
            1024,
            serde_json::json!({"axis": ppm(8.0, 600.13), "sweep_axis": {"axis": ppm(8.0, 600.13)}}),
        );
        assert!(is_2d_spectrum(&bare));
        // F2 named, F1 not: not taken as a 2D shift spectrum
        let half = trace(
            "processed_spectrum",
            1024,
            serde_json::json!({"axis": ppm(8.0, 600.13), "nucleus": "<1H>",
                "sweep_axis": {"axis": ppm(8.0, 600.13)}}),
        );
        assert!(!is_2d_spectrum(&half));
        let full = trace(
            "processed_spectrum",
            1024,
            serde_json::json!({"axis": ppm(8.0, 600.13), "nucleus": "<1H>",
                "sweep_axis": {"axis": ppm(160.0, 150.9), "nucleus": "<13C>"}}),
        );
        assert!(is_2d_spectrum(&full));
        // the acquisition names the nuclei
        let fid = trace(
            "time_domain",
            256,
            serde_json::json!({"nucleus": "<1H>",
                "indirect_dimensions": [{"dimension": 2, "nucleus": "<13C>"}]}),
        );
        assert_eq!(
            acquisition_nuclei(&[fid.clone(), bare.clone()]),
            (Some("1H".into()), Some("13C".into()))
        );
        assert_eq!(acquisition_nuclei(&[bare]), (None, None));
        let off = trace(
            "time_domain",
            256,
            serde_json::json!({"nucleus": "<1H>", "indirect_dimensions": [{"nucleus": "<off>"}]}),
        );
        assert_eq!(acquisition_nuclei(&[off]), (Some("1H".into()), None));
    }

    #[test]
    fn axes_outside_any_spectrum_get_a_note() {
        let dim = |nucleus: &str, first: f64| SpectrumDimension {
            nucleus: Some(nucleus.into()),
            axis: axis(first, -0.1, 2048, 125.77),
        };
        assert!(implausible_axis("F1", &dim("13C", 210.0)).is_none());
        assert!(implausible_axis("F2", &dim("1H", 14.0)).is_none());
        // partly outside is still a spectrum
        assert!(implausible_axis("F1", &dim("13C", 450.0)).is_none());
        let note = implausible_axis("F1", &dim("13C", -199_859.7)).unwrap();
        assert!(note.contains("F1 axis") && note.contains("13C"), "{note}");
        assert!(implausible_axis("F2", &dim("1H", -200_016.8)).is_some());
        // nuclei without an entry and unnamed dimensions get none
        assert!(implausible_axis("F1", &dim("59Co", 12000.0)).is_none());
        let unnamed = SpectrumDimension {
            nucleus: None,
            axis: axis(-200_000.0, -0.1, 2048, 125.77),
        };
        assert!(implausible_axis("F1", &unnamed).is_none());
    }

    #[test]
    fn homonuclear_cross_peaks_have_mirror_partners() {
        let a = axis(10.0, -10.0 / 199.0, 200, 400.0);
        let m = synth(
            200,
            200,
            &[
                (50.0, 50.0, 200.0, 1.5),
                (120.0, 120.0, 200.0, 1.5),
                (50.0, 120.0, 60.0, 1.5),
                (120.0, 50.0, 60.0, 1.5),
                (80.0, 150.0, 60.0, 1.5),
            ],
            1.0,
        );
        let p = pick_cross_peaks(&m, &a, &a, &CrossPeakOptions::default(), true).cross_peaks;
        let off: Vec<&CrossPeak> = p.iter().filter(|x| x.on_diagonal == Some(false)).collect();
        assert_eq!(off.len(), 3, "{p:?}");
        assert_eq!(
            off.iter()
                .filter(|x| x.mirror_partner == Some(true))
                .count(),
            2
        );
        assert!(
            p.iter()
                .filter(|x| x.on_diagonal == Some(true))
                .all(|x| x.mirror_partner.is_none())
        );
    }

    #[test]
    fn kinds_from_pulse_programs() {
        assert_eq!(experiment_kind("hsqcedetgpsisp2.3"), Some(("hsqc", true)));
        assert_eq!(experiment_kind("hsqcetgpsi2.ph"), Some(("hsqc", false)));
        assert_eq!(experiment_kind("hmqcgpqf"), Some(("hmqc", false)));
        assert_eq!(experiment_kind("hmbcgpl2ndqf"), Some(("hmbc", false)));
        assert_eq!(experiment_kind("cosygpppqf"), Some(("cosy", false)));
        assert_eq!(experiment_kind("mlevphpp"), Some(("tocsy", false)));
        assert_eq!(experiment_kind("noesygpphzs"), Some(("noesy", false)));
        assert_eq!(experiment_kind("roesyphpr.2"), Some(("roesy", false)));
        assert_eq!(experiment_kind("gCOSY"), Some(("cosy", false)));
        assert_eq!(experiment_kind("zg30"), None);
    }

    #[test]
    fn jcamp_headers_from_the_reader() {
        // MestReNova: observe frequencies and nuclei per column
        let mest = trace(
            "ntuples",
            512,
            serde_json::json!({"dimensions": [
                {"symbol": "F2", "unit": "HZ", "first": 4241.14, "last": -230.39, "size": 512,
                 "observe_frequency_mhz": 400.20, "nucleus": "1H"},
                {"symbol": "F1", "unit": "HZ", "first": 23146.75, "last": -960.65, "size": 512,
                 "observe_frequency_mhz": 100.64, "nucleus": "13C"},
            ]}),
        );
        let h = jcamp_header(&mest);
        assert_eq!(h.f1, Some((23146.75, -960.65, 512, "HZ".into())));
        assert_eq!(h.f2.as_ref().map(|d| d.2), Some(512));
        assert!(h.f1_mhz.is_some_and(|f| (f - 100.64).abs() < 1e-9));
        assert_eq!(h.f1_nucleus.as_deref(), Some("13C"));
        assert!(h.offset.is_empty());
        // TopSpin: its processing parameters once per dimension
        let topspin = trace(
            "ntuples",
            1024,
            serde_json::json!({
                "dimensions": [
                    {"symbol": "F2", "unit": "HZ", "first": 6403.99, "last": 0.0, "size": 1024},
                    {"symbol": "F1", "unit": "HZ", "first": 19061.27, "last": -4953.71, "size": 1024},
                ],
                "dimension_parameters": {
                    "spectrometer_frequency_mhz": [400.08, 100.600_196],
                    "nucleus": ["1H", "13C"],
                    "offset_ppm": [12.714_22, 189.4755],
                },
            }),
        );
        let h = jcamp_header(&topspin);
        assert_eq!(h.sf, vec![400.08, 100.600_196]);
        assert_eq!(h.offset, vec![12.714_22, 189.4755]);
        assert_eq!(h.axnuc, vec!["1H", "13C"]);
        assert!(h.f1_mhz.is_none() && h.f1_nucleus.is_none());
        // a trace without them (a Bruker spectrum, a 1-D JCAMP-DX block)
        let none = trace("ntuples", 4, serde_json::json!({}));
        assert_eq!(jcamp_header(&none), JcampHeader::default());
    }
}
