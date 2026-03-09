mod calculations;
mod constants;
mod file_parser;
mod hybridization;
mod types;
mod vienna_wrapper;

use calculations::calc_dg_mrna;
use clap::Parser;
use constants::*;
use file_parser::{parse_csv_file, parse_fasta_file, parse_fasta_str, OstirInput};
use hybridization::{calc_dg_mrna_rrna, calc_dg_standby_site};
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use serde::Serialize;
use std::io;
use std::path::Path;
use types::DanglesSetting;

extern crate openmp_sys;

// ── Start codon energies (kcal/mol) ──────────────────────────────────────────
fn start_codon_energy(codon: &str) -> f64 {
    match codon.to_uppercase().replace('U', "T").as_str() {
        "ATG" => -1.194,
        "GTG" => -0.0748,
        "TTG" => -0.0435,
        "CTG" => -0.03406,
        _ => 0.0,
    }
}

// ── Expression level ─────────────────────────────────────────────────────────
fn calc_expression_level(dg_total: f64) -> f64 {
    k() * (-dg_total / rt_eff()).exp()
}

// ── Reverse complement ───────────────────────────────────────────────────────
fn reverse_complement(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            'U' => 'A',
            other => other,
        })
        .collect()
}

// ── Find start codons (0-indexed positions) ───────────────────────────────────
fn find_start_codons(sequence: &str, start_range: (usize, usize)) -> Vec<(usize, String)> {
    let seq_len = sequence.len();
    if seq_len < 3 {
        return vec![];
    }
    let end_0 = (start_range.1.saturating_sub(1)).min(seq_len - 3);
    let begin_0 = (start_range.0.saturating_sub(1)).min(end_0);
    let mut result = Vec::new();
    for i in begin_0..=end_0 {
        let codon = &sequence[i..i + 3];
        match codon.to_uppercase().as_str() {
            "ATG" | "AUG" | "GTG" | "GUG" | "TTG" | "UUG" | "CTG" | "CUG" => {
                result.push((i, codon.to_string()));
            }
            _ => {}
        }
    }
    result
}

// ── Output record ─────────────────────────────────────────────────────────────
#[derive(Debug, Clone, Serialize)]
pub struct OstirResult {
    pub name: String,
    pub start_codon: String,
    pub start_position: usize, // 1-indexed
    pub expression: f64,
    #[serde(rename = "RBS_distance_bp")]
    pub rbs_distance_bp: i64,
    #[serde(rename = "dG_total")]
    pub dg_total: f64,
    #[serde(rename = "dG_rRNA:mRNA")]
    pub dg_rrna_mrna: f64,
    #[serde(rename = "dG_mRNA")]
    pub dg_mrna: f64,
    #[serde(rename = "dG_spacing")]
    pub dg_spacing: f64,
    #[serde(rename = "dG_standby")]
    pub dg_standby: f64,
    #[serde(rename = "dG_start_codon")]
    pub dg_start_codon: f64,
}

const DECIMAL_PLACES: u32 = 4;
const HYBRIDIZATION_PENALTY: f64 = 2.481;

fn round4(x: f64) -> f64 {
    let factor = 10f64.powi(DECIMAL_PLACES as i32);
    (x * factor).round() / factor
}

/// Compute OSTIR results for a single (start_pos, codon) pair.
fn compute_one(
    sequence: &str,
    name: &str,
    start_pos: usize,
    codon: &str,
    asd: &str,
) -> Option<OstirResult> {
    let dangles = if start_pos > CUTOFF {
        DanglesSetting::new("none").unwrap()
    } else {
        DanglesSetting::new("all").unwrap()
    };

    let dg_mrna = calc_dg_mrna(sequence, start_pos, &dangles);
    let mrna_rrna_output = calc_dg_mrna_rrna(sequence, asd, start_pos, &dangles, None);

    let (dg_mrna_rrna_withspacing_raw, structure, spacing_value) = mrna_rrna_output?;

    let dg_mrna_rrna_withspacing = dg_mrna_rrna_withspacing_raw - HYBRIDIZATION_PENALTY;
    let dg_mrna_rrna_nospacing = structure.dg_mrna_rrna - HYBRIDIZATION_PENALTY;

    let dg_standby = calc_dg_standby_site(&structure, asd, &dangles, None);
    let dg_start_codon = start_codon_energy(codon);
    let dg_total = dg_mrna_rrna_withspacing + dg_start_codon - dg_mrna - dg_standby;
    let expression = calc_expression_level(dg_total);

    let rbs_distance_bp = if spacing_value.is_finite() {
        spacing_value.round() as i64
    } else {
        -1
    };

    Some(OstirResult {
        name: name.to_string(),
        start_codon: codon.to_uppercase().replace('T', "U"),
        start_position: start_pos + 1,
        expression: round4(expression),
        rbs_distance_bp,
        dg_total: round4(dg_total),
        dg_rrna_mrna: round4(dg_mrna_rrna_nospacing),
        dg_mrna: round4(dg_mrna),
        dg_spacing: round4(structure.dg_spacing),
        dg_standby: round4(dg_standby),
        dg_start_codon: round4(dg_start_codon),
    })
}

// ── Main OSTIR calculation ────────────────────────────────────────────────────
/// Run OSTIR on a sequence with threading and optional progress bar.
/// Returns results sorted by start_position.
pub fn ostir(
    sequence: &str,
    name: &str,
    start_range: (usize, usize), // 1-indexed
    asd: &str,
    circular: bool,
    bidirectional: bool,
    threads: usize,
    verbosity: u8,
) -> Vec<OstirResult> {
    // Circular: append sequence prefix so ORFs crossing the index are found
    let working_seq = if circular {
        let mrna_len = sequence.len();
        let append_len = if mrna_len > 200 { 200 } else { mrna_len };
        format!("{}{}", sequence, &sequence[..append_len])
    } else {
        sequence.to_string()
    };

    // Build list of (start_pos, codon) for the working sequence
    let seq_end = start_range.1.min(sequence.len()); // search only original length
    let start_codons = find_start_codons(&working_seq, (start_range.0, seq_end));

    // Optional progress bar (verbosity >= 1)
    let total_work = start_codons.len()
        + if bidirectional { start_codons.len() } else { 0 };
    let pb: Option<ProgressBar> = if verbosity >= 1 {
        let pb = ProgressBar::new(total_work as u64);
        pb.set_style(
            ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
        );
        Some(pb)
    } else {
        None
    };

    let pb_ref = &pb;

    // Forward pass — parallel when threads > 1, sequential otherwise
    let mut results: Vec<OstirResult> = if threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .unwrap()
            .install(|| {
                start_codons
                    .par_iter()
                    .filter_map(|(start_pos, codon)| {
                        let r = compute_one(&working_seq, name, *start_pos, codon, asd);
                        if let Some(pb) = pb_ref { pb.inc(1); }
                        r
                    })
                    .collect()
            })
    } else {
        start_codons
            .iter()
            .filter_map(|(start_pos, codon)| {
                let r = compute_one(&working_seq, name, *start_pos, codon, asd);
                if let Some(pb) = pb_ref { pb.inc(1); }
                r
            })
            .collect()
    };

    // Bidirectional: reverse complement pass
    if bidirectional {
        let rc = reverse_complement(sequence);
        let rc_len = rc.len();
        let rc_working = if circular {
            let append_len = if rc_len > 200 { 200 } else { rc_len };
            format!("{}{}", rc, &rc[..append_len])
        } else {
            rc.clone()
        };
        let rc_codons = find_start_codons(&rc_working, (1, rc_len.min(start_range.1)));
        let rc_results: Vec<OstirResult> = if threads > 1 {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| {
                    rc_codons
                        .par_iter()
                        .filter_map(|(start_pos, codon)| {
                            let mut r = compute_one(&rc_working, name, *start_pos, codon, asd)?;
                            r.start_position = sequence.len() - start_pos;
                            if let Some(pb) = pb_ref { pb.inc(1); }
                            Some(r)
                        })
                        .collect()
                })
        } else {
            rc_codons
                .iter()
                .filter_map(|(start_pos, codon)| {
                    let mut r = compute_one(&rc_working, name, *start_pos, codon, asd)?;
                    r.start_position = sequence.len() - start_pos;
                    if let Some(pb) = pb_ref { pb.inc(1); }
                    Some(r)
                })
                .collect()
        };
        results.extend(rc_results);
    }

    if let Some(pb) = pb {
        pb.finish_and_clear();
    }

    results.sort_by_key(|r| r.start_position);
    results
}

// ── CSV output helpers ────────────────────────────────────────────────────────

/// Write results to CSV, optionally including sequence and/or asd columns.
/// Column order matches Python: name, [anti-Shine-Dalgarno,] [sequence,] start_codon, ...
fn write_results<W: io::Write>(
    wtr: &mut csv::Writer<W>,
    results: &[OstirResult],
    sequence_col: Option<&str>,
    asd_col: Option<&str>,
    headers_written: &mut bool,
) {
    // Build column list
    let mut cols: Vec<&str> = vec!["name"];
    if asd_col.is_some() { cols.push("anti-Shine-Dalgarno"); }
    if sequence_col.is_some() { cols.push("sequence"); }
    cols.extend_from_slice(&["start_codon", "start_position", "expression",
        "RBS_distance_bp", "dG_total", "dG_rRNA:mRNA", "dG_mRNA",
        "dG_spacing", "dG_standby", "dG_start_codon"]);

    if !*headers_written {
        wtr.write_record(&cols).expect("Failed to write CSV header");
        *headers_written = true;
    }

    for r in results {
        let mut row: Vec<String> = vec![r.name.clone()];
        if let Some(a) = asd_col { row.push(a.to_string()); }
        if let Some(s) = sequence_col { row.push(s.to_string()); }
        row.push(r.start_codon.clone());
        row.push(r.start_position.to_string());
        row.push(r.expression.to_string());
        row.push(r.rbs_distance_bp.to_string());
        row.push(r.dg_total.to_string());
        row.push(r.dg_rrna_mrna.to_string());
        row.push(r.dg_mrna.to_string());
        row.push(r.dg_spacing.to_string());
        row.push(r.dg_standby.to_string());
        row.push(r.dg_start_codon.to_string());
        wtr.write_record(&row).expect("Failed to write CSV record");
    }
}

#[derive(Parser, Debug)]
#[command(name = "ostir", about = "Open Source Translation Initiation Rates")]
struct Cli {
    /// Input: FASTA file path, CSV file path, or raw sequence string
    #[arg(short = 'i', long)]
    input: String,

    /// Output CSV file (stdout if not provided)
    #[arg(short = 'o', long)]
    output: Option<String>,

    /// Start position of search range, 1-indexed (default: 1)
    #[arg(short = 's', long, default_value = "1")]
    start: usize,

    /// End position of search range, 1-indexed (default: sequence length)
    #[arg(short = 'e', long)]
    end: Option<usize>,

    /// Anti-Shine-Dalgarno sequence (3' end of 16S rRNA)
    #[arg(short = 'a', long, default_value = "ACCTCCTTA")]
    asd: String,

    /// Number of threads for parallel per-codon calculations
    #[arg(short = 'j', long, default_value = "1")]
    threads: usize,

    /// Treat sequence as circular (ORFs may cross the start/end index)
    #[arg(short = 'c', long)]
    circular: bool,

    /// Also search reverse complement strand
    #[arg(short = 'b', long)]
    bidirectional: bool,

    /// Include input mRNA sequence in output CSV
    #[arg(short = 'p', long)]
    print_sequence: bool,

    /// Include anti-Shine-Dalgarno sequence in output CSV
    #[arg(short = 'q', long)]
    print_asd: bool,

    /// Verbosity: 0 = silent, 1 = progress bar
    #[arg(short = 'v', long, default_value = "0")]
    verbosity: u8,
}

fn main() {
    let cli = Cli::parse();

    if cli.verbosity > 1 {
        eprintln!("Warning: verbosity levels above 1 are treated as 1");
    }
    let verbosity = cli.verbosity.min(1);

    // Build a default OstirInput from CLI args
    let defaults = OstirInput {
        sequence: String::new(),
        name: "unnamed".to_string(),
        asd: cli.asd.clone(),
        start: cli.start,
        end: cli.end,
        circular: cli.circular,
        print_sequence: cli.print_sequence,
        print_asd: cli.print_asd,
    };

    // Resolve inputs
    let inputs: Vec<OstirInput> = {
        let p = Path::new(&cli.input);
        if p.exists() {
            let ext = p.extension().and_then(|e| e.to_str()).unwrap_or("");
            if ext == "csv" {
                parse_csv_file(p, &defaults).expect("Failed to parse CSV input")
            } else {
                // FASTA (or any other file)
                parse_fasta_file(p, &defaults).expect("Failed to parse FASTA input")
            }
        } else {
            // Raw sequence string
            vec![OstirInput {
                sequence: cli.input.clone(),
                name: "unnamed".to_string(),
                ..defaults
            }]
        }
    };

    // Set up CSV writer
    let writer: Box<dyn io::Write> = match &cli.output {
        Some(path) => Box::new(
            std::fs::File::create(path).expect("Failed to create output file"),
        ),
        None => Box::new(io::stdout()),
    };
    let mut wtr = csv::WriterBuilder::new().has_headers(false).from_writer(writer);
    let mut headers_written = false;

    for input in inputs {
        let end = input.end.unwrap_or(input.sequence.len());
        let results = ostir(
            &input.sequence,
            &input.name,
            (input.start, end),
            &input.asd,
            input.circular,
            cli.bidirectional,
            cli.threads,
            verbosity,
        );
        let seq_col = if input.print_sequence { Some(input.sequence.as_str()) } else { None };
        let asd_col = if input.print_asd { Some(input.asd.as_str()) } else { None };
        write_results(&mut wtr, &results, seq_col, asd_col, &mut headers_written);
    }
    wtr.flush().expect("Failed to flush CSV writer");
}

// ── Tests ─────────────────────────────────────────────────────────────────────
#[cfg(test)]
mod tests {
    use super::*;

    const TEST_SEQ: &str = "ACUUCUAAUUUAUUCUAUUUAUUCGCGGAUAUGCAUAGGAGUGCUUCGAUGUCAU";
    const DEFAULT_ASD: &str = "ACCTCCTTA";

    fn run(seq: &str, start: usize, end: usize, asd: &str) -> Vec<OstirResult> {
        ostir(seq, "unnamed", (start, end), asd, false, false, 1, 0)
    }

    #[test]
    fn test_unit_one_sequence() {
        let results = run(TEST_SEQ, 1, TEST_SEQ.len(), DEFAULT_ASD);
        assert_eq!(results.len(), 3);

        let r0 = &results[0];
        assert_eq!(r0.start_codon, "AUG");
        assert_eq!(r0.start_position, 31);
        assert_eq!(r0.expression, 2.1121);
        assert_eq!(r0.rbs_distance_bp, 16);
        assert_eq!(r0.dg_total, 16.3277);
        assert_eq!(r0.dg_rrna_mrna, -0.481);
        assert_eq!(r0.dg_mrna, -7.2);
        assert_eq!(r0.dg_spacing, 10.8027);
        assert_eq!(r0.dg_standby, 0.0);
        assert_eq!(r0.dg_start_codon, -1.194);

        let r1 = &results[1];
        assert_eq!(r1.start_codon, "GUG");
        assert_eq!(r1.start_position, 41);
        assert_eq!(r1.expression, 845.1532);
        assert_eq!(r1.rbs_distance_bp, 8);
        assert_eq!(r1.dg_total, 1.3491);
        assert_eq!(r1.dg_rrna_mrna, -3.581);
        assert_eq!(r1.dg_mrna, -3.6);
        assert_eq!(r1.dg_spacing, 1.4049);
        assert_eq!(r1.dg_standby, 0.0);
        assert_eq!(r1.dg_start_codon, -0.0748);

        let r2 = &results[2];
        assert_eq!(r2.start_codon, "AUG");
        assert_eq!(r2.start_position, 49);
        assert_eq!(r2.expression, 12448.1756);
        assert_eq!(r2.rbs_distance_bp, 5);
        assert_eq!(r2.dg_total, -5.375);
        assert_eq!(r2.dg_rrna_mrna, -7.981);
        assert_eq!(r2.dg_mrna, -3.6);
        assert_eq!(r2.dg_spacing, 0.0);
        assert_eq!(r2.dg_standby, -0.2);
        assert_eq!(r2.dg_start_codon, -1.194);
    }

    #[test]
    fn test_unit_one_sequence_new_asd() {
        let results = run(TEST_SEQ, 31, 31, "ACGTCCCTA");
        assert_eq!(results.len(), 1);
        let r = &results[0];
        assert_eq!(r.start_codon, "AUG");
        assert_eq!(r.start_position, 31);
        assert_eq!(r.expression, 367.8105);
        assert_eq!(r.rbs_distance_bp, 4);
        assert_eq!(r.dg_total, 3.4289);
        assert_eq!(r.dg_rrna_mrna, -2.581);
        assert_eq!(r.dg_mrna, -7.2);
        assert_eq!(r.dg_spacing, 0.0039);
        assert_eq!(r.dg_standby, 0.0);
        assert_eq!(r.dg_start_codon, -1.194);
    }

    /// Helper: run CLI and compare CSV output vs expected file.
    fn run_cli_csv(args: &[&str], expected_path: &str) {
        use std::path::PathBuf;
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let expected = PathBuf::from(manifest_dir).join("tests").join("expected").join(expected_path);
        let out_dir = PathBuf::from(manifest_dir).join("tests").join("output");
        std::fs::create_dir_all(&out_dir).ok();
        let out_path = out_dir.join(expected_path);

        // Build args list including output path
        let mut full_args: Vec<String> = args.iter().map(|s| s.to_string()).collect();
        full_args.push("-o".to_string());
        full_args.push(out_path.to_str().unwrap().to_string());

        // Build OstirInput list same way main() does it
        let mut cli_args = vec!["ostir".to_string()];
        cli_args.extend(full_args);
        let cli = Cli::try_parse_from(&cli_args).expect("CLI parse failed");

        let defaults = OstirInput {
            sequence: String::new(),
            name: "unnamed".to_string(),
            asd: cli.asd.clone(),
            start: cli.start,
            end: cli.end,
            circular: cli.circular,
            print_sequence: cli.print_sequence,
            print_asd: cli.print_asd,
        };

        let inputs: Vec<OstirInput> = {
            let p = Path::new(&cli.input);
            if p.exists() {
                let ext = p.extension().and_then(|e| e.to_str()).unwrap_or("");
                if ext == "csv" {
                    parse_csv_file(p, &defaults).expect("Failed to parse CSV")
                } else {
                    parse_fasta_file(p, &defaults).expect("Failed to parse FASTA")
                }
            } else {
                vec![OstirInput { sequence: cli.input.clone(), name: "unnamed".to_string(), ..defaults }]
            }
        };

        let writer = std::fs::File::create(&out_path).expect("create output");
        let mut wtr = csv::WriterBuilder::new().has_headers(false).from_writer(writer);
        let mut headers_written = false;
        for input in inputs {
            let end = input.end.unwrap_or(input.sequence.len());
            let results = ostir(&input.sequence, &input.name, (input.start, end), &input.asd, input.circular, cli.bidirectional, cli.threads, 0);
            let seq_col = if input.print_sequence { Some(input.sequence.as_str()) } else { None };
            let asd_col = if input.print_asd { Some(input.asd.as_str()) } else { None };
            write_results(&mut wtr, &results, seq_col, asd_col, &mut headers_written);
        }
        wtr.flush().unwrap();

        // Compare CSVs
        let out_data = std::fs::read_to_string(&out_path).unwrap();
        let exp_data = std::fs::read_to_string(&expected).unwrap_or_else(|_| panic!("Missing expected file: {}", expected.display()));
        assert_eq!(out_data, exp_data, "CSV mismatch for {}", expected_path);
    }

    #[test]
    fn test_integration_fasta_input() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let input_path = format!("{}/tests/input/command_line_FASTA_input.fa", manifest_dir);
        run_cli_csv(&["-j", "4", "-i", &input_path, "-v", "0"], "command_line_FASTA_input.csv");
    }

    #[test]
    fn test_integration_string_input() {
        run_cli_csv(
            &["-j", "4", "-p", "-i", "TTCTAGAAAAAAAATAAGGAGGTAAAATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATG", "-v", "0"],
            "command_line_string_input.csv",
        );
    }

    #[test]
    fn test_integration_csv_input() {
        let manifest_dir = env!("CARGO_MANIFEST_DIR");
        let input_path = format!("{}/tests/input/command_line_CSV_input.csv", manifest_dir);
        run_cli_csv(&["-j", "4", "-i", &input_path, "-v", "0"], "command_line_CSV_input.csv");
    }

    #[test]
    fn test_vienna_mfe_after_subopt() {
        use crate::vienna_wrapper::{mfe, subopt};
        let dangles = DanglesSetting::new("all").unwrap();
        // Step 1: mfe on full seq (calc_dg_mrna)
        let seqs1 = vec!["ATAAGGAGGTATG"];
        let r1 = mfe(&seqs1, "", 37.0, &dangles);
        assert!(r1.is_ok(), "mfe on full seq should succeed");
        drop(r1);
        // Step 2: subopt on mrna + rrna
        let seqs2 = vec!["ATAAGGAGGT", DEFAULT_ASD];
        let r2 = subopt(&seqs2, "", 3.0, 37.0, &dangles);
        assert!(!r2.is_empty(), "subopt should return results");
        drop(r2);
        // Step 3: mfe on pre-sequence "A" (single char)
        let seqs3 = vec!["A"];
        let r3 = mfe(&seqs3, "", 37.0, &dangles);
        assert!(r3.is_ok(), "mfe on 'A' should succeed");
        // Step 4: mfe on "AA" (two chars)
        let seqs4 = vec!["AA"];
        let r4 = mfe(&seqs4, "", 37.0, &dangles);
        assert!(r4.is_ok(), "mfe on 'AA' should succeed");
    }

    #[test]
    fn test_fasta_seq() {
        let results = run("ATAAGGAGGTATG", 1, 13, DEFAULT_ASD);
        assert!(!results.is_empty(), "Expected at least one result");
        // Previously crashed due to ViennaRNA density_of_states negative index bug
        let r = &results[0];
        assert_eq!(r.start_codon, "AUG");
        assert_eq!(r.start_position, 11);
    }

    #[test]
    fn test_fasta_seq2() {
        let results = run("TTCTAGAAAAAAAATAAGGAGGTATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAA", 1, 80, DEFAULT_ASD);
        assert!(!results.is_empty(), "Expected at least one result");
    }
}
