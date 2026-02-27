mod calculations;
mod constants;
mod file_parser;
mod hybridization;
mod types;
mod vienna_wrapper;

use calculations::calc_dg_mrna;
use clap::Parser;
use constants::*;
use hybridization::{calc_dg_mrna_rrna, calc_dg_standby_site};
use serde::Serialize;
use std::io;
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

// ── Main OSTIR calculation ────────────────────────────────────────────────────
pub fn ostir(
    sequence: &str,
    name: &str,
    start_range: (usize, usize), // 1-indexed
    asd: &str,
) -> Vec<OstirResult> {
    let mut results: Vec<OstirResult> = Vec::new();

    let start_codons = find_start_codons(sequence, start_range);

    for (start_pos, codon) in start_codons {
        // Set dangles based on distance from 5' end (Python: "all" if start_pos <= cutoff)
        let dangles = if start_pos > CUTOFF {
            DanglesSetting::new("none").unwrap()
        } else {
            DanglesSetting::new("all").unwrap()
        };

        // dG of mRNA secondary structure
        let dg_mrna = calc_dg_mrna(sequence, start_pos, &dangles);

        // dG of mRNA:rRNA hybridization
        let mrna_rrna_output =
            calc_dg_mrna_rrna(sequence, asd, start_pos, &dangles, None);

        let (dg_mrna_rrna_withspacing_raw, structure, spacing_value) = match mrna_rrna_output {
            Some(v) => v,
            None => continue,
        };

        // Apply hybridization penalty correction to match NUPACK
        let dg_mrna_rrna_withspacing = dg_mrna_rrna_withspacing_raw - HYBRIDIZATION_PENALTY;
        let dg_mrna_rrna_nospacing = structure.dg_mrna_rrna - HYBRIDIZATION_PENALTY;

        // Standby site penalty
        let dg_standby = calc_dg_standby_site(&structure, asd, &dangles, None);

        // Start codon energy
        let dg_start_codon = start_codon_energy(&codon);

        // Total free energy change
        let dg_total = dg_mrna_rrna_withspacing + dg_start_codon - dg_mrna - dg_standby;

        // Expression level
        let expression = calc_expression_level(dg_total);

        let rbs_distance_bp = if spacing_value.is_finite() {
            spacing_value.round() as i64
        } else {
            -1
        };

        results.push(OstirResult {
            name: name.to_string(),
            start_codon: codon.to_uppercase().replace('T', "U"),
            start_position: start_pos + 1, // 1-indexed output
            expression: round4(expression),
            rbs_distance_bp,
            dg_total: round4(dg_total),
            dg_rrna_mrna: round4(dg_mrna_rrna_nospacing),
            dg_mrna: round4(dg_mrna),
            dg_spacing: round4(structure.dg_spacing),
            dg_standby: round4(dg_standby),
            dg_start_codon: round4(dg_start_codon),
        });
    }

    results.sort_by_key(|r| r.start_position);
    results
}

// ── CLI ───────────────────────────────────────────────────────────────────────
#[derive(Parser, Debug)]
#[command(name = "ostir", about = "Open Source Translation Initiation Rates")]
struct Cli {
    /// Input: FASTA file path or raw sequence string
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

    /// Number of threads (currently unused)
    #[arg(short = 'j', long, default_value = "1")]
    threads: usize,

    /// Circular sequence
    #[arg(short = 'c', long)]
    circular: bool,

    /// Verbosity level
    #[arg(short = 'v', action = clap::ArgAction::Count)]
    verbosity: u8,
}

fn main() {
    let cli = Cli::parse();

    // Determine input sequences: either a FASTA file or a raw sequence string
    let sequences: Vec<(String, String)> = if std::path::Path::new(&cli.input).exists() {
        // Parse as FASTA
        let content =
            std::fs::read_to_string(&cli.input).expect("Failed to read input file");
        parse_fasta_text(&content)
    } else {
        // Treat as raw sequence
        vec![("unnamed".to_string(), cli.input.clone())]
    };

    // Set up CSV writer
    let writer: Box<dyn io::Write> = match &cli.output {
        Some(path) => Box::new(
            std::fs::File::create(path).expect("Failed to create output file"),
        ),
        None => Box::new(io::stdout()),
    };

    let mut wtr = csv::Writer::from_writer(writer);

    for (name, seq) in sequences {
        let end = cli.end.unwrap_or(seq.len());
        let results = ostir(&seq, &name, (cli.start, end), &cli.asd);
        for r in results {
            wtr.serialize(&r).expect("Failed to write CSV record");
        }
    }
    wtr.flush().expect("Failed to flush CSV writer");
}

/// Parse multi-FASTA text, returning (name, sequence) pairs.
fn parse_fasta_text(content: &str) -> Vec<(String, String)> {
    let mut sequences = Vec::new();
    let mut current_name = String::new();
    let mut current_seq = String::new();

    for line in content.lines() {
        let line = line.trim();
        if line.starts_with('>') {
            if !current_seq.is_empty() {
                sequences.push((current_name.clone(), current_seq.clone()));
                current_seq.clear();
            }
            current_name = line[1..].trim().to_string();
            if current_name.is_empty() {
                current_name = "unnamed".to_string();
            }
        } else if !line.is_empty() {
            current_seq.push_str(line);
        }
    }
    if !current_seq.is_empty() {
        sequences.push((current_name, current_seq));
    }
    sequences
}
