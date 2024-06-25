mod calculations;
mod constants;
mod file_parser;
mod types;
pub use file_parser::fileparser;
pub use file_parser::fileparser::DNASequence;
use indicatif::ProgressBar;
mod vienna_wrapper;
use polars::prelude::*;
use vienna_wrapper::{mfe, subopt};
extern crate openmp_sys;

fn main() {
    println!("Hello, world!");
    let sequences = fileparser::parse_file("test.fasta", 20).unwrap();

    for sequence in sequences {
        for _next in sequence {
            //println!("{}", next.sequence);
        }
    }

    let test_seq1 = "TTATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCATGCGTTTCAAAGTTCGTATGGAAGGT";
    let _test_seq2 = "TTATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCAT";
    let _test_seq_vec1 = vec![test_seq1];
    let _test_seq_vec2 = vec![test_seq1, _test_seq2];
    test_mfe(&_test_seq_vec1);
    test_subopt(&_test_seq_vec2);
}

fn ostir(
    sequence: DNASequence,
    start: i64,
    end: i64,
    name: &str,
    asd: &str,
    circular: bool,
    threads: i32,
    bidirectional: bool,
    verbosity: i32,
) {
    // Get start codon positions
    let mut start_codon_positions: Vec<(usize, &str)> = vec![];
    let start_codons: Vec<String> = vec![
        "ATG".to_string(),
        "GTG".to_string(),
        "TTG".to_string(),
        "CTG".to_string(),
        "AUG".to_string(),
        "GUG".to_string(),
        "UUG".to_string(),
        "CUG".to_string(),
    ];
    let tgt_sequence = &sequence.record;
    for start_codon in start_codons {
        let mut result: Vec<_> = tgt_sequence.match_indices(&start_codon).collect();
        start_codon_positions.append(&mut result);
    }

    // Create Dataframe
    let (_position, _codon): (Vec<usize>, Vec<&str>) =
        start_codon_positions.iter().cloned().unzip();
    let _df: DataFrame = df!(
        "start_base" => Vec::<i64>::new(),
        "start_codon" => Vec::<&str>::new(),
    )
    .unwrap();

    // Set up progress bar
    let _bar = ProgressBar::new(start_codon_positions.len().try_into().unwrap());

    // Run calculations

    // Return results
}

fn test_mfe(test_seq_vec: &Vec<&str>) {
    // test mfe
    println!("MFE");

    let mfe = mfe(
        test_seq_vec,
        "placeholder",
        37.0,
        &types::DanglesSetting::new("all").unwrap(),
    )
    .unwrap();
    // print mfe result as f32 and as string
    println!("{}", mfe.get_d_g());
    println!("{}", mfe.get_dots());
    println!("{:?}", mfe.get_bp_x());
    println!("{:?}", mfe.get_bp_y());
}

fn test_subopt(test_seq_vec: &Vec<&str>) {
    // test subopt
    println!("Subopt");

    let subopt = subopt(
        test_seq_vec,
        "placeholder",
        50.0,
        37.0,
        &types::DanglesSetting::new("all").unwrap(),
    );

    // print subopt result
    for result in subopt {
        let dg = result.get_d_g();
        let bp_x = result.get_bp_x();
        let bp_y = result.get_bp_y();
        println!("{}", dg);
        println!("{}", result.get_dots());
        println!("{:?}", bp_x);
        println!("{:?}", bp_y);
    }
}
