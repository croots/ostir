use crate::constants::*;
use crate::types::{CoFoldResult, DanglesSetting, FoldResult, MonoFoldResult};
use crate::vienna_wrapper::*;
use std::cmp::{max, min};
use std::error::Error;

///Calculate a "kinetic score", a heuristic measure of the maximum time required for
///the mRNA secondary structure to form. This is related to the RNA polymer model by David et. al.
///This heuristic should not be used in any way to quantify the folding kinetics of an mRNA sequence
///because it completely ignores cooperative RNA folding mechanisms, such as zipping or strand
///displacement. Here, we use it to eliminate mRNA sequences that MAY fold slowly.
fn calc_kinetic_score(fold: &MonoFoldResult) -> Result<(f64, f64), Box<dyn Error>> {
    let mrnalen = fold.seqs.len() as usize;
    let mut largest_range_helix = 0;

    for (nt_x, nt_y) in fold.bp_x.iter().zip(fold.bp_y.iter()) {
        if *nt_x <= mrnalen && *nt_y <= mrnalen {
            let val = nt_y - nt_x;
            largest_range_helix = largest_range_helix.max(val);
        }
    }

    let kinetic_score = largest_range_helix as f64 / mrnalen as f64;
    let min_bp_prob = if largest_range_helix > 0 {
        (largest_range_helix as f64).powf(-1.44)
    } else {
        1.0
    };

    Ok((kinetic_score, min_bp_prob))
}

fn calc_dg_mrna<'a>(
    trimmed_mrna: &'a Vec<&'a str>,
    start_pos: usize,
    dangles: &'_ DanglesSetting,
    constraints: &'_ str,
) -> Result<FoldResult<'a>, Box<dyn Error>> {
    // Calculates the dG_mRNA given the mRNA sequence

    let constraints = &constraints[std::cmp::max(0, start_pos as isize - CUTOFF as isize) as usize
        ..std::cmp::min(trimmed_mrna.len(), start_pos + CUTOFF)];

    Ok(mfe(trimmed_mrna, constraints, DEFAULT_TEMP, dangles)?)
}

///Calculates the dG_standby given the structure of the mRNA:rRNA complex
///To calculate the mfe structure while disallowing base pairing at the standby site,
///we split the folded mRNA sequence into three parts: (i) a pre-sequence (before the standby
///site) that can fold; (ii) the standby site, which can not fold; (iii) the 16S rRNA binding
///site and downstream sequence, which has been previously folded.
fn calc_dg_standby_site(fold: &CoFoldResult, dangles: &DanglesSetting, constraints: &str) -> f64 {
    //todo
    todo!();
    // Put the sets of base pairs together

    // Calculate its energy
}

///Calculates a dG-like penalty for the ribosome binding away from the optimal start position
fn calc_spacing_penalty(aligned_spacing: usize) -> f64 {
    let ds: f64;
    let spacing_penalty: f64;
    ds = (aligned_spacing - OPTIMAL_SPACING) as f64;
    if aligned_spacing < OPTIMAL_SPACING {
        spacing_penalty = DG_SPACING_PUSH[0]
            / (1.0 + (DG_SPACING_PUSH[1] * (ds + DG_SPACING_PUSH[2])).exp())
                .powf(DG_SPACING_PUSH[3]);
    } else {
        spacing_penalty =
            DG_SPACING_PULL[0] * ds * ds + DG_SPACING_PULL[1] * ds + DG_SPACING_PULL[2];
    }
    return spacing_penalty;
}

///Figure out where exactly the ribosome is binding
fn find_binding_position(start_pos: usize, fold: &CoFoldResult) -> Result<usize, Box<dyn Error>> {
    let seqs = fold.seqs;

    let mut last_bound_rrna: Option<usize> = None;
    let mut last_bound_mrna: Option<usize> = None;
    let len_rrna = seqs.1.as_bytes().len() as usize;
    let len_mrna = seqs.0.as_bytes().len() as usize;
    let num_folded_bases = fold.bp_x.len();
    let mut counter: usize = 1;

    while last_bound_rrna.is_none() {
        match (
            fold.bp_y[num_folded_bases - counter],
            fold.bp_x[num_folded_bases - counter],
        ) {
            (rrna, _) if rrna < len_mrna => (), // This case is mRNA backfolding
            (_, mrna) if mrna >= start_pos as usize => {
                return Err("Ribosome is sitting on the start codon")?
            }
            (rrna, mrna) => (last_bound_rrna, last_bound_mrna) = (Some(rrna), Some(mrna)),
            _ => return Err("Ran out of bases")?,
        };

        counter = counter + 1
    }

    let calibrated_position: usize = (start_pos - last_bound_mrna.unwrap() as usize)
        - (last_bound_rrna.unwrap() as usize - len_rrna as usize);

    Ok(calibrated_position)
}

fn cutoff_mrna<'a>(mrna: &'a str, start_pos: usize) -> &'a str {
    &mrna
        [max(0, start_pos as isize - CUTOFF as isize) as usize..min(mrna.len(), start_pos + CUTOFF)]
}
