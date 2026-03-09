use crate::constants::*;
use crate::types::DanglesSetting;
use crate::vienna_wrapper::mfe;
use std::cmp::{max, min};

/// Calculates the minimum free energy of the mRNA structure in the window around start_pos.
/// Returns the MFE energy in kcal/mol.
pub fn calc_dg_mrna(mrna: &str, start_pos: usize, dangles: &DanglesSetting) -> f64 {
    let begin = if start_pos > CUTOFF { start_pos - CUTOFF } else { 0 };
    let end = min(mrna.len(), start_pos + CUTOFF);
    let trimmed = &mrna[begin..end];
    let seqs = vec![trimmed];
    match mfe(&seqs, "", DEFAULT_TEMP, dangles) {
        Ok(fold) => {
            // Round to 2 decimal places to match Python's ostir ViennaRNA.mfe() behavior
            (*fold.get_d_g() as f64 * 100.0).round() / 100.0
        }
        Err(_) => 0.0,
    }
}

/// Calculates a dG-like penalty for the ribosome binding away from the optimal start position.
pub fn calc_spacing_penalty(aligned_spacing: usize) -> f64 {
    let ds: f64 = (aligned_spacing as isize - OPTIMAL_SPACING as isize) as f64;
    if aligned_spacing < OPTIMAL_SPACING {
        DG_SPACING_PUSH[0]
            / (1.0 + (DG_SPACING_PUSH[1] * (ds + DG_SPACING_PUSH[2])).exp())
                .powf(DG_SPACING_PUSH[3])
    } else {
        DG_SPACING_PULL[0] * ds * ds + DG_SPACING_PULL[1] * ds + DG_SPACING_PULL[2]
    }
}

#[allow(dead_code)]
fn cutoff_mrna<'a>(mrna: &'a str, start_pos: usize) -> &'a str {
    &mrna[max(0, start_pos as isize - CUTOFF as isize) as usize..min(mrna.len(), start_pos + CUTOFF)]
}
