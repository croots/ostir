use crate::constants::*;
use crate::types::DanglesSetting;
use crate::vienna_wrapper::{coordinates_to_dots, eval_structure, mfe, subopt};
use std::cmp::min;

pub struct MRNARRNAStructure {
    pub mrna: String,
    pub dg_mrna_rrna: f64,
    pub dg_mrna_rrna_withspacing: f64,
    pub dg_spacing: f64,
    pub bp_x: Vec<usize>,
    pub bp_y: Vec<usize>,
}

fn calc_dg_spacing(aligned_spacing: f64) -> f64 {
    if aligned_spacing.is_infinite() || aligned_spacing > 1e8 {
        return 1e9;
    }
    let ds = aligned_spacing - OPTIMAL_SPACING as f64;
    if aligned_spacing < OPTIMAL_SPACING as f64 {
        DG_SPACING_PUSH[0]
            / (1.0 + (DG_SPACING_PUSH[1] * (ds + DG_SPACING_PUSH[2])).exp())
                .powf(DG_SPACING_PUSH[3])
    } else {
        DG_SPACING_PULL[0] * ds * ds + DG_SPACING_PULL[1] * ds + DG_SPACING_PULL[2]
    }
}

/// Calculates the aligned spacing between the 16S rRNA binding site and the start codon.
/// Returns infinity if no suitable binding is found.
pub fn calc_aligned_spacing(
    rrna: &str,
    mrna: &str,
    start_pos: usize,
    bp_x: &[usize],
    bp_y: &[usize],
) -> f64 {
    let rrna_len = rrna.len();
    let mrna_len = mrna.len();
    let seq_len = mrna_len + rrna_len;

    // Iterate from farthest 3' rRNA nt (seq_len, 1-indexed) down to first rRNA nt (mrna_len+1)
    for rrna_nt in (mrna_len + 1..=seq_len).rev() {
        if let Some(rrna_pos) = bp_y.iter().position(|&y| y == rrna_nt) {
            let x_start = bp_x[rrna_pos];
            if x_start < start_pos {
                let farthest_3_prime_rrna = (rrna_nt - mrna_len) as f64;
                let mrna_nt = bp_x[rrna_pos];
                let distance_to_start = (start_pos - mrna_nt + 1) as f64;
                return distance_to_start - farthest_3_prime_rrna;
            } else {
                break;
            }
        }
    }

    f64::INFINITY
}

/// Calculates the dG_mRNA_rRNA from the mRNA and rRNA sequence.
/// Considers all feasible 16S rRNA binding sites and includes the effects of non-optimal spacing.
/// Returns (dG_mRNA_rRNA_withspacing, MRNARRNAStructure, spacing_value) or None.
pub fn calc_dg_mrna_rrna(
    mrna_in: &str,
    rrna: &str,
    start_pos: usize,
    dangles: &DanglesSetting,
    _constraints: Option<&str>,
) -> Option<(f64, MRNARRNAStructure, f64)> {
    let begin = if start_pos > CUTOFF { start_pos - CUTOFF } else { 0 };
    let mrna_len = min(mrna_in.len(), start_pos.saturating_add(CUTOFF));
    let start_pos_in_subsequence = min(start_pos, CUTOFF);
    let startpos_to_end_len = mrna_len - start_pos_in_subsequence - begin;

    let mrna = &mrna_in[begin..start_pos];
    if mrna.is_empty() {
        return None; // leaderless start codon
    }

    let seqs = vec![mrna, rrna];
    let subopt_results = subopt(
        &seqs,
        "",
        ENERGY_CUTOFF,
        DEFAULT_TEMP,
        dangles,
    );

    // Calculate aligned spacing and energy+spacing for each result
    let mut aligned_spacings: Vec<f64> = Vec::new();
    let mut dg_mrna_rrna_list: Vec<f64> = Vec::new();
    let mut dg_spacing_list: Vec<f64> = Vec::new();
    let mut dg_with_spacing_list: Vec<f64> = Vec::new();

    for result in &subopt_results {
        let spacing = calc_aligned_spacing(
            rrna,
            mrna,
            start_pos_in_subsequence,
            result.get_bp_x(),
            result.get_bp_y(),
        );
        let energy = *result.get_d_g() as f64;
        let spacing_penalty = calc_dg_spacing(spacing);
        aligned_spacings.push(spacing);
        dg_mrna_rrna_list.push(energy);
        dg_spacing_list.push(spacing_penalty);
        dg_with_spacing_list.push(energy + spacing_penalty);
    }

    // Find minimum energy + spacing
    let (index, _) = dg_with_spacing_list
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())?;

    let dg_spacing_final = dg_spacing_list[index];
    let spacing_value = aligned_spacings[index];

    let best_result = &subopt_results[index];
    let bp_x = best_result.get_bp_x();
    let bp_y = best_result.get_bp_y();

    let mrna_len_local = mrna.len();

    // Identify rRNA-mRNA base pairs and find most 5' mRNA binding position
    let mut most_5p_mrna = usize::MAX;
    let mut bp_x_target: Vec<usize> = Vec::new();
    let mut bp_y_target: Vec<usize> = Vec::new();

    for (&nt_x, &nt_y) in bp_x.iter().zip(bp_y.iter()) {
        if nt_y > mrna_len_local {
            most_5p_mrna = most_5p_mrna.min(nt_x);
            bp_x_target.push(nt_x);
            bp_y_target.push(nt_y);
        }
    }

    if most_5p_mrna == usize::MAX {
        return None; // No rRNA binding found
    }

    // Build pre-sequence (mRNA before rRNA binding site)
    // most_5p_mrna is 1-indexed in mrna; begin + most_5p_mrna - 1 is the position in mrna_in
    let pre_end = begin.saturating_add(most_5p_mrna.saturating_sub(1));
    let mrna_pre = if pre_end > begin { &mrna_in[begin..pre_end] } else { "" };

    // Post-sequence (footprint=1000 means it's always empty in practice)
    let post_window_begin = min(start_pos + FOOTPRINT as usize, mrna_len);
    let mrna_post = if post_window_begin < mrna_len {
        &mrna_in[post_window_begin..mrna_len]
    } else {
        ""
    };

    let mut total_bp_x: Vec<usize> = Vec::new();
    let mut total_bp_y: Vec<usize> = Vec::new();

    // 1. Pre-sequence MFE folding
    if !mrna_pre.is_empty() {
        let seqs_pre = vec![mrna_pre];
        if let Ok(fold) = mfe(&seqs_pre, "", DEFAULT_TEMP, dangles) {
            total_bp_x.extend_from_slice(fold.get_bp_x());
            total_bp_y.extend_from_slice(fold.get_bp_y());
        }
    }

    // 2. rRNA binding site pairs, with rRNA offset applied
    let rrna_offset = startpos_to_end_len;
    total_bp_x.extend_from_slice(&bp_x_target);
    for &nt_y in &bp_y_target {
        total_bp_y.push(nt_y + rrna_offset);
    }

    // 3. Post-sequence MFE folding (typically empty)
    if !mrna_post.is_empty() {
        let seqs_post = vec![mrna_post];
        if let Ok(fold) = mfe(&seqs_post, "", DEFAULT_TEMP, dangles) {
            let offset = post_window_begin - begin;
            for (&nt_x, &nt_y) in fold.get_bp_x().iter().zip(fold.get_bp_y().iter()) {
                total_bp_x.push(nt_x + offset);
                total_bp_y.push(nt_y + offset);
            }
        }
    }

    // Evaluate total energy using the full window mRNA + rRNA
    let mrna_full = &mrna_in[begin..mrna_len];
    let dots = coordinates_to_dots(&vec![mrna_full, rrna], &total_bp_x, &total_bp_y);
    let total_energy =
        eval_structure(&vec![mrna_full, rrna], &dots, DEFAULT_TEMP, dangles) as f64;
    let total_energy_withspacing = total_energy + dg_spacing_final;

    let structure = MRNARRNAStructure {
        mrna: mrna_full.to_string(),
        dg_mrna_rrna: total_energy,
        dg_mrna_rrna_withspacing: total_energy_withspacing,
        dg_spacing: dg_spacing_final,
        bp_x: total_bp_x,
        bp_y: total_bp_y,
    };

    Some((total_energy_withspacing, structure, spacing_value))
}

/// Calculates the dG_standby penalty.
/// Returns dG_standby_site (≤ 0.0).
pub fn calc_dg_standby_site(
    structure: &MRNARRNAStructure,
    rrna: &str,
    dangles: &DanglesSetting,
    _constraints: Option<&str>,
) -> f64 {
    let mrna = &structure.mrna;
    let bp_x = &structure.bp_x;
    let bp_y = &structure.bp_y;
    let energy_before = structure.dg_mrna_rrna;
    let mrna_len = mrna.len();

    // Find most 5' mRNA nt bound to rRNA (first pair where nt_x ≤ mrna_len and nt_y > mrna_len)
    let mut most_5p_mrna: usize = 0;
    for (&nt_x, &nt_y) in bp_x.iter().zip(bp_y.iter()) {
        if nt_x <= mrna_len && nt_y > mrna_len {
            most_5p_mrna = nt_x;
            break;
        }
    }

    if most_5p_mrna == 0 {
        return 0.0;
    }

    // Extract base pairs 3' of most_5p_mrna (keep them fixed)
    let bp_x_3p: Vec<usize> = bp_x
        .iter()
        .zip(bp_y.iter())
        .filter(|(&x, _)| x >= most_5p_mrna)
        .map(|(&x, _)| x)
        .collect();
    let bp_y_3p: Vec<usize> = bp_x
        .iter()
        .zip(bp_y.iter())
        .filter(|(&x, _)| x >= most_5p_mrna)
        .map(|(_, &y)| y)
        .collect();

    // mRNA subsequence before the standby site (standby_site_length = 4 nt before binding)
    let subseq_end = most_5p_mrna.saturating_sub(STANDBY_SITE_LEN as usize).saturating_sub(1);
    let mrna_subsequence = if subseq_end > 0 { &mrna[0..subseq_end] } else { "" };

    // Fold the pre-standby subsequence
    let (bp_x_5p, bp_y_5p) = if !mrna_subsequence.is_empty() {
        let seqs = vec![mrna_subsequence];
        match mfe(&seqs, "", DEFAULT_TEMP, dangles) {
            Ok(fold) => (fold.get_bp_x().clone(), fold.get_bp_y().clone()),
            Err(_) => (vec![], vec![]),
        }
    } else {
        (vec![], vec![])
    };

    // Combine: pre-standby pairs + rRNA-binding/downstream pairs
    let mut bp_x_after: Vec<usize> = bp_x_5p;
    let mut bp_y_after: Vec<usize> = bp_y_5p;
    bp_x_after.extend_from_slice(&bp_x_3p);
    bp_y_after.extend_from_slice(&bp_y_3p);

    // Calculate energy of modified structure
    let dots = coordinates_to_dots(&vec![mrna.as_str(), rrna], &bp_x_after, &bp_y_after);
    let energy_after = eval_structure(
        &vec![mrna.as_str(), rrna],
        &dots,
        DEFAULT_TEMP,
        dangles,
    ) as f64;

    let dg_standby = energy_before - energy_after;
    if dg_standby > 0.0 { 0.0 } else { dg_standby }
}
