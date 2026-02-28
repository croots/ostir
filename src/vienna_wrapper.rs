use crate::types::{DanglesSetting, FoldResult};
use librna_sys::{
    vrna_eval_structure, vrna_fold_compound, vrna_md_set_default, vrna_md_t, vrna_mfe,
    vrna_params_load_defaults, vrna_subopt_cb,
};
use std::error::Error;
use std::ffi::{c_char, c_double, c_float, c_void, CStr, CString};
use std::sync::Mutex;

// ViennaRNA thread-safety notes:
//
// vrna_mfe and vrna_eval_structure are safe to call concurrently when each
// caller creates its own vrna_fold_compound_t — they do not write to shared
// global state under normal operation.
//
// vrna_subopt_cb is NOT thread-safe: it writes to the PUBLIC global array
// `density_of_states[MAXDOS+1]` on every result, and in ViennaRNA ≤ 2.6.4 it
// also has a latent bug where floating-point rounding can produce a negative
// index (e = -1), writing one word before the array and corrupting the adjacent
// `last_param_file` global pointer.  A subsequent call to vrna_fold_compound
// will then crash inside strncpy trying to copy from that invalid address.
//
// Hardening strategy:
//  1. Hold VIENNA_LOCK for the entire subopt call so density_of_states writes
//     and the potential corruption of last_param_file are serialised against
//     any concurrent mfe / eval_structure caller that would read last_param_file
//     inside vrna_fold_compound.
//  2. While still holding the lock, call vrna_params_load_defaults() after
//     vrna_subopt_cb returns.  This resets last_param_file to a valid heap
//     string ("RNA - Turner 2004") so that the corrupted value never escapes
//     the critical section.
//
// mfe and eval_structure also acquire the lock so they are protected against a
// concurrent subopt that might corrupt last_param_file between their lock
// acquisition and the vrna_fold_compound call inside them.
static VIENNA_LOCK: Mutex<()> = Mutex::new(());

fn make_md(dangles: &DanglesSetting, temp: f32) -> Box<vrna_md_t> {
    unsafe {
        // Heap-allocate so the address is stable after vrna_md_set_default.
        let mut md = Box::new(std::mem::zeroed::<vrna_md_t>());
        vrna_md_set_default(md.as_mut());
        md.temperature = temp as c_double;
        md.noLP = 1;
        if let Ok(t) = dangles.as_int() {
            md.dangles = t;
        }
        md
    }
}

// MFE ------------

pub fn mfe<'a>(
    sequences: &'a Vec<&'a str>,
    constraints: &'_ str,
    temp: f32,
    dangles: &'_ DanglesSetting,
) -> Result<FoldResult<'a>, Box<dyn Error>> {
    // Guard: ViennaRNA requires at least 2 nucleotides for meaningful folding.
    // For sequences that are too short, return a trivial result (all unpaired, energy 0).
    let total_len: usize = sequences.iter().map(|s| s.len()).sum();
    if total_len < 2 {
        let dots: String = sequences
            .iter()
            .enumerate()
            .map(|(i, s)| {
                if i > 0 { format!("&{}", ".".repeat(s.len())) } else { ".".repeat(s.len()) }
            })
            .collect();
        let (bp_x, bp_y) = dots_to_coordinates(&dots);
        return Ok(FoldResult::create(Some(sequences), 0.0, dots, bp_x, bp_y));
    }

    let _lock = VIENNA_LOCK.lock().unwrap();
    let joined = sequences.join("&").replace("T", "U").to_uppercase();
    let c_sequence = CString::new(joined.clone())?;
    let md = make_md(dangles, temp);
    let dot_vec = vec![0u8; joined.len() + 1];
    let dot_ptr = dot_vec.as_ptr() as *mut i8;
    let (result, dot_string) = unsafe {
        let fc = vrna_fold_compound(c_sequence.as_ptr(), md.as_ref(), 1);
        let result = vrna_mfe(fc, dot_ptr);
        let dot_cstr = CStr::from_ptr(dot_ptr as *const c_char);
        let dot_string = dot_cstr.to_str().expect("ViennaRNA returned invalid UTF-8").to_owned();
        (result, dot_string)
    };
    let coordinates = dots_to_coordinates(&dot_string);
    let result_rounded = (result as f64 * 100.0).round() as f32 / 100.0;
    Ok(FoldResult::create(
        Some(sequences),
        result_rounded,
        dot_string,
        coordinates.0,
        coordinates.1,
    ))
}

// Subopt ----------------
unsafe extern "C" fn subopt_cb_fun(x: *const c_char, y: c_float, z: *mut c_void) {
    // ViennaRNA calls the callback with NULL to signal the end of subopt results
    if x.is_null() {
        return;
    }

    let data: &mut Vec<FoldResult> = unsafe { &mut *(z as *mut Vec<FoldResult>) };

    let placeholder = unsafe { CStr::from_ptr(x) };
    let dots_string = match placeholder.to_str() {
        Ok(s) => s,
        Err(_) => return,
    };
    let coordinates = dots_to_coordinates(dots_string);
    let result = FoldResult::create(
        None,
        y,
        dots_string.to_string(),
        coordinates.0,
        coordinates.1,
    );

    data.push(result);
}

pub fn subopt<'a>(
    sequences: &'a Vec<&'a str>,
    constraints: &'_ str,
    energy_gap: f32,
    temp: f32,
    dangles: &'_ DanglesSetting,
) -> Vec<FoldResult<'a>> {
    let _lock = VIENNA_LOCK.lock().unwrap();
    let joined = sequences.join("&").replace("T", "U").to_uppercase();
    let c_sequence = CString::new(joined).expect("Sequence contains null byte");
    let md = make_md(dangles, temp);

    let mut resultholder: Vec<FoldResult> = vec![];
    let holder_ptr: *mut c_void = &mut resultholder as *mut _ as *mut c_void;
    let energy_gap_rounded: i32 = ((energy_gap + 2.481) * 100.0).round() as i32;

    unsafe {
        let fc = vrna_fold_compound(c_sequence.as_ptr(), md.as_ref(), 1);
        vrna_subopt_cb(
            fc,
            energy_gap_rounded,
            Some(subopt_cb_fun as _),
            holder_ptr,
        );
        // Harden against the ViennaRNA ≤ 2.6.4 density_of_states[-1] bug:
        // if the OOB write corrupted last_param_file, reset it to a valid
        // pointer by reloading the default energy parameters.  This must
        // happen while the lock is still held so that no concurrent
        // vrna_fold_compound call can observe the corrupted value.
        vrna_params_load_defaults();
    }

    resultholder.sort_by(|b, a| b.get_d_g().partial_cmp(&a.get_d_g()).unwrap());

    // When folding multiple sequences, only keep results where rRNA strand has base pairs
    if sequences.len() > 1 {
        resultholder.retain(|result| {
            let dots = result.get_dots();
            if let Some(rrna_part) = dots.split('&').nth(1) {
                rrna_part.contains('(') || rrna_part.contains(')')
            } else {
                false
            }
        });
    }

    resultholder
}

// Evaluate Fold for Energy ----------------
pub fn eval_structure(
    sequences: &Vec<&str>,
    dots: &str,
    temp: f32,
    dangles: &DanglesSetting,
) -> f64 {
    let _lock = VIENNA_LOCK.lock().unwrap();
    let joined = sequences.join("&").replace("T", "U").to_uppercase();
    let c_sequence = CString::new(joined).expect("Sequence contains null byte");
    let adj_dots = dots.replace("&", "");
    let c_dots = CString::new(adj_dots).expect("Structure string contains null byte");
    let md = make_md(dangles, temp);

    let energy: c_float = unsafe {
        let fc = vrna_fold_compound(c_sequence.as_ptr(), md.as_ref(), 1);
        vrna_eval_structure(fc, c_dots.as_ptr())
    };

    // Round to 2 decimal places to match Python's ostir ViennaRNA.energy() behavior
    (energy as f64 * 100.0).round() / 100.0
}

// Utilities ----------------

pub fn coordinates_to_dots(strands: &Vec<&str>, bp_x: &Vec<usize>, bp_y: &Vec<usize>) -> String {
    let bp_x: Vec<usize> = bp_x.iter().map(|&pos| pos - 1).collect(); // Shift so that 1st position is 0
    let bp_y: Vec<usize> = bp_y.iter().map(|&pos| pos - 1).collect(); // Shift so that 1st position is 0

    let mut bracket_notation = Vec::new();
    let mut counter = 0;

    for (strand_number, &seq) in strands.iter().enumerate() {
        let seq_len = seq.len();
        if strand_number > 0 {
            bracket_notation.push('&');
        }
        for pos in counter..(seq_len + counter) {
            if bp_x.contains(&pos) {
                bracket_notation.push('(');
            } else if bp_y.contains(&pos) {
                bracket_notation.push(')');
            } else {
                bracket_notation.push('.');
            }
        }
        counter += seq_len;
    }

    bracket_notation.iter().collect::<String>()
}

pub fn dots_to_coordinates(dots_string: &str) -> (Vec<usize>, Vec<usize>) {
    let mut bp_x: Vec<usize> = vec![];
    let mut bp_y: Vec<usize> = vec![];

    for _ in 0..dots_string.matches(")").count() {
        bp_y.push(0); // Placeholder value to be replaced later
    }

    let mut last_nt_x_list: Vec<usize> = Vec::new();
    let mut num_strands = 0;

    for (pos, letter) in dots_string.chars().enumerate() {
        match letter {
            '.' => {}
            '(' => {
                bp_x.push((pos - num_strands).try_into().unwrap());
                last_nt_x_list.push(pos - num_strands);
            }
            ')' => {
                let nt_x = last_nt_x_list.pop().unwrap();
                let nt_x_pos = bp_x
                    .iter()
                    .position(|&x| x == nt_x.try_into().unwrap())
                    .unwrap();
                bp_y[nt_x_pos] = (pos - num_strands).try_into().unwrap();
            }
            '&' => {
                num_strands += 1;
            }
            _ => {
                println!("Error! Invalid character in bracket notation.");
            }
        }
    }

    if !last_nt_x_list.is_empty() {
        println!("Error! Leftover unpaired nucleotides when converting from bracket notation to numbered base pairs.");
    }

    if bp_y.len() > 1 {
        bp_x.iter_mut().for_each(|x| *x += 1); // Shift so that 1st position is 1
        bp_y.iter_mut().for_each(|y| *y += 1); // Shift so that 1st position is 1
    }

    return (bp_x, bp_y);
}
