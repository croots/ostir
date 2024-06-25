use crate::types::{DanglesSetting, FoldResult};
use librna_sys::{
    vrna_eval_structure, vrna_fold_compound, vrna_fold_compound_t, vrna_md_set_default, vrna_md_t,
    vrna_mfe, vrna_subopt_cb,
};
use std::error::Error;
use std::ffi::{c_char, c_double, c_float, c_void, CStr};
use std::mem::MaybeUninit;

// Fold compound struct with safeguards --------
struct FoldCompound {
    c: *mut vrna_fold_compound_t,
}

impl FoldCompound {
    fn new(sequences: &Vec<&str>, _constraints: &str, dangles: &DanglesSetting, temp: f32) -> Self {
        let sequence;
        if sequences.len() > 1 {
            sequence = sequences.join("&").replace("T", "U").to_uppercase()
        } else {
            sequence = sequences[0].replace("T", "U").to_uppercase()
        }

        unsafe {
            let mut md = MaybeUninit::<vrna_md_t>::uninit();
            let md_ptr = md.as_mut_ptr();
            vrna_md_set_default(md_ptr);
            let mut initialized_md = md.assume_init();

            initialized_md.temperature = temp as c_double;

            let _dangles_int = dangles.as_int();
            match dangles.as_int() {
                Ok(t) => initialized_md.dangles = t,
                Err(_e) => {}
            }

            let c = vrna_fold_compound(sequence.as_ptr() as *const i8, &initialized_md, 1 as u32);
            // TODO: Add constraints

            FoldCompound { c }
        }
    }
}

#[cfg(test)]
mod tests {
    use pyo3::prelude::*;
    use pyo3::types::IntoPyDict;

    #[test]
    fn python_test() -> PyResult<()> {
        pyo3::prepare_freethreaded_python();
        Python::with_gil(|py| {
            println!("Testing python");
            let sys = py.import_bound("sys")?;
            let path = sys.getattr("path")?;
            let sys = py.import_bound("ostir")?;

            let version: String = sys.getattr("version")?.extract()?;
            let binary: String = sys.getattr("executable")?.extract()?;

            let locals = [("os", py.import_bound("os")?)].into_py_dict_bound(py);
            let code = "os.getenv('USER') or os.getenv('USERNAME') or 'Unknown'";
            let user: String = py.eval_bound(code, None, Some(&locals))?.extract()?;

            println!("Hello {}, I'm Python {} at {}", user, version, binary);
            Ok(())
        })
    }
}

// MFE ------------

pub fn mfe<'a>(
    sequences: &'a Vec<&'a str>,
    constraints: &'_ str,
    temp: f32,
    dangles: &'_ DanglesSetting,
) -> Result<FoldResult<'a>, Box<dyn Error>> {
    // @TODO: Add constraints option

    let dot_vec = vec![0; sequences.join("&").len() + 1];
    let dot_ptr = dot_vec.as_ptr() as *mut i8;
    let fold_compound = FoldCompound::new(sequences, constraints, dangles, temp);
    let result;
    unsafe {
        result = vrna_mfe(fold_compound.c, dot_ptr);
    }
    let dot_string = std::str::from_utf8(&dot_vec).expect("TODO: Handle invalid UTF-8");
    let coordinates = dots_to_coordinates(dot_string);

    return Ok(FoldResult::create(
        Some(sequences),
        result,
        dot_string.to_string(),
        coordinates.0,
        coordinates.1,
    ));
}

// Subopt ----------------
unsafe extern "C" fn subopt_cb_fun(x: *const c_char, y: c_float, z: *mut c_void) {
    let optional_char = x.as_ref(); // Catch a returned null pointer
    match optional_char {
        Some(_t) => {}
        None => return,
    }

    let data: &mut Vec<FoldResult> = unsafe { &mut *(z as *mut Vec<FoldResult>) };

    let placeholder = CStr::from_ptr(x);
    let dots_string = placeholder.to_str().unwrap();
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
    // error if temp < 0
    // error if dangles no 'all', 'some', or 'none'
    // energy_gap in kcal/mol

    let fold_compound = FoldCompound::new(sequences, constraints, dangles, temp);
    let _vienna_temperature = (temp * 100 as f32).round() as i32;

    let mut resultholder: Vec<FoldResult> = vec![];
    let holder_ptr: *mut c_void = &mut resultholder as *mut _ as *mut c_void;

    let hybridization_penalty = 2.481 as f32;

    let _energy_gap_adjusted = (energy_gap + hybridization_penalty) * 100.0;
    let energy_gap_rounded: i32 = (energy_gap as f32).round() as i32;

    unsafe {
        vrna_subopt_cb(
            fold_compound.c,
            energy_gap_rounded,
            Some(subopt_cb_fun as _),
            holder_ptr,
        );
    }

    resultholder.sort_by(|b, a| b.get_d_g().partial_cmp(&a.get_d_g()).unwrap());

    return resultholder;
}

// Evaluate Fold for Energy ----------------
pub fn eval_structure(
    sequences: &Vec<&str>,
    dots: &str,
    temp: f32,
    dangles: &DanglesSetting,
) -> f32 {
    let adj_dots = dots.replace("&", "");
    let fold_compound = FoldCompound::new(sequences, "", dangles, temp);

    let energy: c_float;
    unsafe {
        energy = vrna_eval_structure(fold_compound.c, adj_dots.as_ptr() as *const i8);
    }

    return energy;
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
    let _unpaired_x_index: Vec<usize> = vec![];
    let _unpaired_x_pos: Vec<usize> = vec![];
    let mut bp_y: Vec<usize> = vec![];

    let _i = 1; // Dot positions are 1 indexed
    let _x_counter = 0;
    let _strand_count = 0;

    let mut _last_x_pos: usize;
    let mut _last_x_index: usize;

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
                let nt_x = last_nt_x_list.pop().unwrap(); // nt_x is list of "(" except last entry
                let nt_x_pos = bp_x
                    .iter()
                    .position(|&x| x == nt_x.try_into().unwrap())
                    .unwrap();
                bp_y[nt_x_pos] = (pos - num_strands + 2).try_into().unwrap();
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
