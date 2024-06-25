pub struct DanglesSetting {
    setting: String,
}

impl DanglesSetting {
    pub fn new(setting_str: &str) -> Result<DanglesSetting, &'static str> {
        if setting_str == "all"
            || setting_str == "none"
            || setting_str == "some"
            || setting_str == "default"
        {
            Ok(DanglesSetting {
                setting: setting_str.to_string(),
            })
        } else {
            Err("Invalid dangle setting")
        }
    }

    pub fn as_int(&self) -> Result<i32, &'static str> {
        match self.setting.as_str() {
            "all" => Ok(2),
            "some" => Ok(1),
            "none" => Ok(0),
            _ => Err("Must pull default dangle setting"),
        }
    }
}

// Represents a fold result from ViennaRNA
pub enum FoldResult<'a> {
    Unknown(UnknownFoldResult<'a>),
    Mono(MonoFoldResult<'a>),
    Co(CoFoldResult<'a>),
}

pub struct UnknownFoldResult<'a> {
    pub seqs: Option<&'a Vec<&'a str>>,
    pub d_g: f32,
    pub dots: String,
    pub bp_x: Vec<usize>,
    pub bp_y: Vec<usize>,
}
pub struct MonoFoldResult<'a> {
    pub seqs: &'a str,
    pub d_g: f32,
    pub dots: String,
    pub bp_x: Vec<usize>,
    pub bp_y: Vec<usize>,
}
pub struct CoFoldResult<'a> {
    pub seqs: (&'a str, &'a str),
    pub d_g: f32,
    pub dots: String,
    pub bp_x: Vec<usize>,
    pub bp_y: Vec<usize>,
}

impl<'a> FoldResult<'_> {
    pub fn get_d_g(&self) -> &f32 {
        // Some code here
        match self {
            FoldResult::Mono(x) => return &x.d_g,
            FoldResult::Co(x) => return &x.d_g,
            FoldResult::Unknown(x) => return &x.d_g,
        }
    }

    pub fn get_dots(&self) -> &String {
        // Some code here
        match self {
            FoldResult::Mono(x) => return &x.dots,
            FoldResult::Co(x) => return &x.dots,
            FoldResult::Unknown(x) => return &x.dots,
        }
    }

    pub fn get_bp_x(&self) -> &Vec<usize> {
        // Some code here
        match self {
            FoldResult::Mono(x) => return &x.bp_x,
            FoldResult::Co(x) => return &x.bp_x,
            FoldResult::Unknown(x) => return &x.bp_x,
        }
    }

    pub fn get_bp_y(&self) -> &Vec<usize> {
        // Some code here
        match self {
            FoldResult::Mono(x) => return &x.bp_y,
            FoldResult::Co(x) => return &x.bp_y,
            FoldResult::Unknown(x) => return &x.bp_y,
        }
    }

    pub fn create(
        sequences: Option<&'a Vec<&'a str>>,
        d_g: f32,
        dots: String,
        bp_x: Vec<usize>,
        bp_y: Vec<usize>,
    ) -> FoldResult {
        if sequences.is_none() {
            return FoldResult::Unknown(UnknownFoldResult {
                seqs: None,
                d_g,
                dots,
                bp_x,
                bp_y,
            });
        }
        let inside_sequences = sequences.unwrap();
        match inside_sequences.len() {
            1 => {
                return FoldResult::Mono(MonoFoldResult {
                    seqs: inside_sequences[0],
                    d_g,
                    dots,
                    bp_x,
                    bp_y,
                })
            }
            2 => {
                return FoldResult::Co(CoFoldResult {
                    seqs: (inside_sequences[0], inside_sequences[1]),
                    d_g,
                    dots,
                    bp_x,
                    bp_y,
                })
            }
            _ => {
                return FoldResult::Unknown(UnknownFoldResult {
                    seqs: Some(inside_sequences),
                    d_g,
                    dots,
                    bp_x,
                    bp_y,
                })
            }
        }
    }
}
