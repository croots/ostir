
use std::f64::consts::E;

pub const BETA: f64 = 0.40002512;
pub fn rt_eff() -> f64 {return 1.0/BETA}
const LOGK: f64 = 7.279194329;
pub fn k() -> f64 {return E.powf(LOGK)}
pub const RNAMODEL: &str = "rna2004";
pub const AUTO_DANGLES: bool = true;
pub const DEFAULT_DANGLES: &str = "all";
pub const DEFAULT_TEMP: f32 = 37.0;
pub const OPTIMAL_SPACING: usize = 5;
pub const CUTOFF: usize = 35;

// From OSTIR calibration using Salis2009 data. See calibration directory for procedure
pub const DG_SPACING_PUSH: [f64; 4] = [17.20965071, 3.46341492, 1.790848365, 3.0];
pub const DG_SPACING_PULL: [f64; 3] = [0.06422042, 0.275640836, 0.0];

pub const STANDBY_SITE_LEN: i32 = 4;  // Number of nt before SD sequence that must be unpaired for ribosome binding
pub const FOOTPRINT: i32 = 1000;
pub const ENERGY_CUTOFF: f32 = 3.0;
pub const VERBOSE: bool = false;