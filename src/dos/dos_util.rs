use std::fs;

use crate::parser::bands::Bands;

use super::{
    dos_compute::{compute_energy_range, total_dos, HATREE_TO_EV},
    TDOS,
};

pub fn calculate_total_dos_from_bands(
    bands_file: &str,
    energy_range_limit: Option<(f64, f64)>,
    e_num_per_ev: Option<u32>,
) -> TDOS {
    let bands_file_string = fs::read_to_string(bands_file).expect("Failed to open .bands file!");
    let bands_data = Bands::parse(&bands_file_string)
        .expect("Failed to parse .bands!")
        .1;
    let (e_min, e_max) = energy_range_limit.unwrap_or((
        (bands_data.min_eigen_energy() - bands_data.e_fermi()[0]) * HATREE_TO_EV - 2.00,
        (bands_data.max_eigen_energy() - bands_data.e_fermi()[0]) * HATREE_TO_EV + 2.00,
    ));
    println!("{}, {}", e_min, e_max);
    let num_per_ev = e_num_per_ev.unwrap_or(100);
    let num_points = ((e_max - e_min) * num_per_ev as f64) as u32;
    let energy_range = compute_energy_range(e_min, e_max, num_points);
    let total_dos = total_dos(&bands_data, &energy_range);
    TDOS::new(bands_data.num_spins(), num_points, energy_range, total_dos)
}
