use std::fs;

use crate::parser::pdos_weights::{EigenWeightPerOrb, PDOSWeight, UnitOrbitalWeight};

fn read_pdos_weights(pdos_weights_file: &str) -> PDOSWeight {
    let file = fs::read(pdos_weights_file).expect("Failed to open pdos_weights file!");
    let (_, pdos_weights) = PDOSWeight::parse(&file).unwrap();
    pdos_weights
}

/**
Get the eigenweight at given `orbital_id` (0-indexed) and spin
*/
fn get_orbital_eigenweights(
    weights: &PDOSWeight,
    orbital_id: u32,
    spin: u32,
) -> &UnitOrbitalWeight {
    let eigen_weights = weights.eigen_weights_at_orbital(orbital_id);
    match spin {
        1 => eigen_weights.upspin(),
        2 => eigen_weights.downspin(),
        _ => panic!("Wrong spin number :{spin}"),
    }
}
