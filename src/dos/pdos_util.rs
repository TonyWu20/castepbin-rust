#![allow(soft_unstable)]
use std::{fs, ops::Add};

use crate::parser::pdos_weights::{EigenWeightPerOrb, PDOSWeight, UnitOrbitalWeight};

/**
Struct representing a command to generate PDOS.
This intend to generate merged weights data for 1 line in PDOS plot.
Example: To compute PDOS at Pt-1, p orbitals,
PDOSCommand {
    vec![Pt],
    vec![1],
    vec![p]
}
To compute PDOS of Pt1 and Pt2 and their p,d orbitals
PDOSCommand {
    vec![Pt,Pt, Pt, Pt]
    vec![1, 1, 2, 2]
    vec![p, d, p, d]
}
*/
pub struct PDOSCommand {
    species: Vec<String>,
    ion_in_species: Vec<u32>,
    am_channel: Vec<String>,
}

impl PDOSCommand {
    pub fn new(species: Vec<String>, ion_in_species: Vec<u32>, am_channel: Vec<String>) -> Self {
        Self {
            species,
            ion_in_species,
            am_channel,
        }
    }
}

pub struct PDOSComputeConfig {
    seed: String,
    commands: Vec<PDOSCommand>,
}

impl PDOSComputeConfig {
    pub fn new(seed: String, commands: Vec<PDOSCommand>) -> Self {
        Self { seed, commands }
    }
    // Parse one command, to generate a set of weights for PDOS calculation.
    pub fn parse_command(&self, command: &PDOSCommand) -> UnitOrbitalWeight {
        todo!();
    }
}

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
) -> &[UnitOrbitalWeight] {
    let eigen_weights = weights.eigen_weights_at_orbital(orbital_id);
    match spin {
        1 => eigen_weights.upspin(),
        2 => eigen_weights.downspin().unwrap().as_ref(),
        _ => panic!("Wrong spin number :{spin}"),
    }
}
fn merged_weights(orbital_weights: &[&UnitOrbitalWeight]) -> Vec<f64> {
    let number_of_eigens = orbital_weights[0].weights().len();
    orbital_weights
        .iter()
        .fold(vec![0.0; number_of_eigens], |a, b| {
            a.iter()
                .zip(b.weights().iter())
                .map(|(a, b)| a + b)
                .collect()
        })
}

#[cfg(test)]
mod test {
    use std::ops::Add;

    use nalgebra::{vector, DVector};

    use crate::parser::pdos_weights::UnitOrbitalWeight;

    use super::{get_orbital_eigenweights, read_pdos_weights};
}
