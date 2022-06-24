/// Computation routine of DOS
use crate::parser::bands::{Bands, KPointEnergiesVec};

use super::DOS;

use itertools_num::linspace;
use std::f64::consts::{PI, SQRT_2};
const SMEARING_WIDTH: f64 = 0.2;
const HATREE_TO_EV: f64 = 27.211396641308;

/**
Compute the total Density of States (DOS).
# Arguments:
  * `band_data`: The `Bands` struct containing data and information read from `<seed>.bands`.
  * `energy_range`: Range of energy to compute the DOS.
# Returns:
  * `(Vec<f64>, Option<Vec<f64>>)`: Spin 1 is guaranteed to return, while Spin 2 is optional
  given by the existence of the data.
# Notes:
  * The function will panic if `num_spins` are larger than 2 -- which usually would not happen.
*/
pub fn total_dos(band_data: &Bands, energy_range: &Vec<f64>) -> (Vec<f64>, Option<Vec<f64>>) {
    if band_data.num_spins() == 1 {
        (total_dos_spin(&band_data, energy_range, 1 as u8), None)
    } else {
        assert_eq!(
            2,
            band_data.num_spins(),
            "Spin cannot be > 2! Current value: {}",
            band_data.num_spins()
        );
        (
            total_dos_spin(&band_data, energy_range, 1),
            Some(total_dos_spin(&band_data, energy_range, 2)),
        )
    }
}

/**
Helper function to compute the total Density of States (DOS) at given spin.
*/
fn total_dos_spin(band_data: &Bands, energy_range: &Vec<f64>, spin: u8) -> Vec<f64> {
    let eigen_energies_vec = band_data.kpt_eigen_energies_array().eigen_energies();
    let e_fermi = band_data.e_fermi()[spin as usize - 1];
    let current_spin_shifted_eigens: Vec<Vec<f64>> = eigen_energies_vec
        .iter()
        .map(|entry| -> Vec<f64> {
            let eigen_at_spin = entry
                .eigen_at_spin(spin)
                .expect(&format!("No eigenvalues at given spin {}", spin));
            shift_eigen_by_e_fermi(&eigen_at_spin.to_vec(), e_fermi)
        })
        .collect();
    let k_point_weight_array = band_data.kpt_eigen_energies_array().kpt_weight();
    let mut dos = vec![0.0; energy_range.len()];
    current_spin_shifted_eigens
        .iter()
        // to get idx of eigen values, and obtain idx-th k-point weight
        .enumerate()
        .for_each(|(i, entry)| {
            entry // each eigen value list
                .iter()
                /*
                generate delta energies with each value in the list,
                weighted by their k-point weight
                */
                .map(|de| -> Vec<f64> {
                    let this_de_contribution = e_delta(&energy_range, *de, k_point_weight_array[i]);
                    this_de_contribution
                })
                .collect::<Vec<Vec<f64>>>() // Collect to Vector of contribution
                // Iterate over contributions
                .iter()
                // For each contribution
                .for_each(|list: &Vec<f64>| {
                    // Iterate over the contribution values
                    list.iter()
                        // Turn to enumerate to align with dos vector
                        .enumerate()
                        .for_each(|(j, contrib)| {
                            // Add over the corresponding position
                            dos[j] += contrib;
                        })
                })
            /*
            Collection of each eigenvalue's contribution to DOS
            Length = number of eigenvalues in list
            */
        });
    dos
}
/**
Shift the eigenvalues by fermi energy
*/
fn shift_eigen_by_e_fermi(eigen_energies_entry: &Vec<f64>, e_fermi: f64) -> Vec<f64> {
    eigen_energies_entry
        .iter()
        .map(|e| -> f64 { (e - e_fermi) * HATREE_TO_EV })
        .collect()
}

/**
Calculate energy range from min, max and number
of points per eV.
*/
pub fn compute_energy_range(e_min: f64, e_max: f64, num_points: u32) -> Vec<f64> {
    linspace::<f64>(e_min, e_max, num_points as usize)
        .into_iter()
        .collect()
}

/**
Gaussian smearing
e from given range, de from k-point eigenvalue
*/
fn gaussian_smearing(e: f64, de: f64) -> f64 {
    /*
    The smearing here is the standard deviation
    Same as implemented in CASTEP
    */
    let smear = SMEARING_WIDTH;
    /*
    g(x) = 1 / (σ*sqrt(2PI)) * exp(-1/2 * (e - de)^2/σ^2)
    */
    (-((e - de) / (smear)).powi(2) / 2.0).exp() / ((2.0 * PI).sqrt() * smear)
}
/**
Generate energies centered at the given k-point
eigenvalue of energy, broadened by Gaussian smearing,
weighted by given k-point weight value.
*/
fn e_delta(energies: &Vec<f64>, de: f64, k_point_weight: f64) -> Vec<f64> {
    energies
        .iter()
        .map(|e| -> f64 { gaussian_smearing(*e, de) * k_point_weight })
        .collect()
}
