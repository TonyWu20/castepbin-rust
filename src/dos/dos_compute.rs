/// Computation routine of DOS
use crate::{
    parser::bands::{Bands, KPointEnergiesVec},
    util::ElementWiseAdd,
};

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
pub fn total_dos(band_data: &Bands, energy_range: &[f64]) -> (Vec<f64>, Option<Vec<f64>>) {
    if band_data.num_spins() == 1 {
        (total_dos_spin(band_data, energy_range, 1_u8), None)
    } else {
        assert_eq!(
            2,
            band_data.num_spins(),
            "Spin cannot be > 2! Current value: {}",
            band_data.num_spins()
        );
        (
            total_dos_spin(band_data, energy_range, 1),
            Some(total_dos_spin(band_data, energy_range, 2)),
        )
    }
}

/**
Helper function to compute the total Density of States (DOS) at given spin.
*/
fn total_dos_spin(band_data: &Bands, energy_range: &[f64], spin: u8) -> Vec<f64> {
    /*
    vec of energy eigenvalues for every k-point
    length = number of k-points
    */
    let eigen_energies_vec = band_data.kpt_eigen_energies_array().eigen_energies();
    let e_fermi = band_data.e_fermi()[spin as usize - 1];
    // Shift each k-point's eigen energies by fermi energy
    let current_spin_shifted_eigens: Vec<Vec<f64>> = eigen_energies_vec
        .iter()
        .map(|entry| -> Vec<f64> {
            let eigen_at_spin = entry
                .eigen_at_spin(spin)
                .unwrap_or_else(|| panic!("No eigenvalues at given spin {}", spin));
            shift_eigen_by_e_fermi(eigen_at_spin, e_fermi)
        })
        .collect();
    let k_point_weight_array = band_data.kpt_eigen_energies_array().kpt_weight();
    current_spin_shifted_eigens
        .iter()
        // Pairs with k-point weights
        .zip(k_point_weight_array.iter())
        .map(|(de_list_per_k, kpt_weight)| -> Vec<f64> {
            // E_delta over given energy range
            // Total DOS so orbital_weights are `None`
            e_deltas_per_kpt(energy_range, de_list_per_k, *kpt_weight, None)
        })
        // Fold vector of vectors into one vector
        .reduce(|prev, next| -> Vec<f64> { prev.add(&next) })
        .unwrap()
}
/**
Shift the eigenvalues by fermi energy and convert Hatree to eV
*/
pub fn shift_eigen_by_e_fermi(eigen_energies_entry: &[f64], e_fermi: f64) -> Vec<f64> {
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
fn e_delta(energies: &[f64], de: f64, k_point_weight: f64, orbital_weight: f64) -> Vec<f64> {
    energies
        .iter()
        .map(|e| -> f64 { gaussian_smearing(*e, de) * k_point_weight * orbital_weight })
        .collect()
}

pub fn e_deltas_per_kpt(
    energies: &[f64],
    de_array: &[f64],
    k_point_weight: f64,
    orbital_weight_array: Option<&[f64]>,
) -> Vec<f64> {
    // orbital_weight_array has the same length as de_array
    match orbital_weight_array {
        Some(orbital_weights) => de_array
            .iter()
            .zip(orbital_weights.iter())
            .map(|(de, orb_weight)| e_delta(energies, *de, k_point_weight, *orb_weight))
            // Now we have n Vec<f64> of e_delta
            .reduce(|prev, next| -> Vec<f64> {
                // merge Vec<_> by zip and map
                prev.add(&next)
            })
            .unwrap(),
        None => de_array
            .iter()
            .map(|de| e_delta(energies, *de, k_point_weight, 1.0))
            // Now we have n Vec<f64> of e_delta
            .reduce(|prev, next| -> Vec<f64> {
                // merge Vec<_> by zip and map
                prev.add(&next)
            })
            .unwrap(),
    }
}
