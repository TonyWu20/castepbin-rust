/// Computation routine of DOS
use crate::{
    parser::bands::{Bands, KPointEnergiesVec},
    util::ElementWiseAdd,
};

use super::DOS;

use itertools_num::linspace;
use nalgebra::{
    Const, DMatrix, Dim, Dynamic, Matrix, Matrix1xX, Matrix2x1, MatrixXx1, MatrixXx2, SMatrix,
    VecStorage, Vector2, U32,
};
use ndarray::{Array, Array1, Array2, Axis, ShapeBuilder};
use rayon::iter::{
    IndexedParallelIterator, IntoParallelRefIterator, ParallelBridge, ParallelIterator,
};
use std::{
    f64::consts::{PI, SQRT_2},
    ops::{Mul, Sub},
};
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
// pub fn total_dos(band_data: &Bands, energy_range: &[f64]) -> (Vec<f64>, Option<Vec<f64>>) {
//     if band_data.num_spins() == 1 {
//         (total_dos_spin(band_data, energy_range, 1_u8), None)
//     } else {
//         assert_eq!(
//             2,
//             band_data.num_spins(),
//             "Spin cannot be > 2! Current value: {}",
//             band_data.num_spins()
//         );
//         (
//             total_dos_spin(band_data, energy_range, 1),
//             Some(total_dos_spin(band_data, energy_range, 2)),
//         )
//     }
// }

/**
Helper function to compute the total Density of States (DOS) at given spin.
*/
// fn _total_dos_spin(band_data: &Bands, energy_range: &[f64], spin: u8) -> Vec<f64> {
//     /*
//     vec of energy eigenvalues for every k-point
//     length = number of k-points
//     */
//     let eigen_energies_vec = band_data.kpt_eigen_energies_array().eigen_energies();
//     let e_fermi = band_data.e_fermi()[spin as usize - 1];
//     // Shift each k-point's eigen energies by fermi energy
//     // Generates a 2-dimensional Vec<Vec<f64>>,
//     // eigen_energies_vec[nth_kpoint_data][nth_shifted_eigenvalues]
//     let current_spin_shifted_eigens = eigen_energies_vec
//         .iter()
//         .map(|entry| -> Vec<f64> {
//             let eigen_at_spin = entry.
//                 .unwrap_or_else(|| panic!("No eigenvalues at given spin {}", spin));
//             shift_eigen_by_e_fermi(eigen_at_spin, e_fermi)
//         })
//         .collect();
//     let k_point_weight_array = band_data.kpt_eigen_energies_array().kpt_weight();
//     current_spin_shifted_eigens
//         .par_iter()
//         // Pairs with k-point weights
//         .zip(k_point_weight_array.par_iter())
//         .map(|(de_list_per_k, kpt_weight)| -> Vec<f64> {
//             // E_delta over given energy range
//             // Total DOS so orbital_weights are `None`
//             e_deltas_per_kpt(energy_range, de_list_per_k, *kpt_weight, None)
//         })
//         /*
//         Fold vector of vectors into one vector
//         k1:   [e, ...]
//         k2:   [e, ...]
//         ...     ...
//         kx:   [e, ...]
//                 ||
//                 \/
//         Total: [e, ... ]
//         */
//         .reduce_with(|prev, next| -> Vec<f64> { prev.add(&next) })
//         .unwrap()
// }

fn dos_spin(
    band_data: &Bands,
    energy_range: &[f64],
    spin: u8,
    orbital_weight_array: &[f64],
) -> Vec<f64> {
    let eigen_values_at_spin = band_data
        .kpt_eigen_energies_array()
        .eigen_energies()
        .select(Axis(0), &[spin as usize]);
    let e_fermi = band_data.e_fermi()[0];
    let num_eigens = band_data.num_eigenvalues()[0] as usize;
    let shifted_eigens = (eigen_values_at_spin - e_fermi) * HATREE_TO_EV;
    let shifted_eigens = shifted_eigens
        .to_shape((band_data.num_eigenvalues()[0] as usize, 1 as usize))
        .unwrap();
    let energy_range = Array2::from_shape_vec(
        (400, num_eigens).f(),
        (0..num_eigens)
            .into_iter()
            .map(|_| linspace(-20.0, 20.0, 400).into_iter().collect::<Vec<f64>>())
            .collect::<Vec<Vec<f64>>>()
            .into_iter()
            .flatten()
            .collect::<Vec<f64>>(),
    )
    .unwrap();
    let mut shifted_eigens_mat: Array2<f64> = Array2::zeros((400, num_eigens));
    shifted_eigens_mat
        .axis_iter_mut(Axis(1))
        .zip(shifted_eigens.axis_iter(Axis(0)))
        .for_each(|(mut lcol, rcol)| lcol.fill(rcol[0]));
    println!("{:#.2?}", energy_range);
    println!("{:#.2?}", shifted_eigens_mat);
    let mut energy_matrix = energy_range - shifted_eigens_mat;
    energy_matrix.par_mapv_inplace(|de| gaussian_smearing_delta(de));
    let result = energy_matrix.sum_axis(Axis(1));
    println!("{:#.2?}", result);
    todo!()
}

/**
Shift the eigenvalues by fermi energy and convert Hatree to eV
*/
pub fn shift_eigen_by_e_fermi(eigen_energies_entry: &[f64], e_fermi: f64) -> Vec<f64> {
    eigen_energies_entry
        .par_iter()
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

fn gaussian_smearing_delta(de: f64) -> f64 {
    /*
    The smearing here is the standard deviation
    Same as implemented in CASTEP
    */
    let smear = SMEARING_WIDTH;
    /*
    g(x) = 1 / (σ*sqrt(2PI)) * exp(-1/2 * (e - de)^2/σ^2)
    */
    (-(de / (smear)).powi(2) / 2.0).exp() / ((2.0 * PI).sqrt() * smear)
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
            .par_iter()
            .zip(orbital_weights.par_iter())
            .map(|(de, orb_weight)| e_delta(energies, *de, k_point_weight, *orb_weight))
            // Now we have n Vec<f64> of e_delta
            .reduce_with(|prev, next| -> Vec<f64> {
                // merge Vec<_> by zip and map
                prev.add(&next)
            })
            .unwrap(),
        None => de_array
            .par_iter()
            .map(|de| e_delta(energies, *de, k_point_weight, 1.0))
            // Now we have n Vec<f64> of e_delta
            .reduce_with(|prev, next| -> Vec<f64> {
                // merge Vec<_> by zip and map
                prev.add(&next)
            })
            .unwrap(),
    }
}
#[cfg(test)]
#[test]
fn test_band() {
    use std::fs;

    use ndarray::{Array2, Array3, Axis, Order, ShapeBuilder};
    use rayon::iter::IntoParallelIterator;

    use crate::parser::pdos_weights::PDOSWeight;

    let file = "./Pt_311/Pt_311_12lyr_v20_CO.bands";
    let band_file = fs::read_to_string(file).unwrap();
    let pdos_file =
        fs::read("./Pt_311/Pt_311_12lyr_v20_CO.pdos_weights").expect("Error opening pdos_weights");
    let (_, pdos_weight) = PDOSWeight::parse(&pdos_file).unwrap();
    let (_, band_data) = Bands::parse(&band_file).unwrap();
    let eigen_energies = band_data.kpt_eigen_energies_array().eigen_energies();
    let num_eigens = band_data.num_eigenvalues()[0] as usize;
    let num_kpts = band_data.num_kpts() as usize;
    let num_points = 1000;
    let kpt_weights = band_data.kpt_eigen_energies_array().kpt_weight();

    // PDOS weight
    let species = pdos_weight.species_no();
    let ranks = pdos_weight.rank_in_species();
    let am_channels = pdos_weight.am_channel();

    let orbitals_id: Vec<usize> = species
        .iter()
        .zip(ranks.iter())
        .zip(am_channels.iter())
        .enumerate()
        .filter(|(_, ((&spe, &rank), &am))| spe == 3_u8 && rank == 4_u8 && am == 2_u8)
        .map(|(i, _)| i)
        .collect();
    let projections = pdos_weight
        .orbital_at_eigen_weights()
        .select(Axis(3), &[0])
        .select(Axis(0), &orbitals_id);
    let merged_weights = projections
        .sum_axis(Axis(0))
        .into_shape((num_eigens, num_kpts))
        .unwrap();
    // DOS compute
    let eigen_energies_upspin: Array2<f64> = eigen_energies
        .select(Axis(2), &[0])
        .into_shape((num_eigens, num_kpts))
        .unwrap();
    let eigen_energies_upspin = (eigen_energies_upspin - band_data.e_fermi()[0]) * HATREE_TO_EV;
    let spin_up_dos = eigen_energies_upspin
        .columns()
        .into_iter()
        .zip(merged_weights.columns().into_iter())
        .zip(kpt_weights.iter())
        .par_bridge()
        .map(|((col, orb_weight), kpt_weight)| {
            let e_range: Array2<f64> = Array2::from_shape_vec(
                (num_points, 1).f(),
                linspace(-10.0, 0.0, num_points)
                    .into_iter()
                    .collect::<Vec<f64>>(),
            )
            .unwrap();
            let e_width: Array2<f64> = Array2::ones((1, num_eigens));
            let energy_matrix = e_range.dot(&e_width);
            let col_t = col.into_shape((1, num_eigens)).unwrap();
            let col_height: Array2<f64> = Array2::ones((num_points, 1).f());
            let eigen_mat = col_height.dot(&col_t);
            let mut de_mat = energy_matrix - eigen_mat;
            let orb_weight = orb_weight.to_shape((num_eigens, 1)).unwrap();
            println!("{:?}, {:?}", de_mat.dim(), orb_weight.dim());
            de_mat.par_mapv_inplace(|de| gaussian_smearing_delta(de));
            de_mat.dot(&orb_weight) * *kpt_weight
        })
        .reduce_with(|a, b| a + b)
        .unwrap();
    println!("{:#.5?}", spin_up_dos);
}
