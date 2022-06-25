use crate::{
    parser::{
        bands::Bands,
        pdos_weights::{EigenWeightPerOrb, UnitOrbitalWeight},
    },
    util::ElementWiseAdd,
};

use super::dos_compute::{e_deltas_per_kpt, shift_eigen_by_e_fermi};

pub fn pdos(
    band_data: &Bands,
    energy_range: &[f64],
    weight_matrix: EigenWeightPerOrb,
) -> (Vec<f64>, Option<Vec<f64>>) {
    if band_data.num_spins() == 1 {
        (
            pdos_spin(band_data, energy_range, 1_u8, weight_matrix.upspin()),
            None,
        )
    } else {
        assert_eq!(
            2,
            band_data.num_spins(),
            "Spin cannot be > 2! Current value: {}",
            band_data.num_spins()
        );
        (
            pdos_spin(band_data, energy_range, 1_u8, weight_matrix.upspin()),
            Some(pdos_spin(
                band_data,
                energy_range,
                2_u8,
                weight_matrix.downspin().unwrap(),
            )),
        )
    }
}

fn pdos_spin(
    band_data: &Bands,
    energy_range: &[f64],
    spin: u8,
    orbital_weights: &[UnitOrbitalWeight],
) -> Vec<f64> {
    /*
    vec of energy eigenvalues for every k-point
    length = number of k-points
    */
    let kpt_ne_vec = band_data.kpt_eigen_energies_array().eigen_energies();
    let e_fermi = band_data.e_fermi()[(spin - 1) as usize];
    let current_spin_shifted_eigens: Vec<Vec<f64>> = kpt_ne_vec
        .iter()
        .map(|kpt_ne| -> Vec<f64> {
            let eigen_at_spin = kpt_ne
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
        .zip(orbital_weights.iter())
        .map(|((de_list_per_k, kpt_weight), orb_weight)| -> Vec<f64> {
            // E_delta over given energy range
            // Total DOS so orbital_weights are `None`
            e_deltas_per_kpt(
                energy_range,
                de_list_per_k,
                *kpt_weight,
                Some(orb_weight.weights()),
            )
        })
        // Fold vector of vectors into one vector
        .reduce(|prev, next| -> Vec<f64> { prev.add(&next) })
        .unwrap()
}
