use crate::config::pdos_config::{PDOSTargets, PDOSTask};
use crate::config::TaskProcess;
use crate::dos::dos_compute::{compute_energy_range, dos_spin, HATREE_TO_EV};
use crate::parser::bands::Bands;
use crate::util::ElementWiseAdd;
use std::{fmt, fs};

use ndarray::prelude::*;
use ndarray::{ArrayBase, ArrayViewD, Dim, ShapeBuilder, ViewRepr};
use nom::{error::Error, IResult};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::parser::{
    castep_bin::{FromCastepCheck, UnitCell},
    pdos_weights::PDOSWeight,
};

use super::PDOS;

#[derive(Debug, Clone)]
struct AngularMomentParsingError;

impl fmt::Display for AngularMomentParsingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Invalid angular moment id (>3)")
    }
}

impl TaskProcess for PDOSTask {
    type Output = PDOS;
    type Target = PDOSTargets;
    type Reference = (UnitCell, PDOSWeight, Bands);
    fn compute_target(
        &self,
        target: &PDOSTargets,
        reference: &(UnitCell, PDOSWeight, Bands),
    ) -> PDOS {
        let (cell, pdos_weights, band_data) = reference;
        let num_eigens = band_data.num_eigenvalues()[0] as usize;
        let num_kpts = band_data.num_kpts() as usize;
        let command = target.translate(&cell);
        let species_arr = pdos_weights.species_no();
        let ranks_arr = pdos_weights.rank_in_species();
        let am_channel_arr = pdos_weights.am_channel();
        let mut species_symbols: Vec<String> = vec![];
        let mut ranks: Vec<u8> = vec![];
        let mut am_channels: Vec<String> = vec![];
        // find projector orbitals, grouped by am channel
        let orbital_ids: Vec<Vec<usize>> = command
            .species()
            .iter()
            .zip(command.rank_in_species().iter())
            .zip(command.am_channel().iter())
            .map(|combo| -> Vec<usize> {
                let ((&spe, &rank), &am) = combo;
                let symbol = cell.species_symbol()[spe as usize - 1].to_string();
                let find_am_symbol = |am: u8| -> Result<String, AngularMomentParsingError> {
                    match am {
                        0_u8 => Ok("s".to_string()),
                        1_u8 => Ok("p".to_string()),
                        2_u8 => Ok("d".to_string()),
                        3_u8 => Ok("f".to_string()),
                        _ => Err(AngularMomentParsingError),
                    }
                };
                let am_symbol = find_am_symbol(am).unwrap();
                species_symbols.push(symbol);
                ranks.push(rank);
                am_channels.push(am_symbol);
                species_arr
                    .iter()
                    .zip(ranks_arr.iter())
                    .zip(am_channel_arr.iter())
                    .enumerate()
                    .filter(|(_, ((&spe_id, &rank_id), &am_id))| {
                        spe_id == spe && rank_id == rank && am_id == am
                    })
                    .map(|(i, _)| i)
                    .collect()
            })
            .collect();
        let (e_min, e_max) = (
            (band_data.min_eigen_energy() - band_data.e_fermi()[0]) * HATREE_TO_EV - 2.00,
            (band_data.max_eigen_energy() - band_data.e_fermi()[0]) * HATREE_TO_EV + 2.00,
        );
        let num_per_ev = 100;
        let num_points = ((e_max - e_min) * num_per_ev as f64) as u32;
        let energy_range = compute_energy_range(e_min, e_max, num_points);

        match self.spin() {
            true => {
                let spin_up_res: Vec<f64> = orbital_ids
                    .par_iter()
                    .flat_map(|projector_ids| {
                        let merged_orbital_weights = get_merged_weights(
                            &pdos_weights,
                            0_u8,
                            &projector_ids,
                            num_eigens,
                            num_kpts,
                        );
                        dos_spin(band_data, &energy_range, 0_u8, &merged_orbital_weights)
                    })
                    .collect();
                let spin_down_res: Vec<f64> = orbital_ids
                    .par_iter()
                    .flat_map(|projector_ids| {
                        let merged_orbital_weights = get_merged_weights(
                            &pdos_weights,
                            1_u8,
                            &projector_ids,
                            num_eigens,
                            num_kpts,
                        );
                        dos_spin(band_data, &energy_range, 1_u8, &merged_orbital_weights)
                    })
                    .collect();
                PDOS::new(
                    2_u8,
                    num_points,
                    orbital_ids.len(),
                    energy_range,
                    vec![spin_up_res, spin_down_res],
                    species_symbols,
                    ranks,
                    am_channels,
                )
            }
            false => {
                let spin_up_res: Vec<f64> = orbital_ids
                    .par_iter()
                    .flat_map(|projector_ids| {
                        let merged_orbital_weights = get_merged_weights(
                            &pdos_weights,
                            0_u8,
                            &projector_ids,
                            num_eigens,
                            num_kpts,
                        );
                        dos_spin(band_data, &energy_range, 0_u8, &merged_orbital_weights)
                    })
                    .collect();
                PDOS::new(
                    1_u8,
                    num_points,
                    orbital_ids.len(),
                    energy_range,
                    vec![spin_up_res],
                    species_symbols,
                    ranks,
                    am_channels,
                )
            }
        }
    }
    fn task_execute(&self) -> Vec<Self::Output> {
        let castep_bin_content = fs::read(self.castep_bin_filename())
            .unwrap_or_else(|e| panic!("Reading {}.castep_bin error, {}", self.seed(), e));
        let pdos_weights_content = fs::read(self.pdos_weight_filename())
            .unwrap_or_else(|e| panic!("Reading {}.pdos_weights error, {}", self.seed(), e));
        let bands_content = fs::read_to_string(self.bands_filename())
            .unwrap_or_else(|e| panic!("Reading {}.bands error, {}", self.seed(), e));
        let (_, cell) = UnitCell::parser(&castep_bin_content).unwrap_or_else(|e| {
            panic!(
                "Parsing unit cell from {}.castep_bin error, {}",
                self.seed(),
                e
            )
        });
        let (_, pdos_weights) = PDOSWeight::parse(&pdos_weights_content).unwrap_or_else(|e| {
            panic!(
                "Parsing pdos weights from {}.pdos_weights error, {}",
                self.seed(),
                e
            )
        });
        let (_, band_data) = Bands::parse(&bands_content)
            .unwrap_or_else(|e| panic!("Parsing bands from {}.bands error, {}", self.seed(), e));
        let reference = (cell, pdos_weights, band_data);
        self.targets()
            .iter()
            .map(|target| self.compute_target(target, &reference))
            .collect()
    }
}

fn get_merged_weights(
    pdos_weights: &PDOSWeight,
    spin_channel: u8,
    orbital_ids: &[usize],
    num_eigens: usize,
    num_kpts: usize,
) -> Array2<f64> {
    assert!(spin_channel < 2, "Incorrect spin number! (0 or 1 only)");
    let projections = pdos_weights
        .orbital_at_eigen_weights()
        .select(Axis(3), &[spin_channel as usize])
        .select(Axis(0), &orbital_ids);
    projections
        .sum_axis(Axis(0))
        .into_shape((num_eigens, num_kpts))
        .unwrap()
}

// #[cfg(test)]
// mod test {
//     use std::{fs, ops::Add};
//
//     use ndarray::{Array, Array1, Array2, Axis};
//     use ndarray::{Order, ShapeBuilder};
//
//     use crate::config::{Config, TaskProcess};
//     use crate::{parser::pdos_weights::PDOSWeight, util::ElementWiseAdd};
//
// }
