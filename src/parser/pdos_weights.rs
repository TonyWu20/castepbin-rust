use std::{collections::HashMap, ops::Add, vec};

use ndarray::{s, Array1, Array2, Array3, Array4, ShapeBuilder};
use nom::{
    multi::{count, many_m_n},
    IResult,
};

use crate::util::ElementWiseAdd;

use super::general::{
    parse_multi_f64, parse_multi_f64_from_record, parse_record, parse_u32, parse_u32_from_record,
    parse_u32_vec_from_record, parse_u8_from_record, parse_u8_vec_from_u32_vec_record,
};

type WeightsForEigenAtK = Vec<Vec<f64>>;
/**
Parsed from <seed>.pdos_weight
Information to calculate projected-DOS
* Notes:
    * num_kpoints, num_spins, species_no, am_channel should not exceed 256 in practice,
    so use u8 to save memory.
    * rank_in_species means the Atom ID for a species, currently we should not be
    dealing with systems where any element has over 256 atoms. So u8 should also
    be fine for the moment.
*/
#[derive(Debug)]
pub struct PDOSWeight {
    num_kpoints: u8,
    num_spins: u8,
    num_orbitals: u32,
    num_eigenvalues: u32,
    /// length = num_orbitals
    species_no: Vec<u8>,
    /// length = num_orbitals
    rank_in_species: Vec<u8>,
    /// length = num_orbitals
    am_channel: Vec<u8>,
    /// dimension = (num_kpoints, num_spins, num_eigenvalues, num_orbitals)
    orbital_at_eigen_weights: Array4<f64>,
}

impl PDOSWeight {
    pub fn new(
        num_kpoints: u8,
        num_spins: u8,
        num_orbitals: u32,
        num_eigenvalues: u32,
        species_no: Vec<u8>,
        rank_in_species: Vec<u8>,
        am_channel: Vec<u8>,
        orbital_at_eigen_weights: Array4<f64>,
    ) -> Self {
        Self {
            num_kpoints,
            num_spins,
            num_orbitals,
            num_eigenvalues,
            species_no,
            rank_in_species,
            am_channel,
            orbital_at_eigen_weights,
        }
    }
    pub fn filter_species_rank_am(
        &self,
        species_id: u8,
        rank_id: u8,
        am_channel_id: u8,
    ) -> Option<Vec<usize>> {
        let projections: Vec<usize> = self
            .species_no()
            .iter()
            .zip(self.rank_in_species().iter())
            .zip(self.am_channel().iter())
            .enumerate()
            .filter(|(_, ((&s, &r), &am))| s == species_id && r == rank_id && am == am_channel_id)
            .map(|(i, _)| i)
            .collect();
        match projections.is_empty() {
            true => {
                println!("Did not find projections!");
                None
            }
            false => Some(projections),
        }
    }

    pub fn parse(data: &[u8]) -> IResult<&[u8], Self> {
        let (i, n_kpts) = parse_u8_from_record(data)?; // Line for num of k-points
        let (i, n_spins) = parse_u8_from_record(i)?; // Line for num of spins
        let (i, n_orbitals) = parse_u32_from_record(i)?; // Line for num of orbitals
        let (i, n_eigens) = parse_u32_from_record(i)?; // Line for num of eigenvalues

        // Array marking species id in the system
        // The numbers are be_u32 format, but we turn them to be u8 to maximize
        // memory efficiency
        let (i, arr_species_no) = parse_u8_vec_from_u32_vec_record(i)?;
        // Array marking nth-atom of the species
        // E.g. We have two C and 1 O.
        // It will be [1,...,2,...,1,...] to denote C1, C2 and O1
        // The numbers are be_u32 format, but we turn them to be u8 to maximize
        // memory efficiency
        let (i, arr_rank_in_species) = parse_u8_vec_from_u32_vec_record(i)?;
        // Array marking the orbital angular moment (s, p, d, f,..)
        // Begin from 0 for s. Obviously it could not overflow u8 type.
        // The numbers are be_u32 format, but we turn them to be u8 to maximize
        // memory efficiency
        let (i, arr_am_channel) = parse_u8_vec_from_u32_vec_record(i)?;
        // Collect all weights for each k-point
        // let (i, weights_for_kpoints) = many0(KpointOrbitalWeights::parse)(i)?;
        let (i, weights_for_kpoints) =
            count(KpointOrbitalWeightsAtEigen::parse, n_kpts as usize)(i)?;
        let all_weights_for_kpoints: Vec<Vec<Vec<Vec<f64>>>> = (0..n_spins as usize)
            .into_iter()
            .map(|spin| -> Vec<WeightsForEigenAtK> {
                weights_for_kpoints
                    .iter()
                    .map(|item| item.orbital_weight_items()[spin].to_vec())
                    .collect::<Vec<WeightsForEigenAtK>>()
            })
            .collect();
        let flattened_weights: Vec<f64> = all_weights_for_kpoints
            .into_iter()
            .flatten()
            .flatten()
            .flatten()
            .collect();
        let weights_array = Array4::from_shape_vec(
            (
                n_orbitals as usize,
                n_eigens as usize,
                n_kpts as usize,
                n_spins as usize,
            )
                .f(),
            flattened_weights,
        )
        .unwrap();
        Ok((
            i,
            PDOSWeight::new(
                n_kpts,
                n_spins,
                n_orbitals,
                n_eigens,
                arr_species_no,
                arr_rank_in_species,
                arr_am_channel,
                weights_array,
            ),
        ))
    }

    pub fn num_kpoints(&self) -> u8 {
        self.num_kpoints
    }

    pub fn num_spins(&self) -> u8 {
        self.num_spins
    }

    pub fn num_orbitals(&self) -> u32 {
        self.num_orbitals
    }

    pub fn num_eigenvalues(&self) -> u32 {
        self.num_eigenvalues
    }

    pub fn species_no(&self) -> &[u8] {
        self.species_no.as_ref()
    }

    pub fn rank_in_species(&self) -> &[u8] {
        self.rank_in_species.as_ref()
    }

    pub fn am_channel(&self) -> &[u8] {
        self.am_channel.as_ref()
    }

    pub fn orbital_at_eigen_weights(&self) -> &Array4<f64> {
        &self.orbital_at_eigen_weights
    }
}

/**
 Struct for orbital weights of each k-point
 The orbital weight data in <seed>.pdos_weights are formatted as:
 ...(general information lines)
 nth_k-point, [x, y, z] (k-point vector)
 Spin number (1, 2 if spin-polarised and down-spin)
 weight for 1st eigenvalue, [f64; number of orbitals]
 ... (continue for other eigenvalue)
 Spin number
 weights
 next k-point, [x,y,z],
 ... (continue as above)
 EN
*/
#[derive(Debug)]
pub struct KpointOrbitalWeightsAtEigen {
    nth_kpoint: u32,          // Id for n-th kpoint (start from 1)
    kpoint_vectors: Vec<f64>, // Vector of kpoint, length = 3
    /**!
    Orbital weigthts data for every kpoint and spin
    This struct represents the data serialization in
    .pdos_weights but not optimal for indexing through
    the info-array recorded above.
    Store spin-orbital_weights, we will reorganize it later
    */
    orbital_weight_items: Vec<WeightsForEigenAtK>,
}

impl KpointOrbitalWeightsAtEigen {
    pub fn new(
        nth_kpoint: u32,
        kpoint_vectors: Vec<f64>,
        orbital_weight_items: Vec<WeightsForEigenAtK>,
    ) -> Self {
        Self {
            nth_kpoint,
            kpoint_vectors,
            orbital_weight_items,
        }
    }

    pub fn parse(data: &[u8]) -> IResult<&[u8], Self> {
        let (i, kpoint_info) = parse_record(data)?;
        let (vec_str, nth_kpoint) = parse_u32(kpoint_info)?;
        let (_, kpoint_vectors) = parse_multi_f64(vec_str)?;
        let (i, orbital_weight_items): (&[u8], Vec<WeightsForEigenAtK>) =
            many_m_n(1, 2, Self::orbital_weights_parse)(i)?;
        Ok((
            i,
            KpointOrbitalWeightsAtEigen::new(nth_kpoint, kpoint_vectors, orbital_weight_items),
        ))
    }

    fn orbital_weights_parse(data: &[u8]) -> IResult<&[u8], WeightsForEigenAtK> {
        let (i, _) = parse_u32_from_record(data)?;
        let (i, num_eigen) = parse_u32_from_record(i)?;
        let (i, all_weights) = count(parse_multi_f64_from_record, num_eigen as usize)(i)?;
        Ok((i, all_weights))
    }

    pub fn nth_kpoint(&self) -> u32 {
        self.nth_kpoint
    }

    pub fn kpoint_vectors(&self) -> &[f64] {
        self.kpoint_vectors.as_ref()
    }

    pub fn orbital_weight_items(&self) -> &[WeightsForEigenAtK] {
        self.orbital_weight_items.as_ref()
    }
}
