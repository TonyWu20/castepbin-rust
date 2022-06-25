use std::{collections::HashMap, vec};

use nom::{
    multi::{count, many_m_n},
    IResult,
};

use super::general::{
    parse_multi_f64, parse_multi_f64_from_record, parse_record, parse_u32, parse_u32_from_record,
    parse_u32_vec_from_record, parse_u8_from_record, parse_u8_vec_from_u32_vec_record,
};

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
    species_no: Vec<u8>,                           // length = num_orbitals
    rank_in_species: Vec<u8>,                      // length = num_orbitals
    am_channel: Vec<u8>,                           // length = num_orbitals
    orbital_eigen_weights: Vec<EigenWeightPerOrb>, // length = num_orbitals
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
        orbital_eigen_weights: Vec<EigenWeightPerOrb>,
    ) -> Self {
        Self {
            num_kpoints,
            num_spins,
            num_orbitals,
            num_eigenvalues,
            species_no,
            rank_in_species,
            am_channel,
            orbital_eigen_weights,
        }
    }
    pub fn search_species_rank_am(
        species_id: u32,
        rank_id: u32,
        am_channel_id: u32,
        search_hash: &HashMap<(u8, u8, u8), usize>,
    ) -> Option<usize> {
        let return_idx = search_hash.get(&(species_id as u8, rank_id as u8, am_channel_id as u8));
        match return_idx {
            Some(id) => Some(*id),
            None => None,
        }
    }

    fn hash_identity_arrays(&self) -> HashMap<(u8, u8, u8), usize> {
        let mut hashtable = HashMap::new();
        self.species_no()
            .into_iter()
            .zip(self.rank_in_species().into_iter())
            .zip(self.am_channel().into_iter())
            .into_iter()
            .enumerate()
            .for_each(|(i, record)| {
                let ((spe, rank), am) = record;
                hashtable.insert((*spe, *rank, *am), i);
            });
        hashtable
    }

    pub fn parse<'a>(data: &'a [u8]) -> IResult<&'a [u8], Self> {
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
        let orbital_weight_items: Vec<EigenWeightPerOrb> = (0..n_orbitals)
            .map(|i| -> EigenWeightPerOrb {
                EigenWeightPerOrb::from_array_of_kpoint_orbital_weight(
                    &weights_for_kpoints,
                    i as usize,
                )
            })
            .collect();
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
                orbital_weight_items,
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

    pub fn orbital_eigen_weights(&self) -> &[EigenWeightPerOrb] {
        self.orbital_eigen_weights.as_ref()
    }
    pub fn eigen_weights_at_orbital(&self, orbital_id: u32) -> &EigenWeightPerOrb {
        &self.orbital_eigen_weights[orbital_id as usize]
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
    nth_kpoint: u32,                                 // Id for n-th kpoint (start from 1)
    kpoint_vectors: Vec<f64>,                        // Vector of kpoint, length = 3
    orbital_weight_items: Vec<OrbitalWeightAtEigen>, // At most 2, record the weights of orbital at each E_nk of k-point
}

impl KpointOrbitalWeightsAtEigen {
    pub fn new(
        nth_kpoint: u32,
        kpoint_vectors: Vec<f64>,
        orbital_weight_items: Vec<OrbitalWeightAtEigen>,
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
        let (i, orbital_weight_items) = many_m_n(1, 2, OrbitalWeightAtEigen::parse)(i)?;
        Ok((
            i,
            KpointOrbitalWeightsAtEigen::new(nth_kpoint, kpoint_vectors, orbital_weight_items),
        ))
    }

    pub fn nth_kpoint(&self) -> u32 {
        self.nth_kpoint
    }

    pub fn kpoint_vectors(&self) -> &[f64] {
        self.kpoint_vectors.as_ref()
    }

    pub fn orbital_weight_items(&self) -> &Vec<OrbitalWeightAtEigen> {
        self.orbital_weight_items.as_ref()
    }
}

/**
Orbital weigthts data for every kpoint and spin
This struct represents the data serialization in
.pdos_weights but not optimal for indexing through
the info-array recorded above.
*/
#[derive(Debug)]
pub struct OrbitalWeightAtEigen {
    spin: u32,
    weights_for_eigen: Vec<Vec<f64>>, // 2D array, weights_for_eigen[nth_eigen][nth_weight]
}
impl OrbitalWeightAtEigen {
    fn new(spin: u32, weights_for_eigen: Vec<Vec<f64>>) -> Self {
        Self {
            spin,
            weights_for_eigen,
        }
    }
    pub fn parse(data: &[u8]) -> IResult<&[u8], Self> {
        let (i, spin_number) = parse_u32_from_record(data)?;
        let (i, num_eigen) = parse_u32_from_record(i)?;
        let (i, all_weights) = count(parse_multi_f64_from_record, num_eigen as usize)(i)?;
        Ok((i, Self::new(spin_number, all_weights)))
    }

    pub fn spin(&self) -> u32 {
        self.spin
    }

    pub fn weights_for_eigen(&self) -> &Vec<Vec<f64>> {
        self.weights_for_eigen.as_ref()
    }
}

/// Struct containing the weights for every eigenvalues of this orbital.
#[derive(Debug)]
pub struct EigenWeightPerOrb {
    ///! Length of the vector equals to number of k-points.
    num_spins: u8,
    eigen_weights_for_each_k: Vec<UnitOrbitalWeight>,
}

impl EigenWeightPerOrb {
    pub fn new(num_spins: u8, eigen_weights_for_each_k: Vec<UnitOrbitalWeight>) -> Self {
        Self {
            num_spins,
            eigen_weights_for_each_k,
        }
    }

    pub fn num_spins(&self) -> u8 {
        self.num_spins
    }

    pub fn eigen_weights_for_each_k(&self) -> &[UnitOrbitalWeight] {
        self.eigen_weights_for_each_k.as_ref()
    }
    pub fn upspin(&self) -> &UnitOrbitalWeight {
        &self.eigen_weights_for_each_k[0]
    }
    pub fn downspin(&self) -> &UnitOrbitalWeight {
        &self.eigen_weights_for_each_k[1]
    }
    pub fn from_array_of_kpoint_orbital_weight(
        array_kpow: &Vec<KpointOrbitalWeightsAtEigen>,
        orbital_id: usize,
    ) -> Self {
        if array_kpow[0].orbital_weight_items().len() == 2 {
            let spin_up_weight_total = Self::collect_weight_for_orbital(array_kpow, orbital_id, 1);
            let spin_down_weight_total =
                Self::collect_weight_for_orbital(array_kpow, orbital_id, 2);
            return Self::new(2 as u8, vec![spin_up_weight_total, spin_down_weight_total]);
        } else {
            let total = Self::collect_weight_for_orbital(array_kpow, orbital_id, 1);
            return Self::new(1 as u8, vec![total]);
        }
    }
    /// Helper function to collect merged weights for selected orbital
    fn collect_weight_for_orbital(
        array_kpow: &Vec<KpointOrbitalWeightsAtEigen>,
        orbital_id: usize,
        spin: usize,
    ) -> UnitOrbitalWeight {
        // for every item in array_kpow
        //      item is a KpointOrbitalWeights (index = k-point num)
        // for every item in KpointOrbitalWeights.orbital_weight_items
        //      item is a OrbitalWeight for every eigenvalues, and contain up-spin down-spin
        //      (index = eigen-values)
        // We should: get the data for each spin, make into a matrix, get the orbital_id-th column,
        // repeat across the k-point, get an array of [(UnitOrbitalWeight); num_kpts]
        // For the moment, directly merge with k-point weight to return the
        // k-point-weighted eigenvalue weight for the chosen orbital.
        // spin_up_matrix[nth_kpt][nth_eigen_values][nth_orbitals]
        let spin_matrix: Vec<Vec<f64>> = array_kpow
            .iter()
            .map(|item| -> Vec<f64> {
                item.orbital_weight_items()[spin - 1] // up-spin
                    .weights_for_eigen() // get &Vec<Vec<f64>
                    .iter()
                    .map(|weight_at_eigen| -> f64 { weight_at_eigen[orbital_id] })
                    .collect()
            })
            .collect();
        let mut merged_vector = vec![0.0; spin_matrix[0].len()];
        spin_matrix
            .iter()
            // every k-point
            .for_each(|item| {
                //item is the weight at eigenvalue
                item.iter().enumerate().for_each(|(i, val)| {
                    // Merge weights at all eigenvalue for this orbital
                    merged_vector[i] += val;
                })
            });
        UnitOrbitalWeight {
            weights_for_each_eigen: merged_vector,
        }
    }
}

/**
Unit struct for weights array of eigenvalue for each orbital
*/
#[derive(Debug)]
pub struct UnitOrbitalWeight {
    weights_for_each_eigen: Vec<f64>,
}

impl UnitOrbitalWeight {
    pub fn new(weights_for_each_eigen: Vec<f64>) -> Self {
        Self {
            weights_for_each_eigen,
        }
    }

    pub fn weights(&self) -> &[f64] {
        self.weights_for_each_eigen.as_ref()
    }
}
