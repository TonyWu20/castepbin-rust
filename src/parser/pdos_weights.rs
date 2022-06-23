use nom::{multi::count, IResult};

use super::general::{
    parse_multi_f64, parse_multi_f64_from_record, parse_record, parse_u32, parse_u32_from_record,
    parse_u32_vec_from_record,
};

// Parsed from <seed>.pdos_weight
// Information to calculate projected-DOS
#[derive(Debug)]
pub struct PDOSWeight {
    num_kpoints: u32,
    num_spins: u32,
    num_orbitals: u32,
    num_eigenvalues: u32,
    species_no: Vec<u32>,                   // length = num_orbitals
    rank_in_species: Vec<u32>,              // length = num_orbitals
    am_channel: Vec<u32>,                   // length = num_orbitals
    weight_data: Vec<KpointOrbitalWeights>, // length = num_kpoints
}

impl PDOSWeight {
    pub fn new(
        num_kpoints: u32,
        num_spins: u32,
        num_orbitals: u32,
        num_eigenvalues: u32,
        species_no: Vec<u32>,
        rank_in_species: Vec<u32>,
        am_channel: Vec<u32>,
        weight_data: Vec<KpointOrbitalWeights>,
    ) -> Self {
        Self {
            num_kpoints,
            num_spins,
            num_orbitals,
            num_eigenvalues,
            species_no,
            rank_in_species,
            am_channel,
            weight_data,
        }
    }

    pub fn parse<'a>(data: &'a [u8]) -> IResult<&'a [u8], Self> {
        let (i, n_kpts) = parse_u32_from_record(data)?; // Line for num of k-points
        let (i, n_spins) = parse_u32_from_record(i)?; // Line for num of spins
        let (i, n_orbitals) = parse_u32_from_record(i)?; // Line for num of orbitals
        let (i, n_eigens) = parse_u32_from_record(i)?; // Line for num of eigenvalues

        // Array marking species id in the system
        let (i, arr_species_no) = parse_u32_vec_from_record(i)?;
        // Array marking nth-atom of the species
        // E.g. We have two C and 1 O.
        // It will be [1,...,2,...,1,...] to denote C1, C2 and O1
        let (i, arr_rank_in_species) = parse_u32_vec_from_record(i)?;
        // Array marking the orbital angular moment (s, p, d, f,..)
        // Begin from 0 for s
        let (i, arr_am_channel) = parse_u32_vec_from_record(i)?;
        // Collect all weights for each k-point
        // let (i, weights_for_kpoints) = many0(KpointOrbitalWeights::parse)(i)?;
        let parse_kpoint_orbital_weights =
            |data: &'a [u8]| -> IResult<&[u8], KpointOrbitalWeights> {
                let (i, kpoint_info) = parse_record(data)?;
                let (vec_str, nth_kpoint) = parse_u32(kpoint_info)?;
                let (_, kpoint_vectors) = parse_multi_f64(vec_str)?;
                let (i, orbital_weight_items) = count(OrbitalWeight::parse, n_spins as usize)(i)?;
                Ok((
                    i,
                    KpointOrbitalWeights::new(nth_kpoint, kpoint_vectors, orbital_weight_items),
                ))
            };
        let (i, weights_for_kpoints) = count(parse_kpoint_orbital_weights, n_kpts as usize)(i)?;
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
                weights_for_kpoints,
            ),
        ))
    }

    pub fn num_kpoints(&self) -> u32 {
        self.num_kpoints
    }

    pub fn num_spins(&self) -> u32 {
        self.num_spins
    }

    pub fn num_orbitals(&self) -> u32 {
        self.num_orbitals
    }

    pub fn num_eigenvalues(&self) -> u32 {
        self.num_eigenvalues
    }

    pub fn species_no(&self) -> &[u32] {
        self.species_no.as_ref()
    }

    pub fn rank_in_species(&self) -> &[u32] {
        self.rank_in_species.as_ref()
    }

    pub fn am_channel(&self) -> &[u32] {
        self.am_channel.as_ref()
    }

    pub fn weight_data(&self) -> &[KpointOrbitalWeights] {
        self.weight_data.as_ref()
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
 END
*/
#[derive(Debug)]
pub struct KpointOrbitalWeights {
    nth_kpoint: u32,                          // Id for n-th kpoint (start from 1)
    kpoint_vectors: Vec<f64>,                 // Vector of kpoint, length = 3
    orbital_weight_items: Vec<OrbitalWeight>, // At most 2, record the weights of orbital at each E_nk of k-point
}

impl KpointOrbitalWeights {
    pub fn new(
        nth_kpoint: u32,
        kpoint_vectors: Vec<f64>,
        orbital_weight_items: Vec<OrbitalWeight>,
    ) -> Self {
        Self {
            nth_kpoint,
            kpoint_vectors,
            orbital_weight_items,
        }
    }

    pub fn nth_kpoint(&self) -> u32 {
        self.nth_kpoint
    }

    pub fn kpoint_vectors(&self) -> &[f64] {
        self.kpoint_vectors.as_ref()
    }

    pub fn orbital_weight_items(&self) -> &[OrbitalWeight] {
        self.orbital_weight_items.as_ref()
    }
}

// Orbital weigthts data for every kpoint and spin
#[derive(Debug)]
pub struct OrbitalWeight {
    spin: u32,
    weights_for_eigen: Vec<Vec<f64>>, // 2D array, weights_for_eigen[nth_eigen][nth_weight]
}
impl OrbitalWeight {
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

    pub fn weights_for_eigen(&self) -> &[Vec<f64>] {
        self.weights_for_eigen.as_ref()
    }
}
