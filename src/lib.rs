#![allow(dead_code, unused_imports)]
use std::{fs, str::from_utf8};

use na::{Matrix3, MatrixXx3, OMatrix, Vector3};
use nom::{
    bytes::complete::take_until,
    multi::{count, many0},
    number::complete::{be_f64, be_u32},
    sequence::tuple,
    IResult,
};
use parser::{general::parse_multi_f64_from_record, *};

use crate::parser::{
    castep_bin::skip_to_optimized_cell,
    general::{
        parse_data, parse_f64, parse_multi_f64, parse_record, parse_u32, parse_u32_from_record,
        parse_u32_vec_from_record,
    },
};

extern crate nalgebra as na;
pub mod dos;
pub mod parser;
pub mod test;

pub const BOHR_RADIUS: f64 = 1.8897259886;
pub const AMU: f64 = 1822.8839;
#[derive(Debug, Clone)]
pub struct UnitCell {
    real_lattice_vectors: Matrix3<f64>,
    recip_lattice_vectors: Matrix3<f64>,
    volume: f64,
    num_species: u32,
    num_ions: u32,
    max_ions_in_species: u32,
    num_ions_in_species: u32,
    ionic_positions: Option<MatrixXx3<f64>>,
    ionic_charge_real: f64,
    species_mass: Vec<f64>,
    species_symbol: Vec<String>,
}

impl UnitCell {
    pub fn new(
        real_lattice_vectors: Matrix3<f64>,
        recip_lattice_vectors: Matrix3<f64>,
        volume: f64,
        num_species: u32,
        num_ions: u32,
        max_ions_in_species: u32,
        num_ions_in_species: u32,
        ionic_positions: Option<MatrixXx3<f64>>,
        ionic_charge_real: f64,
        species_mass: Vec<f64>,
        species_symbol: Vec<String>,
    ) -> Self {
        Self {
            real_lattice_vectors,
            recip_lattice_vectors,
            volume,
            num_species,
            num_ions,
            max_ions_in_species,
            num_ions_in_species,
            ionic_positions,
            ionic_charge_real,
            species_mass,
            species_symbol,
        }
    }
}

pub struct CellGlobal {
    n_kpts: u32,
    kpoints: MatrixXx3<f64>,
    kpoint_weights: Vec<f64>,
    num_symmetry_operations: u32,
    num_crystal_symmetry_operations: u32,
    cell_symmetry_ion_constraints: f64,
    bs_kpoints: MatrixXx3<f64>,
    bs_kpoint_weights: Vec<f64>,
}

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
}

/*
** Struct for orbital weights of each k-point
** The orbital weight data in <seed>.pdos_weights are formatted as:
** ...(general information lines)
** nth_k-point, [x, y, z] (k-point vector)
** Spin number (1, 2 if spin-polarised and down-spin)
** weight for 1st eigenvalue, [f64; number of orbitals]
** ... (continue for other eigenvalue)
** Spin number
** weights
** next k-point, [x,y,z],
** ... (continue as above)
** END
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
}

trait FromCastepCheck {
    type ParsedStruct;
    fn parser(data: &[u8]) -> IResult<&[u8], Self::ParsedStruct>;
}

impl FromCastepCheck for UnitCell {
    type ParsedStruct = UnitCell;
    fn parser(data: &[u8]) -> IResult<&[u8], UnitCell> {
        let (left, _) = parse_record(data)?; // Skip "BEGIN_UNIT_CELL"
        let (left, real_lat) = parse_data(left)?;
        let (left, recip_lat) = parse_data(left)?;
        let (left, cell_volume) = parse_data(left)?;
        let (left, num_species) = parse_data(left)?;
        let (left, num_ions) = parse_data(left)?;
        let (left, max_ions_in_species) = parse_data(left)?;
        let (left, _) = parse_data(left)?; // Skip "CELL_VERSION_NUMBER"
        let (left, num_ions_in_species) = parse_data(left)?;
        let (left, ionic_positions) = parse_data(left)?;
        let (left, _) = parse_data(left)?; // Skip "IONIC_VELOCITIES"
        let (left, _) = parse_data(left)?; // Skip "IONIC_CHARGE"
        let (left, ionic_charge_real) = parse_data(left)?;
        let (left, _) = parse_data(left)?; // Skip "ATOM_MOVE"
        let (left, species_mass) = parse_data(left)?;
        let (remains, species_symbol) = parse_data(left)?;
        let real_lat: Vec<f64> = parse_multi_f64(real_lat)
            .unwrap()
            .1
            .iter()
            .map(|i| -> f64 { i / BOHR_RADIUS })
            .collect();
        println!("{}", &real_lat.len());
        let cell_real_lat: Matrix3<f64> = Matrix3::from_vec(real_lat);
        let cell_volume =
            parse_f64(cell_volume).unwrap().1 / BOHR_RADIUS / BOHR_RADIUS / BOHR_RADIUS;
        let num_species = parse_u32(num_species).unwrap().1;
        let num_ions = parse_u32(num_ions).unwrap().1;
        let max_ions_in_species = parse_u32(max_ions_in_species).unwrap().1;
        let num_ions_in_species = parse_u32(num_ions_in_species).unwrap().1;
        let ionic_positions: Option<MatrixXx3<f64>> = Some(MatrixXx3::from_vec(
            parse_multi_f64(ionic_positions).unwrap().1,
        ));
        let ionic_charge_real = parse_f64(ionic_charge_real).unwrap().1;
        let species_mass = parse_multi_f64(species_mass)
            .unwrap()
            .1
            .iter()
            .map(|i| -> f64 { i / AMU })
            .collect();
        let species_symbol: Vec<String> = from_utf8(species_symbol)
            .unwrap()
            .split_whitespace()
            .map(|str| -> String { str.to_string() })
            .collect();
        // let recip_a = 2.0 * PI / cell_volume * cell_real_lat.column(1).cross(column(2));
        let recip_lattice_vectors = Matrix3::from_vec(
            parse_multi_f64(recip_lat)
                .unwrap()
                .1
                .iter()
                .map(|i| -> f64 { i * BOHR_RADIUS })
                .collect(),
        );
        Ok((
            remains,
            UnitCell::new(
                cell_real_lat,
                recip_lattice_vectors,
                cell_volume,
                num_species,
                num_ions,
                max_ions_in_species,
                num_ions_in_species,
                ionic_positions,
                ionic_charge_real,
                species_mass,
                species_symbol,
            ),
        ))
        // let todo!();
    }
}
#[ignore]
#[test]
fn test_unit_cell() {
    let file = fs::read("./Si2.castep_bin").unwrap();
    let (cell_content, _) = skip_to_optimized_cell(&file).unwrap();
    let (_, cell) = UnitCell::parser(&cell_content).unwrap();
    println!("{:#?}", cell);
    let file = fs::read("./example.castep_bin").unwrap();
    let (cell_content, _) = skip_to_optimized_cell(&file).unwrap();
    let (_, cell) = UnitCell::parser(&cell_content).unwrap();
    println!("{:#?}", cell);
}
#[test]
fn test_pdos_weights_struct() {
    let file = fs::read("./Si_2_custom CASTEP GeomOpt/Si_2_custom.pdos_weights").unwrap();
    let (_, pdos_weight) = PDOSWeight::parse(&file).unwrap();
    println!("Si_2");
    println!("{:#?}", pdos_weight);
    let file = fs::read("./Pt_310_12lyr_v20_CO_DOS/Pt_310_12lyr_v20_CO_DOS.pdos_weights").unwrap();
    let (_, pdos_weight) = PDOSWeight::parse(&file).unwrap();
    println!("Pt_310");
    println!("{:#?}", pdos_weight);
}
