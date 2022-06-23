#![allow(clippy::too_many_arguments, dead_code, unused_imports)]
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
