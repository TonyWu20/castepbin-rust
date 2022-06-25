use nalgebra::Matrix3;
use nom::{
    bytes::complete::{tag, take_until},
    multi::many0,
    sequence::preceded,
    IResult,
};
use std::{fs, str::from_utf8};

use crate::{AMU, BOHR_RADIUS};

use super::{
    general::{parse_data, parse_f64, parse_multi_f64, parse_multi_u32, parse_record, parse_u32},
    END_CELL_GLOBAL, END_DUMP,
};

pub fn parse_castep_bin(data: &[u8]) -> IResult<&[u8], Vec<&[u8]>> {
    preceded(
        match_end_of_cell_global,
        preceded(tag(END_CELL_GLOBAL), many0(parse_record)),
    )(data)
}

fn match_end_of_param_dump(data: &[u8]) -> IResult<&[u8], &[u8]> {
    take_until(END_DUMP)(data)
}

pub fn skip_to_optimized_cell(data: &[u8]) -> IResult<&[u8], &[u8]> {
    preceded(match_end_of_param_dump, tag(END_DUMP))(data)
}

pub fn match_end_of_cell_global(data: &[u8]) -> IResult<&[u8], &[u8]> {
    take_until(END_CELL_GLOBAL)(data)
}

/**
Struct representing the data between
"BEGIN_UNIT_CELL" to "END_UNIT_CELL"
The fields here are not complete since
there are too many.
Can be added if needed.
*/
#[derive(Debug, Clone)]
pub struct UnitCell {
    // CASTEP uses Bohr instead of Angstrom. Conversion needed
    real_lattice_vectors: Matrix3<f64>,
    // CASTEP uses Bohr instead of Angstrom. Conversion needed
    recip_lattice_vectors: Matrix3<f64>,
    // CASTEP uses Bohr instead of Angstrom. Conversion needed
    volume: f64,
    // Number of species (elements) in cell
    num_species: u32,
    // Number of total ions (atoms) in cell
    num_ions: u32,
    /*
    Max number of each species' ions.
    Important: It defines the largest length
    of the ionic_positions array for each species.
    */
    max_ions_in_species: u32,
    /*
    An array of u32, ordered by appearance of species, each number
    tells the number of ions for the species.
    */
    num_ions_in_species: Vec<u32>,
    /*
    A 3-dimensional array.
    The first layer is nth-species.
    The second layeris nth-ions-of-species
    The third layer is xyz.
    */
    ionic_positions: Vec<Vec<Vec<f64>>>,
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
        num_ions_in_species: Vec<u32>,
        ionic_positions: Vec<Vec<Vec<f64>>>,
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
    kpoints: Vec<f64>,
    kpoint_weights: Vec<f64>,
    num_symmetry_operations: u32,
    num_crystal_symmetry_operations: u32,
    cell_symmetry_ion_constraints: f64,
    bs_kpoints: Vec<f64>,
    bs_kpoint_weights: Vec<f64>,
}

pub trait FromCastepCheck {
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
        let cell_real_lat: Matrix3<f64> = Matrix3::from_vec(real_lat);
        let cell_volume =
            parse_f64(cell_volume).unwrap().1 / BOHR_RADIUS / BOHR_RADIUS / BOHR_RADIUS;
        let num_species = parse_u32(num_species).unwrap().1;
        let num_ions = parse_u32(num_ions).unwrap().1;
        let max_ions_in_species = parse_u32(max_ions_in_species).unwrap().1;
        let num_ions_in_species = parse_multi_u32(num_ions_in_species).unwrap().1;
        let ionic_positions: Vec<f64> = parse_multi_f64(ionic_positions).unwrap().1;
        let ionic_positions: Vec<Vec<f64>> =
            ionic_positions.chunks(3).map(|xyz| xyz.to_vec()).collect();
        let ionic_positions: Vec<Vec<Vec<f64>>> = ionic_positions
            .chunks(max_ions_in_species as usize)
            .map(|pos_per_species| pos_per_species.to_vec())
            .collect();
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
// #[ignore]
