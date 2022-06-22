use std::fs;

use nom::{
    bytes::complete::tag,
    character::complete::{multispace0, multispace1},
    combinator::recognize,
    multi::{count, many0, many1},
    sequence::{delimited, preceded, terminated, tuple},
    IResult,
};

use super::general::{decimal, float};

/// Struct representing the <seed>.bands file structure from CASTEP.
#[derive(Debug)]
pub struct Bands {
    num_kpts: u32,
    num_spins: u32,
    ///! Length = num_spins
    num_electrons: Vec<f64>,
    ///! Length = num_spins
    num_eigenvalues: Vec<u32>,
    ///! Length = num_spins
    e_fermi: Vec<f64>,
    ///! flatten vector representing a 3x3 column-major matrix
    unit_cell_vectors: Vec<f64>,
    ///! Length = num_kpts
    kpt_eigen_energies: Vec<KPointEnergies>,
}

impl Bands {
    pub fn new(
        num_kpts: u32,
        num_spins: u32,
        num_electrons: Vec<f64>,
        num_eigenvalues: Vec<u32>,
        e_fermi: Vec<f64>,
        unit_cell_vectors: Vec<f64>,
        kpt_eigen_energies: Vec<KPointEnergies>,
    ) -> Self {
        Self {
            num_kpts,
            num_spins,
            num_electrons,
            num_eigenvalues,
            e_fermi,
            unit_cell_vectors,
            kpt_eigen_energies,
        }
    }
    pub fn parse(data: &str) -> IResult<&str, Self> {
        let (i, num_kpts) = Self::parse_nkpts(data)?;
        let (i, num_spins) = Self::parse_nspins(i)?;
        let (i, num_electrons) = Self::parse_nelectrons(i, num_spins)?;
        let (i, num_eigenvalues) = Self::parse_eigen_values(i, num_spins)?;
        let (i, e_fermi) = Self::parse_efermi(i, num_spins)?;
        let (i, cell_vectors) = Self::parse_unit_cell_vectors(i)?;
        let (i, kpt_eigen_energies) = many0(KPointEnergies::parse)(i)?;
        Ok((
            i,
            Self::new(
                num_kpts,
                num_spins,
                num_electrons,
                num_eigenvalues,
                e_fermi,
                cell_vectors,
                kpt_eigen_energies,
            ),
        ))
    }
    fn parse_nkpts(data: &str) -> IResult<&str, u32> {
        preceded(
            tag("Number of k-points"),
            delimited(multispace0, decimal, multispace0),
        )(data)
    }
    fn parse_nspins(data: &str) -> IResult<&str, u32> {
        preceded(
            tag("Number of spin components"),
            delimited(multispace0, decimal, multispace0),
        )(data)
    }
    fn parse_nelectrons(data: &str, num_spins: u32) -> IResult<&str, Vec<f64>> {
        preceded(
            tag("Number of electrons"),
            preceded(
                multispace0,
                count(terminated(float, multispace0), num_spins as usize),
            ),
        )(data)
    }
    fn parse_eigen_values(data: &str, num_spins: u32) -> IResult<&str, Vec<u32>> {
        preceded(
            tag("Number of eigenvalues"),
            preceded(
                multispace0,
                count(terminated(decimal, multispace0), num_spins as usize),
            ),
        )(data)
    }
    fn parse_efermi(data: &str, num_spins: u32) -> IResult<&str, Vec<f64>> {
        preceded(
            tag("Fermi energies (in atomic units)"),
            preceded(
                multispace0,
                count(terminated(float, multispace0), num_spins as usize),
            ),
        )(data)
    }
    fn parse_unit_cell_vectors(data: &str) -> IResult<&str, Vec<f64>> {
        delimited(
            tuple((tag("Unit cell vectors"), multispace0)),
            many0(preceded(multispace0, float)),
            multispace0,
        )(data)
    }

    pub fn num_kpts(&self) -> u32 {
        self.num_kpts
    }

    pub fn num_spins(&self) -> u32 {
        self.num_spins
    }

    pub fn num_electrons(&self) -> &[f64] {
        self.num_electrons.as_ref()
    }

    pub fn num_eigenvalues(&self) -> &[u32] {
        self.num_eigenvalues.as_ref()
    }

    pub fn e_fermi(&self) -> &[f64] {
        self.e_fermi.as_ref()
    }

    pub fn unit_cell_vectors(&self) -> &[f64] {
        self.unit_cell_vectors.as_ref()
    }

    pub fn kpt_eigen_energies(&self) -> &[KPointEnergies] {
        self.kpt_eigen_energies.as_ref()
    }
}

/// Struct representing the k-point and eigenvalues
#[derive(Debug)]
pub struct KPointEnergies {
    nth_kpts: u32,
    ///! ID in .bands file
    kpt_coordinate: Vec<f64>,
    kpt_weight: f64,
    ///! One or two subvectors depending on spin-polarised or not
    eigen_energies: Vec<Vec<f64>>,
}

impl KPointEnergies {
    fn new(
        nth_kpts: u32,
        kpt_coordinate: Vec<f64>,
        kpt_weight: f64,
        eigen_energies: Vec<Vec<f64>>,
    ) -> Self {
        Self {
            nth_kpts,
            kpt_coordinate,
            kpt_weight,
            eigen_energies,
        }
    }
    pub fn parse(data: &str) -> IResult<&str, Self> {
        let (i, (nth_kpt, coordinates, weight)) = Self::parse_kpoint_entry(data)?;
        let (i, eigen_energies) = many1(Self::parse_eigen_values)(i)?;
        Ok((i, Self::new(nth_kpt, coordinates, weight, eigen_energies)))
    }
    fn parse_kpoint_entry(data: &str) -> IResult<&str, (u32, Vec<f64>, f64)> {
        let (i, res) = tuple((
            preceded(tuple((tag("K-point"), multispace0)), decimal),
            preceded(multispace0, count(terminated(float, multispace1), 3)),
            terminated(float, multispace0),
        ))(data)?;
        Ok((i, (res.0, res.1, res.2)))
    }
    fn parse_eigen_values(data: &str) -> IResult<&str, Vec<f64>> {
        let (i, _) = recognize(tuple((tag("Spin component "), decimal)))(data)?;
        many0(delimited(multispace0, float, multispace0))(i)
    }

    pub fn nth_kpts(&self) -> u32 {
        self.nth_kpts
    }

    pub fn kpt_coordinate(&self) -> &[f64] {
        self.kpt_coordinate.as_ref()
    }

    pub fn kpt_weight(&self) -> f64 {
        self.kpt_weight
    }

    pub fn eigen_energies(&self) -> &[Vec<f64>] {
        self.eigen_energies.as_ref()
    }
}

#[test]
fn test_parse_bands() {
    let data =
        fs::read_to_string("./Pt_310_12lyr_v20_CO_DOS/Pt_310_12lyr_v20_CO_DOS.bands").unwrap();
    let (_, bands_data) = Bands::parse(&data).unwrap();
    assert_eq!(4, bands_data.num_kpts());
    assert_eq!(2, bands_data.num_spins());
    assert_eq!(vec![137.0, 113.0], bands_data.num_electrons());
    assert_eq!(vec![209, 209], bands_data.num_eigenvalues());
    assert_eq!(vec![-0.324481, -0.324481], bands_data.e_fermi());
    assert_eq!(4, bands_data.kpt_eigen_energies().len());
}
