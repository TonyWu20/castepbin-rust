use std::fs;

use nom::{
    bytes::complete::tag,
    character::complete::{multispace0, multispace1},
    combinator::{map_res, recognize},
    multi::{count, many0, many1},
    sequence::{delimited, preceded, terminated, tuple},
    IResult,
};

use super::general::{decimal, decimal_u8, float};

/// Struct representing the <seed>.bands file structure from CASTEP.
#[derive(Debug)]
pub struct Bands {
    num_kpts: u8,
    num_spins: u8,
    ///! Length = num_spins
    num_electrons: Vec<f64>,
    ///! Length = num_spins
    num_eigenvalues: Vec<u32>,
    ///! Length = num_spins
    e_fermi: Vec<f64>,
    ///! flatten vector representing a 3x3 column-major matrix
    unit_cell_vectors: Vec<f64>,
    ///! Length = num_kpts
    kpt_eigen_energies_array: KPointEnergiesVec,
}

impl Bands {
    pub fn new(
        num_kpts: u8,
        num_spins: u8,
        num_electrons: Vec<f64>,
        num_eigenvalues: Vec<u32>,
        e_fermi: Vec<f64>,
        unit_cell_vectors: Vec<f64>,
        kpt_eigen_energies_array: KPointEnergiesVec,
    ) -> Self {
        Self {
            num_kpts,
            num_spins,
            num_electrons,
            num_eigenvalues,
            e_fermi,
            unit_cell_vectors,
            kpt_eigen_energies_array,
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
        let kpt_eigen_energies_array = KPointEnergiesVec::from_array_of_struct(&kpt_eigen_energies);
        Ok((
            i,
            Self::new(
                num_kpts,
                num_spins,
                num_electrons,
                num_eigenvalues,
                e_fermi,
                cell_vectors,
                kpt_eigen_energies_array,
            ),
        ))
    }
    fn parse_nkpts(data: &str) -> IResult<&str, u8> {
        preceded(
            tag("Number of k-points"),
            delimited(multispace0, decimal_u8, multispace0),
        )(data)
    }
    fn parse_nspins(data: &str) -> IResult<&str, u8> {
        preceded(
            tag("Number of spin components"),
            delimited(multispace0, decimal_u8, multispace0),
        )(data)
    }
    fn parse_nelectrons(data: &str, num_spins: u8) -> IResult<&str, Vec<f64>> {
        preceded(
            tag("Number of electrons"),
            preceded(
                multispace0,
                count(terminated(float, multispace0), num_spins as usize),
            ),
        )(data)
    }
    fn parse_eigen_values(data: &str, num_spins: u8) -> IResult<&str, Vec<u32>> {
        preceded(
            tag("Number of eigenvalues"),
            preceded(
                multispace0,
                count(terminated(decimal, multispace0), num_spins as usize),
            ),
        )(data)
    }
    fn parse_efermi(data: &str, num_spins: u8) -> IResult<&str, Vec<f64>> {
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

    pub fn num_kpts(&self) -> u8 {
        self.num_kpts
    }

    pub fn num_spins(&self) -> u8 {
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

    pub fn kpt_eigen_energies_array(&self) -> &KPointEnergiesVec {
        &self.kpt_eigen_energies_array
    }
}

/**
Struct representing every entry of k-point and eigenvalues
in <seed>.bands
*/
#[derive(Debug)]
pub struct KPointEnergies {
    nth_kpts: u8,
    ///! ID in .bands file
    kpt_coordinate: [f64; 3],
    kpt_weight: f64,
    ///! One or two subvectors depending on spin-polarised or not
    eigen_energies: KPointEigenValues,
}

impl KPointEnergies {
    fn new(
        nth_kpts: u8,
        kpt_coordinate: [f64; 3],
        kpt_weight: f64,
        eigen_energies: KPointEigenValues,
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
        if eigen_energies.len() == 2 {
            let kpt_eigen_val_data = KPointEigenValues::new(
                eigen_energies[0].to_vec(),
                Some(eigen_energies[1].to_vec()),
            );
            Ok((
                i,
                Self::new(nth_kpt, coordinates, weight, kpt_eigen_val_data),
            ))
        } else {
            let kpt_eigen_val_data = KPointEigenValues::new(eigen_energies[0].to_vec(), None);
            Ok((
                i,
                Self::new(nth_kpt, coordinates, weight, kpt_eigen_val_data),
            ))
        }
    }
    fn parse_kpoint_entry(data: &str) -> IResult<&str, (u8, [f64; 3], f64)> {
        let (i, res) = tuple((
            preceded(tuple((tag("K-point"), multispace0)), decimal_u8),
            preceded(multispace0, count(terminated(float, multispace1), 3)),
            terminated(float, multispace0),
        ))(data)?;
        let coordinate: [f64; 3] = res.1.try_into().unwrap_or_else(|v: Vec<f64>| {
            panic!("Expected a Vec of length 3 but it was {}", v.len())
        });
        Ok((i, (res.0, coordinate, res.2)))
    }
    fn parse_eigen_values(data: &str) -> IResult<&str, Vec<f64>> {
        let (i, _) = recognize(tuple((tag("Spin component "), decimal)))(data)?;
        many0(delimited(multispace0, float, multispace0))(i)
    }

    pub fn nth_kpts(&self) -> u8 {
        self.nth_kpts
    }
    fn kpt_coordinate(&self) -> [f64; 3] {
        self.kpt_coordinate
    }

    pub fn kpt_weight(&self) -> f64 {
        self.kpt_weight
    }

    fn eigen_energies(&self) -> &KPointEigenValues {
        &self.eigen_energies
    }
}

/**
Struct of Array representing the k-point-eigen-energies data
in <seed>.bands
*/
#[derive(Debug)]
pub struct KPointEnergiesVec {
    ///! ID in .bands file, also identifier for the other struct of array field
    nth_kpts: Vec<u8>,
    kpt_coordinate: KPointCoordVec,
    kpt_weight: Vec<f64>,
    ///! Length equal to nth_kpts.len()
    eigen_energies: Vec<KPointEigenValues>,
}

impl KPointEnergiesVec {
    pub fn new(
        nth_kpts: Vec<u8>,
        kpt_coordinate: KPointCoordVec,
        kpt_weight: Vec<f64>,
        eigen_energies: Vec<KPointEigenValues>,
    ) -> Self {
        Self {
            nth_kpts,
            kpt_coordinate,
            kpt_weight,
            eigen_energies,
        }
    }
    pub fn from_array_of_struct(kpe_array: &Vec<KPointEnergies>) -> Self {
        let mut nth_kpts_array: Vec<u8> = vec![];
        let (mut kpt_x, mut kpt_y, mut kpt_z) = (vec![], vec![], vec![]);
        let mut kpt_weight_array: Vec<f64> = vec![];
        let mut eigen_energies_array: Vec<KPointEigenValues> = vec![];
        kpe_array.iter().for_each(|item| {
            nth_kpts_array.push(item.nth_kpts());
            let cd = item.kpt_coordinate();
            kpt_x.push(cd[0]);
            kpt_y.push(cd[1]);
            kpt_z.push(cd[2]);
            kpt_weight_array.push(item.kpt_weight());
            eigen_energies_array.push(item.eigen_energies().to_owned());
        });
        Self::new(
            nth_kpts_array,
            KPointCoordVec {
                x: kpt_x,
                y: kpt_y,
                z: kpt_z,
            },
            kpt_weight_array,
            eigen_energies_array,
        )
    }

    pub fn nth_kpts(&self) -> &[u8] {
        self.nth_kpts.as_ref()
    }

    pub fn kpt_weight(&self) -> &[f64] {
        self.kpt_weight.as_ref()
    }

    pub fn kpt_coordinate_mut(&mut self) -> &mut KPointCoordVec {
        &mut self.kpt_coordinate
    }

    pub fn kpt_coordinate(&self) -> &KPointCoordVec {
        &self.kpt_coordinate
    }

    pub fn eigen_energies(&self) -> &[KPointEigenValues] {
        self.eigen_energies.as_ref()
    }
}

#[derive(Debug)]
pub struct KPointCoordVec {
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
}

impl KPointCoordVec {
    fn new(x: Vec<f64>, y: Vec<f64>, z: Vec<f64>) -> Self {
        Self { x, y, z }
    }

    fn x(&self) -> &[f64] {
        self.x.as_ref()
    }

    fn y(&self) -> &[f64] {
        self.y.as_ref()
    }

    fn z(&self) -> &[f64] {
        self.z.as_ref()
    }
}

/** k-point eigenvalue representation in <seed>.bands */
#[derive(Debug, Clone)]
pub struct KPointEigenValues {
    eigen_spin_1: Vec<f64>,
    eigen_spin_2: Option<Vec<f64>>,
}

impl KPointEigenValues {
    pub fn new(eigen_spin_1: Vec<f64>, eigen_spin_2: Option<Vec<f64>>) -> Self {
        Self {
            eigen_spin_1,
            eigen_spin_2,
        }
    }

    pub fn eigen_at_spin(&self, spin: u8) -> Option<&[f64]> {
        if spin == 1 {
            Some(self.eigen_spin_1())
        } else if spin == 2 {
            self.eigen_spin_2()
        } else {
            None
        }
    }

    pub fn eigen_spin_1(&self) -> &[f64] {
        self.eigen_spin_1.as_ref()
    }

    pub fn eigen_spin_2(&self) -> Option<&[f64]> {
        if let Some(spin2) = &self.eigen_spin_2 {
            Some(spin2.as_ref())
        } else {
            None
        }
    }
}
