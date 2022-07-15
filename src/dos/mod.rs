extern crate itertools_num;

pub mod dos_compute;
pub mod dos_util;
mod pdos_compute;
mod pdos_util;

#[derive(Debug)]
pub struct TDOS {
    spin: u8,
    num_points: u32,
    energy_range: Vec<f64>,
    dos: Vec<Vec<f64>>,
}

impl TDOS {
    pub fn new(spin: u8, num_points: u32, energy_range: Vec<f64>, dos: Vec<Vec<f64>>) -> Self {
        Self {
            spin,
            num_points,
            energy_range,
            dos,
        }
    }

    pub fn num_points(&self) -> u32 {
        self.num_points
    }

    pub fn energy_range(&self) -> &[f64] {
        self.energy_range.as_ref()
    }

    pub fn dos(&self) -> &[Vec<f64>] {
        self.dos.as_ref()
    }

    pub fn spin(&self) -> u8 {
        self.spin
    }
}

#[derive(Debug)]
pub struct PDOS {
    spin: u8,
    num_points: u32,
    num_projectors: usize,
    energy_range: Vec<f64>,
    dos: Vec<Vec<f64>>,
    species: Vec<String>,
    ranks: Vec<u8>,
    am_channels: Vec<String>,
}

impl PDOS {
    pub fn new(
        spin: u8,
        num_points: u32,
        num_projectors: usize,
        energy_range: Vec<f64>,
        dos: Vec<Vec<f64>>,
        species: Vec<String>,
        ranks: Vec<u8>,
        am_channels: Vec<String>,
    ) -> Self {
        Self {
            spin,
            num_points,
            num_projectors,
            energy_range,
            dos,
            species,
            ranks,
            am_channels,
        }
    }

    pub fn spin(&self) -> u8 {
        self.spin
    }

    pub fn num_points(&self) -> u32 {
        self.num_points
    }

    pub fn energy_range(&self) -> &[f64] {
        self.energy_range.as_ref()
    }

    pub fn dos(&self) -> &[Vec<f64>] {
        self.dos.as_ref()
    }
    /**
    Slice vector of dos by projector_id and spin channel.
    # Arguments
        - projector_id - 0th-indexed projector id
        - spin - 0 or 1
    */
    pub fn get_projector_dos(&self, projector_id: usize, spin: u8) -> &[f64] {
        &self.dos()[spin as usize][projector_id * self.num_points() as usize
            ..(projector_id + 1) * self.num_points() as usize]
    }

    pub fn num_projectors(&self) -> usize {
        self.num_projectors
    }
}
