extern crate itertools_num;
extern crate plotlib;

pub mod dos_compute;
pub mod dos_util;
mod pdos_compute;
mod pdos_util;

#[derive(Debug)]
pub struct DOS {
    spin: u8,
    num_points: u32,
    energy_range: Vec<f64>,
    dos: (Vec<f64>, Option<Vec<f64>>),
}

impl DOS {
    pub fn new(
        spin: u8,
        num_points: u32,
        energy_range: Vec<f64>,
        dos: (Vec<f64>, Option<Vec<f64>>),
    ) -> Self {
        Self {
            spin,
            num_points,
            energy_range,
            dos,
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

    pub fn dos(&self) -> &(Vec<f64>, Option<Vec<f64>>) {
        &self.dos
    }
}
