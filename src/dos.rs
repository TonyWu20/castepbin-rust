extern crate itertools_num;
extern crate plotlib;
use std::f64::consts::{PI, SQRT_2};

use itertools_num::linspace;
use plotlib::{page::Page, style::LineStyle, view::ContinuousView};
const SMEARING_WIDTH: f64 = 0.2;
const HATREE_TO_EV: f64 = 27.211396641308;
// Compute total dos by given eigenvalues of kpoints and fermi energy
pub fn compute_dos(eigen_values: Vec<f64>, fermi_energy: f64) -> (Vec<f64>, Vec<f64>) {
    let energy_kn: Vec<f64> = eigen_values
        .clone()
        .iter()
        .map(|e| -> f64 { (e - fermi_energy) * HATREE_TO_EV })
        .collect();
    let e_min = -8.0;
    let e_max = 14.0;
    let n_per_ev = 100.0;
    let num_points = ((e_max - e_min) * n_per_ev) as i32;
    let energy_range: Vec<f64> = linspace::<f64>(e_min, e_max, num_points as usize)
        .into_iter()
        .collect();
    let mut dos = vec![0.0; num_points as usize];
    let all_delta_energies: Vec<Vec<f64>> = energy_kn
        .iter()
        .map(|de| e_delta(&energy_range, *de))
        .collect();
    all_delta_energies.iter().for_each(|v| {
        v.iter().enumerate().for_each(|(i, val)| {
            dos[i] += val;
        })
    });
    (energy_range, dos)
}
// e from given range, de from k-point eigenvalue
fn _e_delta(e: f64, de: f64) -> f64 {
    let _a: f64 = 2.0 * (2.0_f64.ln());
    let smear = SMEARING_WIDTH * 2.0 * _a.sqrt();
    let x = (-((e - de) / smear).powi(2)).exp() / (PI.sqrt() * smear);
    x
}
fn e_delta(energies: &Vec<f64>, de: f64) -> Vec<f64> {
    energies
        .iter()
        .map(|e| -> f64 { _e_delta(*e, de) })
        .collect()
}

#[test]
fn test_dos_compute() {
    let eigenvalues: Vec<f64> = vec![
        -0.36823620,
        -0.36749997,
        -0.14518335,
        -0.14469776,
        -0.13037202,
        -0.12981298,
        -0.07961317,
        -0.07470166,
    ];
    let fermi_energy: f64 = -0.137401;
    let dos_result = compute_dos(eigenvalues, fermi_energy);
    let (e, pe) = dos_result;
    let data: Vec<(f64, f64)> = e.into_iter().zip(pe).collect();
    let dos_plot = plotlib::repr::Plot::new(data).line_style(LineStyle::new().colour("#DD3355"));
    let v = ContinuousView::new()
        .add(dos_plot)
        .x_range(-16., 16.)
        .y_range(0., 3.5);
    Page::single(&v).save("dos.svg").unwrap();
}
