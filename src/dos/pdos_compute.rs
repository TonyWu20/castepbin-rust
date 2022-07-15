use ndarray::Array2;
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{parser::bands::Bands, util::ElementWiseAdd};

use super::dos_compute::{e_deltas_per_kpt, shift_eigen_by_e_fermi};

#[cfg(test)]
#[test]
fn test_pdos() {
    use std::fs;

    use itertools_num::linspace;
    use plotters::prelude::*;

    use crate::config::pdos_config::PDOSCommandGroup;
    use crate::config::{pdos_config, Config, TaskProcess};
    use crate::dos::dos_compute;
    use crate::{
        dos::dos_compute::compute_energy_range,
        parser::{
            castep_bin::{skip_to_optimized_cell, FromCastepCheck, UnitCell},
            pdos_weights::PDOSWeight,
        },
    };
    let config_file = "./Pt_311/Pt_311_12lyr_v20_CO.toml";
    let config: Config = toml::from_str(&fs::read_to_string(config_file).unwrap()).unwrap();
    let pdos_tasks = config.tasks().pdos().unwrap();
    let pdos_results = pdos_tasks.task_execute();
    pdos_results.iter().enumerate().for_each(|(i, dos)| {
        let file_name = format!("pdos_{i}.png");
        let root = BitMapBackend::new(&file_name, (1920, 1080)).into_drawing_area();
        root.fill(&WHITE);
        let root = root.margin(10, 10, 10, 10);
        let mut chart = ChartBuilder::on(&root)
            .caption(
                "Projected Density of States",
                ("sans-serif", 40).into_font(),
            )
            .build_cartesian_2d(-20f64..16f64, -2f64..2f64)
            .unwrap();
        chart.configure_mesh().draw().unwrap();
        (0..dos.num_projectors()).into_iter().for_each(|ith| {
            let (e, spin1) = (
                dos.energy_range().to_vec(),
                dos.get_projector_dos(ith, 0_u8).to_vec(),
            );
            let spin2: Vec<f64> = dos
                .get_projector_dos(ith, 1_u8)
                .iter()
                .map(|v| v * -1.0)
                .collect();
            let up: Vec<(f64, f64)> = e.clone().into_iter().zip(spin1).collect();
            let down: Vec<(f64, f64)> = e.into_iter().zip(spin2).collect();
            chart
                .draw_series(LineSeries::new(up.into_iter(), &RED))
                .unwrap();
            chart
                .draw_series(LineSeries::new(down.into_iter(), &GREEN))
                .unwrap();
        });
        root.present().unwrap();
    });

    // let total_dos = calculate_total_dos_from_bands(band_file, None, None);
    // let (e, spin1) = (
    //     total_dos.energy_range().to_vec(),
    //     total_dos.dos().0.to_vec(),
    // );
    // let spin2: Vec<f64> = total_dos
    //     .dos()
    //     .1
    //     .as_ref()
    //     .unwrap()
    //     .iter()
    //     .map(|x| -> f64 { x * -1.0 })
    //     .collect();
}
