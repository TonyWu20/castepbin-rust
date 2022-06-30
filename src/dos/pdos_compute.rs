use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};

use crate::{parser::bands::Bands, util::ElementWiseAdd};

use super::dos_compute::{e_deltas_per_kpt, shift_eigen_by_e_fermi};

pub fn pdos(
    band_data: &Bands,
    energy_range: &[f64],
    projected_weights: &Vec<Vec<UnitOrbitalWeight>>,
) -> (Vec<f64>, Option<Vec<f64>>) {
    if band_data.num_spins() == 1 {
        (
            pdos_spin(band_data, energy_range, 1_u8, &projected_weights[0]),
            None,
        )
    } else {
        assert_eq!(
            2,
            band_data.num_spins(),
            "Spin cannot be > 2! Current value: {}",
            band_data.num_spins()
        );
        (
            pdos_spin(band_data, energy_range, 1_u8, &projected_weights[0]),
            Some(pdos_spin(
                band_data,
                energy_range,
                2_u8,
                &projected_weights[1],
            )),
        )
    }
}

fn pdos_spin(
    band_data: &Bands,
    energy_range: &[f64],
    spin: u8,
    orbital_weights: &[UnitOrbitalWeight],
) -> Vec<f64> {
    /*
    vec of energy eigenvalues for every k-point
    length = number of k-points
    */
    let kpt_ne_vec = band_data.kpt_eigen_energies_array().eigen_energies();
    let e_fermi = band_data.e_fermi()[(spin - 1) as usize];
    let current_spin_shifted_eigens: Vec<Vec<f64>> = kpt_ne_vec
        .iter()
        .map(|kpt_ne| -> Vec<f64> {
            let eigen_at_spin = kpt_ne
                .eigen_at_spin(spin)
                .unwrap_or_else(|| panic!("No eigenvalues at given spin {}", spin));
            shift_eigen_by_e_fermi(eigen_at_spin, e_fermi)
        })
        .collect();
    let k_point_weight_array = band_data.kpt_eigen_energies_array().kpt_weight();
    current_spin_shifted_eigens
        .par_iter()
        // Pairs with k-point weights
        .zip(k_point_weight_array.par_iter())
        .zip(orbital_weights.par_iter())
        .map(|((de_list_per_k, kpt_weight), orb_weight)| -> Vec<f64> {
            // E_delta over given energy range
            // Total DOS so orbital_weights are `None`
            e_deltas_per_kpt(
                energy_range,
                de_list_per_k,
                *kpt_weight,
                Some(orb_weight.weights()),
            )
        })
        // Fold vector of vectors into one vector
        .reduce_with(|prev, next| -> Vec<f64> { prev.add(&next) })
        .unwrap()
}
#[cfg(test)]
#[test]
fn test_pdos() {
    use std::fs;

    use itertools_num::linspace;
    use plotlib::{page::Page, repr::Plot, style::LineStyle, view::ContinuousView};

    use crate::{
        dos::{
            dos_compute::{self, compute_energy_range},
            pdos_util::{PDOSCommand, PDOSCommandGroup, PDOSComputeConfig},
        },
        parser::{
            castep_bin::{skip_to_optimized_cell, FromCastepCheck, UnitCell},
            pdos_weights::PDOSWeight,
        },
    };

    let seed = "./Pt_311/Pt_311_12lyr_v20_CO";
    let pdos_command = PDOSCommand::new("Pt", 4_u8, "d");
    let commands = PDOSCommandGroup::new(vec![pdos_command]);
    let config = PDOSComputeConfig::new(seed, vec![commands], true);
    let cell_path = vec![seed, ".castep_bin"];
    let weight_path = vec![seed, ".pdos_weights"];
    let band_path = vec![seed, ".bands"];
    let cell_file = fs::read(cell_path.join("")).unwrap();
    let (cell_content, _) = skip_to_optimized_cell(&cell_file).unwrap();
    let cell = UnitCell::parser(cell_content).unwrap().1;
    let weight_file = fs::read(weight_path.join("")).unwrap();
    let pdos_weight = PDOSWeight::parse(&weight_file).unwrap().1;
    let band_file = fs::read_to_string(band_path.join("")).unwrap();
    let band_data = Bands::parse(&band_file).unwrap().1;
    let projected_weights = config.projected_weights(
        &config.command_groups()[0].commands()[0],
        &cell,
        &pdos_weight,
    );
    let e_min = -25.0;
    let e_max = 16.0;
    let num_points = ((e_max - e_min) * 100.0) as u32;
    let energy_range = compute_energy_range(-25.0, 16.0, num_points);
    let pdos_result = pdos(&band_data, &energy_range, &projected_weights);
    let (e, spin1) = (energy_range.clone(), pdos_result.0);
    let spin2: Vec<f64> = pdos_result.1.unwrap().iter().map(|x| *x * -1.0).collect();

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
    let data: Vec<(f64, f64)> = e.clone().into_iter().zip(spin1).collect();
    let down: Vec<(f64, f64)> = e.into_iter().zip(spin2).collect();
    let up_dos_plot = plotlib::repr::Plot::new(data).line_style(LineStyle::new().colour("#DD3355"));
    let down_dos_plot =
        plotlib::repr::Plot::new(down).line_style(LineStyle::new().colour("#aaff55"));
    let zeros = plotlib::repr::Plot::new(
        linspace(-26., 20., num_points as usize)
            .into_iter()
            .collect::<Vec<f64>>()
            .into_iter()
            .zip(vec![0.0; num_points as usize])
            .collect(),
    )
    .line_style(LineStyle::new().colour("#000000"));
    let vertical_zero = Plot::new(
        vec![0.0; num_points as usize]
            .into_iter()
            .zip(linspace(-1.15, 1.15, num_points as usize))
            .collect(),
    )
    .line_style(LineStyle::new().colour("#000000"));
    let v = ContinuousView::new()
        .add(zeros)
        .add(vertical_zero)
        .add(up_dos_plot)
        .add(down_dos_plot)
        .x_range(-26., 20.)
        .x_max_ticks(24)
        .y_max_ticks(31);
    Page::single(&v)
        .dimensions(2160, 1080)
        .save("pdos.svg")
        .unwrap();
}
