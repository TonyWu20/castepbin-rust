use std::fs;

use crate::parser::bands::Bands;

use super::{
    dos_compute::{compute_energy_range, total_dos},
    DOS,
};

pub fn calculate_total_dos_from_bands(
    bands_file: &str,
    energy_range_limit: Option<(f64, f64)>,
    e_num_per_ev: Option<u32>,
) -> DOS {
    let bands_file_string = fs::read_to_string(bands_file).expect("Failed to open .bands file!");
    let bands_data = Bands::parse(&bands_file_string)
        .expect("Failed to parse .bands!")
        .1;
    let (e_min, e_max) = energy_range_limit.unwrap_or((-25.00, 20.00));
    let num_per_ev = e_num_per_ev.unwrap_or(100);
    let num_points = ((e_max - e_min) * num_per_ev as f64) as u32;
    let energy_range = compute_energy_range(e_min, e_max, num_points);
    let total_dos = total_dos(&bands_data, &energy_range);
    DOS::new(bands_data.num_spins(), num_points, energy_range, total_dos)
}

#[cfg(test)]
#[test]
fn test_total_dos() {
    use std::fs;

    use itertools_num::linspace;
    use plotlib::{page::Page, repr::Plot, style::LineStyle, view::ContinuousView};

    use crate::dos::dos_compute;

    let band_file = "/Users/tonywu/Library/Mobile Documents/com~apple~CloudDocs/Programming/castepbin-rust/Pt_310_12lyr_v20_CO_DOS/Pt_310_12lyr_v20_CO_DOS.bands";
    let total_dos = calculate_total_dos_from_bands(band_file, None, None);
    let (e, spin1) = (
        total_dos.energy_range().to_vec(),
        total_dos.dos().0.to_vec(),
    );
    let spin2: Vec<f64> = total_dos
        .dos()
        .1
        .as_ref()
        .unwrap()
        .iter()
        .map(|x| -> f64 { x * -1.0 })
        .collect();
    let data: Vec<(f64, f64)> = e.clone().into_iter().zip(spin1).collect();
    let down: Vec<(f64, f64)> = e.into_iter().zip(spin2).collect();
    let up_dos_plot = plotlib::repr::Plot::new(data).line_style(LineStyle::new().colour("#DD3355"));
    let down_dos_plot =
        plotlib::repr::Plot::new(down).line_style(LineStyle::new().colour("#aaff55"));
    let zeros = plotlib::repr::Plot::new(
        linspace(-20., 16., total_dos.num_points() as usize)
            .into_iter()
            .collect::<Vec<f64>>()
            .into_iter()
            .zip(vec![0.0; total_dos.num_points() as usize])
            .collect(),
    )
    .line_style(LineStyle::new().colour("#000000"));
    let vertical_zero = Plot::new(
        vec![0.0; total_dos.num_points() as usize]
            .into_iter()
            .zip(linspace(-22.0, 22.0, total_dos.num_points() as usize))
            .collect(),
    )
    .line_style(LineStyle::new().colour("#000000"));
    let v = ContinuousView::new()
        .add(zeros)
        .add(vertical_zero)
        .add(up_dos_plot)
        .add(down_dos_plot)
        .x_range(-20., 16.)
        .x_max_ticks(73)
        .y_max_ticks(45);
    Page::single(&v)
        .dimensions(2160, 1080)
        .save("dos.svg")
        .unwrap();
}
