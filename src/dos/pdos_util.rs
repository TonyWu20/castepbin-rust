use crate::util::ElementWiseAdd;
use std::fs;

use ndarray::prelude::*;
use ndarray::{ArrayBase, ArrayViewD, Dim, ShapeBuilder, ViewRepr};
use nom::{error::Error, IResult};

use crate::parser::{
    castep_bin::{FromCastepCheck, UnitCell},
    pdos_weights::PDOSWeight,
};

// pub fn projected_weights(
//     command: &PDOSCommand,
//     cell: &UnitCell,
//     pdos_weight: &PDOSWeight,
// ) -> Vec<Vec<UnitOrbitalWeight>> {
//     let species_lookup_table = cell.symbol_order_hashtable();
//     let species_id = species_lookup_table
//         .get(command.species.as_str())
//         .unwrap_or_else(|| panic!("Species {} is not in the cell!", command.species.as_str()));
//     let ang_moment = parse_ang_moment(command.am_channel.as_str()).expect("Wrong ang_moment");
//     let projections = pdos_weight
//         .filter_species_rank_am(*species_id, command.ion_in_species, ang_moment)
//         .unwrap();
//     todo!();
// if self.spin {
//     let upspin_weights: Vec<&[UnitOrbitalWeight]> = projections
//         .iter()
//         .map(|id| get_orbital_eigenweights(pdos_weight, *id as u32, 1))
//         .collect();
//     let downspin_weights: Vec<&[UnitOrbitalWeight]> = projections
//         .iter()
//         .map(|id| get_orbital_eigenweights(pdos_weight, *id as u32, 2))
//         .collect();
//     let merged_upspins = merge_pdos_weights(&upspin_weights);
//     let merged_downspins = merge_pdos_weights(&downspin_weights);
//     return vec![merged_upspins, merged_downspins];
// } else {
//     let weights: Vec<&[UnitOrbitalWeight]> = projections
//         .iter()
//         .map(|id| get_orbital_eigenweights(pdos_weight, *id as u32, 1))
//         .collect();
//     let merged_weights = merge_pdos_weights(&weights);
//     return vec![merged_weights];
// }
// }

fn read_pdos_weights(pdos_weights_file: &str) -> PDOSWeight {
    let file = fs::read(pdos_weights_file).expect("Failed to open pdos_weights file!");
    let (_, pdos_weights) = PDOSWeight::parse(&file).unwrap();
    pdos_weights
}

/**
Get the eigenweight at given `orbital_id` (0-indexed) and spin
*/
// pub fn get_orbital_eigenweights(
//     weights: &PDOSWeight,
//     orbital_id: u32,
//     spin: u32,
// ) -> ArrayBase<ViewRepr<&f64>, Dim<[usize; 2]>> {
//     let eigen_weights: &EigenWeightPerOrb = weights.eigen_weights_at_orbital(orbital_id);
//     match spin {
//         1 => eigen_weights.weight_matrices().slice(s![0, .., ..]),
//         2 => eigen_weights.weight_matrices().slice(s![1, .., ..]),
//         _ => panic!("Wrong spin number :{spin}"),
//     }
// }
// fn merge_pdos_weights(projections: &Vec<&[UnitOrbitalWeight]>) -> Vec<UnitOrbitalWeight> {
//     let merged_weights: Vec<UnitOrbitalWeight> = projections
//         .into_iter()
//         .map(|entry| entry.to_vec())
//         .reduce(
//             |orb_prev: Vec<UnitOrbitalWeight>,
//              orb_next: Vec<UnitOrbitalWeight>|
//              -> Vec<UnitOrbitalWeight> {
//                 orb_prev
//                     .into_iter()
//                     .zip(orb_next.into_iter())
//                     .map(|(vp, vn)| vp.add(&vn))
//                     .collect()
//             },
//         )
//         .unwrap();
//     merged_weights
// }

#[cfg(test)]
mod test {
    use std::{fs, ops::Add};

    use ndarray::{Array, Array1, Array2, Axis};
    use ndarray::{Order, ShapeBuilder};

    use crate::{parser::pdos_weights::PDOSWeight, util::ElementWiseAdd};

    use super::read_pdos_weights;
    // #[ignore]
    #[test]
    fn test_pdos_weights_struct() {
        let file = fs::read("./Pt_310_12lyr_v20_CO_DOS/Pt_310_12lyr_v20_CO_DOS.pdos_weights")
            .expect("Error opening pdos_weights");
        // let file = fs::read("./Si_2_custom CASTEP GeomOpt/Si_2_custom.pdos_weights")
        // .expec;warnt("Error opening pdos_weights");
        // println!("Pt_310");
        // println!("{:#?}", pdos_weight);
        // let file = fs::read("./Au_310/Au_310_12lyr_v20_CHO.pdos_weights")
        //     .expect("Error opening pdos_weights");
        // let file = fs::read("./Si_2_custom CASTEP GeomOpt/Si_2_custom.pdos_weights")
        //     .expect("Error opening pdos_weights");
        let (_, pdos_weight) = PDOSWeight::parse(&file).unwrap();
        let species = pdos_weight.species_no();
        let ranks = pdos_weight.rank_in_species();
        let am_channels = pdos_weight.am_channel();
        let orbitals_id: Vec<usize> = species
            .iter()
            .zip(ranks.iter())
            .zip(am_channels.iter())
            .enumerate()
            .filter(|(_, ((&spe, &rank), &am))| spe == 2_u8 && rank == 1_u8 && am == 1_u8)
            .map(|(i, _)| i)
            .collect();
        println!("orbitals: {:?}", orbitals_id);
        let projections = pdos_weight
            .orbital_at_eigen_weights()
            .select(Axis(3), &[0])
            .select(Axis(0), &orbitals_id);
        let merged_weights = projections.sum_axis(Axis(0));
        // let mut kpt_weights: Array<f64, _> = Array::zeros((4, 1));
        // kpt_weights.fill(0.25);
        println!(
            "{}, {}",
            pdos_weight.num_kpoints(),
            pdos_weight.num_eigenvalues(),
        );
        // println!("{:#.2?}", pdos_weight.orbital_at_eigen_weights());
        println!(
            "{:#.2?}, {}, {:?}",
            merged_weights,
            merged_weights.is_standard_layout(),
            merged_weights.dim()
        );
        // println!("{:#?}, {}", kpt_weights, kpt_weights.is_standard_layout());
        // println!("{:#?}", merged_weights.dot(&kpt_weights));
    }
}
