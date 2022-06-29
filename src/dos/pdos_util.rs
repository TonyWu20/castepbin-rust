use std::fs;

use nom::{error::Error, IResult};

use crate::parser::{
    castep_bin::{FromCastepCheck, UnitCell},
    pdos_weights::{EigenWeightPerOrb, PDOSWeight, UnitOrbitalWeight},
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
pub fn get_orbital_eigenweights(
    weights: &PDOSWeight,
    orbital_id: u32,
    spin: u32,
) -> &[UnitOrbitalWeight] {
    let eigen_weights = weights.eigen_weights_at_orbital(orbital_id);
    match spin {
        1 => eigen_weights.upspin(),
        2 => eigen_weights.downspin().unwrap().as_ref(),
        _ => panic!("Wrong spin number :{spin}"),
    }
}
fn merge_pdos_weights(projections: &Vec<&[UnitOrbitalWeight]>) -> Vec<UnitOrbitalWeight> {
    let merged_weights: Vec<UnitOrbitalWeight> = projections
        .into_iter()
        .map(|entry| entry.to_vec())
        .reduce(
            |orb_prev: Vec<UnitOrbitalWeight>,
             orb_next: Vec<UnitOrbitalWeight>|
             -> Vec<UnitOrbitalWeight> {
                orb_prev
                    .into_iter()
                    .zip(orb_next.into_iter())
                    .map(|(vp, vn)| vp.add(&vn))
                    .collect()
            },
        )
        .unwrap();
    merged_weights
}

#[cfg(test)]
mod test {
    use std::{fs, ops::Add};

    use crate::{
        parser::pdos_weights::{PDOSWeight, UnitOrbitalWeight},
        util::ElementWiseAdd,
    };

    use super::{get_orbital_eigenweights, read_pdos_weights};
    #[ignore]
    #[test]
    fn test_pdos_weights_struct() {
        let file = fs::read("./Pt_310_12lyr_v20_CO_DOS/Pt_310_12lyr_v20_CO_DOS.pdos_weights")
            .expect("Error opening pdos_weights");
        let (_, pdos_weight) = PDOSWeight::parse(&file).expect("Error Pt");
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
        let projections: Vec<Vec<UnitOrbitalWeight>> = orbitals_id
            .iter()
            .map(|i| get_orbital_eigenweights(&pdos_weight, *i as u32, 1).to_vec())
            .collect();
        let merged_weights: Vec<UnitOrbitalWeight> = projections
            .into_iter()
            .reduce(|orb_prev, orb_next| -> Vec<UnitOrbitalWeight> {
                orb_prev
                    .into_iter()
                    .zip(orb_next.into_iter())
                    .map(|(vp, vn)| vp.add(&vn))
                    .collect()
            })
            .unwrap();
        merged_weights
            .iter()
            .for_each(|i| println!("{:?}", i.weights()));
        println!(
            "{}, {}, {}",
            pdos_weight.num_kpoints(),
            pdos_weight.num_eigenvalues(),
            merged_weights.len()
        );
    }
}
