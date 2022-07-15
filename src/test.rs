use std::fs;

use nom::{bytes::complete::tag, multi::many0, sequence::preceded};

use crate::{
    castep_bin, from_utf8, parse_f64, parse_multi_f64, parse_record, parse_u32,
    parser::{
        castep_bin::match_end_of_cell_global,
        general::{parse_f32, parse_multi_u32},
        pdos_weights::PDOSWeight,
    },
    AMU, BOHR_RADIUS,
};

#[cfg(test)]
#[test]
fn test_unit_cell() {
    use crate::{
        dos::dos_compute::HATREE_TO_EV,
        parser::castep_bin::{parse_castep_bin, skip_to_optimized_cell, FromCastepCheck, UnitCell},
    };

    let file = fs::read("./Pt_311/Pt_311_12lyr_v20_CO.castep_bin").unwrap();
    // let (i, _cell) = UnitCell::parser(&file).unwrap();
    let (_, parse_result) = parse_castep_bin(&file).unwrap();
    parse_result.iter().enumerate().for_each(|(i, entry)| {
        let header = entry;
        if i == 2 {
            let vectors: Vec<f64> = parse_multi_f64(header)
                .unwrap()
                .1
                .iter()
                .map(|i| -> f64 { i / BOHR_RADIUS })
                .collect();
            println!("{:?}", vectors);
        }
        if i == 86 || i == 179 || i == 180 {
            println!("{:x?}", header);
        }
        if [8, 10, 12, 16, 20, 22, 82, 84, 104].contains(&i) {
            println!("{}", parse_u32(header).unwrap().1);
            return;
        }
        if [24, 26, 176, 177, 14].contains(&i) {
            if [176, 177].contains(&i) {
                println!("{i}, {:.12}", parse_f64(header).unwrap().1 * HATREE_TO_EV);
                return;
            }
            println!("{:.12}", parse_f64(header).unwrap().1);
            return;
        }
        if i == 34 || i == 76 {
            println!("{:}", parse_u32(header).unwrap().1);
            return;
        }
        if i == 178 {
            println!("{:#?}", parse_multi_f64(header).unwrap().1);
        }
        if i == 28
            || i == 72
            || i == 78
            || i == 80
            || i == 94
            || i == 174
            || i == 100
            || i == 106
            || i == 4
        {
            parse_multi_f64(header).unwrap().1.iter().for_each(|f| {
                if i == 28 {
                    println!("{}", f / AMU);
                }
                println!("{i}: {}", f);
            });

            return;
        }
        match from_utf8(header) {
            Ok(str) => println!(
                "{i}:{}",
                str.trim_end().trim_matches('\n').trim_matches('\0')
            ),
            Err(_) => println!("{i}:{:x?}", header),
        };
    });
}

#[ignore]
#[cfg(test)]
#[test]
fn test_new_parser() {
    use nalgebra::{Matrix3, MatrixXx3};

    let file = fs::read("./Si2.castep_bin").unwrap();
    let (_, parse_result) = crate::castep_bin::parse_castep_bin(&file).unwrap();
    parse_result.iter().enumerate().for_each(|(i, entry)| {
        let header = entry;
        if i == 2 {
            let vectors = parse_multi_f64(header)
                .unwrap()
                .1
                .iter()
                .map(|i| -> f64 { i / BOHR_RADIUS })
                .collect();
            println!("{:?}", Matrix3::from_vec(vectors));
        }
        if i == 18 {
            let vectors: Option<MatrixXx3<f64>> =
                Some(MatrixXx3::from_vec(parse_multi_f64(header).unwrap().1));
            println!("{:?}", vectors);
        }
        if i == 86 || i == 179 || i == 180 {
            println!("{:x?}", header);
        }
        if i == 22
            || i == 8
            || i == 10
            || i == 12
            || i == 16
            || i == 20
            || i == 82
            || i == 84
            || i == 104
            || i == 176
        {
            println!("{}", parse_u32(header).unwrap().1);
            return;
        }
        if i == 24 || i == 26 || i == 177 || i == 14 || i == 178 {
            println!("{:.12}", parse_f64(header).unwrap().1);
            return;
        }
        if i == 34 || i == 76 {
            println!("{:}", parse_u32(header).unwrap().1);
            return;
        }
        if i > 191 {
            println!("{:#?}", parse_multi_f64(header).unwrap().1);
        }
        if i == 28
            || i == 72
            || i == 78
            || i == 80
            || i == 94
            || i == 174
            || i == 100
            || i == 106
            || i == 4
        {
            parse_multi_f64(header).unwrap().1.iter().for_each(|f| {
                if i == 28 {
                    println!("{}", f / AMU);
                }
                println!("{i}: {}", f);
            });

            return;
        }
        match from_utf8(header) {
            Ok(str) => println!(
                "{i}:{}",
                str.trim_end().trim_matches('\n').trim_matches('\0')
            ),
            Err(_) => println!("{i}:{:x?}", header),
        };
    });
}
