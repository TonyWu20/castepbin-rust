use std::fs;

use nom::{bytes::complete::tag, multi::many0, sequence::preceded};

use crate::{
    castep_bin, from_utf8, parse_f64, parse_multi_f64, parse_record, parse_u32,
    parser::{
        castep_bin::match_end_of_cell_global,
        general::{parse_f32, parse_multi_u32},
        pdos_weights::PDOSWeight,
        END_CELL_GLOBAL,
    },
    Matrix3, MatrixXx3, AMU, BOHR_RADIUS,
};

#[cfg(test)]
#[test]
fn test_unit_cell() {
    use crate::parser::castep_bin::{skip_to_optimized_cell, FromCastepCheck, UnitCell};

    let file = fs::read("./Pt_311/Pt_311_12lyr_v20_CO.castep_bin").unwrap();
    let (cell_content, _) = skip_to_optimized_cell(&file).unwrap();
    let (_, cell) = UnitCell::parser(&cell_content).unwrap();
}
#[cfg(test)]
#[test]
fn test_pdos_weights_struct() {
    // let file = fs::read("./Pt_310_12lyr_v20_CO_DOS/Pt_310_12lyr_v20_CO_DOS.pdos_weights")
    // .expect("Error opening pdos_weights");
    let file = fs::read("./Si_2_custom CASTEP GeomOpt/Si_2_custom.pdos_weights")
        .expect("Error opening pdos_weights");
    let (_, pdos_weight) = PDOSWeight::parse(&file).unwrap();
    println!("{}", pdos_weight.orbital_eigen_weights().len());
    println!(
        "{:#?}",
        pdos_weight.orbital_eigen_weights()[0].eigen_weights_for_each_k()[0].weights()
    )
    // let file = fs::read("./Pt_310_12lyr_v20_CO_DOS/Pt_310_12lyr_v20_CO_DOS.pdos_weights")
    //     .expect("Error opening pdos_weights");
    // let (_, pdos_weight) = PDOSWeight::parse(&file).expect("Error Pt");
    // println!("Pt_310");
    // println!("{:#?}", pdos_weight);
}
#[ignore]
#[cfg(test)]
#[test]
fn test_new_parser() {
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
            println!("{:#?}", parse_multi_f64(&header).unwrap().1);
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
            parse_multi_f64(&header).unwrap().1.iter().for_each(|f| {
                if i == 28 {
                    println!("{}", f / AMU);
                }
                println!("{i}: {}", f);
            });

            return;
        }
        match from_utf8(&header) {
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
fn test_pdos_weights() {
    let file = fs::read("./Si_2_custom CASTEP GeomOpt/Si_2_custom.pdos_weights").unwrap();
    // let file = fs::read("./Pt_310_12lyr_v20_CO_DOS/Pt_310_12lyr_v20_CO_DOS.pdos_weights").unwrap();
    let parse_result = many0(parse_record)(&file).unwrap();
    parse_result.1.iter().enumerate().for_each(|(i, record)| {
        if i <= 6 || i == 8 || i == 9 {
            let int_vec = parse_multi_u32(record).unwrap().1;
            println!("{}, {:?}", int_vec.len(), int_vec);
            return;
        }
        if i == 7 || i == 430 {
            let (next, kpt_nth) = parse_u32(record).unwrap();
            let (_, kpt_vec) = parse_multi_f64(next).unwrap();
            println!("{i}: {}, {:?}", kpt_nth, kpt_vec);
            return;
        }
        if i == 8 || i == 18 {
            println!("{i}: {:x?}", record);
            return;
        }
        if i >= 10 && i < 18 || i >= 20 {
            println!("{i}: {:?}", parse_multi_f64(record).unwrap().1);
            return;
        }
        println!("{i}: {:x?}", record);
    });
}
#[ignore]
#[cfg(test)]
#[test]
fn test_parser() {
    let file = fs::read("./Si2.castep_bin").unwrap();
    let parse_result = preceded(
        match_end_of_cell_global,
        preceded(tag(END_CELL_GLOBAL), many0(parse_record)),
    )(&file)
    .unwrap();
    println!("{:x?}", parse_result.0);
    parse_result.1.iter().enumerate().for_each(|(i, record)| {
        let data = record;
        if i == 2 || i == 4 || i == 6 || i == 176 {
            let mut array = parse_multi_f64(data).unwrap().1;
            array.iter_mut().for_each(|val| {
                if i == 6 {
                    *val = *val / BOHR_RADIUS / BOHR_RADIUS;
                }
                *val = *val / BOHR_RADIUS;
            });
            println!("{:#.7?}", array);
            return;
        }
        if i == 18 || i == 26 || i == 42 || i == 44 || i == 116 {
            let array = parse_multi_f64(data).unwrap().1;
            println!("{:#.7?}", array);
            return;
        }
        if i == 28 {
            println!("{:x?}", data);
            println!("{}", parse_f64(data).unwrap().1 / AMU);
            return;
        }
        if i == 30 {
            println!("{}", from_utf8(data).unwrap());
            return;
        }
        if i == 8 || i == 10 || i == 12 || i == 34 || i == 48 || i == 114 {
            println!("{:x?}", data);
            println!("{}", parse_u32(data).unwrap().1);
            return;
        }
        match from_utf8(data) {
            Ok(str) => {
                println!(
                    "{i}:{}",
                    str.trim_end().trim_matches('\0').trim_end_matches('\n') // str
                );
                return;
            }
            Err(_) => match parse_f32(data) {
                Ok(num) => {
                    println!("f32 {i}:{}", num.1);
                    return;
                }
                Err(_) => match parse_f64(data) {
                    Ok(num) => println!("f64: {i}:{}", num.1),
                    Err(e) => println!("{:?}, {:x?}", e, data),
                },
            },
        }
    })
}
