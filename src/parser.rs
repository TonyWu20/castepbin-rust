#![allow(dead_code, unused_imports)]

use std::{fs, str::from_utf8};

use nalgebra::{Matrix3, MatrixXx3};
use nom::{
    branch::alt,
    bytes::complete::{tag, take, take_until, take_while, take_while1},
    character::{
        complete::{char, multispace0, one_of},
        is_alphabetic, is_alphanumeric,
    },
    combinator::{opt, recognize},
    error::Error,
    multi::{count, length_data, many0, many1},
    number::complete::{
        be_f32, be_f64, be_i16, be_i32, be_i64, be_i8, be_u16, be_u32, be_u8, hex_u32, le_f32,
        le_i32,
    },
    sequence::{delimited, preceded, terminated, tuple},
    HexDisplay, IResult,
};

use crate::{AMU, BOHR_RADIUS};
extern crate nom;
const END_DUMP: &[u8] = &[
    0x45, 0x4e, 0x44, 0x5f, 0x50, 0x41, 0x52, 0x41, 0x4d, 0x45, 0x54, 0x45, 0x52, 0x53, 0x5f, 0x44,
    0x55, 0x4d, 0x50, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x00, 0x00,
    0x00, 0x1e,
];
const END_CELL_GLOBAL: &[u8] = &[
    0x45, 0x4e, 0x44, 0x5f, 0x43, 0x45, 0x4c, 0x4c, 0x5f, 0x47, 0x4c, 0x4f, 0x42, 0x41, 0x4c, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20,
    0x00, 0x00, 0x01, 0x00,
];

pub fn parse_record(data: &[u8]) -> IResult<&[u8], &[u8]> {
    // let (_, record_marker) = take(4 as u32)(data)?;
    // let (_, record_length) = be_u32(record_marker)?;
    // tuple((take(4 as u32), take(record_length), take(4 as u32)))(data)
    terminated(length_data(be_u32), take(4 as u32))(data)
}

fn parse_record_and_data(data: &[u8]) -> IResult<&[u8], Vec<&[u8]>> {
    count(parse_record, 2)(data)
}

pub fn parse_data(data: &[u8]) -> IResult<&[u8], &[u8]> {
    preceded(parse_record, parse_record)(data)
}

fn parse_castep_bin(data: &[u8]) -> IResult<&[u8], Vec<&[u8]>> {
    preceded(
        match_end_of_cell_global,
        preceded(tag(END_CELL_GLOBAL), many0(parse_record)),
    )(data)
}

pub fn parse_f64(data: &[u8]) -> IResult<&[u8], f64> {
    be_f64(data)
}
pub fn parse_f32(data: &[u8]) -> IResult<&[u8], f32> {
    be_f32(data)
}

pub fn parse_u32(data: &[u8]) -> IResult<&[u8], u32> {
    be_u32(data)
}
fn match_end_of_param_dump(data: &[u8]) -> IResult<&[u8], &[u8]> {
    take_until(END_DUMP)(data)
}

pub fn skip_to_optimized_cell(data: &[u8]) -> IResult<&[u8], &[u8]> {
    preceded(match_end_of_cell_global, tag(END_CELL_GLOBAL))(data)
}

fn match_end_of_cell_global(data: &[u8]) -> IResult<&[u8], &[u8]> {
    take_until(END_CELL_GLOBAL)(data)
}

pub fn parse_multi_f64(data: &[u8]) -> IResult<&[u8], Vec<f64>> {
    many0(be_f64)(data)
}
#[ignore]
#[test]
fn test_new_parser() {
    let file = fs::read("./Si2.castep_bin").unwrap();
    let (_, parse_result) = parse_castep_bin(&file).unwrap();
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
        if i == 28
            || i == 72
            || i == 78
            || i == 80
            || i == 94
            || i == 174
            || i == 100
            || i == 106
            || i == 214
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
            Err(e) => println!("{i}:{:x?}", header),
        };
    });
}

#[ignore]
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
            Err(e) => match parse_f32(data) {
                Ok(num) => {
                    println!("f32 {i}:{}", num.1);
                    return;
                }
                Err(e) => match parse_f64(data) {
                    Ok(num) => println!("f64: {i}:{}", num.1),
                    Err(e) => println!("{:?}, {:x?}", e, data),
                },
            },
        }
    })
}

// #[test]
// fn test_float() {
//     let file = fs::read("./Pt_311_12lyr_v20_CO.pdos_weights").unwrap();
//     let parse_result = many0(parse_record)(&file).unwrap();
//     println!("{:x?}", parse_result.0);
//     parse_result.1.iter().enumerate().for_each(|(i, record)| {
//         let (mark_st, data, mark_ed) = record;
//         let length = parse_u32(*mark_st).unwrap().1;
//         println!("{i} Length: {}", length);
//         if i == 7 {
//             println!("{:x?}", data);
//         }
//         match parse_multi_f64(data) {
//             Ok(float) => println!("64 {i}: {:#.7?}", float.1),
//             Err(_) => println!("{i}: {:#.7}", parse_u32(data).unwrap().1),
//         }
//     });
// }
