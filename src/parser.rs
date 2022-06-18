#![allow(dead_code, unused_imports)]

use std::{fs, str::from_utf8};

use nom::{
    branch::alt,
    bytes::complete::{tag, take, take_until, take_while, take_while1},
    character::{
        complete::{char, multispace0, one_of},
        is_alphabetic, is_alphanumeric,
    },
    combinator::{opt, recognize},
    error::Error,
    multi::{length_data, many0, many1},
    number::complete::{
        be_f32, be_f64, be_i16, be_i32, be_i64, be_i8, be_u16, be_u32, be_u8, hex_u32, le_i32,
    },
    sequence::{delimited, preceded, terminated, tuple},
    HexDisplay, IResult,
};
extern crate nom;
const END_DUMP: &[u8] = &[
    0x45, 0x4e, 0x44, 0x5f, 0x50, 0x41, 0x52, 0x41, 0x4d, 0x45, 0x54, 0x45, 0x52, 0x53, 0x5f, 0x44,
    0x55, 0x4d, 0x50, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x00, 0x00,
    0x00, 0x1e,
];

fn find_string_header(data: &[u8]) -> IResult<&[u8], &[u8]> {
    let (i, record_marker) = take(4 as usize)(data)?;
    // let (_, record_length) = be_u32(record_marker)?;
    terminated(take_until(record_marker), take(4 as usize))(i)
}

fn parse_record(data: &[u8]) -> IResult<&[u8], (&[u8], &[u8], &[u8])> {
    let (_, record_marker) = take(4 as u32)(data)?;
    let (_, record_length) = be_u32(record_marker)?;
    tuple((take(4 as u32), take(record_length), take(4 as u32)))(data)
}

fn parse_f64(data: &[u8]) -> IResult<&[u8], f64> {
    be_f64(data)
}
fn parse_f32(data: &[u8]) -> IResult<&[u8], f32> {
    be_f32(data)
}

fn parse_u32(data: &[u8]) -> IResult<&[u8], u32> {
    be_u32(data)
}
fn match_end_of_param_dump(data: &[u8]) -> IResult<&[u8], &[u8]> {
    take_until(END_DUMP)(data)
}

#[test]
fn test_parser() {
    let file = fs::read("./Si2.castep_bin").unwrap();
    let parse_result = preceded(
        match_end_of_param_dump,
        preceded(tag(END_DUMP), many0(parse_record)),
    )(&file)
    .unwrap();
    println!("{:x?}", parse_result.0);
    parse_result.1.iter().enumerate().for_each(|(i, record)| {
        let (mark_st, data, mark_ed) = record;
        let length = parse_u32(*mark_st).unwrap().1;
        if i == 28 {
            println!("{:x?}", data);
            println!(
                "{}, {}",
                parse_f32(data).unwrap().1,
                parse_f64(data).unwrap().1
            );
        }
        if i == 8 || i == 10 || i == 12 {
            println!("{:x?}", data);
            println!("{}", parse_u32(data).unwrap().1);
        }
        if length == 8 {
            match many0(parse_f64)(data) {
                Ok(num) => {
                    num.1.iter().for_each(|f| {
                        print!("{f}");
                    });
                    println!("");
                    return;
                }
                Err(e) => println!("Length:{}, {:?}, {:x?}", length, e, data),
            }
        }
        match from_utf8(data) {
            Ok(str) => {
                println!(
                    "{i}:{}",
                    str.trim_end().trim_matches('\0').trim_end_matches('\n')
                );
                return;
            }
            Err(e) => match parse_f32(data) {
                Ok(num) => {
                    println!("{i}:{}", num.1);
                    return;
                }
                Err(e) => match parse_f64(data) {
                    Ok(num) => println!("{i}:{}", num.1),
                    Err(e) => println!("Length:{}, {:?}, {:x?}", length, e, data),
                },
            },
        }
    })
}

// #[test]
// fn test_float() {
//     let unknown: &[u8] = &[0x03, 0xc2, 0xb8, 0xc3, 0x93, 0xc2, 0x93, 0xc2, 0x96];
//     let flow: &[u8] = &[0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x00, 0x1e];
//     let read = find_i32(flow).unwrap();
//     // let readf64: (&[u8], f64) = find_f64(flow).unwrap();
//     println!("{}", read.1);
//     let unknown_parse = find_float(unknown).unwrap();
//     println!("{}, {:?}", unknown_parse.1, unknown_parse.0);
//     // println!("{}", readf64.1);
// }
