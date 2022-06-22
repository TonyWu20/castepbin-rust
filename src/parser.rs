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
    combinator::{map, opt, recognize},
    error::Error,
    multi::{count, length_data, many0, many1},
    number::complete::{
        be_f32, be_f64, be_i16, be_i32, be_i64, be_i8, be_u16, be_u32, be_u8, hex_u32, le_f32,
        le_i32,
    },
    sequence::{delimited, preceded, terminated, tuple},
    HexDisplay, IResult,
};

use crate::{
    parser::{
        castep_bin::match_end_of_cell_global,
        general::{parse_f32, parse_f64, parse_i32, parse_multi_f64, parse_u32},
    },
    AMU, BOHR_RADIUS,
};

use self::general::parse_record;
extern crate nom;
const END_DUMP: &[u8] = &[
    0x45, 0x4e, 0x44, 0x5f, 0x50, 0x41, 0x52, 0x41, 0x4d, 0x45, 0x54, 0x45, 0x52, 0x53, 0x5f, 0x44,
    0x55, 0x4d, 0x50, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x00, 0x00,
    0x00, 0x1e,
];
pub const END_CELL_GLOBAL: &[u8] = &[
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

pub mod general {

    use nom::{
        bytes::complete::take,
        multi::{count, length_data, many0},
        number::complete::{be_f32, be_f64, be_i32, be_u32},
        sequence::{preceded, terminated},
        IResult,
    };

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
    pub fn parse_u32_from_record(data: &[u8]) -> IResult<&[u8], u32> {
        let (i, u32_data) = parse_record(data)?;
        Ok((i, parse_u32(u32_data).unwrap().1))
    }
    pub fn parse_u32_vec_from_record(data: &[u8]) -> IResult<&[u8], Vec<u32>> {
        let (i, u32_data) = parse_record(data)?;
        Ok((i, many0(parse_u32)(u32_data).unwrap().1))
    }
    pub fn parse_f64(data: &[u8]) -> IResult<&[u8], f64> {
        be_f64(data)
    }
    pub fn parse_multi_f64(data: &[u8]) -> IResult<&[u8], Vec<f64>> {
        many0(be_f64)(data)
    }
    pub fn parse_multi_f64_from_record(data: &[u8]) -> IResult<&[u8], Vec<f64>> {
        let (i, f64_data) = parse_record(data)?;
        Ok((i, parse_multi_f64(f64_data).unwrap().1))
    }
    pub fn parse_f32(data: &[u8]) -> IResult<&[u8], f32> {
        be_f32(data)
    }
    pub fn parse_u32(data: &[u8]) -> IResult<&[u8], u32> {
        be_u32(data)
    }
    pub fn parse_multi_u32(data: &[u8]) -> IResult<&[u8], Vec<u32>> {
        many0(be_u32)(data)
    }
    pub fn parse_i32(data: &[u8]) -> IResult<&[u8], i32> {
        be_i32(data)
    }
}

pub mod castep_bin {
    use nom::{
        bytes::complete::{tag, take_until},
        multi::many0,
        sequence::preceded,
        IResult,
    };

    use super::{general::parse_record, END_CELL_GLOBAL, END_DUMP};

    pub fn parse_castep_bin(data: &[u8]) -> IResult<&[u8], Vec<&[u8]>> {
        preceded(
            match_end_of_cell_global,
            preceded(tag(END_CELL_GLOBAL), many0(parse_record)),
        )(data)
    }

    fn match_end_of_param_dump(data: &[u8]) -> IResult<&[u8], &[u8]> {
        take_until(END_DUMP)(data)
    }

    pub fn skip_to_optimized_cell(data: &[u8]) -> IResult<&[u8], &[u8]> {
        preceded(match_end_of_cell_global, tag(END_CELL_GLOBAL))(data)
    }

    pub fn match_end_of_cell_global(data: &[u8]) -> IResult<&[u8], &[u8]> {
        take_until(END_CELL_GLOBAL)(data)
    }
}
