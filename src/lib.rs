#![allow(clippy::too_many_arguments, dead_code, unused_imports)]
use std::{fs, str::from_utf8};

use na::{Matrix3, MatrixXx3, OMatrix, SMatrix, Vector3};
use nom::{
    bytes::complete::take_until,
    multi::{count, many0, many1},
    number::complete::{be_f64, be_u32},
    sequence::tuple,
    IResult,
};
use parser::{general::parse_multi_f64_from_record, *};

use crate::parser::{
    castep_bin::skip_to_optimized_cell,
    general::{
        parse_data, parse_f64, parse_multi_f64, parse_record, parse_u32, parse_u32_from_record,
        parse_u32_vec_from_record,
    },
};

extern crate nalgebra as na;
pub mod dos;
pub mod parser;
pub mod test;

pub const BOHR_RADIUS: f64 = 1.8897259886;
pub const AMU: f64 = 1822.8839;
