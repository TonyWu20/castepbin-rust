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
