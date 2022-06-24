use nom::{
    branch::alt,
    bytes::complete::take,
    character::complete::{char, multispace1, one_of},
    combinator::{map_res, opt, recognize},
    multi::{count, length_data, many0, many1},
    number::complete::{be_f32, be_f64, be_i32, be_u32, be_u8},
    sequence::{preceded, terminated, tuple},
    IResult,
};

pub fn parse_record(data: &[u8]) -> IResult<&[u8], &[u8]> {
    // let (_, record_marker) = take(4 as u32)(data)?;
    // let (_, record_length) = be_u32(record_marker)?;
    // tuple((take(4 as u32), take(record_length), take(4 as u32)))(data)
    terminated(length_data(be_u32), take(4_u32))(data)
}

fn parse_record_and_data(data: &[u8]) -> IResult<&[u8], Vec<&[u8]>> {
    count(parse_record, 2)(data)
}

pub fn parse_data(data: &[u8]) -> IResult<&[u8], &[u8]> {
    preceded(parse_record, parse_record)(data)
}
pub fn parse_u8_vec_from_u32_vec_record(data: &[u8]) -> IResult<&[u8], Vec<u8>> {
    let (i, u32_data) = parse_record(data)?;
    let (_, u32_array) = many0(parse_u32)(u32_data)?;
    let u8_vec = u32_array
        .into_iter()
        .map(|v| -> u8 {
            v.to_string()
                .parse::<u8>()
                .expect(&format!("Parse to u8 error {}", v))
        })
        .collect();
    Ok((i, u8_vec))
}
pub fn parse_u8_from_record(data: &[u8]) -> IResult<&[u8], u8> {
    let (i, u32_data) = parse_record(data)?;
    let (_, u32_val) = parse_u32(u32_data)?;
    let u8_val: u8 = u32_val.to_string().parse::<u8>().unwrap();
    Ok((i, u8_val))
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
pub fn decimal(input: &str) -> IResult<&str, u32> {
    map_res(recognize(many1(one_of("0123456789"))), |out: &str| {
        out.parse::<u32>()
    })(input)
}
pub fn decimal_u8(input: &str) -> IResult<&str, u8> {
    map_res(recognize(many1(one_of("0123456789"))), |out: &str| {
        out.parse::<u8>()
    })(input)
}
pub fn float(input: &str) -> IResult<&str, f64> {
    map_res(
        alt((
            // Case one: .42
            recognize(tuple((
                char('.'),
                decimal,
                opt(tuple((one_of("eE"), opt(one_of("+-")), decimal))),
            ))), // Case two: 42e42 and 42.42e42
            recognize(tuple((
                decimal,
                opt(preceded(char('.'), decimal)),
                one_of("eE"),
                opt(one_of("+-")),
                decimal,
            ))), // Case three: 42. and 42.42
            recognize(tuple((opt(one_of("+-")), decimal, char('.'), opt(decimal)))),
        )),
        |out: &str| out.parse::<f64>(),
    )(input)
}
