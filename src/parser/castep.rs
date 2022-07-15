use std::{error::Error, fs, path::Path};

use nom::{
    bytes::complete::{tag, take_until},
    character::complete::multispace0,
    multi::many0,
    sequence::preceded,
    IResult,
};

use super::general::float;

pub fn get_last_final_energy<P: AsRef<Path>>(castep_filepath: P) -> Result<f64, Box<dyn Error>> {
    let castep_text = fs::read_to_string(castep_filepath)?;
    // Parse all final energy values
    let (_, energies) = many0(parse_final_energy)(&castep_text).unwrap();
    // Return last one
    Ok(energies[energies.len() - 1])
}

fn parse_final_energy(text: &str) -> IResult<&str, f64> {
    preceded(
        // Skip contents before final energy
        take_until("Final energy, E             ="),
        preceded(
            // Skip text
            tag("Final energy, E             ="),
            // Parse float value
            preceded(multispace0, float),
        ),
    )(text)
}
#[cfg(test)]
#[test]
fn test_castep() {
    use nom::multi::many0;
    let final_energy = get_last_final_energy("Pt_311/Pt_311_12lyr_v20_CO.castep").unwrap();
    println!("{}", final_energy);
}
