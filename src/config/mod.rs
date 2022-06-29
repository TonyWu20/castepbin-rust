/**
Define and control Structs for deserialization and serialization
of the config.toml
*/
use std::{num::ParseIntError, str::FromStr, string::ParseError};

use serde::Deserialize;

use self::pdos_config::PDOSTask;
pub mod pdos_config;

/**
Config file struct for deserialization
# Field:
    * title: String,
    * tasks: Task - Struct of Task
*/
#[derive(Deserialize, Debug)]
struct Config {
    title: String,
    tasks: Task,
}
/**
Task of run.
# Field:
    * pdos: Option<PDOSConfig> - Optional field to store config of PDOS calculation
*/
#[derive(Deserialize, Debug)]
struct Task {
    pdos: Option<PDOSTask>,
}

#[cfg(test)]
#[test]
fn test_toml() {
    use std::fs;

    use crate::parser::castep_bin::{FromCastepCheck, UnitCell};

    let config: Config = toml::from_str(
        r#"title = "castepbin-rust config"

[tasks]
[tasks.pdos]
seed = "./Pt_311/Pt_311_12lyr_v20_CO"
spin = true
smearing_width = 0.2
[[tasks.pdos.targets]]
species = ["Pt", "Pt"]
rank = [4, 5]
am_channel = ["s", "p", "d"]
[[tasks.pdos.targets]]
species = ["C"]
rank = [1]
am_channel = ["s", "p"]
"#,
    )
    .unwrap();
    println!("{:?}", config);
    let cell_file = fs::read(&config.tasks.pdos.as_ref().unwrap().castep_bin_filename()).unwrap();
    let cell = UnitCell::parser(&cell_file).unwrap().1;
    println!("{:?}", config.tasks)
}
