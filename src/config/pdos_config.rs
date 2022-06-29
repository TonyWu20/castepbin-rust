use std::fs;

/// Struct for pdos config in .toml
use serde::Deserialize;

use crate::parser::castep_bin::{FromCastepCheck, UnitCell};
/**
Configs of a PDOS calculation task.
# Field:
    * seed: String - for looking up other files
    * spin: bool - Whether consider spin-polarised
    * smearing_width: f64  - in eV, 0.2 to be consistent with MS
    * targets: Vec<Targets>  - Vector of `Targets`
*/
#[derive(Deserialize, Debug)]
pub struct PDOSTask {
    seed: String,
    spin: bool,
    smearing_width: f64,
    targets: Vec<Targets>,
}

impl PDOSTask {
    pub fn seed(&self) -> &str {
        self.seed.as_ref()
    }

    pub fn castep_bin_filename(&self) -> String {
        let mut castep_seed = self.seed.clone();
        castep_seed.push_str(".castep_bin");
        castep_seed
    }

    pub fn bands_filename(&self) -> String {
        let mut band_seed = self.seed.clone();
        band_seed.push_str(".bands");
        band_seed
    }
    pub fn pdos_weight_filename(&self) -> String {
        let mut pdos_seed = self.seed.clone();
        pdos_seed.push_str(".pdos_weights");
        pdos_seed
    }

    pub fn spin(&self) -> bool {
        self.spin
    }

    pub fn smearing_width(&self) -> f64 {
        self.smearing_width
    }

    pub fn targets(&self) -> &[Targets] {
        self.targets.as_ref()
    }
    pub fn gen_commands(&self) -> PDOSCommandGroup {
        let bin_file = fs::read(self.castep_bin_filename()).expect("Reading .castep_bin error");
        let (_, cell) =
            UnitCell::parser(&bin_file).expect("Parsing .castep_bin for UnitCell error");
        let commands: Vec<PDOSCommand> =
            self.targets().iter().map(|t| t.translate(&cell)).collect();
        PDOSCommandGroup::new(commands)
    }
}

/**
Struct to hold targets in `PDOSConfig`. Fields are designed to be able to be
deserialized directly from `.toml`.
# Field:
    * species: Vec<String>, - vector of species symbols
    * rank: Vec<u8>, - rank in species.
    * am_channel: Vec<String> - "s", "p", "d", "f"
*/
#[derive(Deserialize, Debug)]
pub struct Targets {
    species: Vec<String>,
    rank: Vec<u8>,
    am_channel: Vec<String>,
}

impl Targets {
    pub fn new(species: Vec<String>, rank: Vec<u8>, am_channel: Vec<String>) -> Self {
        Self {
            species,
            rank,
            am_channel,
        }
    }

    pub fn species(&self) -> &[String] {
        self.species.as_ref()
    }

    pub fn rank(&self) -> &[u8] {
        self.rank.as_ref()
    }

    pub fn am_channel(&self) -> &[String] {
        self.am_channel.as_ref()
    }

    /**
    Translate representation in species and am_channel to pdos_weights data form
    # Args:
        - cell: &UnitCell
    */
    pub fn translate(&self, cell: &UnitCell) -> PDOSCommand {
        let species_table = cell.symbol_order_hashtable();
        let species_id: Vec<u8> = self
            .species()
            .iter()
            .map(|elm| -> &u8 {
                species_table
                    .get(elm.as_str())
                    .unwrap_or_else(|| panic!("Species {} is not in the cell!", elm.as_str()))
            })
            .copied()
            .collect();
        let am_channel_id: Vec<u8> = self
            .am_channel
            .iter()
            .map(|am| parse_ang_moment(am).unwrap())
            .collect();
        let mut pdos_command = PDOSCommand::new(vec![], vec![], vec![]);
        species_id
            .iter()
            .zip(self.rank().iter())
            .for_each(|(spe, rank)| {
                am_channel_id.iter().for_each(|id| {
                    pdos_command.species_mut().push(*spe);
                    pdos_command.rank_in_species_mut().push(*rank);
                    pdos_command.am_channel_mut().push(*id);
                });
            });
        pdos_command
    }
}

/**
Function to convert am_channel string to u8
# Arguments:
    * am: &str - String slice of the am_channel
*/
fn parse_ang_moment(am: &str) -> Result<u8, &'static str> {
    match am {
        "s" => Ok(0_u8),
        "p" => Ok(1_u8),
        "d" => Ok(2_u8),
        "f" => Ok(3_u8),
        _ => Err("Angular moment not one of s,p,d,f"),
    }
}

/**
PDOSCommand carries vector of species, ion_in_species and am_channel
# Fields:
  - species: Vec<u8>,
  - rank_in_species: Vec<u8>,
  - am_channel: Vec<u8>
*/
#[derive(Debug)]
pub struct PDOSCommand {
    species: Vec<u8>,
    rank_in_species: Vec<u8>,
    am_channel: Vec<u8>,
}

impl PDOSCommand {
    pub fn new(species: Vec<u8>, rank_in_species: Vec<u8>, am_channel: Vec<u8>) -> Self {
        Self {
            species,
            rank_in_species,
            am_channel,
        }
    }

    pub fn species_mut(&mut self) -> &mut Vec<u8> {
        &mut self.species
    }

    pub fn rank_in_species_mut(&mut self) -> &mut Vec<u8> {
        &mut self.rank_in_species
    }

    pub fn am_channel_mut(&mut self) -> &mut Vec<u8> {
        &mut self.am_channel
    }

    pub fn species(&self) -> &[u8] {
        self.species.as_ref()
    }

    pub fn rank_in_species(&self) -> &[u8] {
        self.rank_in_species.as_ref()
    }

    pub fn am_channel(&self) -> &[u8] {
        self.am_channel.as_ref()
    }
}
/**
Vector of PDOSCommand
# Field:
    * commands: Vec<PDOSCommand>
*/
#[derive(Debug)]
pub struct PDOSCommandGroup {
    commands: Vec<PDOSCommand>,
}

impl PDOSCommandGroup {
    pub fn new(commands: Vec<PDOSCommand>) -> Self {
        Self { commands }
    }

    pub fn commands(&self) -> &[PDOSCommand] {
        self.commands.as_ref()
    }
}
