# `castepbin-rust`
Parse `CASTEP` bin files with rust.

# Why? You have the "Powerful" Materials Studio!
- For the sake of **batch-processing, automation, and ergonomics**. 
- Save your wrists, arms, neck and shoulder from the **laborious mouse wandering** in GUI and the remaining thousands of same tasks for your data folder.
- Personal preference: I love to code instead of clicking mouse to test my health and endurance.
## Target
1. Parse essential informations from `.castep_bin` or `.check`
2. Parse data to calculate DOS (density-of-states) and PDOS (projected-DOS) from `.bands` and `.pdos_weights`
3. Automatic, batch-processing of chunks of `CASTEP` output results for data processing and plotting, e.g., computation of PDOS on target atom and visualization.

## Development stage
### 2022-06-18
Project Started.
### 2022-06-21
1. Primitive implementation of DOS calculation with gaussian smearing method.
2. Begin writing parser for `.pdos_weights`
### 2022-06-22
1. First version of codes relating to parsing `.pdos_weights`
2. First commit of `.bands` parser.
### 2022-06-23
1. First commit of calculation routine of total DOS from `.bands`
### 2022-06-24
1. Corrected used formula and smearing parameter for gaussian smearing. Consistent with CASTEP in MS.
2. Optimized parsed data structures for `.pdos_weights`, for better data acquisition based on orbital.

## Notes
I will share some essential information about the binary files from `CASTEP`.

## Credits
- [castepxbin](https://github.com/zhubonan/castepxbin)
    - Inspiration of this project.
- [Optados](https://github.com/optados-developers/optados)
    - Learn the file structure and data formats of `CASTEP` generated files.
- [DOS code from ASE](https://wiki.fysik.dtu.dk/ase/_modules/ase/dft/dos.html#DOS)
    - Implementation in python. Easy and clean to read and learn.
