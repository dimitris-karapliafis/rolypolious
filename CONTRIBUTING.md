# Contributing to RolyPoly

Contributions welcome! Whether it's bug fixes, new features, documentation improvements, packaging, tests, or reports of your exprience / resource usage for your samples - all help is appreciated. Pull requests or forks are the preferred way to contribute and will be considered, and you can also open issues for discussion or contact one of the developers directly.

## Project Roadmap & TODO List
Check out our [project roadmap and TODO list](https://docs.google.com/spreadsheets/d/1udNbxtK1QMfOhVgxHyhrgw7U1hHFeIazlcLM6VIcbJo/edit?gid=0#gid=0) to see what features and improvements are planned.

## Contribution guidelines
- **Primary Language**: Python >=3.9
- **Secondary Languages**: Some system calls to Bash are allowed.
- **Dependency Management**: via pixi (development)
  - Prefer using existing dependencies over adding new ones.
  - Avoid pandas, and use polars
  - Avoid biopython if possible, check if an existing feature is implemented in `src/rolypoly/utils/bio/*` or use polars-bio.
- **lazy / eager**: the CLI commands are lazy evaluted (see `/src/rolypoly/utils/lazy_group.py`), and need to be explictly added to the `src/rolypoly/rolypoly.py` file. This makes debugging/tracing slightly harder, but it also isolates the commands, so we can break one of them without worrying on it effect it may have on others.

## Code Organization
1. **File Structure**:
   - `src/rolypoly/utils/`: Utility functions and helpers
   - `src/rolypoly/commands/`: Command-line interface modules (using click). 
  - `rolypoly.utils.various` for general-purpose functions that don't fit into other categories (e.g. dataframe operations)
  - `rolypoly.utils.logging` for logging, configuration, output tracking etc

2. **Naming Conventions**:
   - CLI arguments: No positional arguments unless absolutly necceray. Instead, prefer 'decalerd' and explict named arguemnt. Must support both short and long options, e.g. `-s` and `--skip-existing`. Optionally, provide support for json file (`--config config.json`) or json string (`--override-params '{"skip_existing": true}`).
   - Functions and Internal variables: Snake case (e.g., `skip_existing`). Try and reuse variable names from other commands for the same purpose. Long descriptive names are ok.
   - Classes: PascalCase (though use classes sparingly).
   - Environment or Global variables: UPPERCASE or CamelCase.

3. **Temporary Files**:
   - Optionally, create temp directory (hidden argument in some commands `--temp-dir`, if not specified it's within user's output path).  
   - When done, move only final output files to user's output path, or rename the temp-dir if it's easier (same parent path maybe).
   - Try to clean up tmp files unless `--keep-tmp` flag is used.

## Testing & Benchmarking
1. **Testing**:
   - Add code `src/tests/*`
   - You can take a look at `/REDACTED_HPC_PATH/tests/rp_tests/` (on dori) which conntains example data for different commands, seperated by input size.
2. **Benchmarking**:
   - Use `/usr/bin/time` for resource monitoring. Alternatively, hyperfine is great too but. Ideallt - use SLURM and keep track of the job IDs for later analysis with seff/pyseff.

## **Note**
This project is governed under the LBNL IP office. By contributing, you agree that your contributions will be subject to the terms of the GPLv3 license.
