use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use rand::seq::SliceRandom;
use rand::thread_rng;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;

/// Sample a specific number of sequences
#[pyfunction]
pub fn sample_sequences(
    input_file: &str,
    output_file: &str,
    n: usize,
) -> PyResult<usize> {
    // Read all sequences into memory
    let mut sequences = Vec::new();
    let mut current_header = String::new();
    let mut current_sequence = String::new();
    
    // Open the input file
    let input_path = Path::new(input_file);
    let file = File::open(input_path)
        .map_err(|e| PyValueError::new_err(format!("Failed to open input file: {}", e)))?;
    
    let reader = io::BufReader::new(file);
    for line in reader.lines() {
        let line = line.map_err(|e| PyValueError::new_err(format!("Error reading line: {}", e)))?;
        if line.starts_with('>') {
            if !current_header.is_empty() {
                sequences.push((current_header.clone(), current_sequence.clone()));
                current_sequence = String::new();
            }
            current_header = line[1..].trim().to_string();
        } else {
            current_sequence.push_str(line.trim());
        }
    }
    
    // Add the last sequence
    if !current_header.is_empty() {
        sequences.push((current_header, current_sequence));
    }
    
    // Sample sequences
    let mut rng = thread_rng();
    let sample_size = n.min(sequences.len());
    sequences.shuffle(&mut rng);
    let sampled = &sequences[..sample_size];
    
    // Write sampled sequences
    let output_path = Path::new(output_file);
    let output_file = File::create(output_path)
        .map_err(|e| PyValueError::new_err(format!("Failed to create output file: {}", e)))?;
    let mut writer = BufWriter::new(output_file);
    
    for (header, sequence) in sampled {
        writeln!(writer, ">{}", header)
            .map_err(|e| PyValueError::new_err(format!("Error writing header: {}", e)))?;
        
        // Write sequence with line breaks every 60 characters
        for chunk in sequence.as_bytes().chunks(60) {
            writeln!(writer, "{}", std::str::from_utf8(chunk).unwrap())
                .map_err(|e| PyValueError::new_err(format!("Error writing sequence: {}", e)))?;
        }
    }
    
    Ok(sample_size)
}

/// Sample a proportion of sequences
#[pyfunction]
pub fn sample_sequences_by_proportion(
    input_file: &str,
    output_file: &str,
    proportion: f64,
) -> PyResult<usize> {
    if proportion <= 0.0 || proportion > 1.0 {
        return Err(PyValueError::new_err("Proportion must be between 0 and 1"));
    }
    
    // Count total sequences first
    let input_path = Path::new(input_file);
    let file = File::open(input_path)
        .map_err(|e| PyValueError::new_err(format!("Failed to open input file: {}", e)))?;
    
    let reader = io::BufReader::new(file);
    let total_sequences = reader.lines()
        .filter(|line| line.as_ref().map_or(false, |l| l.starts_with('>')))
        .count();
    
    let n = (total_sequences as f64 * proportion).round() as usize;
    sample_sequences(input_file, output_file, n)
}

/// Python module for sequence sampling
#[pymodule]
pub fn seq_sample(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sample_sequences, m)?)?;
    m.add_function(wrap_pyfunction!(sample_sequences_by_proportion, m)?)?;
    Ok(())
} 