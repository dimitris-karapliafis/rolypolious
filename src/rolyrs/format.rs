use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};

static LINE_WIDTH: AtomicUsize = AtomicUsize::new(60);

/// Set the line width for sequence formatting
#[pyfunction]
pub fn set_line_width(width: usize) {
    LINE_WIDTH.store(width, Ordering::SeqCst);
}

/// Format sequences with specified line width
#[pyfunction]
pub fn format_sequence(
    input_file: &str,
    output_file: &str,
) -> PyResult<()> {
    // Get current line width
    let line_width = LINE_WIDTH.load(Ordering::SeqCst);
    
    // Open input file
    let input_path = Path::new(input_file);
    let file = File::open(input_path)
        .map_err(|e| PyValueError::new_err(format!("Failed to open input file: {}", e)))?;
    
    // Open output file
    let output_path = Path::new(output_file);
    let output_file = File::create(output_path)
        .map_err(|e| PyValueError::new_err(format!("Failed to create output file: {}", e)))?;
    let mut writer = BufWriter::new(output_file);
    
    let reader = io::BufReader::new(file);
    let mut current_header = String::new();
    let mut current_sequence = String::new();
    
    for line in reader.lines() {
        let line = line.map_err(|e| PyValueError::new_err(format!("Error reading line: {}", e)))?;
        if line.starts_with('>') {
            // Write previous sequence if exists
            if !current_header.is_empty() {
                writeln!(writer, ">{}", current_header)
                    .map_err(|e| PyValueError::new_err(format!("Error writing header: {}", e)))?;
                
                // Write sequence with specified line width
                for chunk in current_sequence.as_bytes().chunks(line_width) {
                    writeln!(writer, "{}", std::str::from_utf8(chunk).unwrap())
                        .map_err(|e| PyValueError::new_err(format!("Error writing sequence: {}", e)))?;
                }
            }
            
            // Start new sequence
            current_header = line[1..].trim().to_string();
            current_sequence = String::new();
        } else {
            current_sequence.push_str(line.trim());
        }
    }
    
    // Write last sequence
    if !current_header.is_empty() {
        writeln!(writer, ">{}", current_header)
            .map_err(|e| PyValueError::new_err(format!("Error writing header: {}", e)))?;
        
        for chunk in current_sequence.as_bytes().chunks(line_width) {
            writeln!(writer, "{}", std::str::from_utf8(chunk).unwrap())
                .map_err(|e| PyValueError::new_err(format!("Error writing sequence: {}", e)))?;
        }
    }
    
    Ok(())
}

/// Python module for sequence formatting
#[pymodule]
pub fn seq_format(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(format_sequence, m)?)?;
    m.add_function(wrap_pyfunction!(set_line_width, m)?)?;
    Ok(())
} 