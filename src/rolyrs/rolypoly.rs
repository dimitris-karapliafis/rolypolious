use pyo3::prelude::*;

// Import our modules
mod translate;
mod dedup;
mod filter;
mod sample;
mod format;

// Expose the PyO3 modules
#[pymodule]
#[pyo3(name = "rolypoly")]
fn init_module(py: Python, m: &PyModule) -> PyResult<()> {
    // Create submodules
    let _sys = py.import("sys")?;

    let translate_module = PyModule::new(py, "translate")?;
    translate::seq_translate(py, translate_module)?;
    m.add_submodule(translate_module)?;

    let dedup_module = PyModule::new(py, "dedup")?;
    dedup::seq_dedup(py, dedup_module)?;
    m.add_submodule(dedup_module)?;

    let filter_module = PyModule::new(py, "filter")?;
    filter::seq_filter(py, filter_module)?;
    m.add_submodule(filter_module)?;

    let sample_module = PyModule::new(py, "sample")?;
    sample::seq_sample(py, sample_module)?;
    m.add_submodule(sample_module)?;

    let format_module = PyModule::new(py, "format")?;
    format::seq_format(py, format_module)?;
    m.add_submodule(format_module)?;

    // Add direct functions
    m.add_function(wrap_pyfunction!(translate::reverse_complement, m)?)?;
    
    Ok(())
}
