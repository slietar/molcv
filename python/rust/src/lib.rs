use numpy::{IntoPyArray, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use pyo3::{exceptions::PyRuntimeError, prelude::*};
use std::error::Error;

use ::molcv::Engine;


/// Formats the sum of two numbers as string.
#[pyfunction]
#[pyo3(pass_module)]
fn compute_cv<'py>(
    py: &Bound<'py, PyModule>,
    atom_counts_arr: PyReadonlyArray1<'py, u32>,
    atoms_data_arr: PyReadonlyArray2<'py, f32>,
    cutoffs_arr: PyReadonlyArray1<'py, f32>,
) -> PyResult<Bound<'py, PyArray2<f32>>> {
    async fn run<'py>(
        m: &Bound<'py, PyModule>,
        residue_atom_counts: &[u32],
        atoms_data: &[f32],
        cutoffs: &[f32],
    ) -> Result<Bound<'py, PyArray2<f32>>, Box<dyn Error>> {
        let mut engine = Engine::new().await?;

        let result = engine.run(
            residue_atom_counts,
            atoms_data,
            ..,
            cutoffs,
        ).await?;

        Ok(result.into_pyarray_bound(m.py()))
    }

    match pollster::block_on(run(
        py,
        atom_counts_arr.as_slice()?,
        atoms_data_arr.as_slice()?,
        cutoffs_arr.as_slice()?,
    )) {
        Ok(x) => Ok(x),
        Err(e) => Err(PyErr::new::<PyRuntimeError, _>(e.to_string())),
    }
}

#[pymodule]
fn _molcv(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_cv, m)?)?;
    Ok(())
}
