use numpy::{ndarray::{Array2, Axis}, IntoPyArray, PyArray2, PyReadonlyArray1, PyReadonlyArray2};
use std::error::Error;

use ::molcv::Engine;
use pyo3::{exceptions::PyRuntimeError, prelude::*};


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
        engine.set_residues(residue_atom_counts, atoms_data, &..);

        let mut output = Array2::<f32>::zeros((cutoffs.len(), residue_atom_counts.len()));

        for (cutoff_index, &cutoff) in cutoffs.iter().enumerate() {
            let mut output_arr_slice = output.index_axis_mut(Axis(0), cutoff_index);
            let output_slice = output_arr_slice.as_slice_mut().unwrap();
            engine.run(cutoff, output_slice).await?;
        }

        Ok(output.into_pyarray_bound(m.py()))
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
