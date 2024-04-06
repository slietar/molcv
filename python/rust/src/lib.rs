use numpy::{ndarray::{Array2, ArrayView1, Axis}, IntoPyArray, PyArray2, PyArrayDyn, PyReadonlyArray1, PyReadonlyArray2};
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

        let mut output = Array2::<f32>::zeros((residue_atom_counts.len(), cutoffs.len()));

        for (cutoff_index, &cutoff) in cutoffs.iter().enumerate() {
            let cutoff_output_vec = engine.run(cutoff).await?;
            let cutoff_output_arr = ArrayView1::from(&cutoff_output_vec);

            let mut output_slice = output.index_axis_mut(Axis(1), cutoff_index);
            // let x = output_slice.as_slice_mut().unwrap();
            output_slice.assign(&cutoff_output_arr);
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

/// A Python module implemented in Rust.
#[pymodule]
fn _molcv(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(compute_cv, m)?)?;
    Ok(())
}
