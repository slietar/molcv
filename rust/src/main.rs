use std::error::Error;
use clap::Parser;

use molcv::Engine;


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// The path to the input .pdb file.
    pdb_input_path: String,

    /// The path at which to save computed CV as the B factor of a .pdb file.
    #[arg(long)]
    pdb_output_path: Option<String>,

    /// The path at which to save computed CV as a .npy file.
    #[arg(long)]
    data_output_path: Option<String>,

    /// The cutoffs at which to compute the CV.
    #[arg(long)]
    cutoff: Vec<f32>,
}


async fn run() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let (mut structure, _) = pdbtbx::open(&args.pdb_input_path, pdbtbx::StrictnessLevel::Loose)
        .map_err(|errors| format!("Failed to open PDB: {:?}", errors))?;

    let residue_atom_counts = structure
        .residues()
        .map(|residue| residue.atoms().count() as u32)
        .collect::<Vec<_>>();

    let atoms_data = structure
        .atoms()
        .flat_map(|atom| [
            atom.x() as f32,
            atom.y() as f32,
            atom.z() as f32,
            0.0
        ])
        .collect::<Vec<_>>();

    let mut engine = Engine::new().await?;

    engine.set_residues(&residue_atom_counts, &atoms_data, &..);

    let result = engine.run_return(args.cutoff[0]).await?;

    if let Some(pdb_output_path) = &args.pdb_output_path {
        if args.cutoff.len() != 1 {
            return Err("Only one cutoff is supported when saving as a PDB file".into());
        }

        for (residue_index, residue) in structure.residues_mut().enumerate() {
            for atom in residue.atoms_mut() {
                atom.set_b_factor(result[residue_index] as f64)?;
            }
        }

        pdbtbx::save(&structure, pdb_output_path, pdbtbx::StrictnessLevel::Medium)
            .map_err(|errors| format!("Failed to save PDB: {:?}", errors))?;
    }

    Ok(())
}

fn main() {
    if let Err(err) = pollster::block_on(run()) {
        eprintln!("{}", err);
        std::process::exit(1);
    }
}
