use std::{collections::HashSet, error::Error};
use clap::Parser;

use crate::Engine;


#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// The chain to use for computing CV. If empty, all chains are used.
    #[arg(long)]
    chain: Vec<String>,

    /// The path to the input .pdb file.
    pdb_input_path: String,

    /// The path at which to save computed CV as the B factor of a .pdb file.
    #[arg(long)]
    pdb_output_path: Option<String>,

    /// The path at which to save computed CV as a .npy file.
    #[arg(long)]
    npy_output_path: Option<String>,

    /// The cutoffs at which to compute the CV.
    #[arg(long)]
    cutoff: Vec<f32>,
}


async fn run(args: &[String]) -> Result<(), Box<dyn Error>> {
    let args = Args::parse_from(args);

    let (mut structure, _) = pdbtbx::open(&args.pdb_input_path, pdbtbx::StrictnessLevel::Loose)
        .map_err(|errors| format!("Failed to open PDB: {:?}", errors))?;

    let chain_ids = if args.chain.is_empty() {
        None
    } else {
        let structure_chain_ids = structure.chains().map(|chain| chain.id()).collect::<HashSet<_>>();

        for chain_id in &args.chain {
            if !structure_chain_ids.contains(chain_id.as_str()) {
                return Err(format!("Chain {} not found", chain_id).into());
            }
        }

        Some(HashSet::<String>::from_iter(args.chain))
    };

    let is_chain_visible = |chain_id: &str| match &chain_ids {
        Some(chain_ids) => chain_ids.contains(chain_id),
        None => true,
    };

    let residues = structure
        .chains()
        .filter(|&chain| is_chain_visible(chain.id()))
        .flat_map(|chain| chain.residues())
        .collect::<Vec<_>>();

    let residue_atom_counts = residues
        .iter()
        .map(|&residue| residue.atoms().count() as u32)
        .collect::<Vec<_>>();

    let atoms_data = residues
        .iter()
        .flat_map(|&residue| residue.atoms())
        .flat_map(|atom| [
            atom.x() as f32,
            atom.y() as f32,
            atom.z() as f32,
            0.0
        ])
        .collect::<Vec<_>>();

    let mut engine = Engine::new().await?;

    let result = engine.run(
        &residue_atom_counts,
        &atoms_data,
        ..,
        &args.cutoff
    ).await?;

    if let Some(pdb_output_path) = &args.pdb_output_path {
        if args.cutoff.len() != 1 {
            return Err("Exactly one cutoff is supported when saving as a PDB file".into());
        }

        let mut residue_index = 0;

        for chain in structure.chains_mut() {
            let chain_visible = is_chain_visible(chain.id());

            for residue in chain.residues_mut() {
                let b_factor = if chain_visible {
                    let cv = result[[0, residue_index]];
                    residue_index += 1;

                    cv
                } else {
                    0.0
                };

                for atom in residue.atoms_mut() {
                    atom.set_b_factor(b_factor.into())?;
                }
            }
        }

        pdbtbx::save(&structure, pdb_output_path, pdbtbx::StrictnessLevel::Medium)
            .map_err(|errors| format!("Failed to save PDB: {:?}", errors))?;
    }

    if let Some(data_output_path) = &args.npy_output_path {
        ndarray_npy::write_npy(data_output_path, &result)?;
    }

    Ok(())
}


pub fn cli(args: &[String]) {
    if let Err(err) = pollster::block_on(run(args)) {
        eprintln!("Error: {}", err);
        std::process::exit(1);
    }
}
