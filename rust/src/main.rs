use std::{error::Error, fs::File, ops::RangeBounds as _, path::Path};
use serde::Serialize;
use serde_pickle::SerOptions;

use molcv::Engine;


// async fn run_one() -> Result<(), Box<dyn Error>> {
//     let mut engine = Engine::new().await?;

//     for filename in &[
//         // "2h1l.pdb",
//         // "5j7o.pdb",
//         "../drive/FBN1_AlphaFold.pdb"
//     ] {
//         let path = Path::new(filename);

//         let output_dir_path = format!("output/{}", path.file_stem().unwrap().to_str().unwrap());
//         std::fs::create_dir_all(&output_dir_path)?;

//         let (mut structure, _) = pdbtbx::open(filename, pdbtbx::StrictnessLevel::Loose)
//             .map_err(|e| format!("Failed to open PDB: {:?}", e))?;

//         let residues = structure.residues().collect::<Vec<_>>();
//         let residue_range = ..;

//         engine.set_residues(&residues, &residue_range);

//         for &cutoff in &[
//             10.0,
//             20.0,
//             30.0,
//             40.0,
//             50.0,
//             60.0,
//             70.0,
//             80.0,
//             90.0,
//             100.0,
//         ] {
//             eprintln!("Processing {} with cutoff {} Å", filename, cutoff);

//             let cv = engine.run(cutoff).await?;

//             let mut current_relative_residue_index = 0usize;

//             for (residue_index, residue) in structure.residues_mut().enumerate() {
//                 let b_factor = if residue_range.contains(&residue_index) {
//                     current_relative_residue_index += 1;
//                     cv[current_relative_residue_index - 1]
//                 } else {
//                     0.0
//                 };

//                 for atom in residue.atoms_mut() {
//                     atom.set_b_factor(b_factor as f64)?;
//                 }
//             }

//             pdbtbx::save(&structure, &format!("{}/{:03}.pdb", output_dir_path, cutoff as u32), pdbtbx::StrictnessLevel::Medium)
//                 .map_err(|e| format!("Failed to save PDB: {:?}", e))?;
//         }
//     }

//     Ok(())
// }


#[derive(Debug, Serialize)]
struct Output {
    cutoffs: Vec<f64>,
    data: Vec<Vec<Vec<f32>>>,
}


// async fn run() -> Result<(), Box<dyn Error>> {
//     // run_one().await?;
//     // return Ok(());


//     let data = project_preprocessing::deserialize("../structure/output/data.pkl")?;

//     let cutoffs = [
//         10.0,
//         20.0,
//         30.0,
//         40.0,
//         50.0,
//         60.0,
//         70.0,
//         80.0,
//         90.0,
//         100.0,
//     ];

//     let mut engine = Engine::new().await?;
//     let mut cv_all = Vec::with_capacity(data.domains.len());

//     for (domain_index, domain) in data.domains.iter().enumerate() {
//         eprintln!("Processing {:?} {}", domain.kind, domain.number);

//         let (mut structure, _) = pdbtbx::open(&format!("../output/structures/alphafold-contextualized/{:04}.pdb", domain_index), pdbtbx::StrictnessLevel::Loose)
//             .map_err(|e| format!("Failed to open PDB: {:?}", e))?;

//         let residues = structure.residues().collect::<Vec<_>>();

//         let mut cv_domain = Vec::with_capacity(cutoffs.len());

//         let residue_start = domain.start_position - (if domain_index > 0 { data.domains[domain_index - 1].start_position } else { 1 });
//         let residue_end = residue_start + (domain.end_position - domain.start_position + 1);

//         let residue_range = residue_start..residue_end;

//         // eprintln!("{:?}", &residue_range);
//         // eprintln!("{:?}", &domain);
//         // eprintln!("{}", structure.residue_count());

//         engine.set_residues(&residues, &residue_range);

//         let output_dir_path = if domain_index < 5 {
//             let path = format!("output/domains/{:04}", domain_index);
//             std::fs::create_dir_all(&path)?;
//             Some(path)
//         } else {
//             None
//         };

//         for &cutoff in &cutoffs {
//             // eprintln!("Processing {:?} {} with cutoff {} Å", domain.kind, domain.number, cutoff);

//             let cv = engine.run(cutoff).await?;


//             if let Some(output_dir_path) = &output_dir_path {
//                 let mut current_relative_residue_index = 0usize;

//                 for (residue_index, residue) in structure.residues_mut().enumerate() {
//                     let b_factor = if residue_range.contains(&residue_index) {
//                         current_relative_residue_index += 1;
//                         cv[current_relative_residue_index - 1]
//                     } else {
//                         0.0
//                     };

//                     for atom in residue.atoms_mut() {
//                         atom.set_b_factor(b_factor as f64)?;
//                     }
//                 }

//                 pdbtbx::save(&structure, &format!("{}/{:03}.pdb", output_dir_path, cutoff as u32), pdbtbx::StrictnessLevel::Medium)
//                     .map_err(|e| format!("Failed to save PDB: {:?}", e))?;
//             }


//             cv_domain.push(cv);
//         }

//         if let Some(output_dir_path) = output_dir_path {
//             let mut script = String::new();

//             for cutoff in cutoffs {
//                 script += &format!("load {}/{:03}.pdb, C{:03}\n", output_dir_path, cutoff as u32, cutoff as u32);
//             }

//             script += "spectrum b\n";

//             std::fs::write(&format!("{}/script.pml", output_dir_path), script)?;
//         }

//         cv_all.push(cv_domain);
//     }

//     let mut writer = File::create("../output/cv.pkl")?;

//     serde_pickle::to_writer(&mut writer, &Output {
//         cutoffs: cutoffs.to_vec(),
//         data: cv_all,
//     }, SerOptions::new())?;

//     Ok(())
// }

fn main() {
    // pollster::block_on(run()).unwrap();
}
