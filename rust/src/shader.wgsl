struct Residue {
    atom_count: u32,
    atom_offset: u32,
}

struct Settings {
    atom_count: u32,
    residue_count: u32,
    residue_start: u32,
    cutoffs: array<f32>,
}


@group(0)
@binding(0)
var<storage, read> atoms: array<vec3<f32>>;

@group(0)
@binding(1)
var<storage, read> residues: array<Residue>;

@group(0)
@binding(2)
var<storage, read> settings: Settings;

@group(0)
@binding(3)
var<storage, read_write> output: array<f32>;

var<workgroup> workgroup_cv: array<f32, 16>;

@compute
@workgroup_size(32)
fn main(
    @builtin(global_invocation_id) global_id: vec3<u32>,
    @builtin(local_invocation_id) local_id: vec3<u32>,
    @builtin(workgroup_id) workgroup_id: vec3<u32>,
) {
    let cutoff = settings.cutoffs[global_id.y];
    let residue = residues[settings.residue_start + workgroup_id.x];

    if (local_id.x < residue.atom_count) {
        let current_atom = atoms[residue.atom_offset + local_id.x];

        var count = 0u;
        var vec_sum = vec3<f32>();

        for (var atom_index = 0u; atom_index < settings.atom_count; atom_index += 1u) {
            let other_atom = atoms[atom_index];
            let diff = other_atom - current_atom;
            let dist = length(diff);

            if ((dist < cutoff) && (dist > 0.0)) {
                count += 1u;
                vec_sum += diff / dist;
            }
        }

        if (count > 0u) {
            workgroup_cv[local_id.x] = 1.0 - length(vec_sum) / f32(count);
        } else {
            workgroup_cv[local_id.x] = 1.0;
        }
    }

    workgroupBarrier();

    if (local_id.x < 1u) {
        var sum = 0.0;

        for (var atom_index = 0u; atom_index < residue.atom_count; atom_index += 1u) {
            sum += workgroup_cv[atom_index];
        }

        output[settings.residue_count * global_id.y + workgroup_id.x] = sum / f32(residue.atom_count);
    }
}
