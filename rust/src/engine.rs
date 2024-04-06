use std::{error::Error, ops::{Bound, RangeBounds}};


static SHADER_CODE: &'static str = include_str!("./shader.wgsl");

enum BufferInfo<'a> {
    Data(&'a [u8]),
    Size(usize),
}

fn write_buffer(device: &wgpu::Device, queue: &wgpu::Queue, buffer_source: &mut Option<wgpu::Buffer>, info: BufferInfo, usage: wgpu::BufferUsages) {
    // *buffer_source = None;

    let size = match info {
        BufferInfo::Data(data) => data.len(),
        BufferInfo::Size(size) => size,
    };

    let buffer = match buffer_source {
        Some(ref buffer) if buffer.size() >= size as u64 => buffer,
        _ => {
            // eprintln!("ALLOC {}", size);

            let atoms_buffer = device.create_buffer(&wgpu::BufferDescriptor {
                label: None,
                mapped_at_creation: false /* match info {
                    BufferInfo::Data(_) => true,
                    BufferInfo::Size(_) => false,
                } */,
                size: (((size as f32) * 1.2 / 4.0).ceil() as u64) * 4,
                usage,
            });

            *buffer_source = Some(atoms_buffer);
            buffer_source.as_ref().unwrap()
        }
    };

    if let BufferInfo::Data(data) = info {
        queue.write_buffer(&buffer, 0, data);

        // {
        //     let mut buffer_mut = buffer.slice(..).get_mapped_range_mut();
        //     buffer_mut[..data.len()].copy_from_slice(data);
        // }

        // buffer.unmap();
    }

}


/// An object used to calculate circular variance in a molecule.
#[derive(Debug)]
pub struct Engine {
    bind_group_layout: wgpu::BindGroupLayout,
    device: wgpu::Device,
    pipeline: wgpu::ComputePipeline,
    queue: wgpu::Queue,

    atom_count: Option<usize>,
    residue_count: Option<usize>,
    target_residue_range: Option<(usize, usize)>,

    atoms_buffer: Option<wgpu::Buffer>,
    output_buffer: Option<wgpu::Buffer>,
    read_buffer: Option<wgpu::Buffer>,
    residues_buffer: Option<wgpu::Buffer>,
    settings_buffer: wgpu::Buffer,
}

impl Engine {
    /// Create a new instance of the engine.
    pub async fn new() -> Result<Self, Box<dyn Error>> {
        // Request instance and adapter

        let instance = wgpu::Instance::default();

        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions::default())
            .await
            .ok_or("Failed to request adapter")?;

        let (device, queue) = adapter.request_device(
            &wgpu::DeviceDescriptor {
                label: None,
                required_features: wgpu::Features::empty(),
                required_limits: wgpu::Limits::downlevel_defaults(),
            },
            None,
        )
        .await?;

        let _info = adapter.get_info();


        // Create pipeline

        let shader_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: None,
            source: wgpu::ShaderSource::Wgsl(SHADER_CODE.into()),
        });

        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: None,
            layout: None,
            module: &shader_module,
            entry_point: "main",
        });

        let bind_group_layout = pipeline.get_bind_group_layout(0);


        // Create buffers

        let settings_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Settings buffer"),
            mapped_at_creation: false,
            size: 12,
            usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        });

        Ok(Self {
            bind_group_layout,
            device,
            pipeline,
            queue,

            atom_count: None,
            residue_count: None,
            target_residue_range: None,

            read_buffer: None,
            output_buffer: None,
            settings_buffer,
            atoms_buffer: None,
            residues_buffer: None,
        })
    }

    pub fn set_residues<R: RangeBounds<usize>>(&mut self, residue_atom_counts: &[u32], atoms_data: &[f32], target_residue_range: &R) {
        let atom_count = atoms_data.len() / 4;
        let residue_count = residue_atom_counts.len();

        let target_residue_start = match target_residue_range.start_bound() {
            Bound::Excluded(&x) => x + 1,
            Bound::Included(&x) => x,
            Bound::Unbounded => 0,
        };

        let target_residue_end = match target_residue_range.end_bound() {
            Bound::Excluded(&x) => x,
            Bound::Included(&x) => x + 1,
            Bound::Unbounded => residue_count,
        };

        let target_residue_count = target_residue_end - target_residue_start;

        self.atom_count = Some(atom_count);
        self.residue_count = Some(residue_count);
        self.target_residue_range = Some((target_residue_start, target_residue_end));


        // let mut atoms_data = vec![0f32; atom_count * 4];
        let mut residues_data = vec![0u32; residue_count * 2];

        let mut current_atom_offset = 0;

        for (residue_index, &residue_atom_count) in residue_atom_counts.iter().enumerate() {
            residues_data[residue_index * 2 + 0] = residue_atom_count as u32;
            residues_data[residue_index * 2 + 1] = current_atom_offset as u32;

            current_atom_offset += residue_atom_count;

            // for atom in 0..residue_atom_count {
            //     atoms_data[current_atom_offset * 4 + 0] = atom.x() as f32;
            //     atoms_data[current_atom_offset * 4 + 1] = atom.y() as f32;
            //     atoms_data[current_atom_offset * 4 + 2] = atom.z() as f32;
            //     current_atom_offset += 1;
            // }
        }

        write_buffer(&self.device, &self.queue, &mut self.atoms_buffer, BufferInfo::Data(bytemuck::cast_slice(&atoms_data)), wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST);
        write_buffer(&self.device, &self.queue, &mut self.residues_buffer, BufferInfo::Data(bytemuck::cast_slice(&residues_data)), wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST);
        write_buffer(&self.device, &self.queue, &mut self.output_buffer, BufferInfo::Size(target_residue_count * 4), wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC);
        write_buffer(&self.device, &self.queue, &mut self.read_buffer, BufferInfo::Size(target_residue_count * 4), wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST);
    }

    pub async fn run(&mut self, cutoff_distance: f32, output: &mut [f32]) -> Result<(), Box<dyn Error>> {
        // Write settings

        let (target_residue_start, target_residue_end) = self.target_residue_range.unwrap();
        let target_residue_count = target_residue_end - target_residue_start;

        let settings_data = [
            &(self.atom_count.unwrap() as u32).to_le_bytes()[..],
            &cutoff_distance.to_le_bytes(),
            &(target_residue_start as u32).to_le_bytes()[..]
        ].concat();

        self.queue.write_buffer(&self.settings_buffer, 0, &settings_data);


        // Create bind group

        let bind_group = self.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &self.bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: self.atoms_buffer.as_ref().unwrap().as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: self.residues_buffer.as_ref().unwrap().as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: self.settings_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: self.output_buffer.as_ref().unwrap().as_entire_binding(),
                },
            ],
        });


        // Encode commands

        let mut encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });

        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes: None,
            });

            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(
                (target_residue_end - target_residue_start) as u32,
                1,
                1
            );
        }

        let read_buffer = self.read_buffer.as_ref().unwrap();

        encoder.copy_buffer_to_buffer(&self.output_buffer.as_ref().unwrap(), 0, read_buffer, 0, read_buffer.size());

        // queue.write_buffer(&settings_buffer, 4, bytemuck::cast_slice(&[current_y as u32]));
        self.queue.submit(Some(encoder.finish()));

        let buffer_slice = read_buffer.slice(..((target_residue_count * 4) as u64));
        let (sender, receiver) = flume::bounded(1);

        buffer_slice.map_async(wgpu::MapMode::Read, move |v| sender.send(v).unwrap());
        self.device.poll(wgpu::Maintain::wait()).panic_on_timeout();

        if let Ok(Ok(())) = receiver.recv_async().await {
            let data = buffer_slice.get_mapped_range();
            let cast_data = bytemuck::cast_slice::<u8, f32>(&data);

            output.copy_from_slice(cast_data);

            drop(data);

            read_buffer.unmap();

            Ok(())
        } else {
            Err("Failed to run compute on gpu!")?
        }
    }

    pub async fn run_return(&mut self, cutoff_distance: f32) -> Result<Vec<f32>, Box<dyn Error>> {
        let (target_residue_start, target_residue_end) = self.target_residue_range.unwrap();
        let target_residue_count = target_residue_end - target_residue_start;

        let mut output = vec![0f32; target_residue_count];
        self.run(cutoff_distance, &mut output).await?;

        Ok(output)
    }
}
