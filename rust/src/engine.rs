use std::{error::Error, ops::{Bound, RangeBounds}};
use ndarray::Array2;


static SHADER_CODE: &'static str = include_str!("./shader.wgsl");

enum BufferInfo<'a> {
    Data(&'a [u8]),
    Size(usize),
}

fn write_buffer(device: &wgpu::Device, queue: &wgpu::Queue, buffer_source: &mut Option<wgpu::Buffer>, info: BufferInfo, usage: wgpu::BufferUsages) {
    let size = match info {
        BufferInfo::Data(data) => data.len(),
        BufferInfo::Size(size) => size,
    };

    let buffer = match buffer_source {
        Some(ref buffer) if buffer.size() >= size as u64 => buffer,
        _ => {
            let atoms_buffer = device.create_buffer(&wgpu::BufferDescriptor {
                label: None,
                mapped_at_creation: false,
                size: (((size as f32) * 1.2 / 4.0).ceil() as u64) * 4,
                usage,
            });

            *buffer_source = Some(atoms_buffer);
            buffer_source.as_ref().unwrap()
        }
    };

    if let BufferInfo::Data(data) = info {
        queue.write_buffer(&buffer, 0, data);
    }

}


/// An object used to calculate circular variance in a molecule.
#[derive(Debug)]
pub struct Engine {
    bind_group_layout: wgpu::BindGroupLayout,
    device: wgpu::Device,
    pipeline: wgpu::ComputePipeline,
    queue: wgpu::Queue,

    atoms_buffer: Option<wgpu::Buffer>,
    output_buffer: Option<wgpu::Buffer>,
    read_buffer: Option<wgpu::Buffer>,
    residues_buffer: Option<wgpu::Buffer>,
    settings_buffer: Option<wgpu::Buffer>,
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

        Ok(Self {
            bind_group_layout,
            device,
            pipeline,
            queue,

            atoms_buffer: None,
            output_buffer: None,
            read_buffer: None,
            residues_buffer: None,
            settings_buffer: None,
        })
    }

    /// Calculate the circular variance of protein residues.
    pub async fn run<R: RangeBounds<usize>>(
        &mut self,
        residue_atom_counts: &[u32],
        atoms_data: &[f32],
        target_residue_range: R,
        cutoff_distances: &[f32],
    ) -> Result<Array2<f32>, Box<dyn Error>> {
        let atom_count = atoms_data.len() / 4;
        let residue_count = residue_atom_counts.len();


        // Calculate target range

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

        if cutoff_distances.is_empty() {
            return Ok(Array2::zeros((0, target_residue_count)));
        }


        // Prepare residues data

        let mut current_atom_offset = 0;
        let mut residues_data = vec![0u32; residue_count * 2];

        for (residue_index, &residue_atom_count) in residue_atom_counts.iter().enumerate() {
            residues_data[residue_index * 2 + 0] = residue_atom_count as u32;
            residues_data[residue_index * 2 + 1] = current_atom_offset as u32;

            current_atom_offset += residue_atom_count;
        }


        // Prepare settings data

        let settings_data = [
            &(atom_count as u32).to_le_bytes()[..],
            &(target_residue_count as u32).to_le_bytes()[..],
            &(target_residue_start as u32).to_le_bytes()[..],
            bytemuck::cast_slice(&cutoff_distances),
        ].concat();


        // Write buffers

        let output_buffer_size = target_residue_count * cutoff_distances.len() * 4;

        write_buffer(&self.device, &self.queue, &mut self.atoms_buffer, BufferInfo::Data(bytemuck::cast_slice(&atoms_data)), wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST);
        write_buffer(&self.device, &self.queue, &mut self.residues_buffer, BufferInfo::Data(bytemuck::cast_slice(&residues_data)), wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST);
        write_buffer(&self.device, &self.queue, &mut self.output_buffer, BufferInfo::Size(output_buffer_size), wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC);
        write_buffer(&self.device, &self.queue, &mut self.read_buffer, BufferInfo::Size(output_buffer_size), wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST);
        write_buffer(&self.device, &self.queue, &mut self.settings_buffer, BufferInfo::Data(&settings_data), wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST);


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
                    resource: self.settings_buffer.as_ref().unwrap().as_entire_binding(),
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
                cutoff_distances.len() as u32,
                1,
            );
        }

        let read_buffer = self.read_buffer.as_ref().unwrap();

        encoder.copy_buffer_to_buffer(&self.output_buffer.as_ref().unwrap(), 0, read_buffer, 0, read_buffer.size());

        // queue.write_buffer(&settings_buffer, 4, bytemuck::cast_slice(&[current_y as u32]));
        self.queue.submit(Some(encoder.finish()));

        let buffer_slice = read_buffer.slice(..(output_buffer_size as u64));
        let (sender, receiver) = flume::bounded(1);

        buffer_slice.map_async(wgpu::MapMode::Read, move |v| sender.send(v).unwrap());
        self.device.poll(wgpu::Maintain::wait()).panic_on_timeout();

        if let Ok(Ok(())) = receiver.recv_async().await {
            let data = buffer_slice.get_mapped_range();
            let cast_data = bytemuck::cast_slice::<u8, f32>(&data);

            let output = ndarray::Array::from_shape_vec(
                (
                    cutoff_distances.len(),
                    target_residue_count,
                ),
                cast_data.to_vec(),
            )?;

            drop(data);
            read_buffer.unmap();

            Ok(output)
        } else {
            Err("Failed to run compute on gpu!")?
        }
    }
}
