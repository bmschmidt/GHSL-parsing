use gdal::spatial_ref::{CoordTransform, SpatialRef};
use gdal::Dataset;
use rand::Rng;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;
use std::error::Error;

use arrow::array::{ArrayRef, Float64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::ipc::writer::FileWriter;
use arrow::record_batch::RecordBatch;

/// Process a GeoTIFF population density file to generate random points within pixels.
fn process_tiff_to_points(
    output_path: &str,
    row: usize,
    column: usize,
    min_population: f64,
    chunk_size: usize,
    scaling_factor: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    // Open the dataset
    // let dataset = Dataset::open(input_path)?;
    // Get raster size

    let years = [1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025];

    let rasters: Result<Vec<(i32, Dataset)>, gdal::errors::GdalError> = years.iter().map(|&year| {
        let input_file = format!(
            "data/GHS_POP_E{year}_GLOBE_R2023A_54009_100_V1_0_R5_C12.tif",
            year = year
        );
        Dataset::open(&input_file).map(|dataset| (year, dataset))
    }).collect();
    let rasters = rasters?; // Propagate errors if any

    let dataset = rasters.iter().find(|(year, _)| *year == 2015).unwrap().1;
    let raster_size = dataset.raster_size();
    let width = raster_size.0 as usize;
    let height = raster_size.1 as usize;

    // Get geotransform
    let geo_transform = dataset.geo_transform()?;

    // Get the source spatial reference
    let source_srs = dataset.spatial_ref().unwrap();

    // Create the target spatial reference (WGS84)
    let target_srs = SpatialRef::from_epsg(3857)?;

    // Create coordinate transformation
    let coord_transform = CoordTransform::new(&source_srs, &target_srs)?;

    let band = dataset.rasterband(1)?;

    let mut x_coords = Vec::new();
    let mut y_coords = Vec::new();
    let mut total_points = 0;

    let mut y_start = 0;

    while y_start < height {
        let y_end = std::cmp::min(y_start + chunk_size, height);
        let y_size = y_end - y_start;

        // Read chunk of data
        let buffer_size = width * y_size;
        let mut data = vec![0.0f64; buffer_size];

        band.read_into_slice::<f64>(
            (0, y_start as isize),
            (width, y_size),
            (width, y_size),
            &mut data,
            None,
        )?;

        // Process each pixel in chunk
        for y_offset in 0..y_size {
            for x in 0..width {
                let idx = y_offset * width + x;
                let value = data[idx];

                if value <= min_population {
                    continue;
                }

                let y = y_start + y_offset;

                // Calculate pixel bounds in source CRS
                let (x_geo, y_geo) = pixel_to_geo(x as f64, y as f64, &geo_transform);
                let (x_next_geo, y_next_geo) =
                    pixel_to_geo((x + 1) as f64, (y + 1) as f64, &geo_transform);

                // Create rectangle coordinates
                let coords = vec![
                    [x_geo, y_geo],
                    [x_next_geo, y_geo],
                    [x_next_geo, y_next_geo],
                    [x_geo, y_next_geo],
                    [x_geo, y_geo],
                ];

                // Transform coordinates to WGS84
                let transformed_coords = transform_coords(&coords, &coord_transform)?;

                // Divide rectangle into triangles
                let (triangles_array, share) = triangles(&transformed_coords);

                // Determine the number of points to generate
                let exact_points = value * scaling_factor as f64;

                // Number of points in each triangle
                let num_points_a = ((exact_points as f64) * share) as f64;
                let num_points_b = exact_points - num_points_a;

                // Generate random points in each triangle
                let points_a = random_points(&triangles_array[0], num_points_a);
                let points_b = random_points(&triangles_array[1], num_points_b);

                // Collect points and values
                for point in points_a.into_iter().chain(points_b) {
                    x_coords.push(point[0]);
                    y_coords.push(point[1]);
                    total_points += 1;
                }
            }
        }

        println!(
            "Processed {:.1}% ({} points)",
            ((y_start + chunk_size) as f64 / height as f64) * 100.0,
            total_points
        );

        y_start += chunk_size;
    }

    // Write all points at once
    write_points(&x_coords, &y_coords, output_path)?;

    Ok(())
}

/// Converts pixel coordinates to geographical coordinates using the geotransform.
fn pixel_to_geo(x_pix: f64, y_pix: f64, transform: &[f64; 6]) -> (f64, f64) {
    let x_geo = transform[0] + x_pix * transform[1] + y_pix * transform[2];
    let y_geo = transform[3] + x_pix * transform[4] + y_pix * transform[5];
    (x_geo, y_geo)
}

/// Transforms coordinates from source CRS to target CRS.
fn transform_coords(
    coords: &[[f64; 2]],
    coord_transform: &CoordTransform,
) -> gdal::errors::Result<Vec<[f64; 2]>> {
    let mut xs: Vec<f64> = coords.iter().map(|coord| coord[0]).collect();
    let mut ys: Vec<f64> = coords.iter().map(|coord| coord[1]).collect();
    let mut zs = vec![0.0; coords.len()]; // Assuming 2D coordinates

    coord_transform.transform_coords(&mut xs, &mut ys, &mut zs)?;

    let transformed_coords = xs
        .into_iter()
        .zip(ys.into_iter())
        .map(|(x, y)| [x, y])
        .collect();

    Ok(transformed_coords)
}

/// Splits a rectangle into two triangles and calculates the share of the first triangle.
fn triangles(coords: &[[f64; 2]]) -> ([[[f64; 2]; 3]; 2], f64) {
    // Triangle A: coords[0], coords[1], coords[2]
    let triangle_a = [coords[0], coords[1], coords[2]];
    // Triangle B: coords[0], coords[2], coords[3]
    let triangle_b = [coords[0], coords[2], coords[3]];

    // Calculate areas
    let area_a = triangle_area(&triangle_a);
    let area_b = triangle_area(&triangle_b);
    let total_area = area_a + area_b;

    // Share of the area for triangle A
    let share = area_a / total_area;

    ([triangle_a, triangle_b], share)
}

/// Calculates the area of a triangle using the shoelace formula.
fn triangle_area(triangle: &[[f64; 2]; 3]) -> f64 {
    let (x0, y0) = (triangle[0][0], triangle[0][1]);
    let (x1, y1) = (triangle[1][0], triangle[1][1]);
    let (x2, y2) = (triangle[2][0], triangle[2][1]);

    ((x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1)).abs()) / 2.0
}

/// Generates random points within a triangle.
fn random_points(triangle: &[[f64; 2]; 3], num_points: f64) -> Vec<[f64; 2]> {
    let mut rng = rand::thread_rng();
    let last_seed : f64 = rng.gen();
    let mut total_points = num_points.round() as usize;
    if last_seed > num_points % 1.0 {
        total_points += 1;
    }

    let mut points = Vec::with_capacity(total_points);

    for _ in 0..total_points {
        let r1: f64 = rng.gen();
        let r2: f64 = rng.gen();

        let sqrt_r1 = r1.sqrt();

        let u = 1.0 - sqrt_r1;
        let v = sqrt_r1 * (1.0 - r2);
        let w = sqrt_r1 * r2;

        let x = u * triangle[0][0] + v * triangle[1][0] + w * triangle[2][0];
        let y = u * triangle[0][1] + v * triangle[1][1] + w * triangle[2][1];

        points.push([x, y]);
    }


    points
}

/// Write points to an Arrow IPC file.
fn write_points(
    x_coords: &[f64],
    y_coords: &[f64],
    output_path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    // Create Arrow arrays
    let x_array = Float64Array::from(x_coords.to_vec());
    let y_array = Float64Array::from(y_coords.to_vec());

    // Create fields
    let fields = vec![
        Field::new("x", DataType::Float64, false),
        Field::new("y", DataType::Float64, false),
    ];

    let schema = Arc::new(Schema::new(fields));

    // Create RecordBatch
    let columns: Vec<ArrayRef> = vec![
        Arc::new(x_array),
        Arc::new(y_array),
    ];
    let batch = RecordBatch::try_new(schema.clone(), columns)?;

    // Write to file
    let file = File::create(output_path)?;
    let mut writer = FileWriter::try_new(file, &schema)?;

    writer.write(&batch)?;
    writer.finish()?;

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let output_file = "population_points.arrow";

    if Path::new(output_file).exists() {
        std::fs::remove_file(output_file)?;
    }

    // Adjust the scaling_factor as needed
    let scaling_factor = 0.1;

    process_tiff_to_points(output_file, 5, 12, 1e-5, 100, scaling_factor)?;

    Ok(())
}
