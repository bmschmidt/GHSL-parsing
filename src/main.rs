use gdal::spatial_ref::{CoordTransform, SpatialRef};
use gdal::Dataset;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64;
use rayon::prelude::*;
use std::error::Error;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{ArrayRef, BooleanArray, Float64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

/// Process multiple GeoTIFF population density files to generate points with year-based attributes.
fn process_tiffs_to_points(
    datasets: &[(usize, Dataset)],
    years: &[usize],
    output_path: &str,
    min_population: f64,
    chunk_size: usize,
    scaling_factor: f64,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let dataset = &datasets[0].1;
    let raster_size = dataset.raster_size();
    let width = raster_size.0 as usize;
    let height = raster_size.1 as usize;

    // Geotransform and coordinate transformation
    let geo_transform = dataset.geo_transform()?;
    let source_srs: SpatialRef = dataset.spatial_ref().unwrap(); // Mollweide
    let target_srs: SpatialRef = SpatialRef::from_epsg(3857)?; // Web Mercator
    let coord_transform: CoordTransform = CoordTransform::new(&source_srs, &target_srs)?;

    let num_years = datasets.len();

    let mut y_start = 0;
    let tmp_path = format!("{output_path}.tmp", output_path = output_path);
    let mut writer: ArrowWriter<File> = make_writer(&tmp_path)?;

    while y_start < height {
        let y_end = std::cmp::min(y_start + chunk_size, height);
        let y_size = y_end - y_start;

        // Read data for all years
        let mut data_per_year: Vec<Vec<f64>> = Vec::with_capacity(num_years);
        let mut x_coords: Vec<f64> = Vec::new();
        let mut y_coords = Vec::new();
        let mut per_year_bools: Vec<Vec<bool>> = vec![Vec::new(); num_years];

        for (_, dataset) in datasets.iter() {
            let band = dataset.rasterband(1)?;
            let buffer_size = width * y_size;
            let mut data = vec![0.0f64; buffer_size];

            band.read_into_slice::<f64>(
                (0, y_start as isize),
                (width, y_size),
                (width, y_size),
                &mut data,
                None,
            )?;

            data_per_year.push(data);
        }

        // Process each pixel
        for y_offset in 0..y_size {
            for x in 0..width {
                let idx = y_offset * width + x;

                // Get population values for all years
                let mut num_points_per_year = Vec::with_capacity(num_years);
                let mut max_num_points = 0;

                for data in &data_per_year {
                    let value = data[idx];
                    let num_points = if value <= min_population {
                        0
                    } else {
                        let exact_points = value * scaling_factor;
                        exact_points.floor() as usize
                    };
                    num_points_per_year.push(num_points);
                    if num_points > max_num_points {
                        max_num_points = num_points;
                    }
                }

                if max_num_points == 0 {
                    continue;
                }

                // Generate max_num_points points within the pixel using a deterministic RNG
                let mut rng =
                    Pcg64::seed_from_u64(((x + y_start + y_offset) as u64) << 32 | x as u64);

                // Convert pixel to geographical coordinates and transform
                let (x_geo, y_geo) =
                    pixel_to_geo(x as f64, (y_start + y_offset) as f64, &geo_transform);
                let (x_next_geo, y_next_geo) = pixel_to_geo(
                    (x + 1) as f64,
                    (y_start + y_offset + 1) as f64,
                    &geo_transform,
                );
                let coords = vec![
                    [x_geo, y_geo],
                    [x_next_geo, y_geo],
                    [x_next_geo, y_next_geo],
                    [x_geo, y_next_geo],
                    [x_geo, y_geo],
                ];
                let transformed = transform_coords(&coords, &coord_transform)?;
                // on error skip this pixel. This only happens in some weird places right along the equator and poles, right?

                // let transformed_coords = match transformed {
                //     Ok(coords) => coords,
                //     Err(_) => continue,
                // };

                let (triangles_array, share) = triangles(&transformed);

                let mut points = Vec::with_capacity(max_num_points);

                for _ in 0..max_num_points {
                    let triangle_idx = if rng.gen::<f64>() < share { 0 } else { 1 };
                    let triangle = &triangles_array[triangle_idx];

                    let r1: f64 = rng.gen();
                    let r2: f64 = rng.gen();

                    let sqrt_r1 = r1.sqrt();

                    let u = 1.0 - sqrt_r1;
                    let v = sqrt_r1 * (1.0 - r2);
                    let w = sqrt_r1 * r2;

                    let x_coord = u * triangle[0][0] + v * triangle[1][0] + w * triangle[2][0];
                    let y_coord = u * triangle[0][1] + v * triangle[1][1] + w * triangle[2][1];

                    points.push((x_coord, y_coord));
                }

                // For each point, determine existence in each year
                for (point_idx, (x_coord, y_coord)) in points.into_iter().enumerate() {
                    x_coords.push(x_coord);
                    y_coords.push(y_coord);
                    for (year_idx, &num_points) in num_points_per_year.iter().enumerate() {
                        let exists = point_idx < num_points;
                        per_year_bools[year_idx].push(exists);
                    }
                }
            }
        }

        println!(
            "{output_path} {percent:.1}% completed, {n:.1}K points written",
            output_path = output_path,
            percent = ((y_start + chunk_size) as f64 / height as f64) * 100.0,
            n = x_coords.len() as f64 / 1000.0
        );
        write_points_with_years(&x_coords, &y_coords, &per_year_bools, years, &mut writer)?;

        y_start += chunk_size;
    }

    let _ = writer.close();

    // Remove tmp prefix from filename

    std::fs::rename(&tmp_path, output_path)?;

    // Write points and year attributes to Parquet file

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
    let triangle_a = [coords[0], coords[1], coords[2]];
    let triangle_b = [coords[0], coords[2], coords[3]];

    let area_a = triangle_area(&triangle_a);
    let area_b = triangle_area(&triangle_b);
    let total_area = area_a + area_b;

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

fn schema() -> Arc<Schema> {
    // Just writing these out by hand.
    let fields = vec![
        Field::new("x", DataType::Float64, false),
        Field::new("y", DataType::Float64, false),
        Field::new("1975", DataType::Boolean, false),
        Field::new("1980", DataType::Boolean, false),
        Field::new("1985", DataType::Boolean, false),
        Field::new("1990", DataType::Boolean, false),
        Field::new("1995", DataType::Boolean, false),
        Field::new("2000", DataType::Boolean, false),
        Field::new("2005", DataType::Boolean, false),
        Field::new("2010", DataType::Boolean, false),
        Field::new("2015", DataType::Boolean, false),
        Field::new("2020", DataType::Boolean, false),
        Field::new("2025", DataType::Boolean, false),
    ];
    let schema: Arc<Schema> = Arc::new(Schema::new(fields));
    return schema;
}

fn make_writer(output_path: &str) -> Result<ArrowWriter<File>, Box<dyn Error + Send + Sync>> {
    let file = File::create(output_path)?;

    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .build();

    let writer = ArrowWriter::try_new(file, schema(), Some(props)).unwrap();

    Ok(writer)
}

/// Write points and year attributes to an Arrow IPC file.
fn write_points_with_years(
    x_coords: &[f64],
    y_coords: &[f64],
    per_year_bools: &[Vec<bool>],
    years: &[usize],
    writer: &mut ArrowWriter<File>,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    // Create Arrow arrays

    let x_array: arrow::array::PrimitiveArray<arrow::datatypes::Float64Type> =
        Float64Array::from(x_coords.to_vec());
    let y_array = Float64Array::from(y_coords.to_vec());

    let mut columns: Vec<ArrayRef> = vec![Arc::new(x_array), Arc::new(y_array)];

    // Add a BooleanArray for each year
    for (year_idx, _year) in years.iter().enumerate() {
        let bool_array = BooleanArray::from(per_year_bools[year_idx].clone());
        columns.push(Arc::new(bool_array));
    }

    // Create RecordBatch
    let batch = RecordBatch::try_new(schema(), columns)?;

    writer.write(&batch).expect("Writing batch");

    Ok(())
}

fn parse(row: usize, column: usize) -> Result<(), Box<dyn Error + Send + Sync>> {
    let output_file = format!("points/{row}-{column}.parquet", row = row, column = column);
    if Path::new(&output_file).exists() {
        return Ok(());
        //std::fs::remove_file(&output_file)?;
    }

    let years = [
        1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025,
    ];
    let scaling_factor = 1.0;

    // Load all datasets upfront
    let rasters: Result<Vec<(usize, Dataset)>, gdal::errors::GdalError> = years
        .iter()
        .enumerate()
        .map(|(idx, &year)| {
            let input_file = format!(
                "data/GHS_POP_E{year}_GLOBE_R2023A_54009_100_V1_0_R{row}_C{column}.tif",
                year = year
            );

            // Fail if the input file does not exist
            if !Path::new(&input_file).exists() {
                println!("File not found: {}", input_file);
            }

            Dataset::open(&input_file).map(|dataset| (idx, dataset))
        })
        .collect();

    let rasters = rasters?; // Propagate errors if any

    if rasters.len() != years.len() {
        eprintln!(
            "No datasets available for row: {}, column: {}. Skipping.",
            row, column
        );
        return Ok(());
    }

    process_tiffs_to_points(&rasters, &years, &output_file, 1e-5, 500, scaling_factor)?;

    Ok(())
}

fn main() -> Result<(), Box<dyn Error + Send + Sync>> {
    let mut tuples = vec![
        (1, 13),
        (1, 14),
        (1, 15),
        (1, 16),
        (1, 17),
        (1, 18),
        (1, 19),
        (1, 20),
        (1, 21),
        (1, 22),
        (1, 23),
        (1, 24),
        (1, 25),
        (2, 8),
        (2, 9),
        (2, 10),
        (2, 11),
        (2, 12),
        (2, 13),
        (2, 14),
        (2, 15),
        (2, 16),
        (2, 17),
        (2, 18),
        (2, 19),
        (2, 20),
        (2, 21),
        (2, 22),
        (2, 23),
        (2, 24),
        (2, 25),
        (2, 26),
        (2, 27),
        (2, 28),
        (2, 29),
        (2, 30),
        (3, 5),
        (3, 6),
        (3, 7),
        (3, 8),
        (3, 9),
        (3, 10),
        (3, 11),
        (3, 12),
        (3, 13),
        (3, 14),
        (3, 15),
        (3, 16),
        (3, 18),
        (3, 19),
        (3, 20),
        (3, 21),
        (3, 22),
        (3, 23),
        (3, 24),
        (3, 25),
        (3, 26),
        (3, 27),
        (3, 28),
        (3, 29),
        (3, 30),
        (3, 31),
        (3, 32),
        (4, 8),
        (4, 9),
        (4, 10),
        (4, 11),
        (4, 12),
        (4, 13),
        (4, 14),
        (4, 18),
        (4, 19),
        (4, 20),
        (4, 21),
        (4, 22),
        (4, 23),
        (4, 24),
        (4, 25),
        (4, 26),
        (4, 27),
        (4, 28),
        (4, 29),
        (4, 30),
        (4, 31),
        (5, 8),
        (5, 9),
        (5, 10),
        (5, 11),
        (5, 12),
        (5, 13),
        (5, 16),
        (5, 17),
        (5, 18),
        (5, 19),
        (5, 20),
        (5, 21),
        (5, 22),
        (5, 23),
        (5, 24),
        (5, 25),
        (5, 26),
        (5, 27),
        (5, 28),
        (5, 29),
        (5, 30),
        (5, 31),
        (6, 2),
        (6, 3),
        (6, 8),
        (6, 9),
        (6, 10),
        (6, 11),
        (6, 13),
        (6, 17),
        (6, 18),
        (6, 19),
        (6, 20),
        (6, 21),
        (6, 22),
        (6, 23),
        (6, 24),
        (6, 25),
        (6, 26),
        (6, 27),
        (6, 28),
        (6, 29),
        (6, 30),
        (6, 31),
        (6, 32),
        (7, 2),
        (7, 3),
        (7, 4),
        (7, 7),
        (7, 8),
        (7, 9),
        (7, 10),
        (7, 11),
        (7, 12),
        (7, 13),
        (7, 16),
        (7, 17),
        (7, 18),
        (7, 19),
        (7, 20),
        (7, 22),
        (7, 23),
        (7, 24),
        (7, 25),
        (7, 26),
        (7, 27),
        (7, 28),
        (7, 29),
        (7, 30),
        (7, 31),
        (7, 32),
        (7, 33),
        (7, 35),
        (8, 8),
        (8, 9),
        (8, 10),
        (8, 11),
        (8, 12),
        (8, 13),
        (8, 16),
        (8, 17),
        (8, 18),
        (8, 19),
        (8, 20),
        (8, 21),
        (8, 22),
        (8, 23),
        (8, 24),
        (8, 26),
        (8, 27),
        (8, 28),
        (8, 29),
        (8, 30),
        (8, 31),
        (8, 32),
        (8, 33),
        (8, 34),
        (8, 35),
        (8, 36),
        (9, 1),
        (9, 2),
        (9, 3),
        (9, 9),
        (9, 10),
        (9, 11),
        (9, 12),
        (9, 13),
        (9, 14),
        (9, 16),
        (9, 17),
        (9, 18),
        (9, 19),
        (9, 20),
        (9, 21),
        (9, 22),
        (9, 23),
        (9, 24),
        (9, 26),
        (9, 27),
        (9, 28),
        (9, 29),
        (9, 30),
        (9, 31),
        (9, 32),
        (9, 33),
        (9, 34),
        (9, 35),
        (9, 36),
        (10, 1),
        (10, 3),
        (10, 5),
        (10, 9),
        (10, 10),
        (10, 11),
        (10, 12),
        (10, 13),
        (10, 14),
        (10, 15),
        (10, 17),
        (10, 19),
        (10, 20),
        (10, 21),
        (10, 22),
        (10, 23),
        (10, 24),
        (10, 26),
        (10, 28),
        (10, 29),
        (10, 30),
        (10, 31),
        (10, 32),
        (10, 33),
        (10, 34),
        (10, 35),
        (10, 36),
        (11, 1),
        (11, 2),
        (11, 3),
        (11, 4),
        (11, 5),
        (11, 11),
        (11, 12),
        (11, 13),
        (11, 14),
        (11, 15),
        (11, 18),
        (11, 20),
        (11, 21),
        (11, 22),
        (11, 23),
        (11, 24),
        (11, 28),
        (11, 29),
        (11, 30),
        (11, 31),
        (11, 32),
        (11, 33),
        (11, 34),
        (11, 35),
        (11, 36),
        (12, 1),
        (12, 2),
        (12, 3),
        (12, 4),
        (12, 5),
        (12, 6),
        (12, 11),
        (12, 12),
        (12, 13),
        (12, 14),
        (12, 15),
        (12, 16),
        (12, 20),
        (12, 21),
        (12, 22),
        (12, 23),
        (12, 24),
        (12, 25),
        (12, 29),
        (12, 30),
        (12, 31),
        (12, 32),
        (12, 33),
        (12, 34),
        (12, 35),
        (12, 36),
        (13, 2),
        (13, 5),
        (13, 6),
        (13, 7),
        (13, 8),
        (13, 9),
        (13, 11),
        (13, 12),
        (13, 13),
        (13, 14),
        (13, 20),
        (13, 21),
        (13, 22),
        (13, 23),
        (13, 29),
        (13, 30),
        (13, 31),
        (13, 32),
        (13, 33),
        (13, 34),
        (14, 11),
        (14, 12),
        (14, 13),
        (14, 14),
        (14, 17),
        (14, 18),
        (14, 20),
        (14, 21),
        (14, 25),
        (14, 29),
        (14, 30),
        (14, 31),
        (14, 32),
        (14, 33),
        (14, 34),
        (15, 4),
        (15, 12),
        (15, 13),
        (15, 14),
        (15, 22),
        (15, 23),
        (15, 24),
        (15, 31),
        (15, 32),
        (15, 33),
        (16, 13),
        (16, 14),
        (16, 15),
        (16, 16),
        (16, 17),
        (16, 19),
        (16, 24),
        (16, 30),
        (16, 31),
        (17, 9),
        (17, 14),
        (17, 15),
        (17, 16),
        (17, 18),
        (17, 19),
        (17, 20),
        (17, 21),
        (17, 22),
        (17, 23),
        (17, 24),
        (17, 25),
        (17, 26),
        (17, 27),
        (17, 28),
        (18, 12),
        (18, 13),
        (18, 14),
        (18, 15),
        (18, 16),
        (18, 17),
        (18, 18),
        (18, 19),
        (18, 20),
        (18, 21),
        (18, 22),
        (18, 23),
        (18, 24),
        (18, 25),
        (18, 26),
    ];

    // Sort by distance from Tokyo or something. I do this just so I can inspect the output from the start.
    tuples.sort_by(|a, b| {
        let a_dist = (a.0 as f64 - 5.0).hypot(a.1 as f64 - 31.0);
        let b_dist = (b.0 as f64 - 5.0).hypot(b.1 as f64 - 31.0);
        a_dist.partial_cmp(&b_dist).unwrap()
    });
    tuples.par_iter().try_for_each(|(row, column)| {
        println!("Processing row: {}, column: {}", row, column);
        parse(*row, *column)
    })?;
    // let truples = [(3, 21)];
    // for (row, column) in tuples {
    //     println!("Processing row: {}, column: {}", row, column);
    //     parse(row, column)?;
    // }
    Ok(())
}
