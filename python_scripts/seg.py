# mag: 40x
# width: 65735px
# height: 57037px
# dim: 3x1x1
# pixel w/h: 0.2533 micrometer
# uncomp size: 10.5gb
# server type: openslide
# pyramid 1 4 16 32
# stain 1: hematoxylin: 0.651 0.701 0.29
# stain 2: eosin 0.216 0.801 0.558
# stain 3: residual 0.316 -0.598 0.737
# bg 255 255 255

import openslide
import numpy as np
import matplotlib.pyplot as plt
import squidpy as sq
from cellpose import models
from pathlib import Path
import cv2
from tqdm import tqdm

def process_tile(model, tile_array, min_size=30, flow_threshold=0.4):
    tile_array = cv2.normalize(tile_array, None, 0, 255, cv2.NORM_MINMAX)
    return model.eval(
        tile_array,
        channels=[2,0],
        diameter=25,
        min_size=min_size,
        flow_threshold=flow_threshold,
        cellprob_threshold=0,
        do_3D=False
    )[0]

def segment_slide(svs_path, level=1, tile_size=1024, overlap=64, output_dir="output2"):
    Path(output_dir).mkdir(exist_ok=True)
    
    # Setup
    slide = openslide.OpenSlide(svs_path)
    level_dims = slide.level_dimensions[level]
    print(f"Processing slide at level {level}: {level_dims}")
    
    # Initialize Cellpose
    model = models.Cellpose(model_type="cyto", gpu=False)
    
    # Calculate tiles
    x_tiles = int(np.ceil(level_dims[0] / (tile_size - overlap)))
    y_tiles = int(np.ceil(level_dims[1] / (tile_size - overlap)))
    
    results = []
    total_cells = 0
    
    # Process tiles
    for y in tqdm(range(y_tiles), desc="Rows"):
        row_results = []
        for x in range(x_tiles):
            # Calculate coordinates
            x_start = x * (tile_size - overlap)
            y_start = y * (tile_size - overlap)
            
            # Read tile
            tile = slide.read_region(
                (int(x_start * slide.level_downsamples[level]), 
                 int(y_start * slide.level_downsamples[level])),
                level,
                (min(tile_size, level_dims[0] - x_start),
                 min(tile_size, level_dims[1] - y_start))
            )
            
            # Process tile
            tile_array = np.array(tile)[:, :, :3]
            try:
                mask = process_tile(model, tile_array)
                cells_in_tile = len(np.unique(mask)) - 1  # Subtract background
                total_cells += cells_in_tile
                row_results.append(mask)
            except Exception as e:
                print(f"Error processing tile at ({x}, {y}): {e}")
                continue
        
        results.append(row_results)
    
    print(f"Total cells detected: {total_cells}")
    
    # Save example tile
    plt.figure(figsize=(15, 7))
    plt.subplot(121)
    plt.imshow(tile_array)
    plt.title("Sample Tile")
    plt.subplot(122)
    plt.imshow(mask, cmap='tab20b')
    plt.title(f"Segmented ({cells_in_tile} cells)")
    plt.savefig(f"{output_dir}/sample_tile.png")
    plt.close()
    
    return results

results = segment_slide("./H1546147_3_194829.svs", level=1)