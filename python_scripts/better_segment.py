import openslide
import numpy as np
import matplotlib.pyplot as plt
from cellpose import models
from pathlib import Path
import cv2
from tqdm import tqdm
import multiprocessing as mp
import sys
from datetime import datetime
import pickle

def process_tile(args):
    x, y, tile_array, model_type = args
    model = models.Cellpose(model_type=model_type, gpu=False)
    tile_array = cv2.normalize(tile_array, None, 0, 255, cv2.NORM_MINMAX)
    try:
        mask = model.eval(
            tile_array,
            channels=[2,0],
            diameter=25,
            min_size=30,
            flow_threshold=0.4,
            cellprob_threshold=0,
            do_3D=False
        )[0]
        return (x, y, mask, len(np.unique(mask)) - 1)
    except Exception as e:
        print(f"Error at tile ({x},{y}): {e}")
        return (x, y, None, 0)

def segment_slide_parallel(svs_path, level=1, tile_size=1024, overlap=64, n_processes=35):
    slide = openslide.OpenSlide(svs_path)
    level_dims = slide.level_dimensions[level]
    print(f"Processing slide {level_dims} using {n_processes} cores")
    
    x_tiles = int(np.ceil(level_dims[0] / (tile_size - overlap)))
    y_tiles = int(np.ceil(level_dims[1] / (tile_size - overlap)))
    
    tasks = [(x, y, np.array(slide.read_region(
        (int(x * (tile_size - overlap) * slide.level_downsamples[level]),
         int(y * (tile_size - overlap) * slide.level_downsamples[level])),
        level,
        (min(tile_size, level_dims[0] - x * (tile_size - overlap)),
         min(tile_size, level_dims[1] - y * (tile_size - overlap)))
    ))[:,:,:3], "cyto") for x in range(x_tiles) for y in range(y_tiles)]
    
    with mp.Pool(n_processes) as pool:
        results = list(tqdm(pool.imap(process_tile, tasks), total=len(tasks)))
    
    # Save results to pickle
    with open('segmentation_results.pkl', 'wb') as f:
        pickle.dump(results, f)
    
    return results

def visualize_results(results=None, output_dir=None, svs_path=None):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = output_dir or f"output_{timestamp}"
    Path(output_dir).mkdir(exist_ok=True)
    
    try:
        if results is None:
            with open('segmentation_results.pkl', 'rb') as f:
                results = pickle.load(f)
        
        # Get actual tile dimensions from first valid mask
        tile_dims = next((mask.shape for _, _, mask, _ in results if mask is not None), (1024, 1024))
        tile_height, tile_width = tile_dims
        
        # Get grid dimensions
        x_max = max(r[0] for r in results) + 1
        y_max = max(r[1] for r in results) + 1
        
        # Create composite with correct orientation
        composite = np.zeros((y_max * tile_height, x_max * tile_width), dtype=np.int32)
        
        # Place tiles with y-axis inverted
        cell_count = 0
        for x, y, mask, cells in results:
            if mask is not None:
                # Invert y coordinate for correct orientation
                y_inv = (y_max - 1) - y
                y_start = y_inv * tile_height
                x_start = x * tile_width
                
                # Handle potential size mismatches
                mask_height, mask_width = mask.shape
                composite[y_start:y_start + mask_height, 
                         x_start:x_start + mask_width] = mask
                cell_count += cells
        
        # Validate dimensions if SVS provided
        if svs_path:
            slide = openslide.OpenSlide(svs_path)
            level_dims = slide.level_dimensions[1]  # Level 1
            print(f"Original dimensions: {level_dims}")
            print(f"Composite dimensions: {composite.shape}")
        
        plt.figure(figsize=(20, 20))
        plt.imshow(composite, cmap='tab20b')
        plt.title(f"Segmentation Map - {cell_count} cells")
        plt.colorbar()
        save_path = f"{output_dir}/composite.png"
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
        
    except Exception as e:
        print(f"Error during visualization: {e}")

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1] == 'visualize':
            print("Starting visualization...")
            visualize_results()
        else:
            svs_path = sys.argv[1]
            results = segment_slide_parallel(svs_path)
            visualize_results(results, svs_path=svs_path)