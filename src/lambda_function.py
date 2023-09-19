import boto3
import rasterio
from rasterio.enums import Resampling
from rasterio.features import geometry_mask
import geopandas as gpd
import numpy as np
from scipy import ndimage
import tempfile
import matplotlib.pyplot as plt
import json
import os
from matplotlib.colors import ListedColormap, BoundaryNorm


s3 = boto3.client('s3')

# Read bucket names from environment variables
geojson_bucket = os.environ.get('geojson_bucket', 'boundary-plot')  # Default to 'boundary-plot' if not set
png_bucket = os.environ.get('png_bucket', 'gis-colourized-png-data')  # Default to 'gis-colourized-png-data' if not set

def lambda_handler(event, context):
    print(event)
    bucket = event['Records'][0]['s3']['bucket']['name']
    key = event['Records'][0]['s3']['object']['key'].replace('+', ' ')
    print(key)

    first_part, second_part = key.split("/")
    parts = first_part.split("_")
    farm_id = parts[0]
    farm_name = "_".join(parts[1:])
    _, index = second_part.split("_")
    index = index[:-4]
    
    key_file = f"{farm_id}_{farm_name}.geojson"
    print(farm_id, farm_name, index, key_file)

    with tempfile.NamedTemporaryFile() as tmp_file, \
         tempfile.NamedTemporaryFile(suffix='.tif') as tmp_upsampled_file, \
         tempfile.NamedTemporaryFile(suffix='.geojson') as tmp_geojson_file, \
         tempfile.NamedTemporaryFile(suffix='.png') as tmp_png_file:

        s3.download_file(bucket, key, tmp_file.name)

        with rasterio.open(tmp_file.name) as src:
            upscale_factor = 10
            new_width = int(src.width * upscale_factor)
            new_height = int(src.height * upscale_factor)
            transform = src.transform * src.transform.scale(
                (src.width / new_width),
                (src.height / new_height)
            )
            upsampled_ndvi = src.read(
                out_shape=(src.count, new_height, new_width),
                resampling=Resampling.bilinear
            )
            upsampled_meta = src.meta.copy()
            upsampled_meta.update({
                "height": new_height,
                "width": new_width,
                "transform": transform
            })

            ndvi_data = upsampled_ndvi[0, :, :]
            mask = np.isnan(ndvi_data)
            ndvi_data[mask] = ndimage.generic_filter(ndvi_data, np.nanmedian, size=19)[mask]

        with rasterio.open(tmp_upsampled_file.name, 'w', **upsampled_meta) as dest:
            dest.write(upsampled_ndvi)

        s3.download_file(geojson_bucket, key_file, tmp_geojson_file.name)
        gdf = gpd.read_file(tmp_geojson_file.name)

        with rasterio.open(tmp_upsampled_file.name) as src:
            raster_crs = src.crs
            raster_transform = src.transform
            raster_shape = (src.height, src.width)
            raster_data = src.read(1)
            raster_meta = src.meta

        gdf = gdf.to_crs(raster_crs)
        mask = geometry_mask(gdf.geometry, transform=raster_transform, invert=False, out_shape=raster_shape)
        raster_data[mask] = np.nan

        rows, cols = np.where(~np.isnan(raster_data))
        min_row, max_row = np.min(rows), np.max(rows)
        min_col, max_col = np.min(cols), np.max(cols)
        cropped_data = raster_data[min_row:max_row+1, min_col:max_col+1]

        # Custom color map and bounds based on index
        if index == "NDMI":
            colors_list = ['#bbd2f0', '#79aaf8', '#4086e3', '#1e60b1', '#0c468f', '#06408c']
            bounds = [-1, -0.2, 0, 0.2, 0.4, 0.6, 1]
        elif index == "NDVI":
            colors_list = ['#808080', '#94f08d', '#4df267', '#108c07', '#0c6d05', '#074003']
            bounds = [-1, 0.009, 0.1, 0.25, 0.4, 0.6, 1]

        cmap = ListedColormap(colors_list)
        norm = BoundaryNorm(bounds, cmap.N)

        fig, ax = plt.subplots()
        plt.imshow(cropped_data, cmap=cmap, norm=norm)
        ax.set_axis_off()
        fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=10)
        plt.savefig(tmp_png_file.name, dpi=200, bbox_inches='tight', pad_inches=0, transparent=True)

        png_key = key.replace('.tif', '.png')  # Replace the original TIFF file extension with .png
        s3.upload_file(tmp_png_file.name, png_bucket, png_key)

    return {
        "statusCode": 200,
        "body": json.dumps("Success")
    }
