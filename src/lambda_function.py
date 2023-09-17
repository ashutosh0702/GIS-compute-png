import boto3
import rasterio
from rasterio.enums import Resampling
from rasterio.features import geometry_mask
import geopandas as gpd
import numpy as np
from scipy import ndimage
import tempfile
import matplotlib.pyplot as plt

s3 = boto3.client('s3')

def lambda_handler(event, context):

    print(event)
    # Get the bucket and object key from the event
    bucket = event['Records'][0]['s3']['bucket']['name']
    key = event['Records'][0]['s3']['object']['key']

    print(key)

    farm_id = key.split("_")[0]
    farm_name = key.split("_")[1]

    print(farm_id, farm_name)

    # Download the GeoTIFF file from S3 to a temp directory
    with tempfile.NamedTemporaryFile() as tmp_file:
        s3.download_file(bucket, key, tmp_file.name)
        
        # Step 1: Upsampling and Post-processing
        with rasterio.open(tmp_file.name) as src:
            upscale_factor = 20
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
        
            # Save the upsampled and post-processed image to another temp file
            with tempfile.NamedTemporaryFile(suffix='.tif') as tmp_upsampled_file:
                with rasterio.open(tmp_upsampled_file.name, 'w', **upsampled_meta) as dest:
                    dest.write(upsampled_ndvi)
        
                # Download the GeoJSON file from another S3 bucket
                with tempfile.NamedTemporaryFile(suffix='.geojson') as tmp_geojson_file:
                    s3.download_file('boundary-plot', f"{farm_id}_{farm_name}.geojson", tmp_geojson_file.name)
                    gdf = gpd.read_file(tmp_geojson_file.name)
                    
                    # Step 2: Masking
                    with rasterio.open(tmp_upsampled_file.name) as src:
                        raster_crs = src.crs
                        raster_transform = src.transform
                        raster_shape = (src.height, src.width)
                        raster_data = src.read(1)
                        raster_meta = src.meta
                        
                    gdf = gdf.to_crs(raster_crs)
                    mask = geometry_mask(gdf.geometry, transform=raster_transform, invert=False, out_shape=raster_shape)
                    nodata_value = np.nan
                    raster_data[mask] = nodata_value

                    # Find the extent of the valid data
                    rows, cols = np.where(~np.isnan(raster_data))
                    min_row, max_row = np.min(rows), np.max(rows)
                    min_col, max_col = np.min(cols), np.max(cols)

                    # Crop the data to the extent of the valid data
                    cropped_data = raster_data[min_row:max_row+1, min_col:max_col+1]
                    
                    # Save the masked raster to another temp file
                    with tempfile.NamedTemporaryFile(suffix='.tif') as tmp_masked_file:
                        with rasterio.open(tmp_masked_file.name, 'w', **raster_meta) as dest:
                            dest.write(cropped_data,1)
        
                # Colorize and save as PNG
                with tempfile.NamedTemporaryFile(suffix='.png') as tmp_png_file:

                    fig, ax = plt.subplots()
                    plt.imshow(cropped_data, cmap='RdYlGn')
                    ax.set_axis_off()
                    fig.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=10)
                    plt.savefig(tmp_png_file.name, dpi=200,bbox_inches='tight', pad_inches = 0 ,transparent=True)
                    
                    # Upload the colorized raster to a new S3 bucket
                    s3.upload_file(tmp_png_file.name, 'gis-colourized-png-data', f'colorized_{key}.png')
