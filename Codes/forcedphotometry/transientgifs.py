#!/usr/bin/env python
"""forcedphotimages.py: Reads in forced photometry results, crossmatches it with existing list, extracts the relevant image(s), makes a light curve(s), plots the 4-panel ML images, animates the images and saves them as a gif;"""
__author__      = "Sumedha Biswas"

import sys
from astropy.io import ascii
import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.style.use('dark_background')
from astropy.table import QTable, Table, Column
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy import time, coordinates as coord, units as u
from astropy.table import Table, hstack, unique
from astropy.io import fits
from astropy.visualization import ZScaleInterval as zscale
from PIL import Image


# Reading in the forced photometry results, using memory mapping
table = Table.read('/idia/projects/meerlicht/PaulV/4Sumedha/all_gaia_data_mjd_MLforcephot_notimelimit_3sigma_20230525.fits', memmap=True)
print(len(table))

# Find matching rows, based on the GaiaX ID name
def find_rows_by_value(table, column_name, value):
    matching_rows = []
    for row in table:
        if row[column_name] == value:
            matching_rows.append(row)
    return matching_rows

column_name = 'Name'
value = 'GaiaX21-115853'

matching_rows = find_rows_by_value(table, column_name, value)

# Extract the rows that we want, and plot the lightcurve with errorbars
if len(matching_rows) > 0:
    print("Matching rows:")
    
    # Store the data for each row in separate lists
    mjd_obs_list = []
    mag_zogy_list = []
    magerr_zogy_list = []
    
    for row in matching_rows:
        # Access the columns of each row
        mjd_obs = row['MJD-OBS']
        mag_zogy = row['MAG_OPT']
        magerr_zogy = row['MAGERR_OPT']
        
        # Append the values to the respective lists
        mjd_obs_list.append(mjd_obs)
        mag_zogy_list.append(mag_zogy)
        magerr_zogy_list.append(magerr_zogy)

    # Plot all the data at once
    plt.figure()
    plt.errorbar(mjd_obs_list, mag_zogy_list, xerr=magerr_zogy_list, markersize=3, color='#2c7bb6', capsize=5, fmt='o', ecolor='lightgray')
    plt.grid()
    plt.gca().invert_yaxis()
    plt.xlim(58730, 58860)
    plt.ylim(18, 21.5)
    plt.xlabel(r'MJD')
    plt.ylabel(r'Magnitude')
    plt.legend()

    plt.savefig('GaiaX-115853-zoom.png')
    plt.close()
    
else:
    print("No matching rows found.")



# Plot light curves for all GaiaX candidate transients with matching ML detections
def test_plot(matching_rows):
    
    # Iterate over each set of duplicated rows
    for duplicate_value in duplicates:
        # Get the duplicate rows based on the specified column and value
        duplicate_rows = [row for row in table if row[column] == duplicate_value]
        
        # Extract values from duplicated rows
        name = [row['Name'] for row in duplicate_rows]
        g_values = [row['GMag'] for row in duplicate_rows]
        gerr_values = [row['GMagErr'] for row in duplicate_rows]
        m_values = [row['MAG_OPT'] for row in duplicate_rows]
        merr_values = [row['MAGERR_OPT'] for row in duplicate_rows]
        gt_values = [row['ObsTime'] for row in duplicate_rows]
        mt_values = [row['MJD-OBS'] for row in duplicate_rows]
        
        plt.errorbar(gt_values, g_values,
                     xerr = gerr_values,
                     markersize = 5, 
                     color='#d7191c',
                     capsize=5,
                     fmt='o',
                     ecolor='lightgray',
                     label = 'GaiaX')
        plt.errorbar(mt_values, m_values,
                     xerr = merr_values,
                     markersize = 3, 
                     color='#2c7bb6',
                     capsize=5,
                     fmt='o',
                     ecolor='lightgray', 
                     label = 'ML')
        plt.grid()
        plt.gca().invert_yaxis()
        plt.xlabel(r'MJD')
        plt.ylabel(r'Magnitude')
        plt.legend()
        
        filename = f'{duplicate_value}_lightcurve.pdf'
        plt.savefig('/users/sumedhabiswas/projects/gaiax_ml_stuff/lightcurves/' + filename)
        plt.close()
        
    return 

# Plot the 4-panel ML image    
def show_thumbs(data_row):
    for i, row in enumerate(matching_rows):
        fig = plt.figure(figsize=(16, 7))
        ncols, nrows = 4, 1
        fields2show = ['THUMBNAIL_RED', 'THUMBNAIL_REF', 'THUMBNAIL_D', 'THUMBNAIL_SCORR']
        
        for ip in range(len(fields2show)):
            data = row[fields2show[ip]]
            vmin, vmax = zscale().get_limits(data)
            fig.add_subplot(nrows, ncols, ip+1)
            plt.imshow(data, vmin=vmin, vmax=vmax, 
                       interpolation='none', cmap='gist_heat', origin='lower')

        directory_path = '/users/sumedhabiswas/projects/gaiax_ml_stuff/lightcurves'
        directory_name = f"{row['Name']}"
        directory_path_with_name = os.path.join(directory_path, directory_name)

        if not os.path.exists(directory_path_with_name):
            os.makedirs(directory_path_with_name)

        filename = f"{row['MJD-OBS']}_MLimage.png"
        file_path = os.path.join(directory_path_with_name, filename)
        plt.savefig(file_path)
        plt.close(fig)

# for row in matching_rows:
#     show_thumbs(row)
    
# Create gifs, using all the matching values
def animate_images(directory_path):
    image_files = glob.glob(os.path.join(directory_path, '*.png'))

    # Open all the images and create an image sequence
    image_sequence = []
    for file_path in image_files:
        image = Image.open(file_path)
        image_sequence.append(image)

    # Create a new image with the same size as the first image
    stacked_image = Image.new('RGB', image_sequence[0].size)

    # Iterate over each pixel position
    for x in range(stacked_image.width):
        for y in range(stacked_image.height):
            # Create a list of pixel values from the corresponding positions in each image
            pixel_values = [np.array(image.getpixel((x, y)), dtype=np.uint8) for image in image_sequence]

            # Calculate the average pixel value for each channel
            average_pixel = tuple(np.round(np.mean(pixel_values, axis=0)).astype(np.uint8))

            # Assign the average pixel value to the stacked image
            stacked_image.putpixel((x, y), average_pixel)

    # Save the stacked image as a single image file
    stacked_image.save(os.path.join(directory_path, 'stacked_image.png'))

    # Create an animated flipbook
    # flipbook_name = f"{row['MJD-OBS']}_flipbook.gif"
    stacked_image.save(os.path.join(directory_path, 'flipbook.gif'), 
                       save_all=True, append_images=image_sequence[1:], loop=0, duration=200)

    return 

# directory._path = '/users/sumedhabiswas/projects/gaiax_ml_stuff/lightcurves/GaiaX21-115863'
# animate_images(directory_path)
    
# for i in range(0, len(matching_rows)):
#     show_thumbs(matching_rows[i])
# # for i in range(0, len(table)):
    # show_thumbs(table[i])
