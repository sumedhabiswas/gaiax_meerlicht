#!/usr/bin/env python
"""forcedphotplots.py: G_vs_q plot, light curves, histogram of the # of ML detections per GaiaX source"""
__author__  = "Sumedha Biswas"

import seaborn as sns
import numpy as np  
import matplotlib.pyplot as plt 
from astropy.table import Table, unique, QTable
from astropy.io import ascii
import warnings
from collections import Counter, defaultdict
import pandas as pd
#warnings.simplefilter(action="ignore", category=FutureWarning)
#plt.rcParams['text.usetex'] = True


#----------------------------------------------------------------
# Reading in the table, printing columns, applying certain conditions
#----------------------------------------------------------------
table = Table.read('/idia/projects/meerlicht/PaulV/4Sumedha/all_gaia_data_mjd_MLforcephot_notimelimit_3sigma_20230525_noThumbnails.fits', memmap=True)
print("length of forced photomety results", len(table))

print(table.columns)

table = table[table['Dec'] < 30]
# table.sort('MJD-OBS')
# table.reverse() # highest value 60088.95355333019
# print(table['MJD-OBS'])
table.sort('Name')
table.sort('MAG_OPT') # ascending order. to plot the brightest detections
table.reverse() # descending order, to plot the faintest detections
table = table[table['FLAGS_MASK'] == 0]
table = table[table['TQC-FLAG'] != 'red']
table = table[table['QC-FLAG'] != 'red']

user = input("Do you want to consider the SNR_OPT or SNR_MAG values?")
user2 = input("Do you want to consider significant (1) or insignificant (0) detections?")
user3 = input("What is the SNR cut-off you want to set i.e. SNR >= x?")

if user == 'SNR_OPT':
    table = table[table['MAG_OPT'] < 99 ]
    table = table[table['SNR_OPT'] >= user3]
    print("Length of table after MAG_OPT < 99 and SNR_OPT >= {user3}", len(table))

    elif user == 'SNR_ZOGY':
        table = table[table['MAG_ZOGY'] < 99 ]
        table = table[table['SNR_ZOGY'] >= user3]
        print("Length of table after MAG_ZOGY < 99 and SNR_ZOGY >= {user3}", len(table))
    

print("length of table", len(table))

#----------------------------------------------------------------
# Creating a dictionary to map each unique GaiaX detection with matching ML images 

# 3181 GaiaX detections have corresponding ML images
#----------------------------------------------------------------

def create_mapping_dict(data_table, col1_name, col2_name, col3_name, 
                        col4_name, col5_name, col6_name, col7_name, col8_name):
    # Create a defaultdict to store the mapping
    mapping_dict = defaultdict(list)

    # Iterate over each row in the table
    for row in data_table:
        col1_value = row[col1_name]
        col2_value = row[col2_name]
        col3_value = row[col3_name]
        col4_value = row[col4_name]
        col5_value = row[col5_name]
        col6_value = row[col6_name]
        col7_value = row[col7_name]
        col8_value = row[col8_name]
       

        # Append the values to the dictionary
        mapping_dict[col1_value].append(col2_value)
        mapping_dict[col1_value].append(col3_value)
        mapping_dict[col1_value].append(col4_value)
        mapping_dict[col1_value].append(col5_value)
        mapping_dict[col1_value].append(col6_value)
        mapping_dict[col1_value].append(col7_value)
        mapping_dict[col1_value].append(col8_value)


    # Convert the defaultdict to a regular dictionary
    mapping_dict = dict(mapping_dict)

    # Sort the dictionary based on the ascending order of col3 values for each unique col1 value
    # mapping_dict = dict(sorted(mapping_dict.items(), key=lambda x: x[1][1]))
    
    # Count the number of unique values in col1
    unique_col1_values = len(mapping_dict)

    return mapping_dict, unique_col1_values

data_table = QTable(table)

col1_name = 'Name'
col2_name = 'FILENAME'
col3_name = 'GMag'
col4_name = 'MAG_ZOGY'
col5_name = 'GMagErr'
col6_name = 'MAGERR_ZOGY'
col7_name = 'ObsTime'
col8_name = 'MJD-OBS'

mapping_dict, unique_values_count = create_mapping_dict(data_table, 
                                                        col1_name, col2_name, col3_name, col4_name, 
                                                        col5_name, col6_name, col7_name, col8_name)

# Printing the dictionary 
for key, values in mapping_dict.items():
        print(f'{key} --> {", ".join(map(str, values))}')

        print(values)

print(f"Number of unique values in col1: {unique_values_count}")

def save_dict_to_txt(dictionary, filename):
    with open(filename, 'w') as file:
        for key, value in dictionary.items():
            file.write(f'{key}: {value}\n')
    return

# save_dict_to_txt(mapping_dict, 'G_vs_q_dict_snrgte3.txt')

#----------------------------------------------------------------
# Making a histogram of the # of ML detections per GaiaX source
#----------------------------------------------------------------
def count_values(dictionary):
    value_counts = []  # Array to store the count of values for each key

    for key, values in dictionary.items():
        count = len(values)/7  # Count the number of values for the current key
        value_counts.append(count)  # Append the count to the array

    return value_counts

value_counts = count_values(mapping_dict)
print(value_counts)
key_with_most_values = max(mapping_dict, key=lambda k: len(mapping_dict[k]))

print(key_with_most_values)  # Output: 'key2'

plt.figure()
plt.hist(value_counts, bins=1000, histtype='step')  
plt.xscale('log')
plt.yscale('log')
plt.xticks(fontsize=13)
plt.yticks(fontsize=11)
plt.xlabel('# of ML Detections', fontsize=15)
plt.ylabel('# of GaiaX Candidate Transients', fontsize=15)
# plt.xlim(0, 100)
# plt.grid()
plt.savefig('numofMLdetpergaiax.pdf')
plt.show()
plt.close()


#----------------------------------------------------------------
# Plotting G vs q
#----------------------------------------------------------------
def create_errorbar_plot(data_dict):
    fig, ax = plt.subplots()

    toofaint = 0
    brightened = 0
    dimmed = 0
    
    for key, values in data_dict.items():
        col3_values = values[1]  # GMag
        col4_values = values[2]  # MAG_OPT
        col5_values = values[3]  # GMagErr
        col6_values = values[4]  # MAGERR_OPT
        col7_values = values[5]  # Gaia MJD
        col8_values = values[6]  # ML MJD

        ax.errorbar(col3_values, col4_values,
                    xerr=col5_values, yerr=col6_values,
                    markersize=1,
                    color='#d7191c',
                    capsize=5,
                    fmt='o',
                    ecolor='lightgray', label='SNR > 3')
        
    # # print("Number of sources", len(key))
        if col4_values >= 20.68:
            toofaint += 1
        if col4_values > col3_values: # MAG_OPT > GMag 
            brightened += 1
        if col4_values < col3_values: # MAG_OPT < GMag 
            dimmed += 1
            
    print(f"Number of sources above q [MAG] > 20.68: {toofaint}")
    print(f"Number of sources that brightened: {brightened}")
    print(f"Number of sources that dimmed: {dimmed}")


    ax.set_xlim(10, 21.3)
    ax.set_ylim(10, 22.5)
    ax.axhline(y=20.68, color='#feb24c', linestyle='--')
    ax.axvline(x=20.68, color='#feb24c', linestyle='--')
    x = np.arange(10, 24, 1)
    y = np.arange(10, 24, 1)
    ax.plot(x,y, '--', color='cyan', zorder=8, alpha=0.6)
    ax.axhspan(20.68, 22, alpha=0.5, color='#9ecae1')
    ax.set_xlabel(r'$G$ [Mag]', fontsize=15)
    ax.set_ylabel(r'$q$ [Mag]', fontsize=15)
    plt.text(20.7, 10.3, s=r'GaiaX', rotation=-90, fontsize=13)
    # ax.legend()
    ax.grid(True)

    plt.savefig('/users/sumedhabiswas/projects/gaiax_ml_stuff/magplot_g_ml.pdf', bbox_inches='tight')

    plt.show()

    return 

create_errorbar_plot(mapping_dict)

#----------------------------------------------------------------
# Considering multiple detections, and looking at their light curves
#----------------------------------------------------------------
def find_duplicate_rows(table, column):
    # Get the values from the specified column
    values = table[column]
    
    # Initialize a dictionary to store counts of each value
    counts = {}
    
    # Iterate over the values and count occurrences
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    
    # Find duplicate values
    duplicates = [value for value, count in counts.items() if count > 1]
    
    # Get the duplicate rows based on the specified column
    duplicate_rows = [row for row in table if row[column] in duplicates]
    
    return duplicate_rows


def count_duplicate_rows(table, column):
    # Group the table rows based on the specified column
    grouped_table = table.group_by(column)

    # Find duplicate rows
    duplicate_rows = grouped_table.groups.size > 1

    # Get the unique values in the specified column
    unique_values = table[column].unique()

    # Count the number of duplicate rows for each value
    duplicate_row_counts = Counter(grouped_table.groups.size)

    # Create a dictionary to store the duplicate row counts
    counts = {}

    # Populate the dictionary with the duplicate row counts
    for value in unique_values:
        count = duplicate_row_counts.get(value, 0)
        counts[value] = count

    return counts

duplicate_counts = count_duplicate_rows(tab, 'Name')

for value, count in duplicate_counts.items():
    print(f"Value: {value}, Duplicate Rows: {count}")
    

def lightcurves(table, column):
    # Get the values from the specified column
    values = table[column]
    
    # Initialize a dictionary to store counts of each value
    counts = {}
    
    # Iterate over the values and count occurrences
    for value in values:
        counts[value] = counts.get(value, 0) + 1
    
    # Find duplicate values
    duplicates = [value for value, count in counts.items() if count > 1]
    
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
        plt.savefig('/users/sumedhabiswas/projects/gaiax_ml_stuff/lightcurves/ZOGY/' + filename)
        plt.close()
        
    return 

lightcurves(table, 'Name')

def lightcurves_fromdict(dictionary):
    for key, values in dictionary.items():
        # Extract the values from the dictionary
        col3_values = values[1]  # GMag
        col4_values = values[2]  # MAG_OPT
        col5_values = values[3]  # GMagErr
        col6_values = values[4]  # MAGERR_OPT
        col7_values = values[5]  # Gaia MJD
        col8_values = values[6]  # ML MJD


        for i in range(0, unique_values_count):
                
            plt.errorbar(col7_values, col3_values,
                         xerr = col5_values,
                         markersize = 2, 
                         color='#d7191c',
                         capsize=5,
                         fmt='o',
                         ecolor='lightgray',
                         label = 'GaiaX')

            plt.errorbar(col8_values, col4_values,
                         xerr = col6_values,
                         markersize = 2, 
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

            filename = f'{key}_lightcurve.pdf'
            plt.savefig('/users/sumedhabiswas/projects/gaiax_ml_stuff/lightcurves/' + filename)
            plt.close()

    return 

lightcurves_fromdict(mapping_dict)














    
