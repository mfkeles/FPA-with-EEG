# EEG and Fiber Photometry Data Processing and Plotting Script

This repository contains a MATLAB script for handling and combining EEG and Fiber Photometry data. The script processes and combines the data, produces summaries, and optionally generates a variety of plots.

## Script Functionality

This script sets a folder path for EEG and Fiber Photometry data and prompts the user to decide whether or not to generate plots during execution.

An `io` class instance is created to handle file input/output operations, and a `DataProcessor` class instance is created to process the data. The processed data is stored in various tables, which are updated as the script loops through each row of a results table.

If certain processed data files do not exist, the script will load the data for the current row, combine EEG and score data, decimate the fiber data, process the fiber data, get corresponding EEG data chunks for the fiber chunks, add scores to the fiber chunks, and save the processed data.

Then, the script calculates the ZT and mouse values for the current row, and adds these to a score averages table. Transients and transitions are calculated and added to the relevant tables.

If the user has chosen to plot, the `Plotter` class is used to generate plots of the fiber chunks with scores, a summary plot, and a plot of data chunks.

Finally, the script saves the combined data, combined transients, and combined transitions as .csv files, and generates additional plots if the user has chosen to do so.

## How to Use

To use the script, simply run it in MATLAB. When prompted, input 'yes' or 'no' to choose whether to generate plots during execution. The script will take care of the rest.

Please ensure that you have the necessary data files and folders as described in the script comments.

**Note:** This script requires MATLAB and uses features that may require specific MATLAB toolboxes.

## Contribution

Feel free to fork this repository and improve. Here are some improvements that you could make:

- Handle different data formats
- Modify the plotting code to create different plots

Pull requests are welcome.
