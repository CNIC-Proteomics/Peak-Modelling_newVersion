#!/usr/bin/python

# -*- coding: utf-8 -*-

# Module metadata variables
__author__ = "Andrea Laguillo Gómez"
__credits__ = ["Andrea Laguillo Gómez", "Jose Rodriguez", "Samuel Lozano-Juarez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "1.0.0"
__maintainer__ = "Jose Rodriguez"
__email__ = "andrea.laguillo@cnic.es;jmrodriguezc@cnic.es"
__status__ = "Development"

# import modules
import os
import sys
import argparse
import configparser
import glob
import logging
import pandas as pd
import pathlib
import numpy as np
import concurrent.futures
from tqdm import tqdm
from itertools import repeat
from bokeh.plotting import figure, save, output_file
from bokeh import palettes
pd.options.mode.chained_assignment = None  # default='warn'

# TODO if empty space in one column, ignore whole row (all modules)

###################
# Local functions #
###################
def plot_target_decoys(x1,y1,x2,y2,name):

    ### Plot target and decoys counting ###

    p=figure(title='Target and Decoy Histogram',x_axis_label='DM',y_axis_label='count',  width=1600,height=700, tools="pan,xzoom_in,xzoom_out,ywheel_zoom,box_zoom,reset,save,undo,hover",
             tooltips=[('Name', '$name'), ("DM", "$x"), ("Y", "$y")])
    p.line(x1,y1,legend_label='Target',color=palettes.Light3[0],name='Target',line_width=2)
    p.line(x2,y2,legend_label='Decoy',color=palettes.Light3[1],name='Decoy',line_width=2)
    file=os.path.join(name,'target_decoy.html')
    p.legend.click_policy='hide'
    output_file(file)
    save(p)


def concatInfiles(infile):
    '''    
    Load dataframe from path, changing Filename field to category
    '''
   
    df = pd.read_feather(infile)
    df['Filename'] = infile
    df['Filename'] = df['Filename'].astype('category')
    return df

def generate_histogram(df, bin_width, dm0, dm1):
    '''
    Group by DeltaMass into bins of the size specified.
    '''
    
    def _decimal_places(x):
        s = str(x)
        if not '.' in s:
            return 0
        return len(s) - s.index('.') - 1
    
    # sort by deltamass
    df.sort_values(by=['cal_dm_mh'], inplace=True)
    df.reset_index(drop=True, inplace=True)
    
    # make bins
    bins = list(np.arange(dm0,
                          dm1+bin_width,
                          bin_width))
    bins = [round(x, _decimal_places(bin_width)) for x in bins]
    df['bin'] = pd.cut(df['cal_dm_mh'], bins=bins)
    
    # make histogram table
    bins_df = df['bin'].value_counts().to_frame().rename(columns = {'bin':'count'})
    bins_df.insert(0, 'bin', bins_df.index)
    bins_df.insert(1, 'midpoint', bins_df['bin'].apply(lambda x: x.mid))
    bins_df.reset_index(drop=True, inplace=True)
    bins_df.sort_values(by=['bin'], inplace=True)
    bins_df.reset_index(drop=True, inplace=True)
    
    return df, bins_df

def _parallelRegression(list_subset, smoothed, second_derivative):
    '''
    Calculate the linear regression line and return the slope
    (and the intercept, if called with smooth == True).
    '''
    bin_subset = list_subset[0]
    pos = list_subset[1]
    result = linear_regression(bin_subset, smoothed, second_derivative)
    if result is tuple:
        result = linear_regression(bin_subset, smoothed, second_derivative)[0]
    result = [pos, result]
    return result

def linear_regression(bin_subset, smoothed, second_derivative):
    '''
    Calculate the linear regression line and return the slope
    (and the intercept, if called with smooth == True).
    '''
    # ignore special cases at beginning and end and use a
    # linear regression function for the rest
    x_list = bin_subset['midpoint'].tolist()
    if smoothed:
        y_list = bin_subset['smooth_count'].tolist()
    elif not second_derivative:
        y_list = bin_subset['count'].tolist()
    else:
        y_list = bin_subset['slope1'].tolist()
    sum1, sum2 = 0, 0
    for i in range(len(x_list)):
        sum1 += (x_list[i] - np.mean(x_list)) * (y_list[i] - np.mean(y_list))
        sum2 += (x_list[i] - np.mean(x_list)) ** 2
    working_slope = sum1 / sum2
    intercept = np.mean(y_list) - working_slope*np.mean(x_list)
    if smoothed or second_derivative:
        return working_slope
    else:
        return working_slope, intercept

def _parallelSmoothing(list_subset):
    bin_subset = list_subset[0]
    pos = list_subset[1]
    working_slope, intercept = linear_regression(bin_subset, False, False)
    smoothed = [pos, working_slope, intercept]
    return smoothed

def smoothing(bins_df, spoints):
    '''
    Calculate the slope (first derivative) for each bin. Calculate new smoothed
    value for the midpoint using the linear regression line.
    '''
    bins_df['smooth_count'] = None
    bin_subsets = []
    for i in range(spoints, len(bins_df)-spoints):
        bin_subset = bins_df[i-spoints:i+spoints+1]
        bin_subsets.append([bin_subset, i])
    logging.info("\tSmoothing...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:   
        intercepts = list(tqdm(executor.map(_parallelSmoothing, bin_subsets, chunksize=1000),
                               total=len(bin_subsets)))
    intercepts = pd.DataFrame(intercepts, columns=['i', 'working_slope', 'intercept'])
    # TODO: add to bins_df.loc
    bins_df = pd.merge(bins_df, intercepts, left_index=True, right_on='i', how='outer')
    bins_df.reset_index(drop=True, inplace=True)
    bins_df['smooth_count'] = bins_df.apply(lambda x: None if np.isnan(x['intercept'])
                                                           else x['intercept'] + (x['working_slope']*x['midpoint']), axis = 1)
    bins_df[["smooth_count"]] = bins_df[["smooth_count"]].apply(pd.to_numeric)
    bins_df = bins_df.drop(['i', 'working_slope', 'intercept'], axis=1)
    return bins_df

def first_derivative(bins_df, points, spoints):
    '''
    Calculate the slope (first derivative) for each bin.
    Returns the slope of the linear regression line through data points in
    known_y's and known_x's. The slope is the vertical distance divided by
    the horizontal distance between any two points on the line, which is the
    rate of change along the regression line.
    Known_y's  Bins. An array or cell range of numeric dependent data points.
    Known_x's  Count. The set of independent data points.
    '''
    if spoints > 0: #smoothing
        bins_df = smoothing(bins_df, spoints)
        j = 2
    else: #no smoothing
        j = 1
    logging.info("\tFirst derivative...")
    bin_subsets = []
    for i in range(points*j, len(bins_df)-points*j):
        bin_subset = bins_df[i-points:i+points+1]
        bin_subsets.append([bin_subset, i])
    
    if spoints > 0:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:   
            firsts = list(tqdm(executor.map(_parallelRegression,
                                            bin_subsets, repeat(True), repeat(False),
                                            chunksize=1000), total=len(bin_subsets)))
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:   
            firsts = list(tqdm(executor.map(_parallelRegression,
                                            bin_subsets, repeat(False), repeat(False),
                                            chunksize=1000), total=len(bin_subsets)))
    firsts = pd.DataFrame(firsts, columns=['i', 'slope1'])
    bins_df = pd.merge(bins_df, firsts, left_index=True, right_on='i', how='outer').replace({np.nan: None})
    bins_df = bins_df.drop('i', axis=1)
    bins_df.reset_index(drop=True, inplace=True)
    return bins_df

def second_derivative(bins_df, points, spoints):
    '''
    Calculate the second derivative for each bin.
    '''
    if spoints > 0: #smoothed
        j = 3
    else: #not smoothed
        j = 2
    logging.info("\tSecond derivative...")
    bin_subsets = []
    for i in range(points*j, len(bins_df)-points*j):
        bin_subset = bins_df[i-points:i+points+1]
        bin_subsets.append([bin_subset, i])
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:
        seconds = list(tqdm(executor.map(_parallelRegression,
                                        bin_subsets, repeat(False), repeat(True),
                                        chunksize=1000), total=len(bin_subsets)))
    seconds = pd.DataFrame(seconds, columns=['i', 'slope2'])
    bins_df = pd.merge(bins_df, seconds, left_index=True, right_on='i', how='outer').replace({np.nan: None})
    bins_df = bins_df.drop('i', axis=1)
    bins_df.reset_index(drop=True, inplace=True)
    return bins_df




#################
# Main function #
#################

def main(args):
    '''
    Main function
    '''
    
    #Main variables
    bins = float(config._sections['PeakModeller']['bins'])
    slope_points = int(config._sections['PeakModeller']['slope_points'])
    smooth_points = int(config._sections['PeakModeller']['smooth_points'])

    dm0 = int(config._sections['PeakSelector']['dm0'])
    dm1 = int(config._sections['PeakSelector']['dm1'])
    
    model_table = pd.read_csv(args.modelling_table, sep="\t")
    
    logging.info("Reading input file list...")

    if '*' in args.infile: # wildcard
        infiles = glob.glob(args.infile)
        h_outfile = os.path.join(os.path.dirname(args.infile), 'PeakModeller_DMHistogram.tsv')
        t_outfile = os.path.join(os.path.dirname(args.infile), 'PeakModeller_DMTable.feather')
    else:
        infiles = glob.glob(args.infile)
        h_outfile = args.infile[:-4] + '_DMHistogram.tsv'
        t_outfile = args.infile[:-4] + '_DMTable.feather'
    for i in infiles:
        logging.info('\t' + str(os.path.basename(i)))


    # Samuel: the input files are not concatenated. Instead they are stored inside a dictionary
    df_dict = dict()
    # iter over the different group values in the modelling histogram
    for model_group in model_table["Histogram_Modelling_Group"].unique():
        #select the files corresponding to the same modelling group
        filtered_files = list(model_table[model_table["Histogram_Modelling_Group"] == model_group]["Spectrum_File"])
        #from all the infiles list keep only those present in the filtered_files list
        infiles_filtered = [f for f in infiles if any(sub in f for sub in filtered_files)]
        #concat those files
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.n_workers) as executor:            
            df_iterator = executor.map(concatInfiles, infiles_filtered)
        df = pd.concat(df_iterator)
        df.reset_index(drop=True, inplace=True)
        #save the concat df in the dict, using the group name (and path of first infile) as key
        df_dict[os.path.join(os.path.dirname(infiles_filtered[0]),f"Group_{model_group}")] = df
    

    try:
        # Samuel: iter over each dataframe and create the output for each dataframe
        for output_path,df in df_dict.items():
            for j in range(3):

                if j==0:

                    if separate_modelling==True:
                        bools=[True if 'DECOY' not in i else False for i in df['protein']]
                        df_copy=df.loc[bools].reset_index(drop=True)
                        h_outfile = output_path + '_DMHistogram_target.tsv' # Samuel: update the line to create one different output name per df
                        t_outfile = output_path + '_DMTable_target.feather' # Samuel: update the line to create one different output name per df
                    else:
                        continue

                elif j==1:

                    if separate_modelling==True:
                        bools=[True if 'DECOY' in i else False for i in df['protein']]
                        df_copy=df.loc[bools].reset_index(drop=True)
                        h_outfile = output_path + '_DMHistogram_decoy.tsv' # Samuel: update the line to create one different output name per df
                        t_outfile = output_path + '_DMTable_decoy.feather' # Samuel: update the line to create one different output name per df
                    else:
                        continue

                else:
                    h_outfile = output_path + '_DMHistogram.tsv' # Samuel: update the line to create one different output name per df
                    t_outfile = output_path + '_DMTable.feather' # Samuel: update the line to create one different output name per df
                    df_copy=df
                    # pass

            
                logging.info(f"Generating DMHistogram for {output_path.split(os.sep)[-1]}...")
                # make bins
                df_copy, bins_df = generate_histogram(df_copy, bins, dm0, dm1)
                bins_df = first_derivative(bins_df, #does 1st smoothing pass and 2nd normal pass
                                        slope_points//2,
                                        smooth_points//2)  
                bins_df = second_derivative(bins_df,
                                            slope_points//2,
                                            smooth_points//2)
                
                if j==0:
                    midpoint_target=bins_df['midpoint']
                    count_target=bins_df['count']
                elif j==1:
                    midpoint_decoy=bins_df['midpoint']
                    count_decoy=bins_df['count']
                

                #### ESTO SE PUEDE SACAR
                    # check which bins pass
                logging.info("Writing output files...")
                # write DMhistogram
                bins_df.to_csv(h_outfile, index=False, sep='\t', encoding='utf-8')
                # write DMtable (input for PeakSelector)
                df_copy = df_copy.astype(str)
                df_copy.to_feather(t_outfile)

            if separate_modelling == True and j >1:
                plot_target_decoys(midpoint_target,count_target,midpoint_decoy,count_decoy,os.path.dirname(args.infile))
            logging.info("Peak Modelling finished")

    except TypeError:

        print('Check if in your directory there are any .feather files that you do not want to modellate')



if __name__ == '__main__':

    # parse arguments
    parser = argparse.ArgumentParser(
        description='Peak Modeller',
        epilog='''
        Example:
            python PeakModeller.py

        ''')
        
    defaultconfig = os.path.join(os.path.dirname(__file__), "config/SHIFTS.ini")
    
    # Samuel: addition of a parameter to read the modelling table
    parser.add_argument('-t', '--modelling_table', required=True, help='Path to the file with relation among Files and Modelling Groups')
    
    parser.add_argument('-i', '--infile', required=True, help='Input file with the peak file(s) to be filtered')
    parser.add_argument('-c', '--config', default=defaultconfig, help='Path to custom config.ini file')
    # TODO: output file path

    parser.add_argument('-b', '--bins', help='Width of the bins')
    parser.add_argument('-p', '--slope_points', help='Number of points (bins) to use for slope calculation')
    parser.add_argument('-s', '--smooth_points', help='Number of points (bins) to use for smoothing')

    parser.add_argument('-w',  '--n_workers', type=int, default=4, help='Number of threads/n_workers (default: %(default)s)')    
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()
    
    # parse config
    config = configparser.ConfigParser(inline_comment_prefixes='#')
    config.read(args.config)
    if args.bins is not None:
        config.set('PeakModeller', 'bins', str(args.bins))
        config.set('Logging', 'create_ini', '1')
    if args.slope_points is not None:
        config.set('PeakModeller', 'slope_points', str(args.slope_points))
        config.set('Logging', 'create_ini', '1')
    if args.smooth_points is not None:
        config.set('PeakModeller', 'smooth_points', str(args.smooth_points))
        config.set('Logging', 'create_ini', '1')
    # if something is changed, write a copy of ini
    if config.getint('Logging', 'create_ini') == 1:
        with open(os.path.dirname(args.infile) + '/SHIFTS.ini', 'w') as newconfig:
            config.write(newconfig)

    separate_modelling=bool(int(config._sections['PeakModeller']['separate_modelling']))
        
    # logging debug level. By default, info level
    if '*' in args.infile: # wildcard
        log_file = os.path.join(os.path.dirname(args.infile), 'PeakModeller_log.txt')
        log_file_debug = os.path.join(os.path.dirname(args.infile), 'PeakModeller_log_debug.txt')
    else:
        log_file = args.infile[:-4] + '_log.txt'
        log_file_debug = args.infile[:-4] + '_log_debug.txt'
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file_debug),
                                      logging.StreamHandler()])
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p',
                            handlers=[logging.FileHandler(log_file),
                                      logging.StreamHandler()])

    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    # merge all the feather files into one and save it
    logging.info('merging all feather files')
    files_path = pathlib.Path(os.path.dirname(args.infile))
    feather_files = list(files_path.glob("*DMTable.feather"))
    dfs = [pd.read_feather(f) for f in feather_files]
    combined_df = pd.concat(dfs, ignore_index=True)
    combined_df.to_feather(os.path.join(os.path.dirname(args.infile), 'DMTable.feather'))
    logging.info('end script')
