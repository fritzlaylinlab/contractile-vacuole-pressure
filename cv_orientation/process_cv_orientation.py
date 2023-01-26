import numpy as np
import pandas as pd
import os
import argparse

#Manifest of the columns present in the output file
COLUMNS = ['trial', #biological replicate (different days)
           'movie', #index of the movie, corresponding to metadata file
           'cell', #1-based cell index within the movie
           'frame', #1=based frame in the timeseries when this cv pump occurred
           'cell_area', #in pixels
           'cv_area', #in pixels 
           'cell_centroid_r', 'cell_centroid_c', #pixel row,col respectively
           'cv_centroid_r', 'cv_centroid_c', #pixel row,col respectively
           'cell_speed', #in units of pixels/frame
           'vel_dx', 'vel_dy', #components of unit vector of cell velocity
           'cv_dx', 'cv_dy', #components of unit vector pointing to cv from centroid
           'cv_angle'] #angle in radians b/ween cv-centroid and cell velocity

def match_cell_to_cv(cv_data, cell_data):
    '''Match each cv to the cell it is contained in.
    Input:
        - cv_data: dataframe about each cv pump, including its
            position in pixel coordinates (X is column and Y is row), and the 
            1-based Slice in which it occurs
        - cell_data: dataframe about each cell in the movie over time,
            including the 0-based frame index, the 1-based cell index,
            and the pixel-based coordinates of the centroid (centroid-0 is row,
            centroid-1 is column)
    Output: 
        - cv_data: dataframe with added 'cell' column
    '''
    #Identify the cell ID that corresponds to each cv pump in a dict
    matching_cell = dict()
    for cv_idx, cv in cv_data.iterrows():
        cv_centroid = (round(cv['Y']), round(cv['X']))
        #Look only through cells present at the same timepoint as the cv
        for _, cell in cell_data.loc[
                            cell_data['frame'] == (cv['Slice'] - 1)].iterrows():
            #convert list of 2-lists into list of 2-tuples
            #(needs to be immutable for 'in' to work as expected)
            coords = list(map(tuple, cell['coords']))
            if cv_centroid in coords:
                matching_cell[cv_idx] = cell['cell']
    #Add that cell ID to the cv_data dataframe
    cv_data['cell'] = cv_data.index.map(matching_cell)
    return cv_data


def calculate_centroid_velocity(cell_data, step = 5):
    '''Calculate the velocity of the cell centroid.
    Input:
        - cell_data: dataframe of info about the cell
        - step: the step-size between timepoints for computing velocity
            (default 5 frames)
    Output:
        - cell_data: dataframe with added fields
            speed: magnitude of velocity vector (in units of pixels/frame)
            dx, dy: x and y components of velocity, normalized to norm of 1
    '''
    tmp = cell_data.copy()
    #Ensure data is in the proper order for computing differences between frames
    tmp.sort_values(by = ['cell', 'frame'], inplace=True)
    tmp['d_row'] = (tmp[['cell', 'centroid-0']]
                    .groupby('cell')
                    .transform(lambda x: x.diff(periods=step))
                    .astype(float))
    tmp['d_col'] = (tmp[['cell', 'centroid-1']]
                    .groupby('cell')
                    .transform(lambda x: x.diff(periods=step))
                    .astype(float))
    tmp['speed'] = np.sqrt(tmp['d_row']**2 + tmp['d_col']**2)
    #Normalize vector components
    tmp['dx'] =  tmp['d_col'] / tmp['speed']
    tmp['dy'] = -tmp['d_row'] / tmp['speed'] #negative due to flip in image coordinates
    #Put speed in units of pixels/frame
    tmp['speed'] = tmp['speed'] / step
    #Add computed values to cell_data
    cell_data['speed'] = tmp['speed']
    cell_data['dx'] = tmp['dx']
    cell_data['dy'] = tmp['dy']
    
    return cell_data


def calculate_cv_orientation(merged_data):
    '''Get distance and orientation of the cv rel. to centroid.
    Input:
        - merged_data: dataframe containing information about the
            cv and the centroid of each cell at different timepoints 
            (corresponding to cv pump times)
    Output:
        - merged_data: dataframe with columns added for
            cv_dx, cv_dy: components of unit vector pointing from centroid
                to cv
            cv_distance: distance (in pixels) from cell centroid to cv
            cv_angle: angle in radians between centroid-cv vector and 
                centroid velocity vector
    '''
    merged_data['cv_dx'] = pd.to_numeric(merged_data['cv_centroid_c'] 
                                        - merged_data['cell_centroid_c'])
    # Negative to account for flip from row,col to x,y
    merged_data['cv_dy'] = pd.to_numeric(-(merged_data['cv_centroid_r']
                                          - merged_data['cell_centroid_r']))
    merged_data['cv_distance'] = np.sqrt(merged_data['cv_dx']**2
                                        + merged_data['cv_dy']**2)
    #Normalize the displacement components
    merged_data['cv_dx'] = merged_data['cv_dx'] / merged_data['cv_distance']
    merged_data['cv_dy'] = merged_data['cv_dy'] / merged_data['cv_distance']
    
    #Compute the components of a rotation matrix between cv-centroid vector and 
    #centroid velocity vector
    rotx = (merged_data['vel_dx'] * merged_data['cv_dx'] 
            + merged_data['vel_dy'] * merged_data['cv_dy'])
    roty = (-merged_data['vel_dy'] * merged_data['cv_dx'] 
            + merged_data['vel_dx'] * merged_data['cv_dy'])
    merged_data['cv_angle'] = np.arctan2(roty, rotx)
    
    return merged_data


def parse_args():
    parser = argparse.ArgumentParser(
            description=("Find orientation of contractile vacuole relative"
                " to the cell's velocity."))
    parser.add_argument('metadata_file', help=("path to a file of information"
                " about the different datasets to process"))
    parser.add_argument('--outfile', help=("path to save csv output "
                      "(defaults to 'cv_orientation.csv' in same folder "
                      "as metadata file)"))

    args = parser.parse_args()
    if args.outfile is None:
        args.outfile = (os.path.join(os.path.dirname(
                                     os.path.abspath(args.metadata_file)),
                                     "cv_orientation.csv"))

    return args


def main():
    args = parse_args()
    metadata = pd.read_csv(args.metadata_file)
    #Initialize dataframe to hold data from all cells and all movies
    all_data = pd.DataFrame(columns=COLUMNS)

    for idx, row in metadata.iterrows():
        cell_data = pd.read_pickle(os.path.join(row['savedir'],
                       os.path.splitext(row['filename'])[0] + "_cell_data.pkl"))
        cv_data = pd.read_csv(os.path.join(row['filepath'],
                              os.path.splitext(row['filename'])[0] + "_cv.csv"))
        cv_data = match_cell_to_cv(cv_data, cell_data)
        cell_data = calculate_centroid_velocity(cell_data)
        #Adjust frame to be 1-based
        cell_data['frame'] = cell_data['frame'] + 1
        #Rename columns prior to joining
        cv_data = cv_data.rename(columns = {'Slice': 'frame',
                                            'X': 'cv_centroid_c',
                                            'Y': 'cv_centroid_r',
                                            'Area': 'cv_area',})
        cell_data = cell_data.rename(columns = {'centroid-0': 'cell_centroid_r',
                                                'centroid-1': 'cell_centroid_c',
                                                'speed': 'cell_speed',
                                                'dx': 'vel_dx',
                                                'dy': 'vel_dy', 
                                                'area': 'cell_area'})
        cv_data.set_index(['cell', 'frame'], inplace = True)
        cell_data.set_index(['cell', 'frame'], inplace = True)
        merged_data = cv_data.join(cell_data)
        merged_data['movie'] = idx
        merged_data['trial'] = row['trial']
        merged_data.reset_index(inplace=True)
        merged_data = calculate_cv_orientation(merged_data)
        merged_data = merged_data[COLUMNS]
        all_data = pd.concat([all_data, merged_data], ignore_index=True)

    all_data.to_csv(args.outfile)
    print(args.outfile)


if __name__ == "__main__":
    main()
