import numpy as np
import scipy.ndimage as ndi
import skimage as ski
import pandas as pd
import argparse
import sys
import os
sys.path.insert(0, '.')
from segmentation_under_agarose import *
import collections
import csv


#Globals
MIN_AREA_RATIO = 0.75

def parse_metadata():
    parser = argparse.ArgumentParser()
    parser.add_argument('metadata_csv', help = "path to csv file of metadata")
    parser.add_argument('--start', help = ("row of movie to start with "
    "(0-based)"), default = 0, type = int)
    parser.add_argument('--stop', help = ("row of last movie to process "
    "(0-based)"), default = -1, type = int)
    args = parser.parse_args()
    
    data_info = []
    with open(args.metadata_csv, 'r', encoding = 'utf-8-sig') as csv_file:
        reader = csv.DictReader(csv_file, delimiter = ",")
        for row in reader:
            data_info.append(initialize_metadata(row))
    return data_info, args


def initialize_metadata(row):
    grayscale_range = (int(row['range_min']), int(row['range_max']))

    info = dict()
    info['movie'] = int(row['movie'])
    info['grayscale_range'] = grayscale_range
    info['inpath'] = row['filepath']
    info['savedir'] = row['savedir']
    info['savepath'] = setup_save_folder(row['savedir'], row['filename'])
    info['savebase'] = os.path.splitext(row['filename'])[0]
    info['trial'] = int(row['trial'])
    return info


def setup_save_folder(savedir, filename):

    savename = '{}_out'.format(os.path.splitext(filename)[0])
    savepath = os.path.join(savedir, savename)

    if not os.path.isdir(savepath):
        os.makedirs(savepath)
    return savepath


def segment_movie(info):
    
    im = ski.io.imread(os.path.join(info['inpath'],
        "{}.tif".format(info['savebase'])))
    im = ski.exposure.rescale_intensity(im, in_range =
            info['grayscale_range'], out_range = (0, 1)).astype(np.float32)

    seeds = pd.read_csv(os.path.join(info['inpath'], 
                                info['savebase'] + '_initial_points.csv'))
    seed_coords = (seeds['Y'].astype(int), seeds['X'].astype(int))
    seed_labels = seeds[' '].astype(int)
    seed_im = initialize_seed_image(seed_coords, seed_labels, im[0].shape)

    props = ('label', 'area', 'centroid', 'bbox', 'coords')
    #Process first image
    mask = binarize_image(im[0])
    labels = perform_watershed(im[0], mask, seed_im)
    new_labels = fix_segmentation_errors(im[0], labels)
    ski.io.imsave(os.path.join(info['savepath'], 't000.jpg'),
                  ski.color.label2rgb(new_labels))
    cell_data = pd.DataFrame(ski.measure.regionprops_table(new_labels,
        properties = props))
    cell_data.insert(0, 'frame', 0)
    cell_data['cell'] = cell_data['label']

    seed_coords = (cell_data['centroid-0'].astype(int),
                   cell_data['centroid-1'].astype(int))
    seed_labels = cell_data['cell'].astype(int)

    prev_cell_data = cell_data.loc[cell_data['frame'] == 0].set_index('cell')

    #Segement the rest of the movie
    for t in range(1, len(im)):
        print("Processing frame {}...".format(t))
        #Process subsequent images
        mask = binarize_image(im[t])
        sharpened = ski.filters.unsharp_mask(im[t], radius=5, amount =3)
        #Label each connected component in the mask (may contain multiple cells)
        conn_comps = ski.morphology.label(mask)
        conn_comp_props = get_regionprops_table(conn_comps)
        good_seed_coords = find_better_seeds(seed_coords, sharpened)
        #Check if any connected components contain more than one seed
        # -- these are candidates for watershed
        seeds_in_conn_comp = []
        seeds_as_tuple = [(a,b) for a,b in zip(*good_seed_coords)]
        for idx, cc in conn_comp_props.iterrows():
            coords_as_tuple = [(a[0], a[1]) for a in cc['coords']]
            seeds_in_conn_comp.append(sum([s in coords_as_tuple for s in
                seeds_as_tuple]))
        if any([s > 1 for s in seeds_in_conn_comp]): #Merged cc, run watershed
            seed_im = initialize_seed_image(good_seed_coords, seed_labels, im[t].shape)
            labels = perform_watershed(sharpened, mask, seed_im)
        else:
            labels = ski.measure.label(mask)
        candidate_cell_data = get_regionprops_table(labels)
        #For each label figure out which connected component it is in
#        for label, row in candidate_cell_data.iterrows():
#            candidate_cell_data.at[label, 'conn-comp'] = conn_comps[round(row['centroid-0']),
#                                                    round(row['centroid-1'])]
#
#        conn_comp_to_fix = identify_watershed_errors(candidate_cell_data, prev_cell_data)
#        if conn_comp_to_fix:
#            raise Exception("still conn-comps to fix!")

        new_labels = fix_segmentation_errors(im[t], labels, thresh = 0.97)
        ski.io.imsave(os.path.join(info['savepath'], 't{:03d}.jpg'.format(t)),
                      ski.color.label2rgb(new_labels))
        new_cell_data = pd.DataFrame(ski.measure.regionprops_table(new_labels,
                                         properties = props))
        new_cell_data.insert(0, "frame", t)
        new_cell_data.insert(len(new_cell_data.columns), "cell", np.nan)
        #Check for any new objects that were in the mask but did not have a seed in the watershed
        for cell_idx, cell in cell_data.loc[cell_data['frame']==t-1].iterrows():
            for new_cell_idx, new_cell in new_cell_data.iterrows():
                if ((round(cell['centroid-0']), round(cell['centroid-1'])) in
                     [tuple(i) for i in new_cell['coords']]):
                    new_cell_data.at[new_cell_idx, 'cell'] = cell['cell']
        #For any remainin cells add a new cell id
        n_new_cells = new_cell_data.loc[new_cell_data['cell'].isna(),
                                        'cell'].shape[0]
        max_current_cells = max(cell_data.loc[cell_data['frame'] < t, 'cell'])
        new_cell_data.loc[new_cell_data['cell'].isna(),
                          'cell'] = np.arange(max_current_cells + 1,
                                              max_current_cells + 1 +
                                              n_new_cells)

        cell_data = pd.concat([cell_data, new_cell_data], ignore_index = True)
        #Update seed coordinates
        seed_coords = (new_cell_data['centroid-0'].astype(int),
                   new_cell_data['centroid-1'].astype(int))
        seed_labels = new_cell_data['cell'].astype(int)
        prev_cell_data = cell_data.loc[cell_data['frame'] == t].set_index('cell')
    
    #Save the cell_data structure
    cell_data.to_pickle(os.path.join(info['savedir'], 
        '{}_cell_data.pkl'.format(info['savebase'])))


def main():
    data_info, args = parse_metadata()
    if args.stop == -1:
        args.stop = len(data_info)
    for info in data_info[args.start:args.stop]:
        print("Processing {}...".format(info['savebase']))
        segment_movie(info)
        print('done!')
    print('All done!')

if __name__ == '__main__':
    main()
