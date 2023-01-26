import numpy as np
import scipy.ndimage as ndi
import skimage as ski
import pandas as pd
import collections

#Default properties for the regionprops tables
DEFAULT_PROPS = ('label', 'area', 'centroid', 'bbox', 'coords')
MIN_AREA_RATIO = 0.75

def identify_watershed_errors(candidate_cell_data, prev_cell_data, 
        thresh=MIN_AREA_RATIO):
    '''Identify connected components in which watershed failed
    Input:
        - candidate_cell_data: regionprops dataframe of potential
        cells
        - prev_cell_data: regionprops dataframe of cells in previous frame
        - thresh: minimum area ratio that cell can change from one
        image to the next. (default to global MIN_AREA_RATIO)
    Output:
        - conn_comp_to_fix: list of connected component ids to fix
    '''
    #Identify cells present in current and previous frame
    matched_labels = set(prev_cell_data.index).intersection(
            candidate_cell_data['label'])
    labels_to_fix = []
    for label in matched_labels:
        ratio = candidate_cell_data.at[label, 'area'] \
                / prev_cell_data.at[label, 'area']
        if (ratio < thresh):
            print("Watershed error for label {}: ratio {}".format(label,
                ratio))
            labels_to_fix.append(label)
    if labels_to_fix:
        conn_comp_to_fix = list(np.unique(
                candidate_cell_data.loc[labels_to_fix,
                                        'conn-comp'].to_numpy(dtype=int)))
    else:
        conn_comp_to_fix = []
    return(conn_comp_to_fix)


def get_nearest_points(p, p_list):
    '''Sort list of points by their closeness to a point p.
    Input:
        - p: a point (2 coordinates)
        - p_list: a 2-tuple of N points, i.e. 2 lists
        of N points, coordinates for point i given by
        (p_list[0][i],p_list[1][i])
    Output:
        - sorted_list: p_list sorted by distance from p
    '''
    distances = (p_list[0] - p[0])**2 + (p_list[1] - p[1])**2

    return((p_list[0][np.argsort(distances)],
            p_list[1][np.argsort(distances)]))


def different_labels(center_value, neighbor_value, *args):
    '''Helper function for find_neighbors'''
    return (center_value != neighbor_value).astype(float)


def find_neighbors(labels):
    '''Get a dict of neighbor labels.
    From Juan Nunez-Iglesias https://stackoverflow.com/a/72456945
    Return a dict containing a list of neighbors for each label
    in a label matrix (neighbors in an 8-neighborhood).
    Input:
    - labels: a label matrix
    Output:
    - pairs: a dict of lists, keys are labels, lists contain 
    the labels of the neighbors of the key label
    '''
    g, nodes = ski.graph.pixel_graph(labels,
                                    mask = labels.astype(bool),
                                    edge_function = different_labels,
                                    connectivity = 2)
    g.eliminate_zeros()
    coo = g.tocoo()
    center_coords = nodes[coo.row]
    neighbor_coords = nodes[coo.col]
    center_values = labels.ravel()[center_coords]
    neighbor_values = labels.ravel()[neighbor_coords]
    pairs = collections.defaultdict(list)
    #Iterate over center, neighbor pairs (set makes them unique)
    for i,j in set([(a,b) for a,b in zip(center_values,
                                        neighbor_values)]):
        pairs[i].append(j)
    return pairs


def identify_seg_errors(labels, thresh= 0.95):
    '''Identify objects that have a segmentation error
    Input:
        - labels: label matrix
        - thresh: min ratio of area post/pre removing of border near hole
        to qualify as a "problem area" (defualt:0.95)
    Output:
        - need_to_fix: list of labels from the label matrix that must be fixed
    '''
    mask = labels > 0
    old_rp = {x.label: x for x in ski.measure.regionprops(labels)}
    #Remove the around the hole and recalculate region properties
    opened_mask = ski.morphology.opening(mask, ski.morphology.square(15))
    new_labels = np.copy(labels)
    new_labels[np.logical_not(opened_mask)] = 0
    new_rp = {x.label:x for x in ski.measure.regionprops(new_labels)}

    need_to_fix = []
    for label, r in new_rp.items():
        relative_area = r.area_convex / old_rp[label].area_convex
        if relative_area < thresh:
            need_to_fix.append(label)

    return need_to_fix


def fix_outline_active_contour(im, labels, idx):
    '''Fix the outline of a cell using active contours.
    Input:
        - im: the original image
        - labels: labeled image corresponding to im
        - idx: the label in the label image corresponding to the cell to fix
    Returns:
        - new_labels: the fixed labeled image
    '''
    single_blob = labels == idx
    #Initialize snake with a slightly expanded convex hull, with
    # a given number of resampled points
    conv_hull = ski.morphology.convex_hull_object(single_blob)
    region = ski.measure.regionprops(ski.morphology.label(conv_hull))[0]
    contour = ski.measure.find_contours(conv_hull)[0]
    poly = ski.measure.approximate_polygon(contour, 1)
    new_points = sample_polygon_contour(poly[:-1,:], 300)
    #Resize the polygon a bit
    new_points = 1.01 * (new_points - region.centroid) + region.centroid

    snake = ski.segmentation.active_contour(im, new_points,
                                            alpha=0.005,
                                            beta=.1,
                                            gamma=0.01,
                                            w_line=1,
                                            w_edge=1)

    polymask = ski.draw.polygon2mask(labels.shape, snake)
    polymask = ski.morphology.dilation(polymask, ski.morphology.square(11))

    #Update the label matrix (keeping other labels the same)
    new_labels = np.zeros(labels.shape, dtype=int)
    new_obj = np.ravel_multi_index(np.nonzero(polymask), labels.shape)
    new_labels.flat[new_obj] = idx
    #Update all other values (other than idx)
    for label in set(np.unique(labels)) - set((0,idx)):
        new_labels[labels == label] = label
    return new_labels


def sample_polygon_contour(p, n_points):
    '''resample polygon outline into n_points evenly spaced points
    Credit: Stack Overflow code: https://stackoverflow.com/a/42032457

    Input:
        - p: Nx2 array of points on the polygon
        - n_points: number of points to resample
    Output:
        - new_p: n_points x 2 array of points on the polygon
    '''

    vectors = np.roll(p, -1, axis=0) - p
    si = np.linalg.norm(vectors, axis=1)
    di = np.linspace(0, np.sum(si), n_points, endpoint=False)
    new_p = []
    #"Walk" along the polygon by finding which segment the current
    # point d is located on
    for d in di:
        for i,s in enumerate(si):
            if d > s:
                d -= s
            else:
                break
        #Factor of proportionality (how far along the segment are you)
        alpha = d / s
        new_p.append([p[i, 0] + alpha * vectors[i, 0],
                      p[i, 1] + alpha * vectors[i, 1]])

    return np.array(new_p)


def pad_bbox(bbox, imshape, padding=10):
    '''Add padding to a bounding box as much as possible
    Input:
        - bbox: [row_min, col_min, row_max, col_max]
        - imshape: the (row, col) of the original image
        - padding: number of pixels to add to each side of the box
    Output:
        - padded_bbox: the padded bbox (same form)
    '''
    pad_top = min(padding, bbox[0])
    pad_bot = min(padding, imshape[0] - (bbox[2] + 1))
    pad_left = min(padding, bbox[1])
    pad_right = min(padding, imshape[1] - (bbox[3] + 1))
    padded_bbox = [bbox[0] - pad_top, bbox[1] - pad_left,
                   bbox[2] + pad_bot, bbox[3] + pad_right]

    return padded_bbox


def fix_segmentation_errors(im, labels, thresh=0.95):
    '''Fix segmentation errors for any cells in the image using active contours
    Input:
        - im: the raw image
        - labels: the label matrix for this image
        - thresh: min ratio of area post/pre removing of border near hole
        to qualify as a "problem area" (defualt:0.95)
    Output:
        - labels_fixed: the fixed label matrix for this image
    '''
    rp = {x.label:x for x in ski.measure.regionprops(labels)}
    need_to_fix = identify_seg_errors(labels, thresh=thresh)
    labels_fixed = np.copy(labels)
    print("Fixing segmentation errors in {} objects...".format(len(need_to_fix)))
    for label in need_to_fix:
        bbox = pad_bbox(rp[label].bbox, labels.shape)
        fixed = fix_outline_active_contour(im[bbox[0]:bbox[2], bbox[1]:bbox[3]],
                                        labels[bbox[0]:bbox[2], bbox[1]:bbox[3]],
                                        label)
        labels_fixed[bbox[0]:bbox[2], bbox[1]:bbox[3]] = fixed
    return labels_fixed


def binarize_image(im):
    '''Make a binary image.
    Uses edge detection and morphology
    Input:
        - im: a (2D) image
    Output:
        - mask: a binary mask of the image
    '''
    
    #Gaussian blur
    g = ndi.gaussian_filter(im, 3)
    # Edge detection
    edges = ski.filters.scharr(g)
    # Binary mask by Otsu thresholding
    mask = edges >= ski.filters.threshold_otsu(edges)
    #Morphology to close contours, remove small objects, and clear borders
    mask = ski.morphology.binary_closing(mask, ski.morphology.square(11))
    mask = ndi.binary_fill_holes(mask)
    mask = ski.morphology.remove_small_objects(mask, 10000)
    mask = ski.segmentation.clear_border(mask)

    return mask


def find_better_seeds(initial_seeds, im, n = 0):
    '''Pick seeds based on local minima in the image (better for watershed)
    Input: 
        - initial_seeds: 2-tuple of lists of N coordinates
        - im: the image the seeds correspond to
        - n: the nth closest local minima will be returned (default 0)
    Output: 
        - new_seeds: list of better seeds, same format as initial_seeds
    '''
    min_filter = ndi.minimum_filter(im, size=4)
    is_local_min =ski.morphology.erosion(min_filter == im,
                        ski.morphology.disk(2))
    local_min_points = np.nonzero(is_local_min)
    new_seeds = ([],[])
    for p in zip(*initial_seeds):
        nearest_minima = get_nearest_points(p, local_min_points)
        new_seeds[0].append(nearest_minima[0][n])
        new_seeds[1].append(nearest_minima[1][n])
    return(new_seeds)



def initialize_seed_image(seed_coords, seed_labels, im_shape, jitter=False):
    '''Create a seed image from a list of coordinates
    Input:
        - seed_coords: tuple of lists/arrays of the (row_coords, col_coords)
        - seed_labels: list/array of integer indices marking each seed
        (order matches seed_coords)
        - im_shape: [n_rows, n_cols]
        - jitter: boolean if the positions of the seeds should be randomized
    Output:
        - seed_im: image with a number for each seed (1-based indexing)
    '''
    if len(seed_coords[0]) < 256:
        dtype = np.uint8
    elif len(seed_coords[0]) < 65536:
        dtype = np.uint16
    else:
        dtype = np.uint32

    seed_im = np.zeros(im_shape, dtype = dtype)
    if jitter:
        seed_rows = seed_coords[0] + 3*np.random.randint(-7,8,
                size=len(seed_coords[0]))
        seed_cols = seed_coords[1] + 3*np.random.randint(-7,8,
                size=len(seed_coords[1]))
        seed_coords = (seed_rows, seed_cols)
    flat_coords = np.ravel_multi_index(seed_coords, im_shape)
    seed_im.flat[flat_coords] = np.array(seed_labels, dtype=dtype)
    
    return seed_im


def perform_watershed(im, mask, seed_im):
    '''Run the watershed with appropriate seed images
    Input:
        - im: the raw image
        - mask: the binary mask of the image
        - seed_im: seeds marked with integer labels (one pixel per seed)
    Output:
        - expanded: the watershed result (expanded to match the original mask)
    '''
    print("Running watershed...")    
    #Exterior phase halo messes up the watershed, so remove it here
    watershed_mask = ski.morphology.erosion(mask, ski.morphology.square(15))
    print(sum(np.nonzero(watershed_mask)))
    ws = ski.segmentation.watershed(im, seed_im, mask = watershed_mask,
            connectivity=1)
    print(sum(np.nonzero(ws)))
    expanded = ski.segmentation.expand_labels(ws, distance=100)
    print(sum(np.nonzero(expanded)))
    #Restore the original mask regions, now with labels from watershed
    expanded[np.logical_not(mask)] = 0
    print(sum(np.nonzero(expanded)))

    return expanded


def get_regionprops_table(labels, props=DEFAULT_PROPS):
    '''Wrapper to get a pandas table of regionprops
    Input:
        - labels: label matrix
        - props: tuple of strings indicating properties to include (see docs)
    Output:
        - df: a pandas DataFrame of the regionprops.
        The 'label' field is also the index for the DataFrame
    '''
    df = pd.DataFrame(ski.measure.regionprops_table(labels,
                                                    properties=props))
    df.set_index('label', drop=False, inplace=True)
    return df


def bbox_crop(im, bbox):
    '''Crop an image using a bounding box
    Input:
        - im: 2d image
        - bbox: [min_row, min_col, max_row, max_col]
    Output:
        - cropped: cropped image
    '''
    return im[bbox[0]:bbox[2], bbox[1]:bbox[3]]
