# Tracking CV Orientation

The code in this folder was used to identify the orientation of the CV in a cell relative to the cell's direction of motion. Specifically, this code

  1. Segments and tracks phase contrast images of amoeba under agarose (`run_segmentation.py` with helper functions imported from `segmentation_under_agarose.py`)
  2. Grabs the positions of contractile vacuoles from overlays in a TIFF stack manually drawn in Fiji (`overlays_to_csv_.ijm`)
  3. Processes the movement of the cell's centroid and the positions of the contractile vacuole. 
In each frame where there is a contractile vacuole, the orientation of this CV relative to the direction of centroid velocity is computed.
  4. Plots the orientations of each CV in a windrose plot.

This code is provided primarily for transparency rather than reuse, and is not intended to work right out of the box. This documentation should help make
sense of the code. 

For python scripts, the `-h` parameter can be used for minimal reference of the parameters and how to call the script.

## Dependencies

This code was run with Python 3.9.1. The package dependencies (from `pip freeze`) are specified in `required_packages.txt`

## Segmenting and Tracking Cells

**Input**: a metadata file with locations of TIFF stack movies, and those movies

**Output**: Python pickle files containing data about the positions of each cell in the movie over time.

This is accomplished by `python run_segmentation.py <metadata_file>`, where `metadata_file` is the path to a CSV file containing metadata about 
the different movies to process. An example is provided in `20221127_isolated_cells_metadata.csv`. 

### Description of metadata fields 

Each row contains information about a movie (indexed by the `movie` field). 

`filepath` and `filename` refer to the directory and filename, respectively,
where the TIFF files of the movie are stored. `savedir` indicates the directory to save the output `.pkl` files. *Note that this file contains the absolute path names
to these files on a local machine, which will not be the same for your machine!* These raw image files are not included in this repository but are available
upon request. 

`range_min` and `range_max` are 16-bit grayscale intensity values setting the minimum and maximum intensity values for rescaling the data
between 0 and 1. These values were chosen for each movie to maximize contrast without leading to over- or under-saturation. 

`trial` is an index indicates the separate days on which the movies were taken (separate biological replicates)

## Processing Cell Data

**Input**: manually annotated CV positions (annotated the frame before the CV begins a pump) as overlays on the TIFF stack movies described above; 
the same `metadata_file` as above, which will direct the program to the `.pkl` files of cell data produced in the previous step.

**Output**: a CSV file containing orientations of the CV over time, relative to the direction of motion. An example is included as `20221221_cv_orientation.csv`.

The CV positions are extracted using a Fiji macro `overlays_to_csv_.ijm` which can be run in Fiji, and will produce a CSV file of the positions of the
CV centroids. For downstream steps it is important that the scale factor on the image in Fiji be in units of pixels, so that the CV positions will also
be in units of pixels.

Next the script `process_cv_orientation.py` takes in the same metadata file as the previous step and computes the CSV output of CV orientations described below.

### Description of the CV orientation CSV fields

Each row describes a single cell at a single *frame* (timepoint) within a movie.

`trial` is the same as in the metadata file, indicating biological replicates. 

`movie` is the index of the move. Note that this is 0-based while the index in the metadata file is 1-based, so this index will be 1 less than the index in the metadata file.
At this point I am not interested in fixing this and potentially causing bugs elsewhere in the code.

`cell` is an index unique to a tracked cell in a given movie. These indexes may not be complete (i.e. the existence of cell 20 does not imply 
the existence of cells 1 - 19 in the dataset). Note that these indexes start at 1, and they are NOT unique across movies (i.e. there are a lot of cell 1s 
from all the different movies)

`cell_area` and `cv_area` describe the number of pixels in the cell/CV segmented regions (units of pixels)

`cell_centroid_r`, `cell_centroid_c`, `cv_centroid_r`, `cv_centroid_c` refer to the rows and columns (r and c) of the cell and CV in the given frame. 
Positions are in units of pixels, though they have sub-pixel resolution. The coordinate origin is at the top left of the image (as is customary)

`cell_speed` is the magnitude of the velocity vector, in units of pixels/frame. The velocity vector is computed from the x and y displacements with
a temporal separation of 5 frames. *Note that because of this 5-frame separation, for frames that were less than 5 frames from the start of a track 
the velocity, and hence the speed of the cell, will be undefined (NA)*

`vel_dx` and `vel_dy`   are the normalized x and y displacements of the centroid (with a temporal separation of 5 frames). Note these are in standard
Cartesian x and y coordinate space rather than row and column space. They are normalized such that $\mathrm{d}x^2 + \mathrm{d}y^2 = 1$. As above, when
the frame of interest is less than 5 frames from the start of the track for that cell, the velocity will be undefined. These rows were not included
in the plot in the paper.

`cv_dx` and `cv_dy` is the (x,y) position of the CV *relative to* the cell's centroid at in this frame. Again, note x,y rather than row,col. These are in
units of pixels.

`cv_angle` is the angle between the cell centroid velocity and the cell-cv displacement vector, i.e. the orientation of the CV centroid relative to the
direction of motion of the cell's centroid. Units of radians.


## Plotting CV Orientation
**Input**: the CV orientation CSV produced in the previous step

**Output**: a windrose plot with the distribution of orientations of CVs, as well as the per-trial average orientation and the average of the per-trial
orientations. 

This is accomplished by the `plot_cv_orientation.py` script, which requires the path to the CV orientation CSV file from above. The save name of the file
can be optionally specified. Note that the extension specieid there will control the type of file produced (e.g. PDF, SVG, etc.)
