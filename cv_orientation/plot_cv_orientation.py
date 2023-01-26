import numpy as np
import pandas as pd
import os
import argparse
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import date


def parse_args():
    parser = argparse.ArgumentParser(
                description="Plot the orientation of contractile vacuole "
                "relative to the cell's velocity as a windrose plot")
    parser.add_argument('data_file', help = ("path to a csv (from "
        "process_cv_orientation.py) that includes information on the "
        "position and orientation of different contractile vacuoles."))
    parser.add_argument('--plotfile', help=("pathname of the file to save "
    "of the plot. Defaults to {todaysdate}_cv_windrose.pdf in the same folder "
    "as the data_file"))

    args = parser.parse_args()
    if args.plotfile is None:
        args.plotfile = os.path.join(os.path.dirname(
                                     os.path.abspath(args.data_file)),
                                     "{}_cv_windrose.pdf".format(
                                         date.today().strftime("%Y%m%d")))
    return args


def avg_angle(x):
    '''Compute the average angle from a vector of angles (in radians).
    '''
    return np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x)))


def plot_cv_orientation(df, save_name):
    '''Plot orientation of contractile vacuoles as a windrose plot.
    Input: 
        - df: dataframe of contractile vacuole positions and orientations
        - save_name: path to save the windrose plot
    '''
    #Initialize the matplotlib settings
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['font.size'] = 10
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['svg.fonttype'] = 'none'
    fig, ax = plt.subplots(figsize=(3.5,2.75))
    #use polar projection for windrose plot
    ax = plt.subplot(projection = 'polar')
    N = 20 #number of bins
    theta = np.linspace(-np.pi, np.pi, N, endpoint = True)
    radii = np.histogram(df['cv_angle'].to_numpy(), bins = theta)
    theta_centers = 0.5 * (theta[1:] + theta[:-1])
    ax.bar(theta_centers, radii[0], width=(2*np.pi/N), 
            bottom = 0.0, edgecolor = 'k',facecolor='#999999')
    #Add average angle of each movie
    angle = (df[['trial','cv_angle']]
                .groupby('trial')['cv_angle']
                .apply(lambda x: avg_angle(x)))
    for a in angle:
        ax.plot(a, 29.5, 'o', mfc='none', ms = 5, clip_on = False)
    #Add average of averages, same color as the bars
    ax.plot(avg_angle(angle), 29.5, 's', mfc = '#999999', 
            mec='k',ms=7, clip_on = False)
    #Add an arrow pointing left -- direction of centroid velocity
    ax.annotate('', xy=[1.35, .5], xytext=[1.15, .5], 
            xycoords='axes fraction', 
            arrowprops = dict(arrowstyle="->", color='b', lw = 2))
    #Add text to arrow
    ax.text(.77, .415, "centroid\nvelocity", 
            transform = plt.gcf().transFigure)
    ax.set_ylim([0,20])
    ax.set_xlabel(r'$\theta$')
    ax.yaxis.set_major_locator(mpl.ticker.FixedLocator(locs=[5,10,15,20]))
    plt.tight_layout()
    plt.savefig(save_name)
    plt.show()


def main():
    args = parse_args()
    data = pd.read_csv(args.data_file)
    plot_cv_orientation(data, args.plotfile)


if __name__ == "__main__":
    main()


