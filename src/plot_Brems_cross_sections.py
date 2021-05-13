import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse
import os
import sys
import glob
import seaborn as sns

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['legend.fontsize'] = 9.0
mpl.rcParams['figure.figsize'] = 7.2, 5.0
sns.set_palette("rocket_r",6)


Labels = {
    'B1': r"$\pi^{\pm} + \ \pi^{\mp} \ \rightarrow \ \pi^{\pm} \ + \ \pi^{\mp}  \ + \ \gamma$",
    'B2': r"$\pi^{\pm} + \ \pi^{\pm} \ \rightarrow \ \pi^{\pm} \ + \ \pi^{\pm}  \ + \ \gamma$",
    'B3': r"$\pi^{\pm} + \ \pi^{0} \ \rightarrow \ \pi^{\pm} \ + \ \pi^{0}  \ + \ \gamma$",
    'B4': r"$\pi^{\pm} + \ \pi^{\mp} \ \rightarrow \ \pi^{0} \ + \ \pi^{0}  \ + \ \gamma$",
    'B5': r"$\pi^{0} + \ \pi^{0} \ \rightarrow \ \pi^{\pm} \ + \ \pi^{\mp}  \ + \ \gamma$"
    }

Linestyles = {
    'B1': '-.',
    'B2': '--',
    'B3': (0, (1, 2, 1, 2, 5, 2)),
    'B4': (0, (4.5, 1, 1, 1, 4.5, 1)),
    'B5': ':',
    'Total': '-',
    }

def read_cross_section_files(directory):
    ''' Read cross sections files from given directory and return lists of
        the entries.
    '''

    Xsec_files = [os.path.basename(x) for x in glob.glob(directory + 'Brems*.txt')]

    for file in Xsec_files:
        filename = file.split('_')
        if 'Total' in filename:
            Total_Xsections_PiRho_PiGamma = np.loadtxt(directory + file, unpack = True)
        elif 'dSigmadk' in filename:
            dSigmadk = np.loadtxt(directory + file, unpack = True)
        elif 'dSigmadTheta' in filename:
            dSigmadTheta = np.loadtxt(directory + file, unpack = True)
        else:
            print 'WARNING: Cannot match filename of file ' + str(file) + '. Skipping file ...'


    # return Total_Xsections_PiRho_PiGamma
    return Total_Xsections_PiRho_PiGamma, dSigmadk, dSigmadTheta

def plot_diff_cross_section_k(xsec_data, output_dir):
    ''' Plot the differential cross segion dSigma/dk '''

    k = xsec_data[0]

    plt.plot(k, xsec_data[1], ls = Linestyles['B1'], label = Labels['B1'])
    plt.plot(k, xsec_data[2], ls = Linestyles['B2'], label = Labels['B2'])
    plt.plot(k, xsec_data[3], ls = Linestyles['B3'], label = Labels['B3'])
    plt.plot(k, xsec_data[4], ls = Linestyles['B4'], label = Labels['B4'])
    plt.plot(k, xsec_data[5], ls = Linestyles['B5'], label = Labels['B5'])

    plt.xlabel(r"k [GeV]")
    plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}k}$ [mb/GeV]")
    # xticks = [0,45,90,135,180]
    # xticks_labels = ["0","45","90","135","180"]
    # plt.xticks(xticks,xticks_labels)

    plt.yscale('log')
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir + "/Brems_dSigma_dk.pdf")
    plt.close()


def plot_diff_cross_section_theta(xsec_data, output_dir):
    ''' Plot the differential cross segion dSigma/dk '''

    theta = xsec_data[0]

    plt.plot(theta, xsec_data[1], ls = Linestyles['B1'], label = Labels['B1'])
    plt.plot(theta, xsec_data[2], ls = Linestyles['B2'], label = Labels['B2'])
    plt.plot(theta, xsec_data[3], ls = Linestyles['B3'], label = Labels['B3'])
    plt.plot(theta, xsec_data[4], ls = Linestyles['B4'], label = Labels['B4'])
    plt.plot(theta, xsec_data[5], ls = Linestyles['B5'], label = Labels['B5'])

    plt.xlabel(r"$\vartheta$ [GeV]")
    plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\vartheta}$ [mb]")
    plt.ylim(0,4.0)

    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir + "/Brems_dSigma_dtheta.pdf")
    plt.close()


def plot_total_cross_section(xsec_data, output_dir):
    ''' Plot the provided cross section data as a function of sqrt(s) '''

    sqrts = xsec_data[0]

    plt.plot(sqrts, xsec_data[1], ls = Linestyles['B1'], label = Labels['B1'])
    plt.plot(sqrts, xsec_data[2], ls = Linestyles['B2'], label = Labels['B2'])
    plt.plot(sqrts, xsec_data[3], ls = Linestyles['B3'], label = Labels['B3'])
    plt.plot(sqrts, xsec_data[4], ls = Linestyles['B4'], label = Labels['B4'])
    plt.plot(sqrts, xsec_data[5], ls = Linestyles['B5'], label = Labels['B5'])
    plt.plot(sqrts, (xsec_data[1] + xsec_data[2] + xsec_data[3] + xsec_data[4] + xsec_data[5]), ls = Linestyles['Total'], label = 'Total')

    plt.ylim(3e-3,1e2)
    plt.xlim(0.2,5)
    plt.xlabel(r"$\sqrt{s}$ [GeV]")
    plt.ylabel(r"$\sigma$ [mb]")
    plt.yscale('log')
    plt.legend(ncol = 2, loc = 'lower right', frameon = False)

    plt.tight_layout()
    plt.savefig(output_dir + "/Total_Xsec_Brems.pdf")
    plt.close()



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--directory", nargs = '+', required = True,
                        help = "Path to the directory that contains the \
                                cross section files")
    parser.add_argument("-o", "--output_dir", nargs = '+', required = False,
                        default = '.', help = 'Path to the output directory \
                        to store the results')
    args = parser.parse_args()

    cross_section_directory = args.directory[0] + '/'
    output_directory = args.output_dir[0]

    Total, dSigma_dk, dSigma_dTheta = read_cross_section_files(cross_section_directory)

    # plot dSigma/dk
    plot_diff_cross_section_k(dSigma_dk, output_directory)
    plot_diff_cross_section_theta(dSigma_dTheta, output_directory)
    # plot_diff_cross_section(Diff_PiRho, rho_mass, False, output_directory)
    #
    # # plot differential Xsection as a function of theta
    # plot_diff_cross_section(Diff_PiPi, rho_mass, True, output_directory)
    # plot_diff_cross_section(Diff_PiRho, rho_mass, True, output_directory)
    #
    # # plot total cross section as a function of sqrt(s)
    plot_total_cross_section(Total, output_directory)
    # plot_total_cross_section(Total_PiRho, rho_mass, output_directory)
