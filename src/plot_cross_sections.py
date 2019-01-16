import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse
import os
import sys
import glob

mpl.rcParams['font.size'] = 16
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.labelsize'] = 15
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10


Labels = {
    'C11': r"$\pi + \ \rho^0 \ \rightarrow \ (\pi, \rho, \gamma) \ \rightarrow \ \pi \ + \ \gamma$",
    'C12': r"$\pi^0 + \ \rho \ \rightarrow \ (\pi, \rho, \gamma) \ \rightarrow \ \pi \ + \ \gamma$",
    'C13': r"$\pi \ + \ \rho \ \rightarrow \ (\pi, \rho, \gamma) \ \rightarrow \ \pi^0 + \ \gamma$",
    'C14': r"$\pi^0 + \ \rho^0 \ \rightarrow \ \omega \ \rightarrow \ \pi^0 \ + \ \gamma$",
    'C15': r"$\pi \ + \ \rho \rightarrow \ \omega \ \rightarrow \ \pi^0 + \ \gamma$",
    'C16': r"$\pi^0 + \ \rho  \rightarrow \ \omega \ \rightarrow \ \pi \ + \ \gamma$",
    'C21': r"$\pi \ + \ \pi \ \rightarrow \ (\pi, \rho, \gamma) \ \rightarrow \ \rho^0 + \ \gamma$",
    'C22': r"$\pi \ + \ \pi^0 \ \rightarrow \ (\pi, \rho, \gamma) \ \rightarrow \ \rho \ + \ \gamma$"
    }

Colours = {
    'C11': 'C0',
    'C12': 'C1',
    'C13': 'C2',
    'C14': 'C4',
    'C15': 'C2',
    'C16': 'C1',
    'C21': 'C3',
    'C22': 'C7'
    }

Linestyles = {
    'C11': '-',
    'C12': '--',
    'C13': '-.',
    'C14': (0, (4.5, 1, 1, 1, 4.5, 1)),
    'C15': '-.',
    'C16': '--',
    'C21': (0, (1, 2, 1, 2, 5, 2)),
    'C22': ':'
    }

def read_cross_section_files(directory):
    ''' Read cross sections files from given directory and return lists of
        the entries.
    '''

    Xsec_files = [os.path.basename(x) for x in glob.glob(directory + '*Xsec*MeV.txt')]

    for file in Xsec_files:
        filename = file.split('_')
        if 'Diff' in filename and 'PiRho' in filename:
            Diff_Xsections_PiRho_PiGamma = np.loadtxt(directory + file, unpack = True, skiprows = 2)
            mrho = float(filename[-2])/1000        #in GeV
        elif 'Diff' in filename and 'PiPi' in filename:
            Diff_Xsections_PiPi_RhoGamma = np.loadtxt(directory + file, unpack = True, skiprows = 2)
            mrho1 = float(filename[-2])/1000
        elif 'Total' in filename and 'PiRho' in filename:
            Total_Xsections_PiRho_PiGamma = np.loadtxt(directory + file, unpack = True, skiprows = 1)
            mrho2 = float(filename[-2])/1000
        elif 'Total' in filename and 'PiPi' in filename:
            Total_Xsections_PiPi_RhoGamma = np.loadtxt(directory + file, unpack = True, skiprows = 1)
            mrho3 = float(filename[-2])/1000
        else:
            print 'WARNING: Cannot match filename of file ' + str(file) + '. Skipping file ...'

    assert(mrho == mrho1)
    assert(mrho == mrho2)
    assert(mrho == mrho3)

    return Diff_Xsections_PiRho_PiGamma, Diff_Xsections_PiPi_RhoGamma, \
           Total_Xsections_PiRho_PiGamma, Total_Xsections_PiPi_RhoGamma, mrho

def plot_diff_cross_section(xsec_data, mrho, as_theta, output_dir):
    ''' Plot the provided cross section data as a function of mandelstam t
        or the scattering angle theta '''

    to_mb = 0.3894          #conversion factor fm^2 -> mbarn
    t = xsec_data[0]
    theta = xsec_data[1]

    if len(xsec_data) == 8:
        is_PiRho = True
        is_PiPi = False
    elif len(xsec_data) == 4:
        is_PiPi = True
        is_PiRho = False
    else:
        print 'Warning: Number of columns does not match expectation.'

    assert (is_PiPi != is_PiRho)

    if is_PiPi:
        C21 = xsec_data[2]
        C22 = xsec_data[3]

        xticks = [t[0],-0.3,-0.2,-0.1,t[-1]]
        xticks_labels = [r"$t_{\mathrm{min}}$","-0.3","-0.2","-0.1",r"$t_{\mathrm{max}}$"]

        #plt.figure(figsize=(5.3,5.3))
        if as_theta:        # dsigma/dtheta
            plt.plot(theta, C21*to_mb, label = Labels['C21'], color = Colours['C21'], ls = Linestyles['C21'])
            plt.plot(theta, C22*to_mb, label = Labels['C22'], color = Colours['C22'], ls = Linestyles['C22'])
            plt.xlabel(r"$\theta$")
            plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\theta}$ [mb]")
            plt.xlim(0,180)
            xticks = [0,45,90,135,180]
            xticks_labels = ["0","45","90","135","180"]
            plt.xticks(xticks,xticks_labels)
        else:
            plt.plot(t, C21*to_mb, label = Labels['C21'], color = Colours['C21'], ls = Linestyles['C21'])
            plt.plot(t, C22*to_mb, label = Labels['C22'], color = Colours['C22'], ls = Linestyles['C22'])
            plt.axvline(t[0], ls = "-", color = 'grey')
            plt.axvline(t[-1], ls = "-", color = 'grey')
            plt.xlabel(r"$t$ [GeV$^2$]")
            plt.ylabel(r"$\sigma$ [mb]")
            xticks = [t[0],-0.3,-0.2,-0.1,t[-1]]
            xticks_labels = [r"$t_{\mathrm{min}}$","-0.3","-0.2","-0.1",r"$t_{\mathrm{max}}$"]
            plt.xticks(xticks,xticks_labels)
        plt.ylim(0,4.3)
        plt.legend(loc="best")
        plt.tight_layout()
        if as_theta:
            plt.savefig(output_dir + "/Diff_Xsec_PiPi_RhoGamma_theta.pdf")
        else:
            plt.savefig(output_dir + "/Diff_Xsec_PiPi_RhoGamma.pdf")
        plt.close()


    if is_PiRho:
        C11 = xsec_data[2]
        C12 = xsec_data[3]
        C13 = xsec_data[4]
        C14 = xsec_data[5]
        C15 = xsec_data[6]
        C16 = xsec_data[7]


        xticks = [t[0],-0.3,-0.2,-0.1,t[-1]]
        xticks_labels = [r"$t_{\mathrm{min}}$","-0.3","-0.2","-0.1",r"$t_{\mathrm{max}}$"]

        plt.figure(figsize=(6.5,9))
        plt.subplot(211)
        if as_theta:
            plt.plot(theta, C11*to_mb, label = Labels['C11'], color = Colours['C11'], ls = Linestyles['C11'])
            plt.plot(theta, C12*to_mb, label = Labels['C12'], color = Colours['C12'], ls = Linestyles['C12'])
            plt.plot(theta, C13*to_mb, label = Labels['C13'], color = Colours['C13'], ls = Linestyles['C13'])
            plt.xlabel(r"$\theta$")
            plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\theta}$ [mb]")
            plt.xlim(0,180)
            plt.ylim(0,1.4)
        else:
            plt.plot(t, C11*to_mb, label = Labels['C11'], color = Colours['C11'], ls = Linestyles['C11'])
            plt.plot(t, C12*to_mb, label = Labels['C12'], color = Colours['C12'], ls = Linestyles['C12'])
            plt.plot(t, C13*to_mb, label = Labels['C13'], color = Colours['C13'], ls = Linestyles['C13'])
            plt.axvline(t[0], ls = "-", color = 'grey')
            plt.axvline(t[-1], ls = "-", color = 'grey')
            plt.ylim(0,1.4)
            plt.xticks(xticks,xticks_labels)
            plt.xlabel(r"$t$ [GeV$^2$]")
            plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}t}$ [mb/GeV$^2$]")
        plt.legend(loc="best")


        plt.subplot(212)
        if as_theta:
            plt.plot(theta, C14*to_mb, label = Labels['C14'], color = Colours['C14'], ls = Linestyles['C14'])
            plt.plot(theta, C15*to_mb, label = Labels['C15'], color = Colours['C15'], ls = Linestyles['C15'])
            plt.plot(theta, C16*to_mb, label = Labels['C16'], color = Colours['C16'], ls = Linestyles['C16'])
            plt.xlabel(r"$\theta$")
            plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}\theta}$ [mb]")
            plt.xlim(0,180)
        else:
            plt.plot(t, C14*to_mb, label = Labels['C14'], color = Colours['C14'], ls = Linestyles['C14'])
            plt.plot(t, C15*to_mb, label = Labels['C15'], color = Colours['C15'], ls = Linestyles['C15'])
            plt.plot(t, C16*to_mb, label = Labels['C16'], color = Colours['C16'], ls = Linestyles['C16'])
            plt.axvline(t[0], ls = "-", color = 'grey')
            plt.axvline(t[-1], ls = "-", color = 'grey')
            plt.ylim(0,0.25)
            plt.xticks(xticks,xticks_labels)
            plt.xlabel(r"$t$ [GeV$^2$]")
            plt.ylabel(r"$\frac{\mathrm{d}\sigma}{\mathrm{d}t}$ [mb/GeV$^2$]")
        plt.legend(loc="best")

        plt.tight_layout()

        if as_theta: plt.savefig(output_dir + "/Diff_Xsec_PiRho_PiGamma_theta.pdf")
        else: plt.savefig(output_dir + "/Diff_Xsec_PiRho_PiGamma.pdf")
        plt.close()

def plot_total_cross_section(xsec_data, mrho, output_dir):
    ''' Plot the provided cross section data as a function of sqrt(s) '''
    to_mb = 0.3894
    sqrts = xsec_data[0]
    mpion = 0.138

    if len(xsec_data) == 7:
        is_PiRho = True
        is_PiPi = False
    elif len(xsec_data) == 3:
        is_PiPi = True
        is_PiRho = False
    else:
        print 'Warning: Number of columns does not match expectation.'

    assert (is_PiPi != is_PiRho)

    if is_PiPi:
        C21 = xsec_data[1]
        C22 = xsec_data[2]

        plt.plot(sqrts, C21*to_mb, label = Labels['C21'], color = Colours['C21'], ls = Linestyles['C21'])
        plt.plot(sqrts, C22*to_mb, label = Labels['C22'], color = Colours['C22'], ls = Linestyles['C22'])
        plt.axvline(x=mrho, color = 'grey', label = r"$\sqrt{s} = m_{\pi} + m_{\rho}$")
        plt.ylim(0,1)
        plt.xlim(0.6,3)
        plt.xlabel(r"$\sqrt{s}$ [GeV]")
        plt.ylabel(r"$\sigma$ [mb]")
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_dir + "/Total_Xsec_PiPi_RhoGamma.pdf")
        plt.close()

    if is_PiRho:
        C11 = xsec_data[1]
        C12 = xsec_data[2]
        C13 = xsec_data[3]
        C14 = xsec_data[4]
        C15 = xsec_data[5]
        C16 = xsec_data[6]

        plt.figure(figsize=(6.5,9))

        plt.subplot(211)
        plt.plot(sqrts, C11*to_mb, label = Labels['C11'], color = Colours['C11'], ls = Linestyles['C11'])
        plt.plot(sqrts, C12*to_mb, label = Labels['C12'], color = Colours['C12'], ls = Linestyles['C12'])
        plt.plot(sqrts, C13*to_mb, label = Labels['C13'], color = Colours['C13'], ls = Linestyles['C13'])
        plt.axvline(x=mrho + mpion, color = 'grey', label = r"$\sqrt{s} = m_{\pi} + m_{\rho}$")
        plt.ylim(0,0.5)
        plt.xlim(mrho -0.05,3)
        plt.xlabel(r"$\sqrt{s}$ [GeV]")
        plt.ylabel(r"$\sigma$ [mb]")
        plt.legend()

        plt.subplot(212)
        plt.plot(sqrts, C14*to_mb, label = Labels['C14'], color = Colours['C14'], ls = Linestyles['C14'])
        plt.plot(sqrts, C15*to_mb, label = Labels['C15'], color = Colours['C15'], ls = Linestyles['C15'])
        plt.plot(sqrts, C16*to_mb, label = Labels['C16'], color = Colours['C16'], ls = Linestyles['C16'])
        plt.axvline(x=mrho + mpion, color = 'grey', label = r"$\sqrt{s} = m_{\pi} + m_{\rho}$")
        plt.ylim(0,0.3)
        plt.xlim(mrho -0.05,3)
        plt.xlabel(r"$\sqrt{s}$ [GeV]")
        plt.ylabel(r"$\sigma$ [mb]")
        plt.legend()

        plt.tight_layout()
        plt.savefig(output_dir + "/Total_Xsec_PiRho_PiGamma.pdf")
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

    Diff_PiRho, Diff_PiPi, Total_PiRho, Total_PiPi, rho_mass = read_cross_section_files(cross_section_directory)

    # plot differential Xsection as a function of t
    plot_diff_cross_section(Diff_PiPi, rho_mass, False, output_directory)
    plot_diff_cross_section(Diff_PiRho, rho_mass, False, output_directory)

    # plot differential Xsection as a function of theta
    plot_diff_cross_section(Diff_PiPi, rho_mass, True, output_directory)
    plot_diff_cross_section(Diff_PiRho, rho_mass, True, output_directory)

    # plot total cross section as a function of sqrt(s)
    plot_total_cross_section(Total_PiPi, rho_mass, output_directory)
    plot_total_cross_section(Total_PiRho, rho_mass, output_directory)
