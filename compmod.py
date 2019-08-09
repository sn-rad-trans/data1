#! /usr/bin/env python

'''compute magnitudes from spectra
'''

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from pdb import set_trace as stop

### ensure Python2 (2.6 or 2.7) and Python3 compatibility
if sys.version_info.major == 2:
    input = raw_input # input() to mean raw_input() when running Python2

###############################################################################

def read_lbol(file):

    """read in lbol_edep_*.txt file
    """

    print('INFO - reading file '+file)
    
    with open(file, 'r') as f:
        
        ### read in header
        okhdr = 0
        while okhdr == 0:
            line = f.readline()
            if 'NTIMES' in line:
                nt = int(line.split()[1])
                okhdr = 1

    ### read numerical data
    t, lbol, edep = np.loadtxt(file, unpack=True)
                        
    #output
    return t, lbol, edep

###############################################################################

def read_mags(file):

    """read in mags_*.txt file
    """

    print('INFO - reading file '+file)

    with open(file, 'r') as f:
        
        ### read in header
        okhdr = 0
        while okhdr == 0:
            line = f.readline()
            if 'NTIMES' in line:
                nt = int(line.split()[1])
            elif 'NBANDS' in line:
                nband = int(line.split()[1])
            elif 'time' in line:
                split_line = line.split()
                bands = split_line[1:]
                okhdr = 1

    ### read numerical data in one go
    vals = np.loadtxt(file)
    t = vals[:,0]
    mags = {}
    for ib, band in enumerate(bands):
        mags[band] = vals[:,ib+1]
    if 'Ks' in mags.keys():
        mags['K'] = mags['Ks']
        
    #output
    return t, mags

###############################################################################

def read_spectra(file):

    """read in spectra_*.txt file
    """

    print('INFO - reading file '+file)

    with open(file, 'r') as f:
        
        ### read in header
        okhdr = 0
        while okhdr == 0:
            line = f.readline()
            if 'NTIMES' in line:
                nt = int(line.split()[1])
            elif 'NWAVE' in line:
                nw = int(line.split()[1])
            elif 'TIMES' in line:
                split_line = line.split()
                tarrstr = split_line[1:]
                tarr = [float(tt) for tt in tarrstr]
                okhdr = 1

    ### read numerical data in one go
    vals = np.loadtxt(file)
    warr = vals[:,0]
    farr = np.zeros((nt,nw))
    for it in range(nt):
        farr[it,:] = vals[:,it+1]
                        
    #output
    return np.array(tarr), warr, farr

###############################################################################

def read_edep(file):

    """read in edep_*.txt file
    """

    print('INFO - reading file '+file)

    with open(file, 'r') as f:
        
        ### read in header
        okhdr = 0
        while okhdr == 0:
            line = f.readline()
            if 'NTIMES' in line:
                nt = int(line.split()[1])
            elif 'NVEL' in line:
                nvel = int(line.split()[1])
            elif 'TIMES' in line:
                split_line = line.split()
                tarrstr = split_line[1:]
                tarr = [float(tt) for tt in tarrstr]
                okhdr = 1

    ### read numerical data in one go
    vals = np.loadtxt(file)
    velarr = vals[:,0]
    edeparr = np.zeros((nt,nvel))
    for it in range(nt):
        edeparr[it,:] = vals[:,it+1]
                        
    #output
    return np.array(tarr), velarr, edeparr

###############################################################################

def read_tgas(file):

    """read in tgas_*.txt file
    """

    print('INFO - reading file '+file)

    with open(file, 'r') as f:

        ### read in header
        okhdr = 0
        while okhdr == 0:
            line = f.readline()
            if 'NTIMES' in line:
                nt = int(line.split()[1])
            elif 'NVEL' in line:
                nvel = int(line.split()[1])
            elif 'TIMES' in line:
                split_line = line.split()
                tarrstr = split_line[1:]
                tarr = [float(tt) for tt in tarrstr]
                okhdr = 1

    ### read numerical data in one go
    vals = np.loadtxt(file)
    velarr = vals[:,0]
    tgasarr = np.zeros((nt,nvel))
    for it in range(nt):
        tgasarr[it,:] = vals[:,it+1]
                        
    #output
    return np.array(tarr), velarr, tgasarr

###############################################################################

def read_eden(file):

    """read in eden_*.txt file
    """

    print('INFO - reading file '+file)

    with open(file, 'r') as f:

        ### read in header
        okhdr = 0
        while okhdr == 0:
            line = f.readline()
            if 'NTIMES' in line:
                nt = int(line.split()[1])
            elif 'NVEL' in line:
                nvel = int(line.split()[1])
            elif 'TIMES' in line:
                split_line = line.split()
                tarrstr = split_line[1:]
                tarr = [float(tt) for tt in tarrstr]
                okhdr = 1

    ### read numerical data in one go
    vals = np.loadtxt(file)
    velarr = vals[:,0]
    edenarr = np.zeros((nt,nvel))
    for it in range(nt):
        edenarr[it,:] = vals[:,it+1]
                        
    #output
    return np.array(tarr), velarr, edenarr

###############################################################################

if __name__ == '__main__':

    import argparse
    import glob

    #
    # TODO
    #
    # - assign same color to a given model
    # - automatically adjust axis ranges
    # - read_ionfrac() function
    # - options to set axis ranges
    # - report time in legend title when appropriate
    # - help and example runs?
    #
    
    # parse command-line options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument('quantity', default=None, type=str, help="quantity to plot (lbol,edeptot,<band>mag,<band1><band2>col,spec,optspec,nirspec,edep,tgas,eden)")
    parser.add_argument('model', default=None, type=str, help="model name")
    parser.add_argument('--time', default=None, type=float, help="days since explosion")
    parser.add_argument('--nointerp', action='store_true', help="do *not* interpolate at required time")
    parser.add_argument('--maxdt', default=1.0, type=float, help="maximum delta_t")
    args = parser.parse_args()

    # read in files corresponding to quantity
    data_dict = {}
    
    if args.quantity in ['lbol','edeptot']:

        files = glob.glob(args.model + '/lbol_edep_' + args.model + '*.txt')

        for lbolfile in files:
            # extract code name
            codenametmp = lbolfile[lbolfile.rfind('_')+1:] # <code>.txt
            codename = codenametmp[:-4]
            # read in file
            t, lbol, edeptot = read_lbol(lbolfile)
            # update dictionary
            if args.quantity == 'lbol':
                idx = np.where((lbol > 0) & (t > 0))[0]
                data_dict[codename] = {'time':t[idx], 'lbol':np.log10(lbol[idx])}
            elif args.quantity == 'edeptot':
                idx = np.where((edeptot > 0) & (t > 0))[0]
                data_dict[codename] = {'time':t[idx], 'edeptot':np.log10(edeptot[idx])}

    elif args.quantity[-3:] in ['mag','col']: #umag bmag etc. or bvcol

        files = glob.glob(args.model + '/mags_' + args.model + '*.txt')
        
        for magsfile in files:
            # extract code name
            codenametmp = magsfile[magsfile.rfind('_')+1:] # <code>.txt
            codename = codenametmp[:-4]
            # read in file
            t, mags = read_mags(magsfile)
            # update dictionary
            if args.quantity[-3:] == 'mag':
                band = args.quantity[0].upper()
                idx = np.where((mags[band] < 99.) & (t > 0.))[0]
                data_dict[codename] = {'time':t[idx], 'mag':mags[band][idx]}
            if args.quantity[-3:] == 'col':
                band1 = args.quantity[0].upper()
                band2 = args.quantity[1].upper()
                idx = np.where((mags[band1] < 99.) & (mags[band2] < 99.) & (t > 0.))[0]
                data_dict[codename] = {'time':t[idx], 'col':mags[band1][idx]-mags[band2][idx]}

    elif args.quantity[-4:] == 'spec': # spec, optspec, nirspec

        files = glob.glob(args.model + '/spectra_' + args.model + '*.txt')
        
        for specfile in files:
            # extract code name
            codenametmp = specfile[specfile.rfind('_')+1:] # <code>.txt
            codename = codenametmp[:-4]
            # read in file
            tarr, warr, farr = read_spectra(specfile)
            nw = len(warr)
            # find closest time within maxdt
            abs_dt = np.abs(tarr-args.time)
            if args.nointerp:
                idx = np.argmin(abs_dt)
                if abs_dt[idx] < args.maxdt:
                    flux = farr[idx,:]
                    # update dictionary
                    data_dict[codename] = {'wave':warr, 'flux':flux}
            else:
                # interpolate to required time
                idxi = np.where(tarr < args.time)[0]
                idxs = np.where(tarr > args.time)[0]
                if idxi.size > 0 and idxs.size > 0:
                    idxi = idxi.max()
                    idxs = idxs.min()
                    if (tarr[idxi] - args.time < args.maxdt) and (args.time - tarr[idxs] < args.maxdt):
                        fluxi = farr[idxi,:]
                        fluxs = farr[idxs,:]
                        finti = interpolate.interp1d(warr, fluxi, bounds_error=True)
                        fints = interpolate.interp1d(warr, fluxs, bounds_error=True)
                        fluxiinterp = finti(warr)
                        fluxsinterp = fints(warr)
                        tt = np.array([tarr[idxi], tarr[idxs]])
                        flux = np.zeros(nw)
                        for iw in range(nw):
                            ff = np.array([fluxiinterp[iw], fluxsinterp[iw]])
                            g = interpolate.interp1d(tt, ff, bounds_error=True)
                            flux[iw] = g(args.time)
                        # update dictionary
                        data_dict[codename] = {'wave':warr, 'flux':flux}
                        
    elif args.quantity == 'edep':

        files = glob.glob(args.model + '/edep_' + args.model + '*.txt')
        
        for edepfile in files:
            # extract code name
            codenametmp = edepfile[edepfile.rfind('_')+1:] # <code>.txt
            codename = codenametmp[:-4]
            # read in file
            tarr, velarr, edeparr = read_edep(edepfile)
            nvel = len(velarr)
            # find closest time within maxdt
            abs_dt = np.abs(tarr-args.time)
            if args.nointerp:
                idx = np.argmin(abs_dt)
                if abs_dt[idx] < args.maxdt:
                    edep = edeparr[idx,:]
                    # update dictionary
                    idx = np.where(edep > 0)[0]
                    data_dict[codename] = {'vel':velarr[idx], 'edep':np.log10(edep[idx])}
            else:
                # interpolate to required time
                idxi = np.where(tarr < args.time)[0]
                idxs = np.where(tarr > args.time)[0]
                if idxi.size > 0 and idxs.size > 0:
                    idxi = idxi.max()
                    idxs = idxs.min()
                    if (tarr[idxi] - args.time < args.maxdt) and (args.time - tarr[idxs] < args.maxdt):
                        edepi = edeparr[idxi,:]
                        edeps = edeparr[idxs,:]
                        finti = interpolate.interp1d(velarr, edepi, bounds_error=True)
                        fints = interpolate.interp1d(velarr, edeps, bounds_error=True)
                        edepiinterp = finti(velarr)
                        edepsinterp = fints(velarr)
                        tt = np.array([tarr[idxi], tarr[idxs]])
                        edep = np.zeros(nvel)
                        for iw in range(nvel):
                            ff = np.array([edepiinterp[iw], edepsinterp[iw]])
                            g = interpolate.interp1d(tt, ff, bounds_error=True)
                            edep[iw] = g(args.time)
                        # update dictionary
                        idx = np.where(edep > 0)[0]
                        data_dict[codename] = {'vel':velarr[idx], 'edep':np.log10(edep[idx])}
                        
    elif args.quantity == 'tgas':

        files = glob.glob(args.model + '/tgas_' + args.model + '*.txt')
        
        for tgasfile in files:
            # extract code name
            codenametmp = tgasfile[tgasfile.rfind('_')+1:] # <code>.txt
            codename = codenametmp[:-4]
            # read in file
            tarr, velarr, tgasarr = read_tgas(tgasfile)
            nvel = len(velarr)
            # find closest time within maxdt
            abs_dt = np.abs(tarr-args.time)
            if args.nointerp:
                idx = np.argmin(abs_dt)
                if abs_dt[idx] < args.maxdt:
                    tgas = tgasarr[idx,:]
                    # update dictionary
                    idx = np.where(tgas > 0)[0]
                    data_dict[codename] = {'vel':velarr[idx], 'tgas':tgas[idx]}
            else:
                # interpolate to required time
                idxi = np.where(tarr < args.time)[0]
                idxs = np.where(tarr > args.time)[0]
                if idxi.size > 0 and idxs.size > 0:
                    idxi = idxi.max()
                    idxs = idxs.min()
                    if (tarr[idxi] - args.time < args.maxdt) and (args.time - tarr[idxs] < args.maxdt):
                        tgasi = tgasarr[idxi,:]
                        tgass = tgasarr[idxs,:]
                        finti = interpolate.interp1d(velarr, tgasi, bounds_error=True)
                        fints = interpolate.interp1d(velarr, tgass, bounds_error=True)
                        tgasiinterp = finti(velarr)
                        tgassinterp = fints(velarr)
                        tt = np.array([tarr[idxi], tarr[idxs]])
                        tgas = np.zeros(nvel)
                        for iw in range(nvel):
                            ff = np.array([tgasiinterp[iw], tgassinterp[iw]])
                            g = interpolate.interp1d(tt, ff, bounds_error=True)
                            tgas[iw] = g(args.time)
                        # update dictionary
                        idx = np.where(tgas > 0)[0]
                        data_dict[codename] = {'vel':velarr[idx], 'tgas':tgas[idx]}
                        
    elif args.quantity == 'eden':

        files = glob.glob(args.model + '/eden_' + args.model + '*.txt')
        
        for edenfile in files:
            # extract code name
            codenametmp = edenfile[edenfile.rfind('_')+1:] # <code>.txt
            codename = codenametmp[:-4]
            # read in file
            tarr, velarr, edenarr = read_eden(edenfile)
            nvel = len(velarr)
            # find closest time within maxdt
            abs_dt = np.abs(tarr-args.time)
            if args.nointerp:
                idx = np.argmin(abs_dt)
                if abs_dt[idx] < args.maxdt:
                    eden = edenarr[idx,:]
                    # update dictionary
                    idx = np.where(eden > 0)[0]
                    data_dict[codename] = {'vel':velarr[idx], 'eden':np.log10(eden[idx])}
            else:
                # interpolate to required time
                idxi = np.where(tarr < args.time)[0]
                idxs = np.where(tarr > args.time)[0]
                if idxi.size > 0 and idxs.size > 0:
                    idxi = idxi.max()
                    idxs = idxs.min()
                    if (tarr[idxi] - args.time < args.maxdt) and (args.time - tarr[idxs] < args.maxdt):
                        edeni = edenarr[idxi,:]
                        edens = edenarr[idxs,:]
                        finti = interpolate.interp1d(velarr, edeni, bounds_error=True)
                        fints = interpolate.interp1d(velarr, edens, bounds_error=True)
                        edeniinterp = finti(velarr)
                        edensinterp = fints(velarr)
                        tt = np.array([tarr[idxi], tarr[idxs]])
                        eden = np.zeros(nvel)
                        for iw in range(nvel):
                            ff = np.array([edeniinterp[iw], edensinterp[iw]])
                            g = interpolate.interp1d(tt, ff, bounds_error=True)
                            eden[iw] = g(args.time)
                        # update dictionary
                        idx = np.where(eden > 0)[0]
                        data_dict[codename] = {'vel':velarr[idx], 'eden':np.log10(eden[idx])}
                        
    else:
        sys.exit("ERROR - unknown quantity: " + args.quantity)

    # plot
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if args.quantity == 'lbol':
        xx = 'time'
        yy = 'lbol'
        xtit = r'Days since Explosion'
        ytit = r'$\log_{10} L_{\mathrm{bol,uvoir}}$ [erg s$^{-1}$]'
    elif args.quantity == 'edeptot':
        xx = 'time'
        yy = 'edeptot'
        xtit = r'Days since Explosion'
        ytit = r'$\log_{10} E_{\mathrm{dep,tot}}$ [erg s$^{-1}$]'
    elif args.quantity[-3:] == 'mag':
        xx = 'time'
        yy = 'mag'
        xtit = r'Days since Explosion'
        ytit = r'$'+band+'$-band Magnitude'
        plt.gca().invert_yaxis()
    elif args.quantity[-3:] == 'col':
        xx = 'time'
        yy = 'col'
        xtit = r'Days since Explosion'
        ytit = r'$'+band1+'-'+band2+'$ Colour Index'
    elif args.quantity[-4:] == 'spec':
        xx = 'wave'
        yy = 'flux'
        xtit = r'Wavelength [Ang]'
        ytit = r'Flux [erg/s/Ang]'
        if args.quantity[:3] == 'opt':
            ax.set_xlim((2500., 10500.))
        elif args.quantity[:3] == 'nir':
            ax.set_xlim((10500., 25000.))
    elif args.quantity == 'edep':
        xx = 'vel'
        yy = 'edep'
        xtit = r'Velocity [km s$^{-1}$]'
        ytit = r'$\log_{10} E_{\mathrm{dep}}$ [erg s$^{-1}$ cm$^{-3}$]'
    elif args.quantity == 'tgas':
        xx = 'vel'
        yy = 'tgas'
        xtit = r'Velocity [km s$^{-1}$]'
        ytit = r'$T_{\mathrm{gas}}$ [K]'
    elif args.quantity == 'eden':
        xx = 'vel'
        yy = 'eden'
        xtit = r'Velocity [km/s]'
        ytit = r'$n_{\mathrm{e}}$ [cm$^{-3}$]'
        
    ax.set_xlabel(xtit)
    ax.set_ylabel(ytit)

    for code in data_dict.keys():
        ax.plot(data_dict[code][xx], data_dict[code][yy], label=code)
        
    ax.legend(fontsize='small',title='model '+args.model)

    plt.draw()
    zzz = input("===> Hit <return> to quit ")
    ax.clear()

