#! /usr/bin/env python

'''read SN code-comparison output files

To produce example plots based on files in https://github.com/sn-rad-trans/data1 simply run:

python read_outputs.py

in the directory containing the model directories (or set --path2data option)
'''

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace as stop

### ensure Python2 (2.6 or 2.7) and Python3 compatibility
if sys.version_info.major == 2:
    input = raw_input # input() to mean raw_input() when running Python2

###############################################################################

def read_lbol_edep(file):

    '''read in lbol_edep_<model>_<code>.txt file

    #NTIMES: 100
    #time[d] Lbol[erg/s] Edep[erg/s]
    '''

    print('INFO - reading file ' + file)

    with open(file, 'r') as f:
        
        ### read in header
        okhdr = 0
        while okhdr == 0:
            line = f.readline()
            if 'NTIMES' in line:
                nt = int(line.split()[1])
                okhdr = 1

    ### read numerical data
    t, lbol, edep = np.loadtxt(file, comments='#', unpack=True)
                        
    #output
    out = {}
    out['time'] = t
    out['lbol'] = lbol
    out['edep'] = edep
    out['units'] = 'time: days\nlbol: erg/s\nedep: erg/s'

    return out
    
###############################################################################

def read_spectra(file):

    '''read in spectra_<model>_<code>.txt file

    #NTIMES: 100
    #NWAVE: 2000 
    #TIMES[d]: 2.0 3.0 4.0 ... 100.0
    #wavelength[Ang] flux_t0[erg/s/Ang] flux_t1[erg/s/Ang] … flux_tn[erg/s/Ang]
    '''

    print('INFO - reading file ' + file)
    
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

    ### read numerical data
    vals = np.loadtxt(file)
    warr = vals[:,0]
    farr = np.zeros((nt,nw))
    for it in range(nt):
        farr[it,:] = vals[:,it+1]
                        
    #output
    out = {}
    out['time'] = np.array(tarr)
    out['wave'] = warr
    out['flux'] = farr
    out['units'] = 'time: days\nwave: angstroms (A)\nflux: erg/s/A'

    return out
    
###############################################################################

def read_edep(file):

    '''read in edep_<model>_<code>.txt file

    #NTIMES: 100
    #NVEL: 100
    #TIMES[d]: 2.0 3.0 4.0 … 100.0
    #vel_mid[km/s] Edep_t0[erg/s/cm^3] Edep_t1[erg/s/cm^3] … Edep_tn[erg/s/cm^3]
    '''

    print('INFO - reading file ' + file)
    
    with open(file, 'r') as f:
        
        ### read in header
        okhdr = 0
        while okhdr == 0:
            line = f.readline()
            if 'NTIMES' in line:
                nt = int(line.split()[1])
            elif 'NVEL' in line:
                nv = int(line.split()[1])
            elif 'TIMES' in line:
                split_line = line.split()
                tarrstr = split_line[1:]
                tarr = [float(tt) for tt in tarrstr]
                okhdr = 1

    ### read numerical data
    vals = np.loadtxt(file)
    velarr = vals[:,0]
    edeparr = np.zeros((nt,nv))
    for it in range(nt):
        edeparr[it,:] = vals[:,it+1]
                        
    #output
    out = {}
    out['time'] = np.array(tarr)
    out['vel'] = velarr
    out['edep'] = edeparr
    out['units'] = 'time: days\nvel: km/s\nedep: erg/s/cm^3'

    return out
    
###############################################################################

def read_phys(file):

    '''read in phys_<model>_<code>.txt file

    #NTIMES: 100
    #TIMES[d]: 2.0 3.0 4.0 … 100.0
    #
    #TIME: 2.0 
    #NVEL: 100
    #vel_mid[km/s] temp[K] rho[gcc] ne[/cm^3] natom[/cm^3] 
    .
    .
    .
    #TIME: 3.0 
    #NVEL: 100
    #vel_mid[km/s] temp[K] rho[gcc] ne[/cm^3] natom[/cm^3]
    .
    .
    .
    '''

    print('INFO - reading file ' + file)
    
    with open(file, 'r') as f:
        
        ### read in header
        okhdr = 0
        while okhdr == 0:
            line = f.readline()
            if 'NTIMES' in line:
                nt = int(line.split()[1])
            elif 'TIMES' in line:
                split_line = line.split()
                tarrstr = split_line[1:]
                tarr = [float(tt) for tt in tarrstr]
                okhdr = 1

        ### output lists of arrays
        vel, temp, dens, ne, natom = [], [], [], [], []
                
        ### loop over times
        for it in range(nt):
            
            while 1:
                line = f.readline()
                if 'NVEL' in line:
                    split_line = line.split()
                    nvel = int(split_line[1])
                elif 'vel_mid' in line:
                    break

            # temporary arrays
            veltmp = np.zeros(nvel)
            temptmp = np.zeros(nvel)
            denstmp = np.zeros(nvel)
            netmp = np.zeros(nvel)
            natomtmp = np.zeros(nvel)
                
            for iv in range(nvel):
                line = f.readline()
                split_line = line.split()
                veltmp[iv] = float(split_line[0])
                temptmp[iv] = float(split_line[1])
                denstmp[iv] = float(split_line[2])
                netmp[iv] = float(split_line[3])
                natomtmp[iv] = float(split_line[4])

            # append to lists
            vel.append(veltmp)
            temp.append(temptmp)
            dens.append(denstmp)
            ne.append(netmp)
            natom.append(natomtmp)
                
    #output
    out = {}
    out['time'] = np.array(tarr)
    out['vel'] = vel
    out['temp'] = temp
    out['dens'] = dens
    out['ne'] = ne
    out['natom'] = natom
    out['units'] = 'time: days\nvel: km/s\ntemp: K\ndens: g/cm^3\nne: \cm^3\nnatom: /cm^3'

    return out
    
###############################################################################

def read_ionfrac(file):

    '''read in ionfrac_<element>_<model>_<code>.txt file

    #NTIMES: 100
    #NSTAGES: 7
    #TIMES[d]: 2.0 3.0 4.0 … 100.0
    #
    #TIME: 2.0 
    #NVEL: 100
    #vel_mid[km/s] fe0 fe1 ... fe6
    .
    .
    .
    #TIME: 3.0 
    #NVEL: 100
    #vel_mid[km/s] fe0 fe1 ... fe6
    .
    .
    .
    '''

    print('INFO - reading file ' + file)
    
    with open(file, 'r') as f:
        
        ### read in header
        okhdr = 0
        while okhdr == 0:
            line = f.readline()
            if 'NTIMES' in line:
                nt = int(line.split()[1])
            elif 'NSTAGES' in line:
                nstages = int(line.split()[1])
            elif 'TIMES' in line:
                split_line = line.split()
                tarrstr = split_line[1:]
                tarr = [float(tt) for tt in tarrstr]
                okhdr = 1

        ### output lists of arrays
        vel, ionfrac = [], []
                
        ### loop over times
        for it in range(nt):
            
            while 1:
                line = f.readline()
                if 'NVEL' in line:
                    split_line = line.split()
                    nvel = int(split_line[1])
                elif 'vel_mid' in line:
                    break

            # temporary arrays
            veltmp = np.zeros(nvel)
            ionfractmp = np.zeros((nvel, nstages))
                
            for iv in range(nvel):
                line = f.readline()
                split_line = line.split()
                veltmp[iv] = float(split_line[0])
                ionfractmp[iv,:] = [float(x) for x in split_line[1:]]
                
            # append to lists
            vel.append(veltmp)
            ionfrac.append(ionfractmp)
                
    #output
    out = {}
    out['nstages'] = nstages
    out['time'] = np.array(tarr)
    out['vel'] = vel
    out['ionfrac'] = ionfrac
    out['units'] = 'nstages: N/A\ntime: days\nvel: km/s\nionfrac: N/A'

    return out
    
###############################################################################

if __name__ == '__main__':

    import argparse
    
    # parse command-line options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter) 
    parser.add_argument('--path2data', default='.', type=str, help='path to directory containing the model files')
    parser.add_argument('--model', default='toy06', type=str, help='model name')
    parser.add_argument('--code', default='cmfgen', type=str, help='code name')
    parser.add_argument('--element', default='fe', type=str, help='element name')
    args = parser.parse_args()

    path2data = args.path2data
    model = args.model
    code = args.code
    elem = args.element
    
    lbol = read_lbol_edep(path2data+'/'+model+'/lbol_edep_'+model+'_'+code+'.txt')
    spec = read_spectra(path2data+'/'+model+'/spectra_'+model+'_'+code+'.txt')
    edep = read_edep(path2data+'/'+model+'/edep_'+model+'_'+code+'.txt')
    phys = read_phys(path2data+'/'+model+'/phys_'+model+'_'+code+'.txt')
    ionf = read_ionfrac(path2data+'/'+model+'/ionfrac_'+elem+'_'+model+'_'+code+'.txt')
    plt.ion()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Lbol and Edep
    ax.plot(lbol['time'], np.log10(lbol['lbol']), label='lbol')
    ax.plot(lbol['time'], np.log10(lbol['edep']), ls='--', label='edep')
    ax.set_xlabel('Time [days]')
    ax.set_ylabel('Log10(Lbol) [erg/s]')
    ax.set_title('Model: '+model+' ; Code: '+code)
    ax.legend()
    plt.draw()
    zzz = input("===> Hit <return> to quit ")
    ax.clear()

    # spectra
    idxt = np.argmin(abs(lbol['time'] - 20.)) # pick time closest to 20d
    # for i in range(len(spec['time'])):
    ax.plot(spec['wave'], spec['flux'][idxt,:], label=str(spec['time'][idxt])+'d')
    ax.set_xlabel('Wavelength [A]')
    ax.set_ylabel('Flux [erg/s/A]')
    ax.set_title('Model: '+model+' ; Code: '+code)
    ax.legend()
    plt.draw()
    zzz = input("===> Hit <return> to quit ")
    ax.clear()

    # Edep
    idxt = np.argmin(abs(lbol['time'] - 20.)) # pick time closest to 20d
    # for i in range(len(edep['time'])):
    idx = np.where(edep['edep'][idxt,:] > 1e-99)[0]
    ax.plot(edep['vel'][idx], edep['edep'][idxt][idx], label=str(edep['time'][idxt])+'d')
    ax.set_xlabel('Velocity [km/s]')
    ax.set_ylabel('Edep [erg/s/cm^3]')
    ax.set_title('Model: '+model+' ; Code: '+code)
    ax.legend()
    plt.draw()
    zzz = input("===> Hit <return> to quit ")
    ax.clear()
    
    # Temperature
    idxt = np.argmin(abs(lbol['time'] - 20.)) # pick time closest to 20d
    # for i in range(len(phys['time'])):
    idx = np.where(phys['temp'][idxt] > 1e-99)[0]
    ax.plot(phys['vel'][idxt][idx], phys['temp'][idxt][idx], label=str(phys['time'][idxt])+'d')
    ax.set_xlabel('Velocity [km/s]')
    ax.set_ylabel('Temperature [K]')
    ax.set_title('Model: '+model+' ; Code: '+code)
    ax.legend()
    plt.draw()
    zzz = input("===> Hit <return> to quit ")
    ax.clear()
    
    # Ionization fractions
    idxt = np.argmin(abs(lbol['time'] - 20.)) # pick time closest to 20d
    for i in range(ionf['nstages']):
        ax.plot(ionf['vel'][idxt], ionf['ionfrac'][idxt][:,i], label=str(i))
    ax.set_xlabel('Velocity [km/s]')
    ax.set_ylabel('Mass Fraction')
    ax.set_title('Model: '+model+' ; Code: '+code)
    ax.legend(title=elem)
    plt.draw()
    zzz = input("===> Hit <return> to quit ")
    ax.clear()
