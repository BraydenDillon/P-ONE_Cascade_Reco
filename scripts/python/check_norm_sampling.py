import sys
import os

sys.path.insert(0, os.path.abspath('..'))
import photospline

from icecube import dataio, icetray, dataclasses, simclasses
import numpy as np
import matplotlib.pyplot as plt
import numpy.random as random
from scripts.SplineEval import evalLogPdf, evalPdf
from scipy.stats import norm


gcd = dataio.I3File('../gcdfile/PONE_800mGrid.i3.gz')
splinefit_3d = photospline.SplineTable('../fits/splinelog_3D.fits')
for frame in gcd:
    gcdframe = frame
    break
c = 299792458
n = 1.34
def datacollect(frame):
    dt = []
    t = []
    dphilst = []
    dr = []
    xyz = []
    

    Epos = frame['I3MCTree'][1].pos
    doms = gcdframe['I3ModuleGeoMap']
    omkeys = frame['new_photons'].keys()
    photons = frame['new_photons']
    for key in omkeys:
        # if key.string not in string_subset:
        #     continue
        modulekey = dataclasses.ModuleKey(key.string, key.om)
        dompos = doms[modulekey].pos
        for photon in photons[key]:
            photon_pos = dompos + photon.pos
            xyz.append([photon_pos.x, photon_pos.y, photon_pos.z])
            flight = dompos + photon.pos - Epos
            dr.append(flight.magnitude)
            offset = flight.magnitude * n / c
            dt.append(photon.time - offset*10**9)
            t.append(photon.time)
            # phi = photon.dir
            # randx, randy, randz = random.uniform(-1, 1, 3)
            # x = phi.x
            # y = phi.y
            # z = phi.z
            # dx = x - randx
            # dy = y - randy
            # dz = z - randz
            # dphi = dataclasses.I3Direction(dx, dy, dz)
            # dphilst.append(dphi.zenith)
            # Etheta.append(flight.azimuth)
            # Ephi.append(flight.zenith)
            
    return np.column_stack([xyz, t, dt])

t = np.linspace(splinefit_3d.extents[2][0], splinefit_3d.extents[2][1], 1000)

def compare_spline(fixed_coords: np.array, spline, frames)->None:
    # fixed_coords has shape 2, and should be [dr, dphi]
    # should be able to take either single frame or list of frames
    dR, dphi = fixed_coords
    dr_range = [dR - 1, dR + 1]
    dphi_range = [dphi - 0.1, dphi + 0.1]
    vals_in_bin = np.array([])
    #pvalues = np.array([])
    
    i=0
    for frame in frames:
        i+=1
        if len(frame['new_photons']) == 0:
           continue
        truth = frame['I3MCTree'][1]
        Event = datacollect(frame)
            # event_xyz = Event[:,0:3]
        event_t = Event[:,3]
    
    
    
            # Calculate Displacement Magnitude
            # diff = np.array([truth.pos.x, truth.pos.y, truth.pos.z]) - event_xyz
            # dr = np.linalg.norm(diff, axis=1)


            # Calculate Time Residual
            # dt = abs(truth.time - event_t) - (1.34*dr/c * 1e9)
        dt = event_t
   

            # Construct Electron direction unit vector from zenith and azimuth
            # Ex = np.sin(truth.dir.zenith)*np.cos(truth.dir.azimuth)
            # Ey = np.sin(truth.dir.zenith)*np.sin(truth.dir.azimuth)
            # Ez = np.cos(truth.dir.zenith)


            # # Calculate angle between electron travel vector and displacement vector
            # Eangle = np.array([Ex, Ey, Ez])
            # Ephi = np.arccos(np.dot(diff, Eangle) / dr)

            # params = np.column_stack([dr, Ephi, dt])

            
            # for i in range(len(params)):
            #     pdf = evalPdf(spline, params[i,0], params[i,1], t)
            #     cdf = np.sum(pdf[:np.searchsorted(t, params[i,2])]) / np.sum(pdf)
            #     pvalues = np.append(pvalues, cdf)


            # mask = (params[:,0] > dr_range[0]) & (params[:,0] < dr_range[1]) & (params[:,1] > dphi_range[0]) & (params[:,1] < dphi_range[1])
            # params_masked = params[mask]
        


            # vals_in_bin = np.append(vals_in_bin, params_masked[:,2].T)
        vals_in_bin = np.append(vals_in_bin, dt)

            # print(len(vals_in_bin) - len(params))

    
    
    t = np.linspace(-5, 5, 1000)
    # pdf = evalPdf(spline, dR, dphi, t)
    pdf = norm.pdf(t)
    hist, edges = np.histogram(vals_in_bin, np.linspace(-5, 5, 500))
    plt.plot(t, pdf, color='red')
    plt.bar(edges[:-1], hist / (np.sum(hist) * (edges[1] - edges[0])), width=edges[1] - edges[0])
    plt.yscale('log')
    plt.savefig('/mnt/home/dillonb5/cascades/plots/sampling_check_norm.png')



    

def main()->None:
    i=0
    frame_list = []
    for i in range(50):
        if i < 10:
            n = '00'+str(i)
        else:
            n = '0'+str(i)
        sampled_file = dataio.I3File(f'/mnt/scratch/dillonb5/sampled_norm/new_{n}.i3.zst')

        for frame in sampled_file:
        
            if len(frame['new_photons']) != 0:
                frame_list.append(frame)
        

    compare_spline([20, 0.3], splinefit_3d, frame_list)

if __name__ == '__main__':
    main()