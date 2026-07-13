import photospline

import sys
import os

sys.path.insert(0, os.path.abspath("../"))
from icecube import dataio, icetray, dataclasses, simclasses
import numpy.random as random
import numpy as np
import scipy.optimize as optimize
import matplotlib.pyplot as plt

c = 299792458
n = 1.34


def datacollect(frame):
    dt = []
    t = []
    dphilst = []
    dr = []
    xyz = []

    Epos = frame["I3MCTree"][1].pos
    doms = gcdframe["I3ModuleGeoMap"]
    omkeys = frame["new_photons"].keys()
    photons = frame["new_photons"]
    for key in omkeys:
        # Un-comment if you want to use only a string subset
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
            dt.append(photon.time - offset * 10**9)
            t.append(photon.time)
            phi = photon.dir
            randx, randy, randz = random.uniform(-1, 1, 3)
            x = phi.x
            y = phi.y
            z = phi.z
            dx = x - randx
            dy = y - randy
            dz = z - randz
            dphi = dataclasses.I3Direction(dx, dy, dz)
            dphilst.append(dphi.zenith)
            # Etheta.append(flight.azimuth)
            # Ephi.append(flight.zenith)

    return np.column_stack([xyz, t, dt])



def Likelihood_3d(coords: np.array, Event: np.array):
    L = 0
    # coords should have shape [x,y,z,theta,phi,t]
    # Event has shape [N, 6] cols:(x,y,z,t,dr,dt). We only use the first 4 here

    event_xyz = Event[:, 0:3]
    event_t = Event[:, 3]

    # Calculate Displacement Magnitude
    diff = coords[0:3] - event_xyz
    dr = np.linalg.norm(diff, axis=1)
    # Calculate Time Residual
    dt = event_t - coords[5] - (1.34*dr/c * 1e9)
    # dt = event_t

    # Construct Electron direction unit vector from zenith and azimuth
    Ex = np.sin(coords[4]) * np.cos(coords[3])
    Ey = np.sin(coords[4]) * np.sin(coords[3])
    Ez = np.cos(coords[4])

    # Calculate angle between electron travel vector and displacement vector
    Eangle = np.array([Ex, Ey, Ez])
    Ephi = np.arccos(np.dot(diff, Eangle) / dr)

    # print(Ephi)
    # Calculate Likelihood from constructed coordinates
    params = np.array([dr, Ephi, dt])
    vals = splinefit_3d.evaluate_simple([params[0], params[1], params[2]])
    L = np.where(vals == 0, -30, vals)
    return -np.sum(L)	





def Likelihood_3d_calculation(truth, Event) -> None:
    xary = np.linspace(truth[0] - 150, truth[0] + 150, 500)
    yary = np.linspace(truth[1] - 150, truth[1] + 150, 500)
    zary = np.linspace(truth[2] - 150, truth[2] + 150, 500)
    tary = np.linspace(truth[5] - 500, truth[5] + 500, 500)
    # phiary = np.linspace((truth[4] - 0.3), (truth[4] + 0.3), 500)
    # ThetaAry = np.linspace((truth[3] - 0.3), (truth[3] + 0.3), 500)
    phiary = np.linspace(0, np.pi, 500)
    ThetaAry = np.linspace(0, 2 * np.pi, 500)

    best_fit = np.array(
        [
            frame80["LLHFit_step5"].pos.x,
            frame80["LLHFit_step5"].pos.y,
            frame80["LLHFit_step5"].pos.z,
            frame80["LLHFit_step5"].dir.azimuth,
            frame80["LLHFit_step5"].dir.zenith,
            frame80["LLHFit_step5"].time,
        ]
    )

    xplot = []
    yplot = []
    zplot = []
    tplot = []
    phiplot = []
    thetaplot = []
    for i in range(len(xary)):
        xplot.append(
            Likelihood_3d(
                np.array([xary[i], truth[1], truth[2], truth[3], truth[4], truth[5]]),
                Event=Event,
            )
        )
        yplot.append(
            Likelihood_3d(
                np.array([truth[0], yary[i], truth[2], truth[3], truth[4], truth[5]]),
                Event=Event,
            )
        )
        zplot.append(
            Likelihood_3d(
                np.array([truth[0], truth[1], zary[i], truth[3], truth[4], truth[5]]),
                Event=Event,
            )
        )
        tplot.append(
            Likelihood_3d(
                np.array([truth[0], truth[1], truth[2], truth[3], truth[4], tary[i]]),
                Event=Event,
            )
        )
        phiplot.append(
            Likelihood_3d(
                np.array([truth[0], truth[1], truth[2], truth[3], phiary[i], truth[5]]),
                Event=Event,
            )
        )
        thetaplot.append(
            Likelihood_3d(
                np.array(
                    [truth[0], truth[1], truth[2], ThetaAry[i], truth[4], truth[5]]
                ),
                Event=Event,
            )
        )

    return [xplot, yplot, zplot, thetaplot, phiplot, tplot, best_fit]


def likelihood_graphs(truth, data, filename: str):
    best_fit = data[-1]
    xary = np.linspace(truth[0] - 150, truth[0] + 150, 500)
    yary = np.linspace(truth[1] - 150, truth[1] + 150, 500)
    zary = np.linspace(truth[2] - 150, truth[2] + 150, 500)
    tary = np.linspace(truth[5] - 500, truth[5] + 500, 500)
    # phiary = np.linspace(
    #     (truth[4] - 0.3) * 360 / (2 * np.pi), (truth[4] + 0.3) * 360 / (2 * np.pi), 500
    # )
    # ThetaAry = np.linspace(
    #     (truth[3] - 0.3) * 360 / (2 * np.pi) - 10,
    #     (truth[3] + 0.3) * 360 / (2 * np.pi),
    #     500,
    # )
    phiary = np.linspace(0, 180, 500)
    ThetaAry = np.linspace(0, 360, 500)

    fig, axs = plt.subplots(2, 3, figsize=(12, 10))

    axs[0, 0].plot(xary, data[0])
    axs[0, 0].set_title("Likelihood Evolution in X")
    axs[0, 0].set_xlabel("X (m)")
    axs[0, 0].set_ylabel("-ln(L)")
    axs[0, 0].vlines(
        truth[0], min(data[0]), max(data[0]), label="Truth", colors="red", ls="--"
    )
    axs[0, 0].vlines(
        best_fit[0], min(data[0]), max(data[0]), label = "Best Fit", colors = "orange", ls = "--"
    )
    axs[0, 0].grid()
    # axs[0, 0].set_ylim(min(data[0]) - 10, min(data[0]) + 50)
    # axs[0, 0].set_xlim(truth[0] - 10, truth[0] + 10)
    # #axs[0,0].set_yscale('log')
    axs[0, 0].legend()

    axs[0, 1].plot(yary, data[1])
    axs[0, 1].set_title("Likelihood Evolution in Y")
    axs[0, 1].set_xlabel("Y (m)")
    axs[0, 1].set_ylabel("-ln(L)")
    axs[0, 1].grid()
    # axs[0, 1].set_ylim(min(data[1]) - 10, min(data[1]) + 70)
    # axs[0, 1].set_xlim(truth[1] - 10, truth[1] + 10)
    # #axs[0,1].set_yscale('log')
    axs[0, 1].vlines(
        truth[1], min(data[1]), max(data[1]), label="Truth", colors="red", ls="--"
    )
    axs[0, 1].vlines(
        best_fit[1], min(data[1]), max(data[1]), label = "Best Fit", colors = "orange", ls = "--"
    )
    axs[0, 1].legend()

    axs[1, 0].plot(zary, data[2])
    axs[1, 0].set_title("Likelihood Evolution in Z")
    axs[1, 0].set_xlabel("Z (m)")
    axs[1, 0].set_ylabel("-ln(L)")
    # axs[1, 0].set_ylim(min(data[2]) - 10, min(data[2]) + 50)
    # axs[1, 0].set_xlim(truth[2] - 10, truth[2] + 10)
    # #axs[1,0].set_yscale('log')
    axs[1, 0].vlines(
        truth[2], min(data[2]), max(data[2]), label="Truth", colors="red", ls="--"
    )
    axs[1, 0].vlines(
        best_fit[2], min(data[2]), max(data[2]), label = "Best Fit", colors = "orange", ls = "--"
    )
    axs[1, 0].grid()
    axs[1, 0].legend()

    axs[1, 1].plot(tary, data[5])
    axs[1, 1].set_title("Likelihood Evolution in t")
    axs[1, 1].set_xlabel("t (ns)")
    axs[1, 1].set_ylabel("-ln(L)")
    # axs[1, 1].set_xlim(truth[5] - 10, truth[5] + 10)
    # axs[1, 1].set_ylim(min(data[5]) - 50, min(data[5]) + 100)
    # #axs[1,1].set_yscale('log')
    axs[1, 1].vlines(
        truth[5], min(data[5]), max(data[5]), label="Truth", colors="red", ls="--"
    )
    axs[1, 1].vlines(
        best_fit[5], min(data[5]), max(data[5]), label = "Best Fit", colors = "orange", ls = "--"
    )
    axs[1, 1].grid()
    axs[1, 1].legend()

    axs[1, 2].plot(ThetaAry, data[3])
    axs[1, 2].set_title("Likelihood Evolution in Azimuth")
    axs[1, 2].set_xlabel("Azimuth (deg)")
    axs[1, 2].set_ylabel("-ln(L)")
    axs[1, 2].vlines(
        truth[3] * 360 / (2 * np.pi),
        min(data[3]),
        max(data[3]),
        label="Truth",
        colors="red",
        ls="--",
    )
    axs[1, 2].vlines(
        best_fit[3]*360 / (2*np.pi), min(data[3]), max(data[3]), label = "Best Fit", colors = "orange", ls = "--"
    )
    axs[1, 2].grid()
    # axs[1, 2].set_ylim(min(data[3]) - 10, min(data[3]) + 40)
    # axs[1, 2].set_xlim(
    #     truth[3] * 360 / (2 * np.pi) - 0.3 * 360 / (2 * np.pi),
    #     truth[3] * 360 / (2 * np.pi) + 0.3 * 360 / (2 * np.pi),
    # )
    # axs[1,2].set_yscale('log')
    axs[1, 2].legend()

    axs[0, 2].plot(phiary, data[4])
    axs[0, 2].set_title("Likelihood Evolution in Zenith")
    axs[0, 2].set_xlabel("Zenith (deg)")
    axs[0, 2].set_ylabel("-ln(L)")
    axs[0, 2].vlines(
        truth[4] * 360 / (2 * np.pi),
        min(data[4]),
        max(data[4]),
        label="Truth",
        colors="red",
        ls="--",
    )
    axs[0, 2].vlines(
        best_fit[4]*360 / (2*np.pi), min(data[4]), max(data[4]), label = "Best Fit", colors = "orange", ls = "--"
    )
    axs[0, 2].grid()
    # axs[0, 2].set_ylim(min(data[4]) - 10, min(data[4]) + 40)
    # axs[0, 2].set_xlim(
    #     truth[4] * 360 / (2 * np.pi) - 0.05 * 360 / (2 * np.pi),
    #     truth[4] * 360 / (2 * np.pi) + 0.05 * 360 / (2 * np.pi),
    # )
    # # axs[0,2].set_yscale('log')
    axs[0, 2].legend()

    fig.tight_layout()
    fig.savefig(f"../plots/{filename}.png")


def build_plot(frame):
	mctruth = frame['MCTruth']
	truth = np.array([mctruth.pos.x, 
			mctruth.pos.y,
			mctruth.pos.z,
			mctruth.t,
			mctruth.dir.azimuth,
			mctruth.dir.zenith])
	EventData = datacollect(frame)
     

	
