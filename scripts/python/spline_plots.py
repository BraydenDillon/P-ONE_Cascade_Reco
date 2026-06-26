import numpy as np
import photospline

with np.load("500_Event_hist.npz") as f:
    hist = f["hist"]
    bins_t = f["t_edges"]
    bins_r = f["r_edges"]
    bins_theta = f["dtheta_edges"]
    bins_Etheta = f["Etheta_edges"]

with np.load("/mnt/home/dillonb5/cascades/outfiles/ary_012.npz") as f:
    d = f["data"]

t_centers = (bins_t[1:] + bins_t[:-1]) / 2
r_centers = (bins_r[1:] + bins_r[:-1]) / 2
theta_centers = (bins_theta[1:] + bins_theta[:-1]) / 2
E_centers = (bins_Etheta[1:] + bins_Etheta[:-1]) / 2

amplitudes = np.sum(hist, axis=0)
print(
    np.count_nonzero(amplitudes == 0)
    / (np.count_nonzero(amplitudes == 0) + np.count_nonzero(amplitudes))
)
pdfs = hist / amplitudes
print(pdfs.shape)

pdfs[~np.isfinite(pdfs)] = 0
cdfs = np.cumsum(pdfs, axis=0)

errs = np.sqrt(hist)

spline = photospline.SplineTable("../fits/splinelog.fits")
t_knots = np.concatenate(
    (
        [-20, -10, -5, -3, -1],
        np.linspace(0, 40, 20),
        np.geomspace(40.1, max(t_centers), 10),
    )
)
import matplotlib.pyplot as plt

t = np.linspace(0, 1400, 2500)

plt.plot(t, (spline.evaluate_simple([t, 25, np.pi / 2, np.pi / 2])))
plt.plot(t_centers, pdfs[:, 50, 50, 39])
plt.scatter(
    t_knots, spline.evaluate_simple([t_knots, 25, np.pi / 2, np.pi / 2]), marker="o"
)
plt.savefig("PDF_knots_overplotted.png")
plt.clf()

plt.plot(t_centers, (spline.evaluate_simple([t_centers, 210.5, np.pi / 2, np.pi / 2])))
# plt.plot(t_centers,(spline.evaluate_simple([t_centers, 20.05, np.pi, np.pi]) - pdfs[:,99, 99, 39])/np.max(pdfs[:,99, 99, 39]))
plt.savefig("Spline_plot_2.png")
plt.clf()
pvalues = []
pdfs = []
# photon shape: cols=['dt', 'dZenith', 'dAzimuth', 'dr', 'EZenith', 'EAzimuth', 'eNRG', 'Lambda']
import time


def truncate(f, n):
    """Truncates/pads a float f to n decimal places without rounding"""
    s = "{}".format(f)
    if "e" in s or "E" in s:
        return "{0:.{1}f}".format(f, n)
    i, p, d = s.partition(".")
    return ".".join([i, (d + "0" * n)[:n]])


start = time.time()
for photon in d[:10000]:
    #        print(photon[0])
    if photon[0] > t[-1]:
        continue
    pdf = spline.evaluate_simple([t, photon[3], photon[1], photon[4]])
    cdf = np.sum(pdf[np.searchsorted(t, photon[0]) :]) / np.sum(pdf)
    pvalues.append(cdf)
# pdfs.append(np.sum(pdf))
end = time.time()
print(f"time to calculate cdf: {end - start}")
i = 0
for photon in d[:100]:
    if photon[0] > t_centers[-1]:
        continue
    plt.figure()
    pdf = spline.evaluate_simple([t, photon[3], photon[1], photon[4]])
    plt.plot(t, pdf)
    # plt.vlines(photon[0], 0, max(pdf), color = 'red', ls = '--',label=f't={truncate(photon[0], 5)}')
    plt.title(
        f"Time PDF, r={truncate(photon[3], 3)}, dphi={truncate(photon[1],3)}, Ephi={truncate(photon[4], 3)}"
    )
    plt.text(
        0.7,
        0.8,
        f"CDF={truncate(np.sum(pdf[np.searchsorted(t, photon[0]):])/np.sum(pdf),5)}",
        transform=plt.gca().transAxes,
        fontsize=10,
    )
    plt.scatter(
        t_knots,
        spline.evaluate_simple([t_knots, photon[3], photon[1], photon[4]]),
        marker="o",
    )
    plt.xlabel("dt (ns)")
    plt.ylabel("Probability Density")
    # 	plt.legend()
    i += 1
    plt.savefig(f"/mnt/home/dillonb5/cascades/plots/cascade_pdf_plot_{i}.png")
    plt.close()
i = 0
r = np.linspace(0, 460, 1000)
for photon in d[:100]:
    if photon[0] > t_centers[-1]:
        continue
    plt.figure()
    pdf = spline.evaluate_simple([t, photon[3], photon[1], photon[4]])
    dt = t[1] - t[0]
    pdf_norm = pdf / sum(pdf * dt)
    # plt.plot(t, pdf_norm, color = 'red', label = 'PDF from Spline Fit')
    r_idx = np.searchsorted(bins_r, photon[3]) - 1
    dth_idx = np.searchsorted(bins_theta, photon[1]) - 1
    Eth_idx = np.searchsorted(bins_Etheta, photon[4]) - 1
    i += 1
    hist_slice = hist[:, r_idx, dth_idx, Eth_idx]
    dt = np.diff(bins_t)
    # 	dt = np.append(dt, dt[-1])
    normalized_slice = hist_slice / sum(hist_slice * dt)
    plt.bar(bins_t[:-1], normalized_slice, width=np.diff(bins_t))
    plt.plot(t, pdf_norm, color="red", label="PDF from Spline Fit")
    plt.title("dt fit")
    plt.xlabel("dt (ns)")
    # plt.xlim(-1, 400)
    plt.ylabel("Probability Density, Counts")
    plt.legend()
    plt.savefig(f"/mnt/home/dillonb5/cascades/plots/t_hist_comparison_{i}.png")
    plt.close()
# 	rpdf = spline.evaluate_simple([photon[0], r, photon[1], photon[4]])
# 	dr = r[1] - r[0]
# 	rpdf_norm = rpdf / sum(rpdf*dr)
# 	plt.plot(r, rpdf_norm, color = 'red')
# 	t_idx = np.searchsorted(bins_t, photon[0]) - 1
# 	rhist_slice = hist[t_idx, :, dth_idx, Eth_idx]
# 	normalized_slice = rhist_slice / sum(rhist_slice)
# 	plt.bar(bins_r[:-1], normalized_slice, width = np.diff(bins_r))
# 	plt.title('dr fit')
# 	plt.savefig(f'/mnt/home/dillonb5/cascades/plots/r_hist_comparison_{i}.png')
# 	plt.close()
# print(pvalues)
# print(pdfs)
plt.figure()
plt.hist(pvalues, bins=np.linspace(0, 1.01, 500))
plt.title("CDF Distribution")
plt.xlabel("CDF")
plt.ylabel("Counts")
plt.savefig("CDF_Plot.png")
