import numpy as np
import photospline
import resource
with np.load('500_Event_hist.npz') as f:
	hist = f['hist']
	bins_t = f['t_edges']
	bins_r = f['r_edges']
	bins_theta = f['dtheta_edges']
	bins_Etheta = f['Etheta_edges']

t_centers = (bins_t[1:] + bins_t[:-1]) / 2
r_centers = (bins_r[1:] + bins_r[:-1]) / 2
theta_centers = (bins_theta[1:] + bins_theta[:-1]) / 2
E_centers = (bins_Etheta[1:] + bins_Etheta[:-1]) / 2

dt = np.diff(bins_t)
dr = np.diff(bins_r)
dth = np.diff(bins_theta)
dE = np.diff(bins_Etheta)
vol = dt[:,None,None,None] * dr[None,:,None,None] * dth[None,None,:,None] * dE[None,None,None,:]

print('histogram loaded, centers defined')

print(resource.getrusage(resource.RUSAGE_SELF)[2])

amplitudes = np.sum(hist, axis=0)


print(np.count_nonzero(amplitudes==0) / (np.count_nonzero(amplitudes==0) + np.count_nonzero(amplitudes)))
pdfs = hist / (amplitudes*vol)
pdfs = np.hstack((np.zeros_like(pdfs), pdfs))
print('pdfs calculated')

pdfs[~np.isfinite(pdfs)] = 0
cdfs = np.cumsum(pdfs, axis=0)

errs = np.sqrt(hist)

print(resource.getrusage(resource.RUSAGE_SELF)[2])
knots = [
np.concatenate(([-30,-20,-10,-5,-3,-1],np.linspace(0,75,20),np.geomspace(75.1, max(t_centers), 10))),
np.concatenate(([-30,-20,-10,-5,-3,-1],np.linspace(0,60,20),np.geomspace(60.1, max(r_centers), 10))),
np.concatenate(([-np.pi/2, -np.pi/3, -np.pi/4],np.linspace(0, np.pi, 20))),
np.concatenate(([-np.pi/2, -np.pi/3, -np.pi/4],np.linspace(0, np.pi, 20)))
]

print(knots[0])
w = errs**(-2)
w = np.hstack((np.full(w.shape, 1), w))
print('w calculated, knots defined')
print(resource.getrusage(resource.RUSAGE_SELF)[2])
w[~np.isfinite(w)] = np.nanmin(w) # was nanmin(w), changed to 0 to hopefully save on memory
#pdfs = np.log(pdfs)
#w = np.log(w)
sparse = photospline.ndsparse.from_data(
        pdfs, w)

print('ndsparse defined')
print(resource.getrusage(resource.RUSAGE_SELF)[2])

spline = photospline.glam_fit(
	*sparse,
	coordinates = [t_centers, r_centers, theta_centers, E_centers],
	knots = knots,
	order = [3,3,3,3],
	smoothing = [0.1, 0.1, 0.1, 0.1],
	penaltyOrder = [2,2,2,2],
	verbose = True
)


print('spline defined')
print(resource.getrusage(resource.RUSAGE_SELF)[2])
spline.write('/mnt/home/dillonb5/cascades/fits/spline_testing_5-22-v1.fits')

