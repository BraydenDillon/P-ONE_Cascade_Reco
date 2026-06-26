import numpy as np
from pathlib import Path
runnumbers = []
for i in range(0,100):
	if i < 10:
		runnumber = '00'+str(i)
	elif i < 100:
		runnumber = '0'+str(i)
	else:
		runnumber = str(i)
	runnumbers.append(runnumber)


t_min = 0.01
t_max = []
r_min = 0.01
r_max = []
dtheta_min = []
dtheta_max = []
Etheta_min = []
Etheta_max = []
for run in runnumbers:
	filename = Path("/mnt/home/dillonb5/cascades/outfiles/ary_"+run+".npz")
	if not filename.exists():
		continue
	with np.load("/mnt/home/dillonb5/cascades/outfiles/ary_"+run+".npz") as f:
		data = f['data']
		t_max.append(data[:,0].max())
		r_max.append(data[:,3].max())
		dtheta_min.append(data[:,1].min())
		dtheta_max.append(data[:,1].max())
		Etheta_min.append(data[:,4].min())
		Etheta_max.append(data[:,4].max())

bins = [100, 100, 100, 40]
t_edges = np.geomspace(t_min, max(t_max), bins[0]+1)
r_edges = np.geomspace(r_min, max(r_max), bins[1]+1)
dtheta_edges = np.geomspace(min(dtheta_min), max(dtheta_max), bins[2]+1)
Etheta_edges = np.geomspace(min(Etheta_min), max(Etheta_max), bins[3]+1)
all_edges = [t_edges, r_edges, dtheta_edges, Etheta_edges]
hist = np.zeros(bins, dtype = np.int64)
for run in runnumbers:
	filename = Path("/mnt/home/dillonb5/cascades/outfiles/ary_"+run+".npz")
	if not filename.exists():
		continue
	with np.load("/mnt/home/dillonb5/cascades/outfiles/ary_"+run+".npz") as f:
		data = f['data']
		arrays = np.column_stack([data[:, 0], data[:,3], data[:,1], data[:,4]])
		histNd, _ = np.histogramdd(arrays, bins = all_edges)
		hist = hist + histNd

np.savez('/mnt/home/dillonb5/cascades/scripts/100_file_hist_v5_courser_bins.npz', hist = hist, 
t_edges = t_edges, 
r_edges = r_edges, 
dtheta_edges = dtheta_edges,
Etheta_edges = Etheta_edges,
labels = ["dt", "dr", "dtheta", "Etheta"])
