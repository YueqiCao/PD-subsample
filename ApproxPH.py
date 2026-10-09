# author: Yueqi Cao
# contact: y.cao21@imperial.ac.uk

import numpy as np
import gudhi as gd
import matplotlib.pyplot as plt
import math
import ot
from scipy.optimize import least_squares


# setting parameters
grid_width = 1
nb_units   = 50
float_error = 1e-8
unit = grid_width/nb_units
mat_size = int(nb_units*(nb_units+1)/2)+1

def sample_torus(n, r1, r2):
	"""
	sample n points from a torus 
	with major radius r1 and minor radius r2
	"""
	theta1 = 2 * np.pi * np.random.rand(n)
	theta2 = 2 * np.pi * np.random.rand(n)
	x = (r1 + r2 * np.cos(theta2)) * np.cos(theta1)
	y = (r1 + r2 * np.cos(theta2)) * np.sin(theta1)
	z = r2 * np.sin(theta2)
	X = np.array([x, y, z]).T
	return X

def sample_annulus(n, r1, r2):
	"""
	sample n points from an anulus 
	with major radius r1 and minor radius r2
	"""
	theta = 2 * np.pi * np.random.rand(n)
	rho = np.random.rand(n)

	x = ((1-rho) * r1 + rho * r2) * np.cos(theta)
	y = ((1-rho) * r1 + rho * r2) * np.sin(theta)

	X = np.array([x, y]).T
	return X

def rescale_points(points):
	"""
	rescale point cloud into a unit box
	"""
	points = np.array(points, dtype=float, copy=True)
	for i in range(points.shape[1]):
		max_scale = np.max(points[:,i])
		min_scale = np.min(points[:,i])
		if max_scale == min_scale:
			points[:,i] = 0
			continue
		a = 2/(max_scale-min_scale)
		b = 1 - a * max_scale
		points[:,i] = a * points[:,i] + b
	return points

def get_subsample(large_set, nb_sub_size, nb_sub):
	"""
	sample multiple subsets from original point set
	each subset contains nb_sub_size points 
	""" 
	row_total = large_set.shape[0]
	subsample_set = []
	for i in range(nb_sub):
		row_sample = np.random.choice(row_total, nb_sub_size, replace = False, p = None)
		subsample_set.append(large_set[row_sample,:])
	return subsample_set

def get_PD(points, max_edge_length, min_persistence, sparse=0.3):
	"""
	compute the persistence diagram for VR filtration of point clouds
	"""
	rips = gd.RipsComplex(points=points,\
						 max_edge_length=max_edge_length,\
						 sparse=sparse)
	rips_st = rips.create_simplex_tree(max_dimension=2)
	pers = rips_st.persistence(min_persistence=min_persistence)
	diag = rips_st.persistence_intervals_in_dimension(1)
	return diag

def mesh_gen():
	"""
	generate mesh on the triangular area
	"""
	grid = []
	for i in range(nb_units):
		for j in range(i+1, nb_units+1):
			grid.append([unit*i, unit*j])
	return grid

def dist_mat(grid, power_index):
	"""
	construct the distance matrix of the grid points.
	The underlying distance is L_infty.
	power_index specifies the type of Wasserstein distance.
	"""
	M = np.zeros([mat_size, mat_size])
	for i in range(mat_size-1):
		for j in range(mat_size-1):
			M[i,j] = max([abs(grid[i][0]-grid[j][0]),\
						abs(grid[i][1]-grid[j][1])])
	# append the diagnal
	for k in range(mat_size-1):
		M[k,mat_size-1] = (grid[k][1]-grid[k][0])/2 
		M[mat_size-1,k] = (grid[k][1]-grid[k][0])/2
	Mp = np.power(M,power_index)
	return Mp

def diag_to_mesr(diag, unit_mass):
	"""
	transform persistence diagrams into persistence measures
	unit_mass specifies the weight of each point in a  
	persistence diagram 
	"""
	mesr = np.zeros(mat_size) + float_error
	mesr_vis = np.zeros([nb_units,nb_units])
	for point in diag:
		if not np.all(np.isfinite(point)) or not (0 <= point[0] <= point[1] <= grid_width):
			raise ValueError("Expected finite 0 <= birth <= death <= grid_width; cap infinite deaths at the filtration cutoff first.")
		i = math.floor(point[0]/unit + float_error)
		j = math.floor(point[1]/unit + float_error)
		if i == j:
			continue
		mesr_vis[i,min(j,nb_units-1)] += unit_mass
		mesr[nb_units*i+j-int(i*(i+1)/2)-1] += unit_mass
	return mesr, mesr_vis

def wass_dist(a, b, Mp, reg=1):
	"""
	Return a powered partial-transport cost with diagonal matching.
	reg=1 preserves the original Sinkhorn calculation (a regularized proxy).
	Use reg=None for exact transport; only its pth root is a distance.
	""" 
	a_mesr = np.asarray(a, dtype=float).copy()
	b_mesr = np.asarray(b, dtype=float).copy()
	a_mesr[-1], b_mesr[-1] = b_mesr[:-1].sum(), a_mesr[:-1].sum()
	if a_mesr.sum() == 0:
		return 0.0
	active_a, active_b = a_mesr > 0, b_mesr > 0
	cost = Mp[np.ix_(active_a, active_b)]
	if reg is None:
		# Scaling keeps small p=8 costs above the solver's absolute tolerance.
		scale = cost.max()
		if scale == 0:
			return 0.0
		return float(ot.emd2(a_mesr[active_a], b_mesr[active_b], cost / scale,
							numItermax=1000000) * scale)
	return float(ot.sinkhorn2(a_mesr[active_a], b_mesr[active_b], cost, reg))

def func(x, a, b, c):
	"""
	function model to fit
	"""
	return a * x ** (-b) + c

def fit_power_law(point_list, dist_list):
	"""Fit a*n**(-b)+c with a free exponent and a nonnegative offset."""
	n, y = np.asarray(point_list, dtype=float), np.asarray(dist_list, dtype=float)
	if n.ndim != 1 or y.shape != n.shape or len(n) < 4 or len(np.unique(n)) < 4:
		raise ValueError("Provide at least four distinct sample sizes and their losses.")
	if not np.all(np.isfinite(n)) or not np.all(np.isfinite(y)) or np.any(n <= 0) or np.any(y < 0) or np.ptp(y) == 0:
		raise ValueError("Sample sizes must be positive and losses finite, nonnegative and nonconstant.")
	x, scale = n / n.min(), y.max()
	# Normalize both axes: the old fixed parameter bounds distorted small losses.
	fits = [least_squares(lambda t: func(x, *t) - y / scale,
						  [1, exponent, 0], bounds=([0, 0.01, 0], [np.inf, 5, np.inf]),
						  max_nfev=5000, ftol=1e-12, xtol=1e-12, gtol=1e-12)
			for exponent in (0.25, 0.5, 1, 2)]
	best = min(fits, key=lambda result: result.cost)
	if not best.success:
		raise RuntimeError(best.message)
	a, b, c = best.x
	return np.array([a * scale * n.min() ** b, b, c * scale])

def plot_diag(diag):
	"""
	plot persistence diagrams
	"""
	fig = plt.figure(figsize=(8,8))
	plt.rcParams.update({'font.family':'Times New Roman', 'font.size':22})
	#plt.rc('text', usetex=True)
	ax = fig.add_subplot(111)
	ax.scatter(diag[:,0], diag[:,1], s=75, marker='o', c='red', alpha=0.8)
	ax.plot([0,1], [0,1], linewidth=0.5)
	plt.fill_between([0,1], [0,1], [0,0], facecolor='green', alpha=0.2)	
	ax.set_xlim((0,1))
	ax.set_ylim((0,1))
	ax.set_xlabel("Births")
	ax.set_ylabel("Deaths")
	ax.set_title("Persistence Diagram")


def plot_mesr(mesr_vis):
	"""
	plot persistence measure
	"""
	fig = plt.figure(figsize=(8,8))
	plt.rcParams.update({'font.family':'Times New Roman', 'font.size':22})
	#plt.rc('text', usetex=True)
	ax = fig.add_subplot(111)
	plt.imshow(mesr_vis.T, origin='lower', cmap='hot_r', interpolation='bilinear')
	L = mesr_vis.shape[0]
	ax.set_xlim((0,L))
	ax.set_ylim((0,L))
	loc, labels = plt.xticks()
	plt.xticks(loc, loc/nb_units)
	plt.yticks(loc, loc/nb_units)
	ax.set_xlabel("Births")
	ax.set_ylabel("Deaths")
	plt.fill_between([0,L], [0,L], [0,0], facecolor='green', alpha=0.2)	
	ax.set_title("Mean Persistence Diagram")
	plt.show()

def plot_fitting_curve(point_list, dist_list):
	point_list = np.asarray(point_list, dtype=float)
	pram = fit_power_law(point_list, dist_list)
	fig = plt.figure(figsize=(8,8))
	plt.rcParams.update({'font.family':'Times New Roman', 'font.size':24})
	plt.rc('text', usetex=False)
	plt.plot(point_list, func(point_list, *pram), linewidth=2,\
			label='Fit: $y = %5.2f x^{-%5.2f}+%5.2f$' % tuple(pram))
	plt.scatter(point_list, dist_list, s=75, color='red', marker='o',\
			label="True Loss")
	plt.plot(point_list, dist_list, 'g--', label='Empirical Curve')
	plt.legend()
	plt.xlabel(r"$m$: Number of Points")
	plt.ylabel(r"$\mathrm{OT}(\hat{D},D_S)$")
	plt.title(r"Convergence Rate")
	plt.show()
