import numpy as np


#returns matrix of size P X K where P is the number of covariates and K is the number of junctions
def correct_betas(beta_raw_object,beta_scale,K_input,P_input):
	#deal with 1 covariate case seperatelyseperately
	if P_input == 1:
		corrected_betas = (beta_raw_object  - (1.0/K_input))*np.asmatrix(beta_scale)[0,0]
	#More than 1 covariate.
	elif P_input > 1:
		corrected_betas = []
		P,K = beta_raw_object.shape
		if len(beta_scale) != P_input or K != K_input:
			print('error')
			pdb.set_trace()
		#loop through covariates
		for p in range(P):
			new_beta = (beta_raw_object[p,:] - (1.0)/K_input)*beta_scale[p]
			corrected_betas.append(new_beta)
		corrected_betas = np.asmatrix(corrected_betas)
	return corrected_betas

def compute_alphas_intercept_only_multi_conc(betas,conc_param):
	term_a = (np.exp(betas)/np.sum(np.exp(betas)))
	alphas = []
	for i,ele in enumerate(conc_param):
		alphas.append(ele*term_a[0,i])
	return alphas


def dirichlet_multinomial_fit(y, seed, DM_GLM):
	#fixed parameters (provide relaxed priors to add optimization)
	concShape=1.0001
	concRate=1e-4
	#Make intercept term (in covariate matrix)
	N,K = y.shape
	x = np.ones((N,1))
	N,P = x.shape
	# Put data in dictionary (required input for pystan)
	data = dict(N=N, K=K, P = P,y = y, x = x, concShape = concShape,concRate = concRate)
	# Optimize GLM using pystan
	op = DM_GLM.optimizing(data = data,verbose=False,iter=5000,seed=seed)
	#Convert betas from simplex space to real space
	betas = correct_betas(op['beta_raw'],op['beta_scale'],K,P)
	#compute actual alpha that defines DM
	alphas = compute_alphas_intercept_only_multi_conc(betas,op['conc'])
	return np.asarray(alphas)

def compute_dm_covariance_matrix(n,alpha):
	alpha_0 = np.sum(alpha)
	p = alpha/alpha_0
	cov = n*((n+alpha_0)/(1+alpha_0))*(np.diag(p) - np.dot(np.transpose(np.asmatrix(p)),np.asmatrix(p)))
	return cov

# Compute mahalanobis distance for one sample
def compute_mahalanobis_distance(x, alpha):
	nn = np.sum(x)
	alpha_0 = np.sum(alpha)
	# Compute covariance matrix (dimension num_jxnsXnum_jxns)
	cov = compute_dm_covariance_matrix(nn,alpha)
	mu = nn*alpha/alpha_0
	diff_mat = np.asmatrix(x-mu)
	distance = np.dot(np.dot(diff_mat,np.linalg.pinv(cov)),np.transpose(diff_mat))[0,0]
	return np.sqrt(distance)

def generate_background_mahalanobis_distances(alpha, num_samples, num_reads):
	##################
	# Draw random sample according to fitted dirichlet-multinomial
	###################
	x = []
	sample_alphas = np.random.dirichlet(alpha,size=num_samples)
	for i in range(num_samples):
		x.append(np.random.multinomial(num_reads, sample_alphas[i,:]))
	x = np.asmatrix(x)

	##################
	# Compute mahalanobis distances for all of these samples 
	###################
	# initialize vector to keep track of MD's
	sample_mahalanobis_distances_old = []
	# Compute quantities required for MD (that are shared across all samples)
	nn = np.sum(x[0,:])
	alpha_0 = np.sum(alpha)
	cov = compute_dm_covariance_matrix(nn, alpha)

	mu = nn*alpha/alpha_0
	inv_cov = np.linalg.pinv(cov)
	diff_mat = x - mu

	# Compute MD for each sample
	# for sample_num in range(num_samples):
		#sample_mahalanobis_distances_old.append(compute_mahalanobis_distance_with_precomputed_quantities(x[sample_num,:], mu, inv_cov))

	# Compute MD for each sample using vector math
	sample_mahalanobis_distances = np.sqrt(np.sum(np.multiply(np.dot(diff_mat,inv_cov),diff_mat),axis=1))
	# Convert from matrix to array
	sample_mahalanobis_distances = np.squeeze(np.asarray(sample_mahalanobis_distances))

	return np.asarray(sample_mahalanobis_distances)

def run_dm_outlier_analysis(X, num_background_samples, num_simulated_reads, seed, DM_GLM):
	#########################################################
	# Fit Dirichlet multinomial based on samples
	#########################################################
	# Fit dirichlet multinomial to X (return dm parameter alpha)
	alpha = dirichlet_multinomial_fit(X, seed, DM_GLM)

	#########################################################
	# Compute mahalanobis distance (MD) for all observed samples
	#########################################################
	# Initialize array to keep track of mahalanobis distances from each samples
	mahalanobis_distances = []
	# Compute mahalanobis distance for each sample. 
	num_samples, num_jxns = X.shape
	# Loop through samples
	for sample_num in range(num_samples):
		mahalanobis_distances.append(compute_mahalanobis_distance(X[sample_num,:], alpha))
	# Convert to numpy array
	mahalanobis_distances = np.asarray(mahalanobis_distances)

	#########################################################
	# Estimate emperical distribution
	#########################################################
	#Take num_samples for the fitted DM. For each draw, compute mahalanobis distance
	sample_mahalanobis_distances = generate_background_mahalanobis_distances(alpha, num_background_samples, num_simulated_reads)

	#########################################################
	# Compute pvalues for observed samples based on emperical distribution
	#########################################################
	# Initialize output vector
	pvalues = []
	# Loop through observed samples' distances
	for distance in mahalanobis_distances:
		# Compute pvalue for this sample
		pvalue = len(np.where(distance <= sample_mahalanobis_distances)[0])/float(len(sample_mahalanobis_distances))
		pvalues.append(pvalue)
	# Convert to numpy array
	pvalues = np.asarray(pvalues)


	return mahalanobis_distances, pvalues, alpha
