import numpy as np
import sys
import pystan
import dirichlet_multinomial
import gzip

# Extract data structure and sample anmes
def extract_raw_cluster_jxn_data_structure(jxn_file):
	# Used to skip header
	head_count = 0
	# Initialize cluster_jxn_data_structure
	cluster_jxn_data_structure = {}
	# Stream input file
	if jxn_file.endswith('.gz'):
		f = gzip.open(jxn_file)
	else:
		f = open(jxn_file)
	for line in f:
		line = line.rstrip()
		data = line.split()
		# Header
		if head_count == 0:
			head_count = head_count + 1
			samples = np.asarray(data[1:])
			continue
		# Standard line
		# Extract Jxn info from current line
		key = data[0]
		key_info = key.split(':')
		chrom_num = key_info[0]
		start = key_info[1]  # 5' ss
		end = key_info[2]  # 3' ss
		cluster_id = key_info[3]  # Name of cluster
		jxn_read_counts = np.asarray(data[1:]).astype(float)  # Vector of raw read counts for this junction

		# Add jxn to  cluster_jxn_data_structure
		if cluster_id not in cluster_jxn_data_structure:  # If we've never seen this cluster before
			cluster_jxn_data_structure[cluster_id] = []
		cluster_jxn_data_structure[cluster_id].append(jxn_read_counts)  # Add read counts from current junction
	f.close()
	# Convert from list of arrays to matrix (in each cluster)
	# Loop through all clusters
	all_clusters = cluster_jxn_data_structure.keys()
	for cluster_id in all_clusters:
		# Create matrix from list of arrays
		jxn_matrix = np.transpose(np.asmatrix(cluster_jxn_data_structure[cluster_id]))
		
		num_jxns = jxn_matrix.shape[1]
		if num_jxns < 2:
			sys.stderr.write('Error: ' + cluster_id + ' only has ' + str(num_jxns) + ' junction mapped to it. 2 junctions are needed for all clusters.\n')
			exit(0)
		# Add this new matrix to the data structure
		cluster_jxn_data_structure[cluster_id] = {}
		cluster_jxn_data_structure[cluster_id]['jxn_matrix'] = jxn_matrix
	return cluster_jxn_data_structure, samples


# Remove clusters with more than $max_number_of_junctions_per_cluster
def max_number_of_jxns_filter_ignore_genes(cluster_jxn_data_structure, max_number_of_junctions_per_cluster):
	# Initialize new cluster_jxn_data_structure
	new_jxn_structure = {}
	# Loop through clusters
	for cluster_id in cluster_jxn_data_structure.keys():
		jxn_mat = cluster_jxn_data_structure[cluster_id]['jxn_matrix']
		# Get dimensionality of this cluster
		N, K = jxn_mat.shape
		# If there are more than max_number_of_junctions_per_cluster
		if K > max_number_of_junctions_per_cluster:
			pass
		else: # Less than or equal
			new_jxn_structure[cluster_id] = {}
			new_jxn_structure[cluster_id]['jxn_matrix'] = jxn_mat
	return new_jxn_structure

# Load in data into dictionary structure where keys are cluster ids and values are junction count matrices for that cluster
def load_in_junction_count_data(junc_file, max_junctions):
	# Extract data structure and sample anmes
	cluster_jxn_data_structure, sample_names = extract_raw_cluster_jxn_data_structure(junc_file)

	# Remove clusters with more than $max_junctions junctions per cluster
	cluster_jxn_data_structure = max_number_of_jxns_filter_ignore_genes(cluster_jxn_data_structure, max_junctions)

	return cluster_jxn_data_structure, sample_names

# Print outlier calling dm results to output file
def outlier_calling_print_helper(arr, all_samples, t, cluster_id):
	cluster_samples = all_samples
	counter = 0
	# Print Row id
	t.write(cluster_id)
	# Loop through all samples
	for sample in all_samples:
		# If sample in cluster specific samples
		if sample in cluster_samples:
			ele = str(arr[counter])
			counter = counter + 1
		# Samples was filtered out of this cluster
		else:
			ele = 'NaN'
		t.write('\t' + ele)
	t.write('\n')
	t.flush()
	return t



def call_splicing_outliers_shell(out_prefix, cluster_jxn_data_structure, sample_names, num_background_samples, num_simulated_reads, seed):
	# Seed for random number generator used in simulating samples for mahalanobis distance emperical distribution
	np.random.seed(seed)
	
	# Load in pystan optimizizer
	DM_GLM = pystan.StanModel(file = "dirichlet_multinomial.stan")


	#Initialize output files
	t_MD = open(out_prefix + '_md.txt','w')  # Filehandle for matrix of mahalanobis distances
	t_pvalue = open(out_prefix + '_emperical_pvalue.txt', 'w')  # Filehandle for matrix of pvalues
	
	# Write headers for output files
	t_MD.write('CLUSTER_ID\t' + '\t'.join(sample_names) + '\n')
	t_pvalue.write('CLUSTER_ID\t' + '\t'.join(sample_names) + '\n')
	
	# Loop through clusters
	for counter, cluster_id in enumerate(sorted(cluster_jxn_data_structure.keys())):
		####################################################################
		# Outlier calling for one cluster
		####################################################################
		# Extract jxn matrix for this gene
		X = cluster_jxn_data_structure[cluster_id]['jxn_matrix']
		# Run outlier analysis:
		# Return:
		#   1: mahalanobis_distances: vector length num_samples where each element is the mahalanobis distance for that sample
		#   2. pvalues: vector of length num_samples where each element is the pvalue for that sample
		#   3. alpha: vector of length num_jxns which defines the fitted dirichlet multinomial distribution
		try:
			mahalanobis_distances, pvalues, alpha = dirichlet_multinomial.run_dm_outlier_analysis(X, num_background_samples, num_simulated_reads, seed, DM_GLM)
			####################################################################
			# Print results to output file
			####################################################################
			# Print Mahalanobis distance results to output file
			t_MD = outlier_calling_print_helper(mahalanobis_distances, sample_names, t_MD, cluster_id)
			# Print emperical pvalue resutls to output file
			t_pvalue = outlier_calling_print_helper(pvalues, sample_names, t_pvalue, cluster_id)
		except:
			print('miss: ' + cluster_id)
	t_MD.close()
	t_pvalue.close()




def main(options):
	# Extract command line arguments
	junc_file = options.juncfiles
	out_prefix = options.outprefix
	max_junctions = options.maxjunctions
	num_background_samples = options.numbackgroundsamples
	num_simulated_reads = options.numsimulatedreads
	seed = options.seed

	# Load in data into dictionary structure where keys are cluster ids and values are junction count matrices for that cluster
	cluster_jxn_data_structure, sample_names = load_in_junction_count_data(junc_file, max_junctions)

	# Shell to call splicing outliers
	call_splicing_outliers_shell(out_prefix, cluster_jxn_data_structure, sample_names, num_background_samples, num_simulated_reads, seed)



if __name__ == "__main__":
	from optparse import OptionParser

	parser = OptionParser()

	parser.add_option("-j", "--juncfile", dest="juncfiles",
		help="text file with all junction files to be processed")

	parser.add_option("-o", "--outprefix", dest="outprefix", default = 'spot',
		help="output prefix (default spot)")

	parser.add_option("-n", "--maxjunctions",type="int", dest="maxjunctions", default = 20,
		help="maximum number of junctions per LeafCutter cluster")

	parser.add_option("-e", "--numbackgroundsamples",type="int", dest="numbackgroundsamples", default = 1000000,
		help="Number of randomly drawn samples per cluster used to generate an emperical mahalanobis distance distribution")

	parser.add_option("-r", "--numsimulatedreads",type="int", dest="numsimulatedreads", default = 20000,
		help="Number of reads per simulated sample")

	parser.add_option("-s", "--seed",type="int", dest="seed", default = 1,
		help="Seed used for random number generator in both optimization and generating random samples for a mahalanobis distance emperical distribution")

	(options, args) = parser.parse_args()
	if options.juncfiles == None:
		sys.stderr.write("Error: no junction file provided...\n")
		exit(0)

	main(options)

