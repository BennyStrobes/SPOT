# SPOT

SPOT (SPlicing Outlier deTection) is a probabilistic framework to detect Splicing Outilers from  RNA-seq data. Briefly, SPOT fits a Dirichlet-Multinomial distribution directly to counts of reads split across alternatively spliced exon-exon junctions for each gene. SPOT then identifies individuals that deviate significantly from the expectation based on this fitted distribution. Please see [ref bioRxiv preprint] for more details.

## Input data

SPOT identifies outliers at the level of a LeafCutter Cluster (see https://github.com/davidaknowles/leafcutter for more details). Therefor, to generate input data for SPOT you must follow the pre-processing steps described in LeafCutter:
1. Align reads and generate exon-exon junction files. (described in step 1 here: http://davidaknowles.github.io/leafcutter/articles/Usage.html)
2. Run LeafCutter "Intron clustering" (described in step 2 here: http://davidaknowles.github.io/leafcutter/articles/Usage.html). This will generate a file with the suffix "perind_numers.counts.gz". This file is the input file for SPOT. SPOT will accept this file either zipped or un-zipped. Briefly, each column in this file corresponds to an RNA-seq sample and each row corresponds to an intron, which are identified as chromosome:intron_start:intron_end:cluster_id. 

An example of a previously generated SPOT input file can be found under 'example_data/exon_exon_junction_file.txt'

It is important to note that in the following paper [ref bioRxiv preprint], we applied the following set of custom filters to the LeafCutter files (before running SPOT) in order to remove exon-exon junctions with low expression while retaining rare exon-exon junctions:
1. Removed exon-exon junctions where no sample has >= 15 split reads
2. Re-defined LeafCutter cluster assignments after removal of exon-exon junctions (according to the above filter) and removed exon-exon junctions that no longer shared a splice site with any other exon-exon junction.
3. Removed all exon-exon junctions belonging to a LeafCutter cluster where more than 10% of the samples had less than 3 reads summed across all exon-exon junctions assigned to that LeafCutter cluster.


## Running SPOT

Once you have generated the SPOT input file (with help from LeafCutter), SPOT can be easily run using the following command:

python spot.py --juncfile $junction_file_name --outprefix $output_root


## Dependencies
Python packages:
1. numpy
2. sys
3. pystan
4. gzip

## Testing environment
SPOT was generated and tested using the following versions:
1. python 2.7.15
2. numpy 1.15.4
3. pystan 2.17.1.0




## Authors

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes)