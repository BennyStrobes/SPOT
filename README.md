# SPOT

SPOT (SPlicing Outlier deTection) is a probabilistic framework to detect Splicing Outilers from  RNA-seq data. Briefly, SPOT fits a Dirichlet-Multinomial distribution directly to counts of reads split across alternatively spliced exon-exon junctions for each gene. SPOT then identifies individuals that deviate significantly from the expectation based on this fitted distribution. Please see [ref bioRxiv preprint] for more details.

## Input data

SPOT identifies outliers at the level of a LeafCutter Cluster (see https://github.com/davidaknowles/leafcutter for more details). Therefor, to generate input data for SPOT you must follow the pre-processing steps described in LeafCutter:
1. Align reads and generate exon-exon junction files. (described in step 1 here: http://davidaknowles.github.io/leafcutter/articles/Usage.html)
2. Run LeafCutter "Intron clustering" (described in step 2 here: http://davidaknowles.github.io/leafcutter/articles/Usage.html). This will generate a file with the suffix "perind_numers.counts.gz". This file is the input file for SPOT. SPOT will accept this file either zipped or un-zipped. Briefly, each column in this file corresponds to an RNA-seq sample and each row corresponds to an intron, which are identified as chromosome:intron_start:intron_end:cluster_id. 

An example SPOT input file can be found under 'example_data/exon_exon_junction_file.txt'

## Running SPOT

Once you have generated the SPOT input file (with help from LeafCutter), SPOT can be easily run using the following command:

python spot.py --juncfile $junction_file_name 


## Testing environment and dependencies
python



## Authors

* **Ben Strober** -- [BennyStrobes](https://github.com/BennyStrobes)