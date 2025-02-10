#! /bin/bash

cd $1/phylo_tree
qiime tools import --input-path fasttree_aa --output-path unrooted_tree --type Phylogeny[Unrooted]	# convert tree Newick format to Phylogeny[Unrooted] qiime2 artifact
qiime phylogeny midpoint-root --i-tree unrooted_tree.qza --o-rooted-tree rooted_tree.qza    # make rooted tree
qiime tools import --type FeatureData[Taxonomy] --input-path $1/drep/mags_metadata.tsv --output-path mags_metadata.qza	# create feature-metadata-file 
qiime empress tree-plot --i-tree rooted_tree.qza --m-feature-metadata-file mags_metadata.qza --o-visualization phylo_tree_viz.qzv