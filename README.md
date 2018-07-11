# matlab_general

Various pieces of matlab code, to perform a range of operations:

Thresholding
consist_thr_length        construct a group consensus structural (tractography-based) connectome, by retaining most constistent edges *based on length* (as in Misic, Betzel et al. Neuron 2015, and Seidlitz et al. Neuron 2018)
mst_threshold             minimum-spanning tree threshold

Statistics / basic operations
ranksum_stats             code to extract the "rank-biserial r" from a "ranksum" (or Wilcoxon rank-sum, or Mann-Whitney U) test; see Kerby, Compr. Psychol 2014
ranksum_effect_size       extended version of the above
conf_int                  confidence interval

Other
mod_relabel               relabel a modular partition to a reference partition, based on maximum overlap (currenly works only for partitions with the same number of modules)        
