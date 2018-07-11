# matlab_general

Various pieces of matlab code, to perform a range of operations:

### Thresholding
  * *consist_thr_length*        &nbsp;&nbsp;&nbsp;  construct a group consensus structural (tractography-based) connectome, by retaining most constistent edges **based on length** (as in Misic, Betzel et al. Neuron 2015, and Seidlitz et al. Neuron 2018)  
  * *mst_threshold*             &nbsp;&nbsp;&nbsp;  minimum-spanning tree threshold  

### Statistics / basic operations
  * *ranksum_stats*             &nbsp;&nbsp;&nbsp;  code to extract the "rank-biserial r" from a "ranksum" (or Wilcoxon rank-sum, or Mann-Whitney U) test; see Kerby, Compr. Psychol 2014  
  * *ranksum_effect_size*       &nbsp;&nbsp;&nbsp;  extended version of the above  
  * *conf_int*                  &nbsp;&nbsp;&nbsp;  confidence interval  

### Other
  * *mod_relabel*               &nbsp;&nbsp;&nbsp;  relabel a modular partition to a reference partition, based on maximum overlap (currenly works only for partitions with the same number of modules)  
