# Rich-clubness test

## References

* *Cannistraci-Muscoloni null-model and rich-clubness test*  
  Muscoloni, A. and Cannistraci, C.V. (2017) "Rich-clubness test: how to determine whether a complex network has or doesn't have a rich-club?". arXiv:1704.03526
* *Maslov-Sneppen null-model*  
  Maslov, S. and Sneppen, K. (2002) "Specificity and Stability in Topology of Protein Networks". Science 296:910

## Folder description

* *randomize_network.m*: it performs a randomization of the network in input preserving the node degrees. The randomization consists of an iterative random reshuffle of link pairs. At each iteration two links are randomly sampled (uniformly or nonuniformly, according to the Maslov-Sneppen or Cannistraci-Muscoloni null-model respectively) and one endpoint of the first is randomly exchanged with one endpoint of the second. If the link that would be created already exists, the attempt is rejected.
* *richclub_test.m*: it performs the statistical test for rich-clubness giving in output a pvalue that indicates whether the network contains a significant rich-club, and a degree-cut that allows to extract the rich-club subnetwork.
* *RUN.m*: usage example that loads an example network (*network.mat*), generates 1000 randomized networks using the CM null-model and performs the rich-clubness test.

## Contact

For any problem, please contact:
* Alessandro Muscoloni: alessandro.muscoloni@gmail.com
* Carlo Vittorio Cannistraci: kalokagathos.agon@gmail.com