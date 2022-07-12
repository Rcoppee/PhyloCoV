# Scripts used for the study "Phylodynamic of SARS-CoV-2 transmission in France at local and international scales during 2020"
<p>xxx<br>
 xxx</p> <br>
<h3>1. Prerequisites</h3>
 <p>To use these scripts, you need the following programs:</p>
 <p>- TreeTime<br>
 - R<br>
 - R packages</p>
 <br>
 <h3>Dataset generation (ANTOINE)</h3>
 <p>xxx</p>
 <br>
 <h3>Dating the phylogenetic trees using TreeTime</h3>
 <p> We used TreeTime to date each phylogenetic tree previously produced. TreeTime requires a phylogenetic tree (newick format), a multiple sequence alignment and text file that indicates the date of sample collection for each tip in the phylogeny. Prior to TreeTime, we developed a R script, named <i><b>get_dates_from_tips.R</b></i> that produces such a file (one output per replicate with the <i>name X_date.txt</i>, where <i>X</i> is the number of the replicate). Then, we produced a shell script, named <i><b>launch_treetime.sh</b></i> that allows to date all the trees in a dataset. For each replicate, one repository is created containing some files, including the dated trees named <i>X_tree_dated.nexus</i> where <i>X</i> is the number of the replicate. </p>
 <br>
 <h3>Check if a dataset positively correlates with epidemiological data</h3>
 <p>Before performing phylodynamics analysis, it is essential to confirm that the dataset is posively correlated with epidemiological data. In this work, we used the number of SARS-CoV-2-related deaths for each country. Epidemiological data for each country or French regions are provided in the Data/ repository. To perform the correlation, we developped R scripts, named <i><b>dataset_correlation_World.R</b></i>, <i><b>dataset_correlation_Europe.R</b></i> and <i><b>dataset_correlation_France.R</b></i> that first compare the cumulative number of deaths and sequences at the end of the time period studied, then directly compare per week the number of deaths and sequences included.
</p>
 <br>
 <h3>Count the number of distinct sequences across the replicates</h3>
 <p>Assuming we focus on one time period and one geographic scale (for example, worldwide scale, first half of 2020), we developed a short R script, named <i><b>count_distinct_sequences_in_set.R</b></i> that count the number of unique identifiers across the replicates. A sequence that appeared multiple times across the replicates is counted once.</p>
  <br>
   <h3>Compare the genetic distances between locations</h3>
 <p>Although optional, we developped R scripts, named <i><b>genetic_distances_World.R</b></i>, <i><b>genetic_distances_Europe.R</b></i> and <i><b>genetic_distances_France.R</b></i> to check whether some locations were associated with higher or lower genetic distances across the phylogenetic trees. Since we perform pairwise comparison of all tips of a same location, it is possible to perform the analysis only on a few replicates. As previously reported in some papers, we observed that all locations have almost systematically similar genetic distances.</p>
  <br>
     <h3>Calculating transmission events (Phylodynamic analysis)</h3>
 <p>This step allows to calculate both inter- and intra-territory transmission events across the replicates for a given time period and geographic scale. This step can be done using the R scripts <i><b>calcul_transmission_World.R</b></i>, <i><b>calcul_transmission_Europe.R</b></i> and <i><b>calcul_transmission_France.R</b></i>. To calculate transmissions, we first determine ancestral state at nodes using the ace function of the R ape package. We recommand the use of SYM or ER models, since the ARD model leads to some problems, possibly because of the outliers. From the ancestral reconstruction, we checked for each node of the phylogenetic tree whether children nodes corresponded to the same geographic region as the current node. If they correspond to different geographic locations, we then assumed the transmission from the parent node-state to the given child node-state (inter-territory transmission). The midpoint of the parent-child branch was chosen as the date of transmission in all cases.</p>
  <br>
       <h3>Check if intra-territory transmissions correlates with epidemiological data</h3>
 <p>To confirm that ancestral state reconstruction using the ace function of ape is efficient, we used the number of intra-territory transmission events as a proxy of the number of SARS-CoV-2-related deaths. XXXXXXXXXX</p>
  <br>
 <h3>5. Citation</h3>
 <p>If you use this pipeline for your own work, please cite:</p>
 <p><i>In preparation.</i></p>
