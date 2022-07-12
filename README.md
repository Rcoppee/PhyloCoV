# Scripts used for the study "Phylodynamic of SARS-CoV-2 transmission in France at local and international scales during 2020"
<p>xxx<br>
 xxx</p> <br>
<h3>1. Prerequisites</h3>
 <p>To use these scripts, you need the following programs:</p>
 <p>- TreeTime<br>
 - R<br>
 - R packages</p>
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
 <h3>4. Producing the .fasta sequence</h3>
 <p>From the simplified simplified_data.tsv table, we then just concatenated positions to form a fasta sequence. At the end of the program, the total length of the sequence and the number of lacking positions and nucleotides are indicated.</p>
  <p><code>python3 extract_fasta_pileup.py -f simplified_data.tsv -o out.fasta</code></p>
  <p>where <code>-f</code> is the input file (simplified_data.tsv) and <code>-o</code> corresponds to the output fasta file (out.fasta).</p>
  <br>
 <h3>5. Citation</h3>
 <p>If you use this pipeline for your own work, please cite:</p>
 <p><i>In preparation.</i></p>
