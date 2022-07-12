# Scripts used for the study "Phylodynamic of SARS-CoV-2 transmission in France at local and international scales during 2020"
<p>xxx<br>
 xxx</p> <br>
<h3>1. Prerequisites</h3>
 <p>To use these scripts, you need the following programs:</p>
 <p>- TreeTime<br>
 - xxx<br>
 - xxx</p>
 <h3>Dataset generation (ANTOINE)</h3>
 <p>xxx</p>
 <p>xxx</p>
 <br>
 <h3>Dating the trees using TreeTime</h3>
 <p> We used TreeTime to date each phylogenetic tree previously produced. TreeTime requires a phylogenetic tree (newick format), a multiple sequence alignment and text file that indicates the date of sample collection for each tip in the phylogeny. Prior to TreeTime, we developed a R script, named <i>get_dates_from_tips.R</i> that produces such a file (one output per replicate). Then, we produced a shell script, named <i>launch_treetime.sh</i> that allows to date all the trees in a dataset. For each replicate, one repository is created containing some files, including the dated trees named <i>X_tree_dated.nexus</i> where <i>X</i> is the number of the replicate. </p>
 <br>
 <h3>3. Generating a table of A, T, G, C and indels content</h3>
 <p>The second step allows to formate a table that indicate the number of A, T, G, C and indels for each position of the genome. For that, we developed a perl script, named extract_data_pileup.pl, that takes the previously produced pileup file as an input (option <code>-p</code>).</p>
  <p><code> perl extract_data_pileup.pl -p file_sorted.pileup -o ./</code></p>
<p>The algorithm will produce a table, named data_from_pileup.tsv. The table is then simplified to retain only the chrosomosomes, the genomic positions and the major bases (that we named here simplified_data.tsv), using the <i>cut</i> function in Linux. Importantly, we used the IUPAC code to produce the major allele, so other letters than A, T, G and C may appear as the result of mixed alleles.</p>
 <p><code>cut data_from_pileup.tsv -f 1,2,8 > simplified_data.tsv</code></p>
 <p>where data_from_pileup.tsv is the file previously produced, and the <code>-f</code> option allows to retain columns 1, 2 and 8 of the file.</p>
 <br>
 <h3>4. Producing the .fasta sequence</h3>
 <p>From the simplified simplified_data.tsv table, we then just concatenated positions to form a fasta sequence. At the end of the program, the total length of the sequence and the number of lacking positions and nucleotides are indicated.</p>
  <p><code>python3 extract_fasta_pileup.py -f simplified_data.tsv -o out.fasta</code></p>
  <p>where <code>-f</code> is the input file (simplified_data.tsv) and <code>-o</code> corresponds to the output fasta file (out.fasta).</p>
  <br>
 <h3>5. Citation</h3>
 <p>If you use this pipeline for your own work, please cite:</p>
 <p><i>In preparation.</i></p>
