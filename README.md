# Scripts used for the study "Phylodynamic of SARS-CoV-2 transmission in France at local and international scales during 2020"
<p>xxx<br>
 xxx</p> <br>
<h3>1. Prerequisites</h3>
 <p>To use this pipeline, you need the following programs:</p>
 <p>- Perl 5<br>
 - Python 3.x<br>
 - Samtools 1.4</p>
 <br>
 <h3>2. Preparing a pileup file</h3>
 <p> The first step consists in the production of a pileup file using Samtools <i>mpileup</i> that calculates for each position across the genome the depth coverage, the quality and content of the reads:</p>
 <p><code> samtools mpileup -a -f reference_genome.fasta file_sorted.bam > file_sorted.pileup</code></p>
 <p>where reference_genome.fasta is the genome of interest in .fasta format (provided with the parameter <code>-f</code>, and file_sorted.bam is the sorted bam file previously produced with samtools <i> sort</i> function. The <code>-a</code> parameter allows to consider positions that were not covered.</p>
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
