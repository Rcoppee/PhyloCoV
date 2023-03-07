#! /bin/bash

# antoine.bridier-nahmias@inserm.fr
# Pairwise alignement of multiple genomes to the same reference
# in parallel

fasta_dir="./${3}/02-genomes/"
gisaid_sub=${fasta_dir}'gisaid_sub.fasta'
split_dir=${fasta_dir}'split_fasta/'

msa_dir="./${3}/03-msa/"
temp_msa=${msa_dir}'temp.msa'
temp_name=${msa_dir}'temp.name'
full_msa=${msa_dir}'full.msa'
ref_msa=${msa_dir}'ref.msa'
final_msa=${msa_dir}"${2}.msa"

ref='./data/Hu-1_OneLine.fasta'

outgroup='>EPI_ISL_402125'

mafft='/tools/mafft/mafft.bat'
seqtk='/tools/seqtk/seqtk'

n_threadus=14

### Split
mkdir ${split_dir}
gawk -v split_dir="${split_dir}"  \
  '/^>/ {close(OUTfile) ; OUTname = gensub(/^>/, "", 1, $0); 
  OUTname = gensub(/\//, "_", "G", OUTname); 
  OUTfile = split_dir OUTname".fasta"} {print > OUTfile} ' < ${1}
  
### Parallel mafft
mkdir ${msa_dir}
rm ${temp_msa}
find ${split_dir} -name '*.fasta' | \
  parallel -j ${n_threadus} --bar \
  "${mafft} --quiet --thread 1 --keeplength --add {} \
    ${ref} >> ${temp_msa}"

# One genome out of two is the ref, we get rid of it here with seqtk subseq 
# (a better solution could be found)

rm ${split_dir}/*.fasta
rm ${temp_name}
grep ">" ${temp_msa} | \
gawk -v ref=${outgroup} '{if (! match($0, "^" ref)) {print gensub(/^>/, "", "1",$0)} }' \
   >> ${temp_name}
${seqtk} subseq ${temp_msa} ${temp_name} >> ${full_msa}
cat ${ref} ${full_msa} > ${ref_msa}

# trim 100bp at each end
# ${seqtk} trimfq -b 100 -e 100 ${ref_msa} > ${final_msa} # not necessary when using masking list from deMaio et al.
cat ${ref_msa} > ${final_msa} # temp fix for masking
rm ${full_msa} ${ref_msa} ${temp_msa} ${temp_name}
