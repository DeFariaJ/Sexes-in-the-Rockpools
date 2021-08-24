#On the fasta nucleotide alignments we run gblocks to remove poor alignment parts:
for i in *.nucl.fasta ; do   Gblocks $i -t=c -b4=21  ; done
-t protein (p), nuclotide (n) or codon (c)
-b4 minimum number of nucleotide in the block
#Remove space from gb files:
for i in *-gb ; do  sed 's/\s//g' $i >$i.fasta; done
#Convert fasta to phylip format:
for f in *.fasta; do python fasta2phylip_new.py $f $f.2 ; done
#Increase number of spaces after the species name:
for f in *.fasta.2; do sed 's/ /   /g' $f >$f.phy ; done
#Now I removed unnecessary files and renamed the files to aligned_nucl_OGXXX.phy:
for f in *:.fa.short-gb.fasta.2.phy; do  mv -- "$f" "${f%:.fa.short-gb.fasta.2.phy}.phy"; done
#Than using a dummy .ctl file I created an individual .ctl file for each aligned .phy:
for f in *.phy; do     var1="$f";     d="$f.ctl";     sed "s/z.phy/$var1/g" "codeml.ctl" > $d; done
#Than run Codeml on the cluster using a script:
conda activate species_phylogeny_venv
for f in *.ctl;
do
  var1="$f"; d="$f.mlc"; codeml $f > $d; mv mlc "$var1.mlc"; mv rst "$var1.rst"; mv rst1 "$var1.rst1"; mv rub "$var1.rub"; mv 2NG.dN "$var1.2NG.dN"; mv 2NG.dS "$var1.2NG.dS"; mv 2NG.t "$var1.2NG.t";
done
