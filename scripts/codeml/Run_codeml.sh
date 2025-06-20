#!/bin/bash

# Create directories for constrained OGs
for i in $(ls /home/PERSONALE/mariangela.iannello2/DUI_project/00_MIX_GENOMES_TRANSCRIPTOME/cds_nt_per_analisi_selezione/06BIS_RELAX_TERMINAL_BRANCHES_NOGTREE/relax/nt_aln/ | grep OG); do mkdir "${i/.fa/}"; cp /home/PERSONALE/mariangela.iannello2/DUI_project/00_MIX_GENOMES_TRANSCRIPTOME/cds_nt_per_analisi_selezione/06BIS_RELAX_TERMINAL_BRANCHES_NOGTREE/relax/nt_aln/"$i" "${i/.fa/}"/; cp /home/PERSONALE/mariangela.iannello2/DUI_project/00_MIX_GENOMES_TRANSCRIPTOME/cds_nt_per_analisi_selezione/05_GENE_TREES/snakemake/model_finder/raxml/"${i/.fa/}".raxml.bestTree "${i/.fa/}"/"${i/.fa/}"_blen.nwk; done

for i in $(grep accelerated /home/PERSONALE/mariangela.iannello2/DUI_project/00_MIX_GENOMES_TRANSCRIPTOME/06_traccer/traccer_MixGT/traccer_MixGT.TRACCER.txt | awk '{print $1}'); do rm -r "$i"; done


# Prepare files for the CODEML run (test trees, alternative trees, phylip alignments, control files)
for i in $(ls | grep OG); do java -jar /usr/local/phyutility/phyutility.jar -pr -in sp_tree_no_bl.nwk -out "$i"/"$i"_backbone_tagged.nwk -names  $(for s in $(cut -f1 species_list_tags.txt); do if grep -q "${s/\#1/}" "$i"/"$i"_blen.nwk; then : ; else echo "$s"; fi ; done | sep " "); if [[ $(grep -oh '\w[a-z]\w*' "$i"/"$i"_blen.nwk | wc -l) == $(wc -l species_list.txt | awk -F " " '{print $1}') ]] ; then cp sp_tree_no_bl.nwk "$i"/"$i"_backbone_tagged.nwk; fi ;  done

for i in $(ls | grep OG); do cd "$i"; ../transfer_blen.R "$i"; cd ../ ; done
for i in $(ls | grep OG); do sed -i 's/\#1\(:[0-9]\.[0-9]\+\)/\1 \#1/g' "$i"/"$i"_alternative.nwk; done


for i in $(ls | grep OG); do java -jar /usr/local/BMGE-1.12/BMGE.jar -i "$i"/"$i".fa -t DNA -h 1 -g 1 -o "$i"/"$i".phy; done
for i in $(ls | grep OG); do sed -i 's/ /  /g' "$i"/"$i".phy; done
for i in $(ls | grep OG); do cp codeml_base.ctl "$i"/"$i"_null.ctl; sed -i 's/INPUTSEQ/'"$i"'.phy/g;s/OUTPUT/'"$i"'_codeml_null.txt/g;s/INPUTTREE/'"$i"'_null.nwk/g;s/MODEL/0/g' "$i"/"$i"_null.ctl; done
for i in $(ls | grep OG); do cp codeml_base.ctl "$i"/"$i"_alternative.ctl; sed -i 's/INPUTSEQ/'"$i"'.phy/g;s/OUTPUT/'"$i"'_codeml_alternative.txt/g;s/INPUTTREE/'"$i"'_alternative.nwk/g;s/MODEL/2/g' "$i"/"$i"_alternative.ctl; done
for i in $(ls | grep OG); do mkdir "$i"/codeml_files_null; mkdir "$i"/codeml_files_alternative; done



# Run CODEML
for i in $(ls | grep OG); do cd "$i"; yes \n | codeml "$i"_null.ctl; mv 2N* codeml_files_null; mv 4fold.nuc codeml_files_null; mv lnf codeml_files_null; mv r* codeml_files_null; yes \n | codeml "$i"_alternative.ctl; mv 2N* codeml_files_alternative; mv 4fold.nuc codeml_files_alternative; mv lnf codeml_files_alternative; mv r* codeml_files_alternative; cd ../; done

# Create final summary file
for i in $(ls | grep OG ); do Ln=$(grep lnL "$i"/"$i"_codeml_null.txt | awk -F " " '{print $5}'); La=$(grep lnL "$i"/"$i"_codeml_alternative.txt | awk -F " " '{print $5}'); Dw=$(grep "w (dN/dS) for branches" "$i"/"$i"_codeml_alternative.txt  | awk -F " " '{print $7-$6}'); echo -e $i"\t"$Ln"\t"$La"\t"$Dw; done > Results_summary.tsv

#Perform likelihood ratio tests
./likelihood_ratio_test.R
