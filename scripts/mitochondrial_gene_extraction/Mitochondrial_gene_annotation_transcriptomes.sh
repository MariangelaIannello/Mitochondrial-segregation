#!/bin/bash 

### $1 = Nome specie
### $2 = Numero processori

### ...Previously: mitochondrial gene database building
## ncbi search string: "mitochondrial OR mitochondrion OR co1 OR co2 OR co3 OR cox OR cytochrome OR cytb OR nadh OR atp6 OR atp8 OR atp AND bivalvia NOT hypothetical NOT orf NOT HORF NOT h open reading frame NOT F-specific ORF NOT female-specific orf NOT MORF NOT F-ORF" => downloading file named ncbi_search_raw_results.fasta
## commands:
# sed -e '/partial/,+1d' ncbi_search_raw_result.fasta > ncbi_search_fullseqs.fasta
# sed -i 's/\[\([^]]*\)\]/\(\1)/g' ncbi_search_fullseqs.fasta
# IFS=$'\n'; for i in $(awk -F "\t" '{print $2}' names.txt); do grep -A1 -i ''"$i"' \[\|'"$i"','  ncbi_search_fullseqs.fasta; done | grep -v "\-\-" > ncbi_search_onlymito.fasta
# IFS=$'\n'; for i in $(awk -F "\t" '{print $2}' names.txt); do n=$(grep -w -i "$i"$ names.txt | awk -F "\t" '{print $1}'); for j in $(grep ">" ncbi_search_onlymito.fasta | grep -w -i $i); do s=$(echo $j | awk -F "(" '{print $NF}' | sed 's/)//g;s/ /_/g' ); echo -e `echo $j | awk -F " " '{print $1}'`"|""$n""|""$s"; grep -A1 "$j" ncbi_search_onlymito.fasta | tail -n1; done; done > ncbi_search_onlymito_simpnames.fasta
# diamond makedb --in ncbi_search_onlymito_simpnames.fasta -d ncbi_search_onlymito_simpnames.dmnd

## Define functions:
function 1linefasta() {	awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $1 | tail -n +2 ; }

########

cd "$1"

### First BLAST of filtered transcriptome against all bivalve mitochondrial genes
diamond blastx --quiet --query "$1"_trinity*fasta --db /home/PERSONALE/giovanni.piccinini5/DUI_project/Mito_annotation/ncbi_search_onlymito_simpnames.dmnd --ultra-sensitive --evalue 1e-3 --max-target-seqs 1 --threads "$2" --outfmt 6 --out "$1"_trinity.fasta_diamond_mito
for i in $(cut -f1 "$1"_trinity.fasta_diamond_mito); do grep -A1 "$i" "$1"_trinity*fasta; done > "$1"_mitohit.fasta

### Transdecoder ORFs extraction of positive hits (all orfs longer than 20 aa - not singelbestonly to account for polycistronic transcripts)
/usr/local/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -m 20 -G Mitochondrial-Invertebrates -t "$1"_mitohit.fasta
diamond blastp --quiet --query "$1"_mitohit.fasta.transdecoder_dir/longest_orfs.pep --db /home/PERSONALE/giovanni.piccinini5/DUI_project/Mito_annotation/ncbi_search_onlymito_simpnames.dmnd --ultra-sensitive --evalue 1e-3 --max-target-seqs 1 --threads "$2" --outfmt 6 --out blastp.outfmt6
/usr/local/TransDecoder-TransDecoder-v5.5.0/TransDecoder.Predict -t "$1"_mitohit.fasta -G Mitochondrial-Invertebrates --retain_blastp_hits blastp.outfmt6

# Directory cleaning
mkdir transdecoder
mv * transdecoder/
mv transdecoder/"$1"_mitohit.fasta.transdecoder.pep ./
mv transdecoder/"$1"_mitohit.fasta.transdecoder.cds ./
mv transdecoder/"$1"_mitohit.fasta ./
mv transdecoder/"$1"_trinity*fasta ./
mv transdecoder/"$1"_trinity*diamond_mito ./
1linefasta "$1"_mitohit.fasta.transdecoder.pep > bla; mv bla "$1"_mitohit.fasta.transdecoder.pep
1linefasta "$1"_mitohit.fasta.transdecoder.cds > bla; mv bla "$1"_mitohit.fasta.transdecoder.cds

### Keeping only orfs with positive hits against all bivalve mitochondrial genes
diamond blastp --quiet --query "$1"_mitohit.fasta.transdecoder.pep --db /home/PERSONALE/giovanni.piccinini5/DUI_project/Mito_annotation/ncbi_search_onlymito_simpnames.dmnd --evalue 1e-5 --ultra-sensitive --max-target-seqs 10 --threads "$2" --outfmt 6 qseqid sseqid bitscore evalue pident length --out "$1"_allorfs_diamond_mito.txt
for i in $(awk '{print $1}' "$1"_allorfs_diamond_mito.txt | sort -u); do grep -w -A1 "$i" "$1"_mitohit.fasta.transdecoder.pep; done > "$1"_allorfs_diamond_mito.fasta
for i in $(awk '{print $1}' "$1"_allorfs_diamond_mito.txt | sort -u); do grep -w -A1 "$i" "$1"_mitohit.fasta.transdecoder.cds; done > "$1"_allorfs_diamond_mito_nc.fasta

### Collapsing fragments (nucleotide sequences)
#~/CAP3/cap3 "$1"_allorfs_diamond_mito_nc.fasta
#1linefasta "$1"_allorfs_diamond_mito_nc.fasta.cap.contigs > bla; mv bla "$1"_allorfs_diamond_mito_nc.fasta.cap.contigs
#1linefasta "$1"_allorfs_diamond_mito_nc.fasta.cap.singlets > bla; mv bla "$1"_allorfs_diamond_mito_nc.fasta.cap.singlets
#cat "$1"_allorfs_diamond_mito_nc.fasta.cap.contigs "$1"_allorfs_diamond_mito_nc.fasta.cap.singlets > "$1"_allorfs_diamond_mito_nc_collapsed.fasta
#mkdir cap3
#mv "$1"_allorfs_diamond_mito_nc.fasta.cap.* cap3/
/usr/local/cdhit-master/cd-hit-est -i "$1"_allorfs_diamond_mito_nc.fasta -o "$1"_allorfs_diamond_mito_nc_collapsed.fasta -t 1 -g 1 -G 0 -aL 0.01 -aS 0.01
sed -i 's/ '"$1"'.*//g' "$1"_allorfs_diamond_mito_nc_collapsed.fasta


### Translating collapsed fragments
mkdir translatorx
perl /usr/local/translatorx_vLocal.pl -i "$1"_allorfs_diamond_mito_nc_collapsed.fasta -c 5 -o translatorx/"$1"_allorfs_diamond_mito_collapsed.fasta
mv translatorx/"$1"_allorfs_diamond_mito_collapsed.fasta.aaseqs.fasta ./"$1"_allorfs_diamond_mito_collapsed.fasta
sed -i 's/_TRINITY.*//g' ./"$1"_allorfs_diamond_mito_collapsed.fasta

### Final blast against all bivalve mitochondrial genes with "clean" names
if grep -q -i "$1" /home/PERSONALE/giovanni.piccinini5/DUI_project/Mito_annotation/ncbi_search_onlymito_simpnames.fasta
	then grep -A1 -i "$1" /home/PERSONALE/giovanni.piccinini5/DUI_project/Mito_annotation/ncbi_search_onlymito_simpnames.fasta | grep -v "\-\-" > "$1"_ncbi_mito.fasta
	diamond makedb --in "$1"_ncbi_mito.fasta -d "$1"_ncbi_mito.dmnd
	diamond blastp --quiet --query "$1"_allorfs_diamond_mito_collapsed.fasta --db "$1"_ncbi_mito.dmnd --evalue 1e-5 --id 90 --ultra-sensitive --max-target-seqs 10 --threads "$2" --outfmt 6 qseqid sseqid bitscore evalue pident length qstart qend sstart send --out "$1"_allorfs_collapsed_mito.txt
	else c=$(echo "$1" | awk -F_ '{print $1}')
	if grep -q -i "$c" /home/PERSONALE/giovanni.piccinini5/DUI_project/Mito_annotation/ncbi_search_onlymito_simpnames.fasta
		then grep -A1 -i "$c" /home/PERSONALE/giovanni.piccinini5/DUI_project/Mito_annotation/ncbi_search_onlymito_simpnames.fasta | grep -v "\-\-" > "$c"_ncbi_mito.fasta
		diamond makedb --in "$c"_ncbi_mito.fasta -d "$c"_ncbi_mito.dmnd
		diamond blastp --quiet --query "$1"_allorfs_diamond_mito_collapsed.fasta --db "$c"_ncbi_mito.dmnd --evalue 1e-5 --ultra-sensitive --max-target-seqs 10 --threads "$2" --outfmt 6 qseqid sseqid bitscore evalue pident length qstart qend sstart send --out "$1"_allorfs_collapsed_mito.txt
		else diamond blastp --quiet --query "$1"_allorfs_diamond_mito_collapsed.fasta --db /home/PERSONALE/giovanni.piccinini5/DUI_project/Mito_annotation/ncbi_search_onlymito_simpnames.dmnd --evalue 1e-5 --ultra-sensitive --max-target-seqs 10 --threads "$2" --outfmt 6 qseqid sseqid bitscore evalue pident length qstart qend sstart send --out "$1"_allorfs_collapsed_mito.txt
	fi
fi

sed -i 's/_'"$1"'.*[-+]\t/\t/g' "$1"_allorfs_collapsed_mito.txt

### Keeping the SINGLE BEST HIT for each mitochondrial gene (based on bistscore, column 3 of custom outfmt 6)
for i in $(sed 's/|/\t/g' "$1"_allorfs_collapsed_mito.txt | sort -k 3,3 -k 5,5nr | awk '!x[$3]++' | awk '{print $1}' | sed 's/_TRINITY.*//g;s/_'"$1"'.*//g')
	do echo -e ">""$1"_`sed 's/|/\t/g' "$1"_allorfs_collapsed_mito.txt | sort -k 3,3 -k 5,5nr | awk '!x[$3]++' | grep -P $i'(?![0-9])' | awk '{print $3}'`
	grep -A1 -P $i'(?![0-9])' "$1"_allorfs_diamond_mito_nc_collapsed.fasta | tail -n1
done > "$1"_mito_genes_nc.fasta

for i in $(sed 's/|/\t/g' "$1"_allorfs_collapsed_mito.txt | sort -k 3,3 -k 5,5nr | awk '!x[$3]++' | awk '{print $1}' | sed 's/_TRINITY.*//g;s/_'"$1"'.*//g')
	do echo -e ">""$1"_`sed 's/|/\t/g' "$1"_allorfs_collapsed_mito.txt | sort -k 3,3 -k 5,5nr | awk '!x[$3]++' | grep -P $i'(?![0-9])' | awk '{print $3}'`
	grep -A1 -P $i'(?![0-9])' "$1"_allorfs_diamond_mito_collapsed.fasta | tail -n1
done > "$1"_mito_genes_aa.fasta


cd ../
