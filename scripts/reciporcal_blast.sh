# in the obiroi_genome folder under data/raw

blastn -query cds_from_genomic_GCA.fna -db genomeGCF_db -out blast_results_GCA2GCF.txt -outfmt 6
grep ">" cds_from_genomic_GCF.fna| tr -d ">" > genomeGCF_ids.txt
seqtk subseq cds_from_genomic_GCF.fna hit_ids.txt > hits_GCA2GCF.fasta
blastn -query hits_GCA2GCF.fasta -db genomeGCA_db -out reciprocal_blast_results_GCF2GCA.txt -outfmt 6