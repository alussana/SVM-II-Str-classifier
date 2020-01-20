SHELL := /bin/bash

AMINO_ACIDS = G A V P L I M F W Y S T C N Q H D E K R
SCOP_CLASSES = a b c d e f g h i j k l

##############################
### Preliminary Statistics ###
##############################

## Count the total number of residues appearing in dssp files
.PHONY: raw_training_count_residues
raw_training_count_residues: 
	for file in $$(ls training_set/dssp); do cat training_set/dssp/$$file | sed '1d; s/./\0\n/g'; done \
	| sed '/^[[:space:]]*$$/d;s/[[:space:]]*$$//' | wc -l

.PHONY: filtered_training_count_residues
filtered_training_count_residues:
	for file in $$(ls training_set/dssp/ | sed 's/.dssp$$//' | grep -f <(cat training_set/training_ids.txt)); do \
		cat training_set/dssp/$${file}.dssp | sed '1d; s/./\0\n/g'; \
	done \
	| sed '/^[[:space:]]*$$/d;s/[[:space:]]*$$//' | wc -l

## .META:
## 1	pdb_identifier
## 2	organism
training_set/pdb2organism: refs/jpred4.pdblist.txt
	pdb2organism $< > $@

DSSP_PIE_SIZE = 10
## .META:
## 1 dssp_type
## 2 residue_count
## 3 residue_percentage
training_set/raw_residue_by_dssp_count.txt:
	for file in $$(ls training_set/dssp); do \
		cat training_set/dssp/$$file | sed '1d; s/./\0\n/g'; done \
	| sed '/^[[:space:]]*$$/d;s/[[:space:]]*$$//' \
	| awk 'BEGIN{H=0;E=0;c=0}; {\
	if($$0=="H"){H=H+1} else if($$0=="E"){E=E+1} \
	else if($$0=="-"){c=c+1}}; END{tot=H+E+c; pH=H*100/tot;\
	pE=E*100/tot; pc=c*100/tot; \
	printf("H\t%i\t%f\nE\t%i\t%f\nC\t%i\t%f\n",H,pH,E,pE,c,pc)}'\
	> $@; \
	cut -f1 training_set/raw_residue_by_dssp_count.txt > training_set/raw_labels.tmp; \
	cut -f3 training_set/raw_residue_by_dssp_count.txt > training_set/raw_data.tmp; \
	pie_chart training_set/raw_data.tmp training_set/raw_labels.tmp $(DSSP_PIE_SIZE) training_set/raw_residue_by_dssp_count.png; \
	rm training_set/raw_labels.tmp training_set/raw_data.tmp

training_set/filtered_residue_by_dssp_count.txt:
	for file in $$(ls training_set/dssp/ | sed 's/.dssp$$//' | grep -f <(cat training_set/training_ids.txt)); do \
		cat training_set/dssp/$${file}.dssp | sed '1d; s/./\0\n/g'; done \
	| sed '/^[[:space:]]*$$/d;s/[[:space:]]*$$//' \
	| awk 'BEGIN{H=0;E=0;c=0}; {\
	if($$0=="H"){H=H+1} else if($$0=="E"){E=E+1} \
	else if($$0=="-"){c=c+1}}; END{tot=H+E+c; pH=H*100/tot;\
	pE=E*100/tot; pc=c*100/tot; \
	printf("H\t%i\t%f\nE\t%i\t%f\nC\t%i\t%f\n",H,pH,E,pE,c,pc)}'\
	> $@; \
	cut -f1 training_set/filtered_residue_by_dssp_count.txt > training_set/filtered_labels.tmp; \
	cut -f3 training_set/filtered_residue_by_dssp_count.txt > training_set/filtered_data.tmp; \
	pie_chart training_set/filtered_data.tmp training_set/filtered_labels.tmp $(DSSP_PIE_SIZE) training_set/filtered_residue_by_dssp_count.png; \
	rm training_set/filtered_labels.tmp training_set/filtered_data.tmp

## .META:
## 1	scop_domain_id
## 2	amino_acid_sequence
training_set/merged_fasta_files.gz:
	for file in $$(ls training_set/fasta/); do \
		cat training_set/fasta/$$file | sed ':a;N;$$!ba;s/\n/\t/1' | sed 's/>//1'; \
	done | sort | gzip > $@

## .META:
## 1    scop_domain_id
## 2    dssp_sequence
training_set/merged_dssp_files.gz:
	for file in $$(ls training_set/dssp/); do \
		cat training_set/dssp/$$file | sed ':a;N;$$!ba;s/\n/\t/1' | sed 's/>//1'; \
	done | sort | gzip > $@

## META:
## 1	amino_acid
## 2	H_count
## 3	E_count
## 4	-_count
training_set/raw_dssp_by_residue_count.txt: training_set/merged_fasta_files.gz training_set/merged_dssp_files.gz refs/jpred4.list.txt
	for aa in $(AMINO_ACIDS); do \
		H=0; E=0; c=0; \
		for protein in $$(cat $$(echo $^ | cut -d ' ' -f3)); do \
			zcat $$(echo $^ | cut -d ' ' -f1) | fgrep $$protein | cut -f2 \
				| tr -d "\n" | sed 's/./&\n/g' > training_set/raw_fasta.tmp;\
			zcat $$(echo $^ | cut -d ' ' -f2) | fgrep $$protein | cut -f2 \
				| tr -d "\n" | sed 's/./&\n/g' > training_set/raw_dssp.tmp; \
			paste training_set/raw_fasta.tmp training_set/raw_dssp.tmp | awk -v aa=$$aa '$$1==aa{print $$2}' \
				 > training_set/raw_subset_dssp.tmp; \
			rm training_set/raw_fasta.tmp training_set/raw_dssp.tmp; \
			for i in $$(cat training_set/raw_subset_dssp.tmp); do \
				if [ $$i = 'H' ]; then H=$$(echo $$H+1 | bc); fi; \
				if [ $$i = 'E' ]; then E=$$(echo $$E+1 | bc); fi; \
				if [ $$i = '-' ]; then c=$$(echo $$c+1 | bc); fi; \
			done; \
			rm training_set/raw_subset_dssp.tmp; \
		done; \
		echo -e "$$aa\t$$H\t$$E\t$$c"; \
	done > $@;

training_set/raw_dssp_by_residue_count.png: training_set/raw_dssp_by_residue_count.txt
	hist_dssp_by_residue $< $@

training_set/filtered_dssp_by_residue_count.txt: training_set/merged_fasta_files.gz training_set/merged_dssp_files.gz training_set/training_ids.txt
	for aa in $(AMINO_ACIDS); do \
		H=0; E=0; c=0; \
		for protein in $$(cat $$(echo $^ | cut -d ' ' -f3)); do \
			zcat $$(echo $^ | cut -d ' ' -f1) | fgrep $$protein | cut -f2 \
				| tr -d "\n" | sed 's/./&\n/g' > training_set/filtered_fasta.tmp;\
			zcat $$(echo $^ | cut -d ' ' -f2) | fgrep $$protein | cut -f2 \
				| tr -d "\n" | sed 's/./&\n/g' > training_set/filtered_dssp.tmp; \
			paste training_set/filtered_fasta.tmp training_set/filtered_dssp.tmp | awk -v aa=$$aa '$$1==aa{print $$2}' \
				> training_set/filtered_subset_dssp.tmp; \
			rm training_set/filtered_fasta.tmp training_set/filtered_dssp.tmp; \
			for i in $$(cat training_set/filtered_subset_dssp.tmp); do \
				if [ $$i = 'H' ]; then H=$$(echo $$H+1 | bc); fi; \
				if [ $$i = 'E' ]; then E=$$(echo $$E+1 | bc); fi; \
				if [ $$i = '-' ]; then c=$$(echo $$c+1 | bc); fi; \
			done; \
			rm training_set/filtered_subset_dssp.tmp; \
		done; \
		echo -e "$$aa\t$$H\t$$E\t$$c"; \
	done > $@;

training_set/filtered_dssp_by_residue_count.png: training_set/filtered_dssp_by_residue_count.txt
	hist_dssp_by_residue $< $@

SCOP_PIE_SIZE = 10
## .META:
## 1    scop_class
## 2    pdb_count
training_set/raw_domain_by_scop_count.txt: refs/dir.cla.scope.2.06-stable.txt refs/jpred4.list.txt
	cat $$(echo $^ | cut -d ' ' -f1) | cut -f1,4 | sort | uniq \
	| fgrep -f $$(echo $^ | cut -d ' ' -f2) > training_set/raw_dir.cla.scope.filtered.tmp; \
	for class in $(SCOP_CLASSES); do \
		count=$$(cat training_set/raw_dir.cla.scope.filtered.tmp | cut -f2 | cut -d '.' -f1 | fgrep $$class | wc -l); \
		echo -e "$$class\t$$count"; \
	done > $@; \
	rm training_set/raw_dir.cla.scope.filtered.tmp; \
	cat training_set/raw_domain_by_scop_count.txt | sed '/0$$/d' > training_set/raw_domain_by_scop_count_nozero.tmp; \
	cut -f1 training_set/raw_domain_by_scop_count_nozero.tmp > training_set/raw_labels.tmp; \
	cut -f2 training_set/raw_domain_by_scop_count_nozero.tmp > training_set/raw_data.tmp; \
	pie_chart training_set/raw_data.tmp training_set/raw_labels.tmp $(SCOP_PIE_SIZE) training_set/raw_domain_by_scop_count.png; \
	rm training_set/raw_data.tmp training_set/raw_labels.tmp training_set/raw_domain_by_scop_count_nozero.tmp

training_set/filtered_domain_by_scop_count.txt: refs/dir.cla.scope.2.06-stable.txt training_set/training_ids.txt
	cat $$(echo $^ | cut -d ' ' -f1) | cut -f1,4 | sort | uniq \
	| fgrep -f $$(echo $^ | cut -d ' ' -f2) > training_set/filtered_dir.cla.scope.filtered.tmp; \
	for class in $(SCOP_CLASSES); do \
		count=$$(cat training_set/filtered_dir.cla.scope.filtered.tmp | cut -f2 | cut -d '.' -f1 | fgrep $$class | wc -l); \
		echo -e "$$class\t$$count"; \
	done > $@; \
	rm training_set/filtered_dir.cla.scope.filtered.tmp; \
	cat training_set/filtered_domain_by_scop_count.txt | sed '/0$$/d' > training_set/filtered_domain_by_scop_count_nozero.tmp; \
	cut -f1 training_set/filtered_domain_by_scop_count_nozero.tmp > training_set/filtered_labels.tmp; \
	cut -f2 training_set/filtered_domain_by_scop_count_nozero.tmp > training_set/filtered_data.tmp; \
	pie_chart training_set/filtered_data.tmp training_set/filtered_labels.tmp $(SCOP_PIE_SIZE) training_set/filtered_domain_by_scop_count.png; \
	rm training_set/filtered_data.tmp training_set/filtered_labels.tmp training_set/filtered_domain_by_scop_count_nozero.tmp

training_set/raw_propensity_heatmap_helix.png: training_set/merged_fasta_files.gz training_set/merged_dssp_files.gz
	propensity_heatmap $$(echo $^ | cut -d ' ' -f1) $$(echo $^ | cut -d ' ' -f2) $@

training_set/raw_propensity_heatmap_strand.png: training_set/merged_fasta_files.gz training_set/merged_dssp_files.gz
	propensity_heatmap $$(echo $^ | cut -d ' ' -f1) $$(echo $^ | cut -d ' ' -f2) $@

training_set/filtered_merged_fasta_files.gz: training_set/merged_fasta_files.gz training_set/training_ids.txt
	zcat $< | grep -f <(cat training_set/training_ids.txt) | gzip > $@

training_set/filtered_merged_dssp_files.gz: training_set/merged_dssp_files.gz training_set/training_ids.txt
	zcat $< | grep -f <(cat training_set/training_ids.txt) | gzip > $@

training_set/filtered_propensity_heatmap_helix.png: training_set/filtered_merged_fasta_files.gz training_set/filtered_merged_dssp_files.gz
	propensity_heatmap $$(echo $^ | cut -d ' ' -f1) $$(echo $^ | cut -d ' ' -f2) $@

training_set/filtered_propensity_heatmap_strand.png: training_set/filtered_merged_fasta_files.gz training_set/filtered_merged_dssp_files.gz
	propensity_heatmap $$(echo $^ | cut -d ' ' -f1) $$(echo $^ | cut -d ' ' -f2) $@

#######################
### Handmade Charts ###
#######################

training_set/raw_training_taxonomy.png: ../local/handmade/raw_training_taxonomy_handmade.txt
	cat $< | cut -f1 > raw_training_sizes.tmp; \
	cat $< | cut -f2 > raw_training_labels.tmp; \
	pie_chart raw_training_sizes.tmp raw_training_labels.tmp $(SCOP_PIE_SIZE) $@; \
	rm raw_training_sizes.tmp raw_training_labels.tmp

training_set/raw_training_organism.png: ../local/handmade/raw_training_organism_handmade.txt
	cat $< | cut -f1 > raw_training_sizes.tmp; \
        cat $< | cut -f2 > raw_training_labels.tmp; \
	pie_chart raw_training_sizes.tmp raw_training_labels.tmp $(SCOP_PIE_SIZE) $@; \
	rm raw_training_sizes.tmp raw_training_labels.tmp

training_set/filtered_training_taxonomy.png: ../local/handmade/filtered_training_taxonomy_handmade.txt
	cat $< | cut -f1 > filtered_training_sizes.tmp; \
	cat $< | cut -f2 > filtered_training_labels.tmp; \
	pie_chart filtered_training_sizes.tmp filtered_training_labels.tmp $(SCOP_PIE_SIZE) $@; \
	rm filtered_training_sizes.tmp filtered_training_labels.tmp

training_set/filtered_training_organism.png: ../local/handmade/filtered_training_organism_handmade.txt
	cat $< | cut -f1 > filtered_training_sizes.tmp; \
        cat $< | cut -f2 > filtered_training_labels.tmp; \
	pie_chart filtered_training_sizes.tmp filtered_training_labels.tmp $(SCOP_PIE_SIZE) $@; \
	rm filtered_training_sizes.tmp filtered_training_labels.tmp

blind_test_set/blindtest_taxonomy.png: ../local/handmade/blindtest_taxonomy_handmade.txt
	cat $< | cut -f1 > blindtest_sizes.tmp; \
	cat $< | cut -f2 > blindtest_labels.tmp; \
	pie_chart blindtest_sizes.tmp blindtest_labels.tmp $(SCOP_PIE_SIZE) $@; \
	rm blindtest_sizes.tmp blindtest_labels.tmp

blind_test_set/blindtest_organism.png: ../local/handmade/blindtest_organism_handmade.txt
	cat $< | cut -f1 > blindtest_sizes.tmp; \
        cat $< | cut -f2 > blindtest_labels.tmp; \
	pie_chart blindtest_sizes.tmp blindtest_labels.tmp $(SCOP_PIE_SIZE) $@; \
	rm blindtest_sizes.tmp blindtest_labels.tmp

#############################
### Select Blind Test Set ###
#############################

# Percent Sequence Alignment Search : PDB can contain Expression Tag sequence = No , and 
# Released between 2015-01-01 and 2019-09-12 and 
# Chain Type: there is a Protein chain and 
# Resolution is between 0.0 and 2.5 and 
# Sequence Length is between 50 and 300 and 
# Representative Structures at 30% Sequence Identity
blind_test_set/fasta.txt.gz: blind_test_set/tabularResults.csv 
	cat $< | grep "Protein" | tr -d '"' | tr "," "\t" | cut -f1,2,6 | sed 's/\t/:/1' | awk '{print ">"$$0}' | sed 's/\t/\n/1' | gzip > $@

## .META:
## 1	pdb_id:chain_id
## 2 	x-ray_resolution
## 3	chain_length
blind_test_set/attributes.gz: blind_test_set/tabularResults.csv
	cat $< | sed '1d' | tr -d '"' | tr "," "\t" | cut -f1,2,5,10 | sed 's/\t/:/1' | gzip > $@	

## 1. filter chains length
MIN_LENGTH = 50
blind_test_set/blind_test.fasta.gz: blind_test_set/fasta.txt.gz
	fasta_filter_length $< $(MIN_LENGTH) | gzip > $@

## 2. clusterize to detect sequence similarity
LENGTH = 0.0
SEQUENCE_IDENTITY = 30
blind_test_set/blind_test_clust.gz: blind_test_set/blind_test.fasta.gz
	zcat $< > blind_test_set/blind_test.fasta.tmp; \
	blastclust -i blind_test_set/blind_test.fasta.tmp -L $(LENGTH) -S $(SEQUENCE_IDENTITY) | gzip > $@; \
	rm blind_test_set/blind_test.fasta.tmp

## 3. choose representatives according to some criteria (at the moment, the first one listed)
blind_test_set/blind_test_clust_representatives.gz: blind_test_set/blind_test_clust.gz
	zcat $< | grep -v "of" | cut -d ' ' -f1 | awk '{print ">"$$0}' | gzip > $@

blind_test_set/blind_test_clust_representatives.fasta: blind_test_set/blind_test_clust_representatives.gz blind_test_set/fasta.txt.gz
	zcat $$(echo $^ | cut -d ' ' -f2) | grep -A 1 -f <(zcat $$(echo $^ | cut -d ' ' -f1)) | sed '/-/d' > $@

## 4. index the blind set database
blind_test_set/blind_test_clust_representatives.fasta.phr: blind_test_set/blind_test_clust_representatives.fasta
	makeblastdb -in $< -out $< -dbtype prot

BLASTP_E_VALUE = 1
## 5. blastp against the JPred4 dataset
blind_test_set/blind_representatives_blastp_vs_jpred.tab: blind_test_set/blind_test_clust_representatives.fasta training_set/training.fasta
	blastp -query blind_test_set/blind_test_clust_representatives.fasta -db training_set/training.fasta -evalue $(BLASTP_E_VALUE) -out $@ -outfmt 6

# 6. Filter out for >30% s.i.
blind_test_set/blastp_negative_filter: blind_test_set/blind_representatives_blastp_vs_jpred.tab
	cut -f1,3 $< | awk '$$2>30 {print $$1}'| sort | uniq | awk '{print ">"$$0}' > $@

blind_test_set/blind_refined.fasta: blind_test_set/blind_test_clust_representatives.fasta blind_test_set/blastp_negative_filter
	cat blind_test_set/blind_test_clust_representatives.fasta | tr '\n' '\t' | sed 's/\t>/\n>/g' \
	| grep -w -v -f blind_test_set/blastp_negative_filter | sed 's/\t/\n/1' > $@

## 7. compute dssp on the set
blind_test_set/blind_refined.dssp.gz: blind_test_set/blind_refined.fasta
	for chain in $$(cat $< | grep ">" | tr -d '>'); do \
		id=$$(echo $$chain | cut -d ':' -f1); ch=$$(echo $$chain | cut -d ':' -f2); \
		wget -P blind_test_set/ http://www.rcsb.org/pdb/files/$$id.pdb; \
		cat blind_test_set/$$id.pdb | awk -v cha=$$ch '{split($$0,line); \
		if(line[1]=="ATOM" && line[5]==cha) {print  $$0}}' > blind_test_set/$$id:$$ch.pdb; \
		mkdssp -i blind_test_set/$$id:$$ch.pdb -o blind_test_set/$$id:$$ch.dssp.tmp; \
		cat <(echo ">$$id:$$ch") <(cat blind_test_set/$$id:$$ch.dssp.tmp) >> blind_test_set/blind_refined.dssp; \
		rm blind_test_set/$$id:$$ch.pdb blind_test_set/$$id:$$ch.dssp.tmp blind_test_set/$$id.pdb; \
	done; gzip blind_test_set/blind_refined.dssp > $@; \

## this also generates blind_test_set/test_set.dssp.gz:
blind_test_set/test_set.fasta.gz: blind_test_set/blind_refined.dssp.gz
	get_dssp $< blind_test_set/test_set.fasta  blind_test_set/test_set.dssp; \
	gzip blind_test_set/test_set.fasta > $@; gzip blind_test_set/test_set.dssp > blind_test_set/test_set.dssp.gz

## 8. filter for dssp length
## blind_test_set/blind_refined_filtered.dssp
blind_test_set/blind_refined_filtered.fasta: blind_test_set/test_set.fasta.gz
	fasta_filter_length $< 50 > $@.tmp; \
	fasta_filter_length blind_test_set/test_set.dssp.gz 50 > blind_test_set/blind_refined_filtered.dssp.tmp; \
	cat $@.tmp | tr '\n' '\t' | sed 's/\t>/\n>/g' | sed '/!/d' | sed 's/\t/\n/1' > $@;
	cat blind_test_set/blind_refined_filtered.dssp.tmp | tr '\n' '\t' | sed 's/\t>/\n>/g' \
	| grep -f <(cat $@ | tr '\n' '\t' | sed 's/\t>/\n>/g' | cut -f1) | sed 's/\t/\n/1' > blind_test_set/blind_refined_filtered.dssp; \
	rm blind_test_set/blind_refined_filtered.dssp.tmp $@.tmp

## 9. generate sequence profiles
blind_test_set/testing_set_ids.txt: blind_test_set/blind_refined_filtered.fasta
	cat $< | grep ">" | sed 's/>//' > $@; \
	for id in $$(cat $@); do \
		cat $< | grep -A 1 $$id > blind_test_set/fasta/$${id}.fasta; \
	done

.PHONY: get_pssm_from_testing_set
get_pssm_from_testing_set: blind_test_set/testing_set_ids.txt
	for blind_id in $$(cat $<); do \
		psiblast -query blind_test_set/fasta/$${blind_id}.fasta -db refs/uniprot_sprot.fasta -evalue $(PSIBLAST_E_VALUE) \
		-out_ascii_pssm blind_test_set/pssm/$${blind_id}.pssm -out blind_test_set/blast/$${blind_id}.blast \
		-num_descriptions 10000 -num_alignments 10000 -num_iterations 3; \
	done;
	rm blind_test_set/pssm/*%*

## 10. filter for successfully generated pssm
## W: needs get_pssm_from_testing_set to be executed first
blind_test_set/testing_set_pssm_pass_ids.txt:
	ls -l blind_test_set/pssm/ | grep pssm | sed -r 's/.*(....:.).pssm$$/\1/' > $@

## 11. filter for no all-zero pssm
blind_test_set/pssm_pass_ids.txt: blind_test_set/testing_set_pssm_pass_ids.txt
	check_pssm --pssm_dir blind_test_set/pssm/ --id_list $< > all_zero_list.tmp;
	cat $< | grep -v -f all_zero_list.tmp > $@;
	rm all_zero_list.tmp

## 12. random sample 150 chains
blind_test_set/blind_selected_150.fasta: blind_test_set/pssm_pass_ids.txt blind_test_set/blind_refined_filtered.fasta
	cat blind_test_set/blind_refined_filtered.fasta | tr '\n' '\t' | sed 's/\t>/\n>/g' | grep -f $< \
	| shuf -n 150 | sed 's/\t/\n/1' > blind_test_set/blind_refined_150.fasta; \
	cat blind_test_set/blind_refined_filtered.dssp | tr '\n' '\t' | sed 's/\t>/\n>/g'\
	| grep -f <(cat blind_test_set/blind_refined_150.fasta | grep ">") | sed 's/\t/\n/1' > blind_test_set/blind_refined_150.dssp; \
	cat blind_test_set/blind_refined_150.fasta | tr '\n' '\t' | sed 's/\t>/\n>/g' | sort -k1,1 | sed 's/\t/\n/1' > $@; \
	cat blind_test_set/blind_refined_150.dssp | tr '\n' '\t' | sed 's/\t>/\n>/g' | sort -k1,1 | sed 's/\t/\n/1' > blind_test_set/blind_selected_150.dssp; \
	rm blind_test_set/blind_refined_150.fasta blind_test_set/blind_refined_150.dssp

## 13. generate single-sequence dssp and fasta files for selected ids
write_blind_test_set_dssp: blind_test_set/blind_selected_150.fasta
	for prot_id in $$(cat blind_test_set/blind_selected_150.dssp | grep ">" | tr -d ">"); do \
		cat blind_test_set/blind_selected_150.fasta | grep -A 1 $${prot_id} > blind_test_set/selected_fasta/$${prot_id}.fasta; \
		cat blind_test_set/blind_selected_150.dssp| grep -A 1 $${prot_id} > blind_test_set/selected_dssp/$${prot_id}.dssp; \
		cp blind_test_set/pssm/$${prot_id}.pssm blind_test_set/selected_pssm/; \
	done

## 14: generate final selected id list
blind_test_set/selected_ids: blind_test_set/blind_selected_150.fasta
	ls blind_test_set/selected_dssp/ | sed 's/.....$$//g' > $@

#################################
### Blind Test Set Statistics ###
#################################

## .META:
## 1 dssp_type
## 2 residue_count
## 3 residue_percentage
blind_test_set/residue_by_dssp_count.txt:
	for file in $$(ls blind_test_set/selected_dssp); do \
		cat blind_test_set/selected_dssp/$$file | sed '1d; s/./\0\n/g'; done \
	| sed '/^[[:space:]]*$$/d;s/[[:space:]]*$$//' \
	| awk 'BEGIN{H=0;E=0;c=0}; {\
	if($$0=="H"){H=H+1} else if($$0=="E"){E=E+1} \
	else if($$0=="-"){c=c+1}}; END{tot=H+E+c; pH=H*100/tot;\
	pE=E*100/tot; pc=c*100/tot; \
	printf("H\t%i\t%f\nE\t%i\t%f\nC\t%i\t%f\n",H,pH,E,pE,c,pc)}'\
	> $@; \
	cut -f1 blind_test_set/residue_by_dssp_count.txt > blind_test_set/labels.tmp; \
	cut -f3 blind_test_set/residue_by_dssp_count.txt > blind_test_set/data.tmp; \
	pie_chart blind_test_set/data.tmp blind_test_set/labels.tmp $(DSSP_PIE_SIZE) blind_test_set/residue_by_dssp_count.png; \
	rm blind_test_set/labels.tmp blind_test_set/data.tmp

## .META:
## 1	pdb_id
## 2	amino_acid_sequence
blind_test_set/merged_fasta_files.gz:
	for file in $$(ls blind_test_set/selected_fasta/); do \
		cat blind_test_set/selected_fasta/$$file | sed ':a;N;$$!ba;s/\n/\t/1' | sed 's/>//1'; \
	done | sort | gzip > $@

## .META:
## 1    pdb_id
## 2    dssp_sequence
blind_test_set/merged_dssp_files.gz:
	for file in $$(ls blind_test_set/selected_dssp/); do \
		cat blind_test_set/selected_dssp/$$file | sed ':a;N;$$!ba;s/\n/\t/1' | sed 's/>//1'; \
	done | sort | gzip > $@

## META:
## 1	amino_acid
## 2	H_count
## 3	E_count
## 4	-_count
blind_test_set/dssp_by_residue_count.txt: blind_test_set/merged_fasta_files.gz blind_test_set/merged_dssp_files.gz blind_test_set/selected_ids
	for aa in $(AMINO_ACIDS); do \
		H=0; E=0; c=0; \
		for protein in $$(cat $$(echo $^ | cut -d ' ' -f3)); do \
			zcat $$(echo $^ | cut -d ' ' -f1) | fgrep $$protein | cut -f2 \
				| tr -d "\n" | sed 's/./&\n/g' > blind_test_set/fasta.tmp;\
			zcat $$(echo $^ | cut -d ' ' -f2) | fgrep $$protein | cut -f2 \
				| tr -d "\n" | sed 's/./&\n/g' > blind_test_set/dssp.tmp; \
			paste blind_test_set/fasta.tmp blind_test_set/dssp.tmp | awk -v aa=$$aa '$$1==aa{print $$2}' > blind_test_set/subset_dssp.tmp; \
			rm blind_test_set/fasta.tmp blind_test_set/dssp.tmp; \
			for i in $$(cat blind_test_set/subset_dssp.tmp); do \
				if [ $$i = 'H' ]; then H=$$(echo $$H+1 | bc); fi; \
				if [ $$i = 'E' ]; then E=$$(echo $$E+1 | bc); fi; \
				if [ $$i = '-' ]; then c=$$(echo $$c+1 | bc); fi; \
			done; \
			rm blind_test_set/subset_dssp.tmp; \
		done; \
		echo -e "$$aa\t$$H\t$$E\t$$c"; \
	done > $@

blind_test_set/dssp_by_residue_count.png: blind_test_set/dssp_by_residue_count.txt
	hist_dssp_by_residue $< $@ 

### TODO ###
## .META:
## 1    scop_class
## 2    pdb_count
blind_test_set/domain_by_scop_count.txt: refs/dir.cla.scope.2.06-stable.txt blind_test_set/selected_ids
	cat $$(echo $^ | cut -d ' ' -f1) | cut -f2,4 | awk '{print toupper($$1), $$2}' | sort | uniq \
	| fgrep -f <(cat $$(echo $^ | cut -d ' ' -f2) | cut -d ':' -f1) > blind_test_set/dir.cla.scope.filtered.tmp; \
	for class in $(SCOP_CLASSES); do \
		count=$$(cat blind_test_set/dir.cla.scope.filtered.tmp | cut -f2 | cut -d '.' -f1 | fgrep $$class | wc -l); \
		echo -e "$$class\t$$count"; \
	done > $@; \
	rm blind_test_set/dir.cla.scope.filtered.tmp; \
	cat blind_test_set/domain_by_scop_count.txt | sed '/0$$/d' > blind_test_set/domain_by_scop_count_nozero.tmp; \
	cut -f1 blind_test_set/domain_by_scop_count_nozero.tmp > blind_test_set/labels.tmp; \
	cut -f2 blind_test_set/domain_by_scop_count_nozero.tmp > blind_test_set/data.tmp; \
	pie_chart blind_test_set/data.tmp blind_test_set/labels.tmp $(SCOP_PIE_SIZE) blind_test_set/domain_by_scop_count.png; \
	rm blind_test_set/data.tmp blind_test_set/labels.tmp blind_test_set/domain_by_scop_count_nozero.tmp

####################
### GOR training ###
####################

## 1. create a list of training proteins ids according to the profile availability
training_set/training_ids.txt: 
	echo $$(ls training_set/pssm/) | sed 's/ /\n/g; s/.pssm//g' > training_set/ids.tmp; \
	check_pssm --pssm_dir training_set/pssm --id_list training_set/ids.tmp > training_set/all_zero_ids.tmp; \
	cat training_set/ids.tmp | grep -v -f training_set/all_zero_ids.tmp > $@; \
	rm training_set/all_zero_ids.tmp training_set/ids.tmp

WINDOW_SIZE = 17
## single sequence gor training:
#gor/single_training.gor: training_set/merged_fasta_files.gz training_set/merged_dssp_files.gz
#	old_gor_training --fasta training_set/merged_fasta_files.gz --dssp training_set/merged_dssp_files.gz --out $@ --window_size $(WINDOW_SIZE)
## profile-based gor training:
gor/profile_training.gor: training_set/training_ids.txt
	gor_training \
		--pssm_dir training_set/pssm \
		--dssp_dir training_set/dssp \
		--id_list training_set/training_ids.txt \
		--window_size $(WINDOW_SIZE) \
		--out $@ 

###################
### GOR testing ###
###################

## single sequence gor testing (not completed)
#impute_dssp_with_single_gor: gor/training.gor
#	for id in $$(cat blind_test_set/blind_selected_150.fasta | grep ">" | tr -d '>'); do \
#		cat blind_test_set/blind_selected_150.fasta | grep -A 1 $${id} > gor/$${id}.fasta.tmp; \
#		old_gor_testing --gor_model gor/training.gor --input gor/$${id}.fasta.tmp --out gor/test_single/$${id}.dssp; \
#		rm gor/$${id}.fasta.tmp; \
#	done 

## profile-based gor testing:
impute_dssp_with_profile_gor: gor/profile_training.gor blind_test_set/blind_selected_150.fasta
	gor_testing \
		--pssm_dir blind_test_set/selected_pssm \
		--id_list blind_test_set/selected_ids \
		--model $< \
		--out_dir gor/profile_testing \

## single_sov_performance or profile_sov_performance
## impute_dssp_with_gor has to be executed before this
gor/%_gor_sov_performance.txt: 
	sov_performance \
		--real_dssp_dir blind_test_set/selected_dssp \
		--imputed_dssp_dir gor/$*_testing \
		--window_size $(WINDOW_SIZE) \
		--id_list blind_test_set/selected_ids \
		--out $@

##################################
### Generate Sequence Profiles ###
##################################

## 1. index databases
refs/uniprot_sprot.fasta.phr: refs/uniprot_sprot.fasta
	makeblastdb -in $< -out $< -dbtype prot

training_set/training.fasta:
	for file in $$(ls fasta/); do cat fasta/$$file; done  > $@

training_set/training.fasta.phr: training_set/training.fasta
	makeblastdb -in $< -out $< -dbtype prot

blind_test_set/blind_selected_150.fasta.phr: blind_test_set/blind_selected_150.fasta
	makeblastdb -in $< -out $< -dbtype prot

## 2. run psiblast on training set
# TODO psiblast ascii_pssm outputs have to be created for each input
# sequence --> split the input into individual sequences before running psiblast

PSIBLAST_E_VALUE = 0.01

.PHONY: get_pssm_from_training_set
get_pssm_from_training_set:
	for jpred_id in $$(cat refs/jpred4.list.txt); do \
		psiblast -query fasta/$${jpred_id}.fasta -db refs/uniprot_sprot.fasta -evalue $(PSIBLAST_E_VALUE) \
		-out_ascii_pssm training_set/pssm/$${jpred_id}.pssm -out training_set/blast/$${jpred_id}.blast \
		-num_descriptions 10000 -num_alignments 10000 -num_iterations 3; \
	done

## ll blind_test_set/blast/ | sed -r 's/^.*(....:.).blast$/\1/g'

############################
### GOR Cross-Validation ###
############################

K_FOLD = $$(seq 0 4)

## K-fold ids definition
split_in_k_folds: 
	mkdir cross_validation/profile_gor/; \
	for i in $(K_FOLD); do \
		make cross_validation/testing$${i}_ids; \
		make cross_validation/training$${i}_ids; \
		mkdir cross_validation/profile_gor/testing$${i}; \
	done

## generate testing folds according to the existence of pssm files
cross_validation/testing%_ids: training_set/training_ids.txt
	cat $< | grep -f cross_validation/test$* > $@

## generate training folds according to the existence of pssm files
cross_validation/training%_ids: training_set/training_ids.txt
	cat $< | grep -v -f cross_validation/test$* > $@

## k-fold profile-based gor training definition
cross_val_train_profile_gor:
	for i in $(K_FOLD); do \
		make cross_validation/profile_gor/training$${i}.gor; \
	done

cross_validation/profile_gor/training%.gor: cross_validation/training%_ids
	gor_training \
		--pssm_dir training_set/pssm \
		--dssp_dir training_set/dssp \
		--id_list $< \
		--window_size $(WINDOW_SIZE) \
		--out $@

## k-fold profile-based gor testing definition:
cross_val_test_profile_gor:
	for i in $(K_FOLD); do \
		mkdir cross_validation/profile_gor/testing$${i}; \
		make cross_val_test$${i}_profile_gor; \
	done

cross_val_test%_profile_gor: cross_validation/profile_gor/training%.gor cross_validation/testing%_ids
	gor_testing \
		--pssm_dir training_set/pssm \
		--id_list cross_validation/testing$*_ids \
		--model $< \
		--out_dir cross_validation/profile_gor/testing$*

cross_val_clear_test_profile_gor:
	rm -fr cross_validation/profile_gor/testing?

## k-fold profile-based gor sov performance definition
cross_val_performance_profile_gor:
	for i in $(K_FOLD); do \
		make cross_validation/profile_gor/testing$${i}_gor_sov_performance.txt; \
	done

cross_validation/profile_gor/testing%_gor_sov_performance.txt: 
	sov_performance \
		--real_dssp_dir training_set/dssp \
		--imputed_dssp_dir cross_validation/profile_gor/testing$* \
		--window_size $(WINDOW_SIZE) \
		--id_list cross_validation/testing$*_ids \
		--out $@

## k-fold profile-based gor identity performance definition
cross_val_conf_matrix_gor:
	for i in $(K_FOLD); do \
		make cross_validation/profile_gor/testing$${i}_gor_conf_matrix.pdf; \
	done

cross_validation/profile_gor/testing%_gor_conf_matrix.pdf:
	confusion_matrix \
		--real_dssp_dir training_set/dssp \
		--imputed_dssp_dir cross_validation/profile_gor/testing$* \
		--window_size $(WINDOW_SIZE) \
		--id_list cross_validation/testing$*_ids \
		--out $@

cross_validation/profile_gor/overall_cv_gor_conf_matrix.pdf:
	overall_cv_performance_indexes \
		--real_dssp_dirs "training_set/dssp training_set/dssp training_set/dssp training_set/dssp training_set/dssp" \
		--imputed_dssp_dirs "cross_validation/profile_gor/testing0 cross_validation/profile_gor/testing1 cross_validation/profile_gor/testing2 \
			cross_validation/profile_gor/testing3 cross_validation/profile_gor/testing4" \
		--window_size $(WINDOW_SIZE) \
		--id_lists "cross_validation/testing0_ids cross_validation/testing1_ids cross_validation/testing2_ids \
			cross_validation/testing3_ids cross_validation/testing4_ids" \
		--out $@

####################
### SVM Training ###
####################

## e.g. svm/rbf_2.0_0.5_model.pkl.gz
svm/%_model.pkl.gz: training_set/training_ids.txt
	kernel=$$(echo "$*" | cut -d '_' -f1); \
	C=$$(echo "$*" | cut -d '_' -f2); \
	gamma=$$(echo "$*" | cut -d '_' -f3); \
	svm_training \
		--pssm_dir training_set/pssm \
		--dssp_dir training_set/dssp \
		--id_list training_set/training_ids.txt \
		--window_size $(WINDOW_SIZE) \
		--out $@ \
		--kernel $${kernel} \
		--C $${C} \
		--gamma $${gamma} > svm/$*_log.txt

###################
### SVM Testing ###
###################

## e.g. impute_dssp_with_SVM_rbf_2.0_0.5_model
impute_dssp_with_SVM_%_model: svm/%_model.pkl.gz blind_test_set/selected_ids
	svm_testing \
		--pssm_dir blind_test_set/selected_pssm \
		--id_list blind_test_set/selected_ids \
		--model $< \
		--out_dir svm/$*_testing \
		--windowSize $(WINDOW_SIZE)

## impute_dssp_with_SVM_%_model has to be executed before this
## e.g. svm/rbf_2.0_0.5_model_performance.txt
svm/%_model_performance.txt: 
	sov_performance \
		--real_dssp_dir blind_test_set/selected_dssp \
		--imputed_dssp_dir svm/$*_testing \
		--window_size $(WINDOW_SIZE) \
		--id_list blind_test_set/selected_ids \
		--out $@

############################
### SVM Cross-Validation ###
############################

## try with 	C = 2 and 4
##		gamma = 0.5 and 2
KERNEL=rbf
C=2.0
GAMMA=2.0
SVM_CV_DIR=cross_validation/$(KERNEL)_$(C)_$(GAMMA)_svm

## k-fold profile-based svm training definition
cross_val_train_svm:
	mkdir cross_validation/$(KERNEL)_$(C)_$(GAMMA)_svm_training; \
	for i in $(K_FOLD); do \
		make cross_validation_svm_training_fold$${i}; \
	done

cross_validation_svm_training_fold%: cross_validation/training%_ids
	svm_training \
		--pssm_dir training_set/pssm \
		--dssp_dir training_set/dssp \
		--id_list cross_validation/training$*_ids \
		--window_size $(WINDOW_SIZE) \
		--out $(SVM_CV_DIR)/$*_model.pkl.gz \
		--kernel $(KERNEL) \
		--C $(C) \
		--gamma $(GAMMA) > $(SVM_CV_DIR)/$*_log.txt

## k-fold svm testing definition:
cross_val_test_svm:
	for i in $(K_FOLD); do \
		mkdir cross_validation/$(KERNEL)_$(C)_$(GAMMA)_svm/testing$${i}; \
		make cross_val_test$${i}_svm; \
	done

cross_val_test%_svm: $(SVM_CV_DIR)/%_model.pkl.gz cross_validation/testing%_ids
	svm_testing \
		--pssm_dir training_set/pssm/ \
		--id_list cross_validation/testing$*_ids \
		--model $< \
		--out_dir $(SVM_CV_DIR)/testing$*/ \
		--windowSize $(WINDOW_SIZE); \
	
## k-fold svm sov performance definition
cross_val_performance_svm:
	for i in $(K_FOLD); do \
		make cross_validation_fold_$${i}_svm_sov_performance; \
	done

cross_validation_fold_%_svm_sov_performance: cross_validation/testing%_ids
	sov_performance \
		--real_dssp_dir training_set/dssp \
		--imputed_dssp_dir $(SVM_CV_DIR)/testing$* \
		--window_size $(WINDOW_SIZE) \
		--id_list $< \
		--out $(SVM_CV_DIR)/testing$*_sov_performance.txt

## k-fold svm identity performance definition
cross_val_conf_matrix_svm:
	for i in $(K_FOLD); do \
		make cross_validation_testing$${i}_svm_conf_matrix.pdf; \
	done

cross_validation_testing%_svm_conf_matrix.pdf:
	confusion_matrix \
		--real_dssp_dir training_set/dssp \
		--imputed_dssp_dir $(SVM_CV_DIR)/testing$* \
		--window_size $(WINDOW_SIZE) \
		--id_list cross_validation/testing$*_ids \
		--out $(SVM_CV_DIR)/testing$*_conf_matrix.pdf

cross_validation_svm_overall_cv_conf_matrix.pdf:
	overall_cv_performance_indexes \
		--real_dssp_dirs "training_set/dssp training_set/dssp training_set/dssp training_set/dssp training_set/dssp" \
		--imputed_dssp_dirs "$(SVM_CV_DIR)/testing0 $(SVM_CV_DIR)/testing1 $(SVM_CV_DIR)/profile_gor/testing2 \
			$(SVM_CV_DIR)/testing3 $(SVM_CV_DIR)/testing4" \
		--window_size $(WINDOW_SIZE) \
		--id_lists "cross_validation/testing0_ids cross_validation/testing1_ids cross_validation/testing2_ids \
			cross_validation/testing3_ids cross_validation/testing4_ids" \
		--out $(SVM_CV_DIR)/overall_cv_conf_matrix.pdf
