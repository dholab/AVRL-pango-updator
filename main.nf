#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channel for consensus sequences 
	ch_consensus_seqs = Channel
		.fromPath( "${params.data_dir}/DHO*/gisaid/*.fasta" )
		.filter { !it.contains(" ") }
		
	
	// Workflow steps for reclassifying pango lineages
	UPDATE_PANGO_DOCKER ( )
	
	UPDATE_PANGO_CONDA ( )
	
	if( params.update_pango == true ){
		println "Pangolin updated to version:"
		UPDATE_PANGO_DOCKER.out ? UPDATE_PANGO_CONDA.out.cue.view() : UPDATE_PANGO_DOCKER.out.cue.view()
	}
	
	IDENTIFY_LINEAGES (
		UPDATE_PANGO_DOCKER.out.cue
			.mix (
				UPDATE_PANGO_CONDA.out.cue
			),
		ch_consensus_seqs
			.map { fasta -> tuple( file(fasta), fasta.getParent(), fasta.getSimpleName() ) }
	)
	
	CONCAT_CSVS (
		IDENTIFY_LINEAGES.out
			.map { report, run_name, parentdir, experiment_number, experiment_date -> report }
			.collect()
	)
	
	// Workflow steps for identifying putative prolonged infections
	GET_DESIGNATION_DATES ( )
	
	FIND_LONG_INFECTIONS (
		GET_DESIGNATION_DATES.out,
		IDENTIFY_LINEAGES.out
	)
	
	CONCAT_LONG_INFECTIONS (
		FIND_LONG_INFECTIONS.out
			.collect()
	)
	
	// Workflow steps for classifying RBD mutation levels relative to BA.2
	GET_LINEAGE_SEQS ( )
	
	UNZIP_LINEAGE_SEQS ( 
		GET_LINEAGE_SEQS.out
	)
	
	ISOLATE_BA_2 ( 
		UNZIP_LINEAGE_SEQS.out
	)
	
	MAP_TO_BA_2 (
		ISOLATE_BA_2.out,
		ch_consensus_seqs
			.filter { it.lastModified() > (new Date("07/12/2021").getTime()) }
			.splitFasta( file: true )
	)
	
	PROCESS_WITH_SAMTOOLS (
		ISOLATE_BA_2.out,
		MAP_TO_BA_2.out
	)
	
	CALL_RBD_VARIANTS (
		ISOLATE_BA_2.out,
		PROCESS_WITH_SAMTOOLS.out
	)
	
	CLASSIFY_LEVELS (
		CALL_RBD_VARIANTS.out.collect(),
		CONCAT_CSVS.out
	)
	
	
}
// --------------------------------------------------------------- //



// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

process UPDATE_PANGO_DOCKER {
	
	// This process builds a new docker image with the latest available pangolin version
	
	when:
	(workflow.profile == 'standard' || workflow.profile == 'docker') && params.update_pango == true
	
	output:
	env version, emit: cue
	
	script:
	"""
	# docker build -t pangolin_updated:${params.date} ${params.dockerfile_path}
	# docker tag pangolin_updated:${params.date} ${params.docker_reg}/pangolin_updated:${params.date}
	# docker push ${params.docker_reg}/pangolin_updated:${params.date}
	pangolin --update --update-data
	version=`pangolin --version | sed 's/pangolin//g' | xargs`
	"""
}

process UPDATE_PANGO_CONDA {
	
	// This process creates a conda environment, where it updates pangolin
	
	when:
	workflow.profile == 'conda' && params.update_pango == true
	
	output:
	env version, emit: cue
	
	script:
	"""
	pangolin --update --update-data
	version=`pangolin --version | sed 's/pangolin//g' | xargs`
	"""
}

process IDENTIFY_LINEAGES {
	
	// This process identifies lineages for all samples in every DHO Lab
	// sequencing run. It pulls FASTAs directly from Google Drive to do so.
	// NOTE: You must have Google Drive installed and mounted to access 
	// the input files for this step.
	
	tag "${experiment_number}"
	// publishDir "${parentdir}", pattern: '*.csv', mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	each cue
	tuple path(fasta), val(parentdir), val(run_name)
	
	output:
	tuple path("*.csv"), val(run_name), val(parentdir), val(experiment_number), env(experiment_date)
	
	script:
	experiment_number = "DHO_" + parentdir.toString().replaceAll('/gisaid','').split("DHO_")[1]
	
	"""
	experiment_date=`date -r ${fasta} "+%Y-%m-%d"`
	
	pangolin \
	--threads ${task.cpus} \
	--outfile ${experiment_number}_lineages_updated_${params.date}.csv \
	"${fasta}"
	"""
	
}

process CONCAT_CSVS {
	
	// This process takes the new pangolin lineage reports for all sequencing runs
	// and combines them into a single, large, CSV file. Users can then search that
	// CSV for lineages of interest 
	
	publishDir params.results, pattern: 'all_lineage_reports*.csv', mode: 'copy'
	
	input:
	path report_list, stageAs: 'report??.csv'
	
	output:
	path "all_lineage_reports*.csv"
	
	script:
	"""
	find . -name "*.csv" > lineage_reports.txt
	
	echo "taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,scorpio_notes,version,pangolin_version,scorpio_version,constellation_version,is_designated,qc_status,qc_notes,note" > all_lineage_reports_${params.date}.csv
	
	for i in `cat lineage_reports.txt`;
	do
		tail -n +2 \$i >> all_lineage_reports_${params.date}.csv
	done
	"""
}


process GET_DESIGNATION_DATES {
	
	// This process downloads a table of pangolin lineage designation dates
	// from Cornelius Roemer's GitHub. These dates represent when each lineage was
	// added to pangolin, after which point sequences could be classified as such 
	
	publishDir params.results, mode: 'copy'
	
	when:
	params.identify_long_infections == true
	
	output:
	path "*.csv"
	
	script:
	"""
	curl -fsSL https://raw.githubusercontent.com/corneliusroemer/pango-designation-dates/main/data/lineage_designation_date.csv > lineage_designation_dates.csv
	"""
}

process FIND_LONG_INFECTIONS {
	
	// This process calls a script that identifies all samples that classify as "old"
	// lineages, i.e. lineages that were first designated in pangolin long before a 
	// sequence classified as it in a sequencing run. By default, this amount of time is
	// eight months, but this is subject to change in the future. In such cases, where
	// a lineage appears long after it arose and subsided, it is more likely that the 
	// infection individual has sustained a prolonged infection since that lineage 
	// was prevalent, and less likely that the old lineage re-appeared despite competition
	// from newer, more fit lineages.
	
	tag "${experiment_number}"
	// publishDir "${parentdir}", pattern: '*.csv', mode: 'copy'
	
	when:
	params.identify_long_infections == true
	
	input:
	each path(lineage_dates)
	tuple path(lineage_csv), val(run_name), val(parentdir), val(experiment_number), val(experiment_date)
	
	output:
	path "*putative_long_infections*.csv"
	
	script:
	"""
	long_infection_finder.R ${run_name} ${experiment_number} ${experiment_date} ${lineage_csv} ${lineage_dates} ${params.days_of_infection}
	"""
}

process CONCAT_LONG_INFECTIONS {
	
	// This step simply concatenates any putative prolonged infection samples into on
	// table, which is then exported into the results directory.
	
	publishDir params.results, mode: 'copy'
	
	input:
	path file_list, stageAs: 'infections??.csv'
	
	output:
	path "*.csv"
	
	script:
	"""
	concat_long_infections.R "${params.days_of_infection}"
	"""
	
}

process GET_LINEAGE_SEQS {
	
	// This process downloads a FASTA with consensus sequences for all pango lineages.
	// As state on Cornelius Roemer's GitHub repo, "[t]hese sequences are not real sequences
	// in databases but are algorithmically constructed consensus sequences that try to 
	// represent the common ancestor sequence of that lineage." It then uses a short R
	// script to pull out the BA.2 sequence, which will serve as a reference sequence
	// downstream 
	
	cpus 1
	
	when:
	params.classify_mutation_levels == true
	
	output:
	path "pango_consensus_sequences.fasta.zstd"
	
	script:
	"""
	curl -fsSL 'https://github.com/corneliusroemer/pango-sequences/blob/main/data/pango_consensus_sequences.fasta.zstd?raw=true' > pango_consensus_sequences.fasta.zstd
	"""
}

process UNZIP_LINEAGE_SEQS {
	
	cpus 1
	
	input:
	path zstd
	
	output:
	path "pango_consensus_sequences.fasta"
	
	script:
	"""
	mv `realpath ${zstd}` pango_consensus_sequences.fasta.zst
	zstd -d pango_consensus_sequences.fasta.zst
	"""
}

process ISOLATE_BA_2 {
	
	publishDir params.results, mode: 'copy'
	
	cpus 1
	
	input:
	path fasta
	
	output:
	path "ba_2_ref.fasta"
	
	script:
	"""
	ba_2_isolator.R
	"""
}

process MAP_TO_BA_2 {
	
	// This process maps each consensus sequence from each run after the designation
	// of BA.2 to the BA.2 consensus sequences. It then sorts the alignment, converts it
	// to a BAM, and then constructs a pile-up that will be used as input for variant-
	// calling downstream.
	
	errorStrategy 'retry'
	maxRetries 4
	
	cpus 1
	
	input:
	each path(refseq)
	path fasta
	
	output:
	tuple path("*.sam"), val(sample)
	
	script:
	sample = fasta.getBaseName()
	"""
	minimap2 -t 1 -a ${refseq} -o ${sample}.sam ${fasta}
	"""
	
}

process PROCESS_WITH_SAMTOOLS {
	
	// This process takes the sequence alignment map (SAM) output from
	// process MAP_TO_BA_2, converts it to binary format (BAM), sorts it (which
	// is only a formality here, as there's just one consensus sequence being 
	// processed), and then "converts" it to the mpileup format iVar requires.
	
	input:
	each path(refseq)
	tuple path(sam), val(sample)
	
	output:
	tuple path("*.mpileup"), val(sample)
	
	script:
	"""
	cat ${sam} \
		| samtools view -Sb - \
		| samtools sort - > tempfile
	samtools mpileup -aa -f ${refseq} --output ${sample}.mpileup tempfile
	"""
	
}

process CALL_RBD_VARIANTS {
	
	// This process creates a simple table of mutations from BA.2 and annotates them
	// with gene, codon, and amino acid information.
	
	input:
	each path(refseq)
	tuple path(mpileup), val(sample)
	
	output:
	path "*.tsv"
	
	script:
	"""
	cat ${mpileup} \
	  | ivar variants -p ${sample}_consensus_variant_table \
	  -t 0 -m 1 -q 1 -r ${refseq} -g ${params.refgff}
	"""
	
}

process CLASSIFY_LEVELS {
	
	// This process uses an R script to trim down all observed variants to only those in
	// the Spike protein receptor binding domain, i.e. Spike amino acid residues 319–541.
	// It then adds the count of RBD mutations——an RBD mutation "level", a la Cornelius 
	// Roemer's method for tracking convergent evolution among SARS-CoV-2 lineages——to
	// the new pango lineage classifications.
	
	input:
	path variant_table_list
	path lineages
	
	output:
	
	script:
	"""
	rbd_lineage_classifer.R
	"""
}
// --------------------------------------------------------------- //
