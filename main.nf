#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channels
	ch_consensus_seqs = Channel
		.fromPath( "${params.data_dir}/DHO*/gisaid/*.fasta" )
		.map { fasta -> tuple( file(fasta), fasta.getParent(), fasta.getSimpleName() ) }
		.filter { !it[2].contains(" copy") }
	
	// Workflow steps for reclassifying pango lineages
	UPDATE_PANGO_DOCKER ( )
	
	UPDATE_PANGO_CONDA ( )
	
	println "Pangolin updated to version:"
	UPDATE_PANGO_DOCKER.out ? UPDATE_PANGO_CONDA.out.view() : UPDATE_PANGO_DOCKER.out.view()
	
	IDENTIFY_LINEAGES (
		UPDATE_PANGO_DOCKER.out.cue
			.mix (
				UPDATE_PANGO_CONDA.out.cue
			),
		ch_consensus_seqs
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
	GET_BA_2_SEQ ( )
	
	MAP_TO_BA_2 (
		GET_BA_2_SEQ.out,
		ch_consensus_seqs
	)
	
	CALL_RBD_VARIANTS (
		GET_BA_2_SEQ.out,
		MAP_TO_BA_2.out
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
	
	// This process builds a new docker image with 
	
	when:
	workflow.profile == 'docker'
	
	output:
	env version, emit: cue
	
	script:
	"""
	docker build -t pangolin_updated:${params.date} ${params.dockerfile_path}
	docker tag pangolin_updated:${params.date} ${params.docker_reg}/pangolin_updated:${params.date}
	docker push ${params.docker_reg}/pangolin_updated:${params.date}
	version=`pangolin --version | sed 's/pangolin//g' | xargs`
	"""
}

process UPDATE_PANGO_CONDA {
	
	when:
	workflow.profile == 'conda'
	
	output:
	env version, emit: cue
	
	script:
	"""
	pangolin --update --update-data
	version=`pangolin --version | sed 's/pangolin//g' | xargs`
	"""
}

process IDENTIFY_LINEAGES {
	
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
	
	// This process double checks that a pango csv exists in the results
	// directory. If it doesn't exist, it creates the file with the correct
	// column names. Then, it appends each row from each lineage report
	// produced during this run.
	
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

process GET_BA_2_SEQ {
	
	publishDir params.results, mode: 'copy'
	
	when:
	params.classify_mutation_levels == true
	
	output:
	path "*.fasta"
	
	script:
	"""
	
	curl -fsSL https://github.com/corneliusroemer/pango-sequences/blob/main/data/pango_consensus_sequences.fasta.zstd?raw=true > pango_consensus_sequences.fasta.zstd
	unzstd pango_consensus_sequences.fasta.zstd
	
	ba_2_isolator.R
	
	"""
}

process MAP_TO_BA_2 {
	
	tag "${experiment_number}"
	
	when:
	file(fasta).lastModified() >> Date.parseToStringDate("2021-12-07").format('yyyy-M-d')
	
	input:
	each path(refseq)
	tuple path(fasta), val(parentdir), val(run_name)
	
	output:
	tuple path("*.mpileup"), val(sample)
	
	script:
	"""
	minimap2 -a ${refseq} ${fasta} \
	  | samtools view -Sb - \
	  | samtools sort - > tempfile
	  samtools mpileup -aa -f ${refseq} --output ${sample}.mpileup tempfile
	"""
	
}

process CALL_RBD_VARIANTS {
	
	tag "${experiment_number}"
	
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
