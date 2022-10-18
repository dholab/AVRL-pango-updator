#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channels
	ch_consensus_seqs = Channel
		.fromPath( "${params.data_dir}/DHO*/gisaid/*.fasta" )
		.map { fasta -> tuple( file(fasta), fasta.getParent(), fasta.getSimpleName() ) }
	
	// Workflow steps
	UPDATE_PANGO_DOCKER ( )
	
	UPDATE_PANGO_CONDA ( )
	
	println "pangolin updated to version:"
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
			.map { report, parentdir, experiment_number, experiment_date -> report }
			.collect()
	)
	
	GET_DESIGNATION_DATES ( )
	
	FIND_LONG_INFECTIONS (
		GET_DESIGNATION_DATES.out,
		IDENTIFY_LINEAGES.out
	)
	
	CONCAT_LONG_INFECTIONS (
		FIND_LONG_INFECTIONS.out.collect()
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
	// publishDir 'parentdir', pattern: '*.csv', mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	each cue
	tuple path(fasta), val(parentdir), val(run_name)
	
	output:
	tuple path("*.csv"), val(parentdir), val(experiment_number), env(experiment_date)
	
	script:
	experiment_number = "DHO_" + parentdir.toString().replaceAll('/gisaid','').split("DHO_")[1]
	
	"""
	experiment_date=`date -r ${fasta} "+%Y-%m-%d"`
	
	pangolin \
	--threads ${task.cpus} \
	--outfile ${experiment_number}_lineage_report_${params.date}.csv \
	${fasta}
	"""
	
}

process CONCAT_CSVS {
	
	// This process double checks that a pango csv exists in the results
	// directory. If it doesn't exist, it creates the file with the correct
	// column names. Then, it appends each row from each lineage report
	// produced during this run.
	
	publishDir params.results, pattern: 'all_lineage_reports*.csv', mode: 'copy'
	
	input:
	path report_list
	
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
	
	when:
	params.identify_long_infections == true
	
	output:
	path "*.csv"
	
	script:
	"""
	curl -fsSL https://raw.githubusercontent.com/corneliusroemer/pango-designation-dates/main/data/lineage_designation_date.csv > lineage_designation_date.csv
	"""
}

process FIND_LONG_INFECTIONS {
	
	tag "${experiment_number}"
	// publishDir parentdir, pattern: '*.csv', mode: 'copy'
	
	when:
	params.identify_long_infections == true
	
	input:
	each path(lineage_dates)
	tuple path(lineage_csv), val(parentdir), val(experiment_number), val(experiment_date)
	
	output:
	path "*putative_long_infections*.csv"
	
	script:
	"""
	long_infection_finder.R ${experiment_number} ${experiment_date} ${lineage_csv} ${lineage_dates}
	"""
}

process CONCAT_LONG_INFECTIONS {
	
	publishDir params.results, mode: 'copy'
	
	input:
	path file_list
	
	output:
	path "*.csv"
	
	script:
	"""
	concat_long_infections.R
	"""
	
}
// --------------------------------------------------------------- //
