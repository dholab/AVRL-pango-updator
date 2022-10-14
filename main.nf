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
	UPDATE_PANGO ( )
	
	IDENTIFY_LINEAGES (
		UPDATE_PANGO.out.cue,
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



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Derivative parameters, mostly for making specific results folders


// --------------------------------------------------------------- //



// PROCESS SPECIFICATION 
// --------------------------------------------------------------- //

// process UPDATE_PANGO {
// 	
// 	// This process builds a new docker image with 
// 	
// 	when:
// 	params.update_pango == true
// 	
// 	output:
// 	val "pango updated", emit: cue
// 	
// 	script:
// 	"""
// 	docker build -t pangolin_updated:${params.date} ${params.dockerfile_path}
// 	docker tag pangolin_updated:${params.date} ${params.docker_reg}/pangolin_updated:${params.date}
// 	docker push ${params.docker_reg}/pangolin_updated:${params.date}
// 	"""
// }

process UPDATE_PANGO {
	
	when:
	params.update_pango == true
	
	conda 'pangolin'
	
	output:
	env env_path, emit: cue
	
	script:
	"""
	env_path=`conda env list | awk 'NR==3' | xargs | sed 's/base//g' | sed 's/*//g' | xargs`
	pangolin --update --update-data 
	"""
}

process IDENTIFY_LINEAGES {
	
	tag "${experiment_number}"
	// publishDir 'parentdir', pattern: '*.csv', mode: 'copy'
	
	when:
	params.update_pango == true && cue == "pango updated" || params.update_pango == false
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	conda "env_path"
	
	input:
	val env_path
	tuple path(fasta), val(parentdir), val(run_name)
	
	output:
	tuple path("*.csv"), val(parentdir), val(experiment_number), env(experiment_date)
	
	script:
	date = new java.util.Date().format('yyyyMMdd')
	experiment_number = "DHO_" + parentdir.toString().replaceAll('/gisaid','').split("DHO_")[1]
	
	"""
	experiment_date=`date -r ${fasta} "+%Y-%m-%d"`
	pangolin --threads ${task.cpus} --outfile lineage_report_${date}.csv ${fasta}
	"""
	
}

process CONCAT_CSVS {
	
	// This process double checks that a pango csv exists in the results
	// directory. If it doesn't exist, it creates the file with the correct
	// column names. Then, it appends each row from each lineage report
	// produced during this run.
	
	tag "DHO_${experiment_number}"
	publishDir params.results
	
	input:
	path report_list
	
	script:
	date = new java.util.Date().format('yyyyMMdd')
	
	"""
	echo "taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,scorpio_notes,version,pangolin_version,scorpio_version,constellation_version,is_designated,qc_status,qc_notes,note" > all_lineage_reports_${date}.csv
	
	find . -type f -name "lineage_report*.csv" > lineage_reports.txt
	
	for i in `cat lineage_reports.txt`;
	do
		tail -n +2 \$i >> all_lineage_reports_${date}.csv
	done
	"""
}

process GET_DESIGNATION_DATES {
	
	output:
	path "*.csv"
	
	script:
	"""
	curl -fsSL https://raw.githubusercontent.com/corneliusroemer/pango-designation-dates/main/data/lineage_designation_date.csv > lineage_designation_date.csv
	"""
}

process FIND_LONG_INFECTIONS {
	
	tag "${experiment_number}"
	// publishDir 'parentdir', pattern: '*.csv', mode: 'copy'
	
	input:
	each path(lineage_dates)
	tuple path(lineage_csv), val(parentdir), val(experiment_number), val(experiment_date)
	
	output:
	path "*.csv"
	
	script:
	"""
	long_infection_finder.R \
	${experiment_number} ${experiment_date} "${parentdir}" \
	${lineage_csv} \
	${lineage_dates}
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
