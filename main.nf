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
			.splitCsv( header: true )
	)
	
	// FIND_LONG_INFECTIONS (
	// 	IDENTIFY_LINEAGES.out.collect()
	// )
	
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
	
	// This process builds a new docker image with 
	
	when:
	params.update_pango == true
	
	output:
	val "pango updated", emit: cue
	
	script:
	"""
	pangolin --update --update-data 
	"""
}

process IDENTIFY_LINEAGES {
	
	tag "DHO_${experiment_number}"
	// publishDir 'parentdir', pattern: '*.csv', mode: 'copy'
	
	when:
	params.update_pango == true && cue == "pango updated" || params.update_pango == false
	
	input:
	val cue
	tuple path(fasta), val(parentdir), val(run_name)
	
	output:
	path csv
	
	script:
	date = new java.util.Date().format('yyyyMMdd')
	experiment_number = parentdir.toString().replaceAll('/gisaid','').split("DHO_")[1]
	
	"""
	pangolin --outfile lineage_report_${date}.csv ${fasta}
	"""
	
}

process CONCAT_CSVS {
	
	// This process double checks that a pango csv exists in the results
	// directory. If it doesn't exist, it creates the file with the correct
	// column names. Then, it appends each row from each lineage report
	// produced during this run.
	
	input:
	val row
	
	script:
	date = new java.util.Date().format('yyyyMMdd')
	
	"""
	if test -f !{params.results}/all_lineage_reports_${date}.csv; then
		touch "${params.results}/all_lineage_reports_${date}.csv"
		echo "taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,scorpio_notes,version,pangolin_version,scorpio_version,constellation_version,is_designated,qc_status,qc_notes,note" > "${params.results}/all_lineage_reports_${date}.csv"
	fi
	
	cat row >> "${params.results}/all_lineage_reports_${date}.csv"
	"""
}

// process FIND_LONG_INFECTIONS {}

// --------------------------------------------------------------- //
