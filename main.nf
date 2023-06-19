#!/usr/bin/env nextflow

nextflow.enable.dsl = 2



// DERIVATIVE PARAMETER SPECIFICATION
// --------------------------------------------------------------- //
// Derivative parameters, mostly for making specific results folders
params.runs = file( params.run_file )
params.runs.text = params.runs_of_interest
params.lineage_reports = params.results + "/" + params.date + "_pangolin_reports"
params.sample_vcfs = params.results + "/" + params.date + "_sample_vcfs"
if( params.data_dir.contains("2019-nCoV open research team/Sequencing Data")){
	params.input_path = params.data_dir + "/DHO*/gisaid/*.fasta"
} else {
	params.input_path = params.data_dir + "/**/*.fasta"
}
// --------------------------------------------------------------- //



// WORKFLOW SPECIFICATION
// --------------------------------------------------------------- //
workflow {
	
	
	// input channel for consensus sequences 
	ch_consensus_seqs = Channel
		.fromPath( params.input_path )
		.filter { !it.toString().contains(" copy") }
	
	ch_runs_of_interest = Channel
		.fromPath( params.runs )
		.splitCsv( header: false, sep: "," )
		.flatten()
	
	
	// Workflow steps for reclassifying all pango lineages
	UPDATE_PANGO_DOCKER ( )
	
	UPDATE_PANGO_CONDA ( )
	
	println "Pangolin updated to version:"
	UPDATE_PANGO_CONDA.out.cue.view()
	UPDATE_PANGO_DOCKER.out.cue.view()
	
	RECLASSIFY_ALL_LINEAGES (
		UPDATE_PANGO_DOCKER.out.cue
			.mix (
				UPDATE_PANGO_CONDA.out.cue
			),
		ch_consensus_seqs
			.map { fasta -> tuple( file(fasta), fasta.getParent(), fasta.getSimpleName() ) }
	)
	
	FIND_TARGET_SEQS (
		UPDATE_PANGO_DOCKER.out.cue
		.mix (
			UPDATE_PANGO_CONDA.out.cue
		),
		ch_runs_of_interest
	)
	
	CLASSIFY_TARGET_SEQS (
		FIND_TARGET_SEQS.out
	)
	
	CONCAT_CSVS (
		RECLASSIFY_ALL_LINEAGES.out
			.map { report, run_name, parentdir, experiment_number, experiment_date -> report }
			.mix (
				CLASSIFY_TARGET_SEQS.out
					.map { report, run_name, parentdir, experiment_number, experiment_date -> report }
			)
			.collect()
	)
	
	// Workflow steps for identifying all putative prolonged infections
	GET_DESIGNATION_DATES ( )
	
	FIND_LONG_INFECTIONS (
		GET_DESIGNATION_DATES.out,
		RECLASSIFY_ALL_LINEAGES.out
			.mix ( CLASSIFY_TARGET_SEQS.out )
	)
	
	CONCAT_LONG_INFECTIONS (
		FIND_LONG_INFECTIONS.out.collect()
	)
	
	// Workflow steps for classifying RBD mutation levels relative to BA.2
	GET_LINEAGE_SEQS ( )
	
	UNZIP_LINEAGE_SEQS ( 
		GET_LINEAGE_SEQS.out
	)
	
	ISOLATE_BA_2 ( 
		UNZIP_LINEAGE_SEQS.out
	)
	
	ALIGN_ALL_TO_BA_2 ( 
		ISOLATE_BA_2.out,
		ch_consensus_seqs
			.filter { it.lastModified() > (new Date("07/12/2021").getTime()) }
			.splitFasta( file: true )
	)
	
	EXTRACT_SAMPLE ( 
		ALIGN_ALL_TO_BA_2.out
	)
	
	MAP_ALL_TO_BA_2 (
		ISOLATE_BA_2.out,
		EXTRACT_SAMPLE.out
	)
	
	ALIGN_TARGETS_TO_BA_2 ( 
		ISOLATE_BA_2.out,
		FIND_TARGET_SEQS.out
			.map { fasta, parentdir, experiment -> file(fasta) }
			.filter { file(it).lastModified() > (new Date("07/12/2021").getTime()) }
			.splitFasta( file: true )
	)
	
	EXTRACT_TARGET_SAMPLE ( 
		ALIGN_TARGETS_TO_BA_2.out
	)
	
	MAP_TARGETS_TO_BA_2 (
		ISOLATE_BA_2.out,
		EXTRACT_TARGET_SAMPLE.out
	)
	
	CALL_RBD_VARIANTS (
		ISOLATE_BA_2.out,
		MAP_ALL_TO_BA_2.out
			.mix ( MAP_TARGETS_TO_BA_2.out )
	)
	
	// BUILD_SNPEFF_DATABASE ( )
	// 
	// ANNOTATE_VCFS (
	// 	BUILD_SNPEFF_DATABASE.out,
	// 	CALL_RBD_VARIANTS.out
	// )
	
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
	workflow.profile == 'standard' || workflow.profile == 'docker' || workflow.profile == 'singularity'
	
	output:
	env version, emit: cue
	
	script:
	if (params.update_pango == true)
		"""
		# docker build -t pangolin_updated:${params.date} ${params.dockerfile_path}
		# docker tag pangolin_updated:${params.date} ${params.docker_reg}/pangolin_updated:${params.date}
		# docker push ${params.docker_reg}/pangolin_updated:${params.date}
		pangolin --update --update-data
		version=`pangolin --version | sed 's/pangolin//g' | xargs`
		"""
	else
		"""
		version=`pangolin --version | sed 's/pangolin//g' | xargs`
		"""
}

process UPDATE_PANGO_CONDA {
	
	// This process creates a conda environment, where it updates pangolin
	
	when:
	workflow.profile == 'conda'
	
	output:
	env version, emit: cue
	
	script:
	if (params.update_pango == true)
		"""
		pangolin --update --update-data
		version=`pangolin --version | sed 's/pangolin//g' | xargs`
		"""
	else
		"""
		version=`pangolin --version | sed 's/pangolin//g' | xargs`
		"""
}

process RECLASSIFY_ALL_LINEAGES {
	
	// This process identifies lineages for all samples in every DHO Lab
	// sequencing run. It pulls FASTAs directly from Google Drive to do so.
	// NOTE: You must have Google Drive installed and mounted to access 
	// the input files for this step.
	
	tag "${experiment_number}"
	publishDir "${run_dir}", pattern: '*.csv', mode: 'copy'
	
	errorStrategy { task.attempt < 3 ? 'retry' : 'ignore' }
	maxRetries 3
	
	input:
	each cue
	tuple path(fasta), val(parentdir), val(run_name)
	
	output:
	tuple path("*.csv"), val(run_name), val(run_dir), val(experiment_number), env(experiment_date)
	
	when:
	params.runs_of_interest.isEmpty()
	
	script:
	if( params.distribute_results == true )
		run_dir = parentdir.toString().replaceAll('/gisaid','')
	else 
		run_dir = params.lineage_reports
	experiment_number = 'DHO_' + parentdir.toString().split("DHO_")[1].replaceAll('/gisaid','')
	
	"""
	experiment_date=`date -r ${fasta} "+%Y-%m-%d"`
	
	pangolin \
	--threads ${task.cpus} \
	--outfile ${experiment_number}_lineages_updated_${params.date}.csv \
	"${fasta}"
	"""
	
}

process FIND_TARGET_SEQS {
	
	tag "${experiment}"
	
	input:
	each cue
	val experiment
	
	output:
	tuple env(fasta), val(parentdir), val(experiment)
	
	script:
	if( params.data_dir.contains("2019-nCoV open research team/Sequencing Data") )
		parentdir = params.data_dir + "/" + experiment + "/gisaid/"
	else 
		parentdir = params.data_dir + "/" + experiment + "/"
	"""
	fasta=`find "${parentdir}" -maxdepth 2 -type f -name "*.fasta"`
	"""
	
}

process CLASSIFY_TARGET_SEQS {
	
	tag "${experiment_number}"
	publishDir "${run_dir}", pattern: '*.csv', mode: 'copy'
	
	cpus 1
	errorStrategy 'retry'
	maxRetries 4
	
	input:
	tuple path(fasta), val(parentdir), val(experiment_number)
	
	output:
	tuple path("*.csv"), val(experiment_number), val(run_dir), val(experiment_number), env(experiment_date)
	
	script:
	if( params.distribute_results == true )
		run_dir = parentdir.toString().replaceAll('/gisaid','')
	else 
		run_dir = params.lineage_reports
	
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
	
	publishDir params.resources, mode: 'copy'
	
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
	// publishDir "${run_dir}", pattern: '*.csv', mode: 'copy'
	
	input:
	each path(lineage_dates)
	tuple path(lineage_csv), val(run_name), val(run_dir), val(experiment_number), val(experiment_date)
	
	output:
	path "*putative_long_infections*.csv"
	
	when:
	params.identify_long_infections == true
	
	script:
	run = run_name.toString().replaceAll(" copy", "")
	
	"""
	long_infection_finder.R ${run} ${experiment_number} ${experiment_date} ${lineage_csv} ${lineage_dates} ${params.days_of_infection}
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
	
	output:
	path "pango_consensus_sequences.fasta.zstd"
	
	when:
	params.classify_mutation_levels == true
	
	script:
	"""
	curl -fsSL 'https://github.com/corneliusroemer/pango-sequences/raw/main/data/pango-consensus-sequences_genome-nuc.fasta.zst' > pango_consensus_sequences.fasta.zstd
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
	cp `realpath ${zstd}` pango_consensus_sequences.fasta.zst
	zstd -d pango_consensus_sequences.fasta.zst
	"""
}

process ISOLATE_BA_2 {
	
	publishDir params.resources, mode: 'copy'
	
	cpus 1
	
	input:
	path fasta
	
	output:
	path "ba_2_ref.fasta"
	
	script:
	"""
	ba_2_isolator.R
	
	head -n 1 ba_2_ref.fasta > ba_2_ref.tmp
	tail -n +2 ba_2_ref.fasta | fold -w 80 >> ba_2_ref.tmp
	rm -f ba_2_ref.fasta
	mv ba_2_ref.tmp ba_2_ref.fasta
	"""
}

process ALIGN_ALL_TO_BA_2 {
	
	errorStrategy 'retry'
	maxRetries 4
	
	time '5minutes'
	
	input:
	each path(refseq)
	path fasta
	
	output:
	tuple path("*.fasta"), val(sample), env(strain)
	
	when:
	params.runs_of_interest.isEmpty()
	
	script:
	sample = fasta.getBaseName() 
	"""
	strain=`head -n 1 ${fasta}`
	cat *.fasta | muscle -diags -maxIters 1 -out ${sample}_BA2_MSA.fasta
	"""
	
}

process EXTRACT_SAMPLE {
	
	input:
	tuple path(fasta), val(sample), val(strain)
	
	output:
	tuple path("*.fasta"), val(sample), val(strain)
	
	script:
	"""
	extract_sample_seq.R ${fasta} ${sample}
	"""
	
}

process MAP_ALL_TO_BA_2 {
	
	// This process maps each consensus sequence from each run after the designation
	// of BA.2 to the BA.2 consensus sequences. It then sorts the alignment, converts it
	// to a BAM, and then constructs a pile-up that will be used as input for variant-
	// calling downstream.
	
	errorStrategy 'retry'
	maxRetries 4
	
	cpus 1
	
	input:
	each path(refseq)
	tuple path(fasta), val(sample), val(strain)
	
	output:
	tuple path("*.sam"), val(sample), val(strain)
	
	when:
	params.runs_of_interest.isEmpty()
	
	script:
	"""
	minimap2 -x asm10 --frag=yes --secondary=yes -N 5 -p 0.8 -t 1 -a \
	${refseq} "${fasta}" > ${sample}.sam
	"""
	
}

process ALIGN_TARGETS_TO_BA_2 {
	
	errorStrategy 'retry'
	maxRetries 4
	
	time '5minutes'
	
	input:
	each path(refseq)
	path fasta
	
	output:
	tuple path("*.fasta"), val(sample), env(strain)
	
	script:
	sample = fasta.getBaseName() 
	"""
	strain=`head -n 1 ${fasta}`
	cat *.fasta | muscle -diags -maxIters 1 -out ${sample}_BA2_MSA.fasta
	"""
	
}

process EXTRACT_TARGET_SAMPLE {
	
	input:
	tuple path(fasta), val(sample), val(strain)
	
	output:
	tuple path("*.fasta"), val(sample), val(strain)
	
	script:
	"""
	extract_sample_seq.R ${fasta} ${sample}
	"""
	
}

process MAP_TARGETS_TO_BA_2 {
	
	errorStrategy 'retry'
	maxRetries 4
	
	cpus 1
	time '5minutes'
	
	input:
	each path(refseq)
	tuple path(fasta), val(sample), val(strain)
	
	output:
	tuple path("*.sam"), val(sample), val(strain)
	
	script:
	"""
	minimap2 -x asm10 --frag=no --secondary=no -t 1 -a \
	${refseq} "${fasta}" > ${sample}.sam
	"""
	
}

process CALL_RBD_VARIANTS {
	
	// This process creates a simple table of mutations from BA.2
	
	publishDir params.sample_vcfs, mode: 'copy', overwrite: true
	
	errorStrategy 'retry'
	maxRetries 4
	
	memory '1 GB'
	time '5minutes'
	
	input:
	each path(refseq)
	tuple path(sam), val(sample), val(strain)
	
	output:
	path "*.vcf"
	
	script:
	strain_name = strain.replaceAll("/","_").replaceAll(">", "")
	"""
	callvariants.sh -Xmx1g \
	in=${sam} out=${strain_name}.vcf \
	ref=${refseq} samstreamer=t clearfilters \
	fixjunk realign=t secondary=t qtrim=f unclip=g ploidy=1 mincov=0 \
	callsub=t calldel=t callins=t overwrite=t
	"""
	
}

process BUILD_SNPEFF_DATABASE {
	
	output:
	path "*.config"
	
	script:
	"""
	snpEff build -gff3 ${params.refgff}
	"""
}

process ANNOTATE_VCFS {
	
	tag "${strain}"
	
	input:
	each path(config)
	tuple path(vcf), val(strain_name)
	
	output:
	path "*annotated.vcf"
	
	script:
	"""
	
	"""
	
}

process CLASSIFY_LEVELS {
	
	// This process uses an R script to trim down all observed variants to only those in
	// the Spike protein receptor binding domain, i.e. Spike amino acid residues 319–541.
	// It then adds the count of RBD mutations——an RBD mutation "level", a la Cornelius 
	// Roemer's method for tracking convergent evolution among SARS-CoV-2 lineages——to
	// the new pango lineage classifications.
	
	publishDir params.results, mode: 'copy'
	
	input:
	path variant_table_list
	path lineages
	
	output:
	path "rbd_classified_lineage_reports*.csv"
	
	script:
	"""
	rbd_lineage_classifer.R && \
	rm -f ${params.results}/all_lineage_reports_${params.date}.csv
	"""
}
// --------------------------------------------------------------- //
