# AVRL-pango-updator

To invoke, run:
```
nextflow run dholab/AVRL-pango-updator \
--update_pango true \
--identify_long_infections true \
--classify_mutation_levels true \
--days_of_infection 240 \
--data_dir "/Volumes/GoogleDrive/Shared drives/2019-nCoV open research team/Sequencing Data" \
--results results/ \
-profile conda
```

Or, if you would prefer to clone the workflow first:
```
git clone https://github.com/dholab/AVRL-pango-updator.git .

nextflow run main.nf \
--update_pango true \
--identify_long_infections true \
--classify_mutation_levels true \
--days_of_infection 240 \
--data_dir "/Volumes/GoogleDrive/Shared drives/2019-nCoV open research team/Sequencing Data" \
--results results/ \
-profile conda
```
