# AVRL-pango-updator

To invoke, run:
```
nextflow run dholab/AVRL-pango-updator \
-profile conda \
--update_pango true \
--data_dir "/Volumes/GoogleDrive/Shared drives/2019-nCoV open research team/Sequencing Data" \
--results results/ \
--identify_long_infections true \
--days_of_infection 240 \
--classify_mutation_levels true \
```


```
nextflow run main.nf \
-profile conda \
--update_pango true \
--data_dir "/Volumes/GoogleDrive/Shared drives/2019-nCoV open research team/Sequencing Data" \
--results results/ \
--identify_long_infections true \
--days_of_infection 240 \
--classify_mutation_levels true \
```
