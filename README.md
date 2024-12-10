# Analysis of modifications in cohort of 1kGP


The base and modification calling, as well as other read processing where done
using in nextflow. A minimal version of `nextflow 23.04.1` is required.

## Dependencies

The software dependencies below should be already available at the command line in $PATH,
as well as the python script in `src/python` directory

- Dorado version 0.7.2 
- Minimap2 version 2.24-r1122
- Samtools version 1.17
- Bedtools version 2.30.0
- Modkit version 0.2.2


## Read calling and processing

The nextflow pipeline should be executed in the next order and you should provide your `yaml` file
containing your own parameters. The template of the config parameters can be found inside every
nextflow pipeline directory in the file `nextflow.config`.

```
# Base and modification calling
nextflow run src/nextflow/basecall/main.nf \
-params-file your_custom_params.yaml

# Filter reads 
nextflow run src/nextflow/filter_reads/main.nf \
-params-file your_custom_params.yaml

# Map reads
nextflow run src/nextflow/map_reads/main.nf \
-params-file your_custom_params.yaml

# Pileup reads per haplotype
nextflow run src/nextflow/pileup_hap_mods/main.nf \
-params-file your_custom_params.yaml

```

