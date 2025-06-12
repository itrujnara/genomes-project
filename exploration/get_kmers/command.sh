nextflow run itrujnara/get_kmers -r dev -profile crg,singularity -c crg_igor.config -c tower.config -params-file params.yml --input samplesheet_all.csv --outdir out -bg > nextflow.log 2> nextflow.err
