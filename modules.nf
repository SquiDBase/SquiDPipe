
/*
 * Concatenate FASTQ files per barcode folder
 */

process CONCATENATE_FASTQ {
    label 'process_low'
    conda params.conda_envs.default_env

    input:
    path fqDir
    val meta

    output:
    tuple path("${meta.folder}_allreads.fastq.gz"), val(meta)

    script:
    """
    cat ${fqDir}/${meta.folder}/*.fastq.gz > ${meta.folder}_allreads.fastq.gz
    """
}

process EXTRACT_READIDS_FROM_FASTQ {
    label 'process_low'
    conda params.conda_envs.default_env
    
    input:
    tuple path(reads_combined), val(meta)
    
    output:
    tuple path("${meta.folder}_readids.txt"), val(meta), emit: tuple_ids_meta
    
    script:
    """
    # Extract read IDs from FASTQ file
    # Remove the '@' symbol and take only the read ID (before first space)
    seqkit seq -n ${reads_combined} | sed 's/^@//' | cut -d' ' -f1 > ${meta.folder}_readids.txt
    """
    
    stub:
    """
    touch ${meta.folder}_readids.txt
    """
}



/*
 * Run Kraken2 on fastq files within barcode files
 */

process RUNKRAKEN2 {

    label 'process_low'
    conda params.conda_envs.default_env

    publishDir "results/kraken2", pattern: "*_krakenReport.txt", mode: 'copy', overwrite: false
    
    tag "Kraken on $meta.species_name"

    input:
    tuple path(reads_combined), val(meta)
    path kraken_database

    output:
    path "*_krakenReport.txt", emit: 'kraken_report'
    tuple path(reads_combined), path("*_ids.txt"), val(meta), optional: true, emit: 'reads_ids_meta'

    script:
    """

    kraken2 \
        --db ${kraken_database} \
        --threads ${task.cpus} \
        --gzip-compressed \
        --minimum-hit-groups 1 \
        --output ${meta.species_name}_krakenReport.txt \
        ${reads_combined}

        awk -F'\\t' -v taxID="${meta.species_taxid}" '\$1 == "C" && \$3 == taxID {print \$2}' ${meta.species_name}_krakenReport.txt > ${meta.species_name}_ids.txt.temp

    # Only move the tmp file to the final location if it's not empty
    if [ -s ${meta.species_name}_ids.txt.temp ]; then
        mv ${meta.species_name}_ids.txt.temp ${meta.species_name}_ids.txt
    else
        echo "${meta.species_name}_ids.txt.temp is empty. No file will be outputted."
        # We simply don't create the output file, leading to "nothing" being emitted.
    fi
    """
}

/*
 * Extract a subset fastq containing reads of interest
 */

process SUBSET_FASTQ {

    label 'process_low'

    conda params.conda_envs.default_env

    input:
    tuple path(reads), path(virus_ids), val(meta)

    output:
    tuple path("${meta.species_name}.fastq.gz"), val(meta)

    script:
    """
    seqtk subseq ${reads} ${virus_ids} > ${meta.species_name}.fastq
    pigz -p ${task.cpus} ${meta.species_name}.fastq
    """

}

/*
 * Retrieve reference genomes based on Taxon IDs
 */

process RETRIEVE_GENOMEREFS {

    label 'process_low'

    conda params.conda_envs.default_env

    publishDir "results/references", pattern: "*.fasta", mode: 'copy', overwrite: false
    publishDir "results/references", pattern: "*.log", mode: 'copy', overwrite: false

    input:
    tuple path(subsetted_fastqs), val(meta)

    output:
    tuple path(subsetted_fastqs), path("*.fasta"), val(meta)

    script:
    """
    genome_retriever.py --taxon_id ${meta.species_taxid} --name ${meta.species_name} -o .
    
    awk -v add="${meta.species_taxid}_" '/^>/ {print ">" add substr(\$0, 2); next} {print}' *.fna > ${meta.species_name}.fasta
    """
}



/*
 * Deduplicate reference
 */


process DEDUPLICATE_REFERENCE {

    label 'process_low'

    conda params.conda_envs.default_env

    input:
    path fasta

    output:
    path 'unique.fasta'

    script:
    """
    seqkit rmdup -n ${fasta} > unique.fasta
    """

}

/*
 * Map reads and perform samtools sort on the resulting bam file
 */

process MAP_READS {

    label 'process_low'

    conda params.conda_envs.default_env

    input:
    tuple path(subsetted_fastq), val(meta)
    path fasta

    output:
    tuple path("*.sorted.bam"), val(meta)

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${fasta} ${subsetted_fastq} > mapped.sam
    samtools view -@ ${task.cpus} -b mapped.sam -o mapped.bam
    samtools sort -l 0 -m 1G -O BAM -@ ${task.cpus} mapped.bam -o ${meta.species_name}_${meta.species_taxid}.sorted.bam
    """

}

/*
 * Mark duplicates and index the bam files
 */
process BAM_PROCESSING {

    label 'process_low'

    conda params.conda_envs.default_env

    input:
    tuple path(bam_sorted), val(meta)

    output:
    tuple path ("*.sorted.dedup.bam"), path ("*.sorted.dedup.bam.bai"), val(meta)

    script:
    """
    samtools markdup -r -@ ${task.cpus} ${bam_sorted} ${meta.species_name}.sorted.dedup.bam
    samtools index ${meta.species_name}.sorted.dedup.bam
    """

}

/*
 * Obtain mapping stats on the maping to all genomes
 */

process MAPPING_STATS_ALL {

    label 'process_low'

    conda params.conda_envs.default_env

    publishDir "results/mapping", pattern: "*.flagstat.txt", mode: 'copy', overwrite: false
    publishDir "results/mapping", pattern: "*.idxstats.txt", mode: 'copy', overwrite: false
    
    input:
    tuple path (bam_dedup), path (bam_dedup_index), val(meta)

    output:
    path "*.flagstat.txt", emit: 'flags'
    path "*.idxstats.txt", emit: 'stats'
    tuple path(bam_dedup), path(bam_dedup_index), val(meta), emit: 'tuple_bam_meta'

    script:
    """
    samtools flagstat ${bam_dedup} > ${meta.species_name}.flagstat.txt
    samtools idxstats ${bam_dedup} > ${meta.species_name}.idxstats.txt
    """
}

/*
 * Extract reads that map to the species of interest
 */

process EXTRACT_READIDS {
    
    label 'process_low'

    conda params.conda_envs.default_env

    publishDir "results/mapping", pattern: "*_filtered.bam", mode: 'copy', overwrite: false
    publishDir "results/mapping", pattern: "*_final_read_ids.txt", mode: 'copy', overwrite: false
    publishDir "results/reads", pattern: "*.fastq.gz", mode: 'copy', overwrite: false

    input:
    tuple path (bam_dedup), path (bam_dedup_index), val(meta)

    output:
    tuple path("*_filtered.bam"), val(meta), emit: 'tuple_bam_meta'
    tuple path("*_final_read_ids.txt"), val(meta), emit: 'tuple_ids_meta'
    tuple path("*.fastq.gz"), val(meta), emit: 'filtered_tp_reads'

    script:
    """
    # extract the correct region name. 

    REGION=\$(samtools view -h ${bam_dedup} | \
    grep "^@SQ" | grep "${meta.species_taxid}" | cut -f2 | sed 's/SN://g')

    # use REGION to filter bam wile retaining header
    
    samtools view -b -h -F 256 -F 4 ${bam_dedup} -o ${meta.species_taxid}_${meta.species_name}_filtered.bam \$REGION
    samtools view ${meta.species_taxid}_${meta.species_name}_filtered.bam | cut -f 1 > ${meta.species_taxid}_${meta.species_name}_final_read_ids.txt
    
    # make final fastq with true positives from bam to be pulished 
    samtools fastq -@ ${task.cpus} ${meta.species_taxid}_${meta.species_name}_filtered.bam > ${meta.species_taxid}_${meta.species_name}.fastq
    pigz -p ${task.cpus} ${meta.species_taxid}_${meta.species_name}.fastq
    """
}

/*
 * Run samtools depth on the subsetted bam files
 */

process SAMTOOLS_DEPTH {

    label 'process_low'

    conda params.conda_envs.default_env

    publishDir "results/mapping", pattern: "*.depth", mode: 'copy', overwrite: false

    input:
    tuple path(filtered_bam), val(meta)

    output:
    tuple path("*.depth"), val(meta)
    
    script:
    """
    samtools depth -aa ${filtered_bam} | grep "${meta.species_taxid}_" | sed 's/^/${meta.species_name}_/' > ${meta.species_taxid}_${meta.species_name}.depth
    """

}

/*
 * Calculate depth of coverage statistics
 */

process DEPTH_OF_COVERAGE {
    label 'process_low'

    conda params.conda_envs.default_env

    publishDir "results", pattern: "*.csv", mode: 'copy', overwrite: false

    input:
    path(depthFile)

    output:
    path "coverage_metrics.csv"
    
    script:
    """
    covstats.py ${depthFile} coverage_metrics.csv
    """

}

/*
 * Subset POD5 files based on the target read IDs extracted from the bam file
 */

process SUBSET_POD5 {

    label 'process_low'

    conda params.conda_envs.default_env

    publishDir "results/pod5", pattern: "*.filtered.pod5", mode: 'copy', overwrite: false

    input:
    path pod5_reads
    tuple path(virus_ids), val(meta)

    output:
    tuple path("*.filtered.pod5"), val(meta)

    script:
    """
    # Check if there are any subfolders 
    if find ${pod5_reads} -mindepth 1 -type d | read; then
        # Subfolders are present, run this command
        pod5 filter ${pod5_reads}/${meta.folder} --output ${meta.species_taxid}_${meta.species_name}.filtered.pod5 --ids ${virus_ids} --missing-ok
    else
        # No subfolders present, just query all files
        pod5 filter ${pod5_reads} --output ${meta.species_taxid}_${meta.species_name}.filtered.pod5 --ids ${virus_ids} --missing-ok
    fi

    """

}