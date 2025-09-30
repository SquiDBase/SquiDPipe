#!/usr/bin/env nextflow


/*

    S Q U I D P I P E - N F 
    ===================================

    ABOUT
    -----
    SQUIDPIPE-NF is a Nextflow pipeline designed for processing nanopore sequencing data, 
    focusing on extracting raw signals and reads from specific taxa in barcoded runs. 
    The pipeline utilizes Kraken2 for taxonomic classification, extracts reads corresponding 
    to specified taxa, and performs subsequent mapping and depth-of-coverage analysis. 
    The pipeline is modular and scalable, supporting the integration of custom reference genomes 
    and POD5 file processing.


*/


nextflow.enable.dsl=2

// include 
include { CONCATENATE_FASTQ                                                 } from '../modules.nf'
include { EXTRACT_READIDS_FROM_FASTQ                                        } from '../modules.nf'
include { RUNKRAKEN2                                                        } from '../modules.nf'
include { SUBSET_FASTQ                                                      } from '../modules.nf'
include { RETRIEVE_GENOMEREFS                                               } from '../modules.nf'
include { DEDUPLICATE_REFERENCE                                             } from '../modules.nf'
include { MAP_READS                                                         } from '../modules.nf'
include { BAM_PROCESSING                                                    } from '../modules.nf'
include { MAPPING_STATS_ALL                                                 } from '../modules.nf'
include { EXTRACT_READIDS                                                   } from '../modules.nf'
include { SAMTOOLS_DEPTH                                                    } from '../modules.nf'
include { DEPTH_OF_COVERAGE                                                 } from '../modules.nf'
include { SUBSET_POD5                                                       } from '../modules.nf'



workflow SQUIDPIPE {

    if( params.csvMeta )
        ch_csv_lines = Channel.fromPath(params.csv_file)
                        .splitCsv(header: true, sep: ";")
                        .map { row ->
                            meta = [
                                folder: "${row.filename}",
                                species_name: "${row.species_name}",
                                species_taxid: "${row.species_taxid}",
                                year_of_isolation: "${row.year_of_isolation}",
                                country_of_isolation: "${row.country_of_isolation}",
                                geographic_origin: "${row.geographic_origin}",
                                strain_lineage: "${row.strain_lineage}",
                                source_id: "${row.source_id}",
                                host_taxid: "${row.host_taxid}",
                                internal_lab_id: "${row.internal_lab_id}",
                                diagnostic_method_id: "${row.diagnostic_method_id}",
                                remarks: "${row.remarks}"
                            ]
                        }
    else
        ch_csv_lines = Channel.fromPath(params.csv_file)
                        .splitCsv(header: false)
                        .map { row ->
                            meta = [
                                folder: row[0],
                                species_name: row[1],
                                species_taxid: row[2]
                            ]
                        }
                        // produces:
                            // val(meta)
                            // meta.folder
                            // meta.species_name
                            // meta.taxid

    concatenate_fastq_results = CONCATENATE_FASTQ(params.fastq_dir, ch_csv_lines)
    
    // Branch workflow based on execution mode


    if (params.execution_mode == 'pure_isolate') {
        // PURE ISOLATE MODE: Skip most steps, directly extract read IDs from FASTQ
        
        // Extract read IDs directly from concatenated FASTQ files
        readids_from_fastq = EXTRACT_READIDS_FROM_FASTQ(concatenate_fastq_results)
        
        // If POD5 splitting is enabled, subset POD5 files using extracted read IDs
        if (params.pod5_split) {
            split_pod5_results = SUBSET_POD5(params.pod5_dir, readids_from_fastq)
            
            // Generate final metadata CSV
            split_pod5_results
                .map { file, meta -> 
                    def headers = "filename," + meta.keySet().drop(1).join(",")
                    def values = "${file.getName()}," + meta.values().drop(1).join(",")
                    return [headers, values]
                }
                .flatten()
                .unique()
                .collectFile(name: 'final_metadata.csv', newLine: true, storeDir: 'results/', sort: false)
        }
        
    } else {
        
        // FULL OR MAPPING_ONLY MODES: Original pipeline logic

        // KrakenDB needs to be provided via a channel to ensure container support
        // ch_krakenDB = channel.fromPath(params.databases.kraken_db)

        kraken_results = RUNKRAKEN2(concatenate_fastq_results, params.databases.kraken_db)

        fastq_subset_results = SUBSET_FASTQ(kraken_results.reads_ids_meta)
        // subset fastq emits a tuple of subsetted reads, and virus ids

        // retrieve genome references
    
        ch_retrieved_refs = RETRIEVE_GENOMEREFS(fastq_subset_results)

        // Split the results into two channels: fna_files and fastq_meta

        ch_mapping_input = ch_retrieved_refs.multiMap{ it ->
            fna_files: it[1]
            reads_meta: it[0, 2]
        }

        // combine all downloaded fna files with the user-provided fna file (typically the "decoy" fasta(s) to avoid false-positives)
        ch_extraRef = channel.fromPath(params.references)

        ch_all_refs = ch_mapping_input.fna_files.concat(ch_extraRef)
    
        ch_refGenome = ch_all_refs.collectFile(name: 'combined_reference.fna')

        // deduplicate and  to a value channel convert ch_refGenome to a value channel

        ch_refGenome_value = DEDUPLICATE_REFERENCE(ch_refGenome).first()

        // map reads per sample to the conbined reference

        mapping_output = MAP_READS(ch_mapping_input.reads_meta, ch_refGenome_value)

        bamProcessing_output = BAM_PROCESSING(mapping_output)

        // todo, how to get correct read ids? 
        // todo, how to get mapsats per correct reference?
        // todo - compare to kraken2 outputs
    
        // extracts the readids of those reads that map to the correct reference genome
    
        mapping_stats_output = MAPPING_STATS_ALL(bamProcessing_output)

        extract_readids_output = EXTRACT_READIDS(mapping_stats_output.tuple_bam_meta)

        // we need only the filtered bam and meta

        samtools_depth_output = SAMTOOLS_DEPTH(extract_readids_output.tuple_bam_meta)

        ch_depthFiles = samtools_depth_output.map{ it -> it[0] }
        ch_concat_deptFile = ch_depthFiles.collectFile(name: 'combined.depth')
    
        depth_of_coverage_output = DEPTH_OF_COVERAGE(ch_concat_deptFile)

        if (params.pod5_split) {
        split_pod5_results = SUBSET_POD5(params.pod5_dir, extract_readids_output.tuple_ids_meta)
        split_pod5_results
        .map { file, meta -> 
            def headers = "filename," + meta.keySet().drop(1).join(",")  // Add "filename" as the first header
            def values = "${file.getName()}," + meta.values().drop(1).join(",")
            return [headers, values]
        }
        .flatten()
        .unique()
        .collectFile(name: 'final_metadata.csv', newLine: true, storeDir: 'results/', sort: false)

        }

    }
}