#!/usr/bin/env nextflow

// File locations - default for testing
params.reads = "$baseDir/data/*_R{1,2}.fastq.gz"
// params.ctrl = "$baseDir/data/control/*_R{1,2}.fastq.gz"
params.ctrl = "none"


// Result directory
params.outdir = "$baseDir/results"

// Script and configurations
params.multiqc_conf="$baseDir/.multiqc_config.yaml"
scriptdir = "$baseDir/bin/"
py_specie = scriptdir+"process_mapping.py"

// Databases locations
params.btindex = "$baseDir/data/db/bowtie/organellome"
params.hgindex = "/mnt/ntfs/databases/hg19/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"
params.nrdb = "/mnt/ntfs/databases/nr_diamond/nr"

// BASTA (LCA) parameters
params.bastamode = "majority"
params.bastaid = "97"
params.bastanum = "5"


nthreads = 24

Channel
    .fromFilePairs( params.reads, size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
	.into { raw_reads_fastqc; raw_reads_trimming }

Channel
    .fromFilePairs(params.ctrl, size: 2)
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.ctrl}"}
    .into { raw_ctrl_fastqc; raw_ctrl_trimming }


/*
* STEP 1 - FastQC
*/
process fastqc {
    // cache false
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
       saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
        set val(name), file(reads) from raw_reads_fastqc

    output:
        file '*_fastqc.{zip,html}' into fastqc_results
        file '.command.out' into fastqc_stdout

    script:
        """
        fastqc -q $reads
        """
}


process fastqc_control {
    // cache false
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
       saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
        set val(name), file(reads) from raw_ctrl_fastqc

    script:
        """
        fastqc -q $reads
        """
}

/*
 * STEP 2 - AdapterRemoval
 */

process adapter_removal {

    cpus = 6
    // cache false
    tag "$name"
    publishDir "${params.outdir}/trimmed", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf(".settings") > 0) "logs/$filename"
        }

    input:
        set val(name), file(reads) from raw_reads_trimming

    output:
        set val(name), file('*.collapsed.fastq') into trimmed_reads
        file '*_fastqc.{zip,html}' into fastqc_results_after_trim
        file '*.settings' into adapter_removal_results

    script:
        """
        AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --collapse --threads ${task.cpus}
        fastqc -q *.collapsed
        rename 's/(.collapsed)/\$1.fastq/' *
        """
}

if (params.ctrl != "none"){
    process adapter_removal_ctrl {

        cpus = 8
        // cache false
        tag "$name"
        publishDir "${params.outdir}/trimmed", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
                else if (filename.indexOf(".settings") > 0) "logs/$filename"
            }

        input:
            set val(name), file(reads) from raw_ctrl_trimming

        output:
            set val(name), file('*.collapsed.fastq') into trimmed_ctrl

        script:
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --collapse --threads ${task.cpus}
            fastqc -q *.collapsed
            rename 's/(.collapsed)/\$1.fastq/' *
            """
    }
}


/*
* STEP 3 - Build Bowtie DB of control
*/

if (params.ctrl != "none"){
    process ctr_bowtie_db {
        cpus = 12
        // cache false


        input:
            file(read) from trimmed_ctrl

        output:
            file "ctrl_index*" into ctrl_index
            """
            sed '/^@/!d;s//>/;N' $read > ctrl.fa
            bowtie2-build --threads ${task.cpus} ctrl.fa ctrl_index
            """
    }
}



/*
* STEP 4 - Align on control, output unaligned reads
*/

if (params.ctrl != "none"){
    process bowtie_align_to_ctrl {

        cpus = 18
        // cache false
        tag "$name"
        publishDir "${params.outdir}/control_removed", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".fastq") > 0)  "./$filename"
            }

        input:
            set val(name), file(reads) from trimmed_reads
            file bt_index from ctrl_index.collect()

        output:
            set val(name), file('*.fastq') into fq_unaligned_ctrl_reads
            file("*.metrics") into ctrl_aln_metrics

        script:
            sam_out = name+".sam"
            fq_out = name+"_unal.fastq"
            metrics = name+".metrics"
            """
            bowtie2 -x ctrl_index -U $reads --no-sq --threads ${task.cpus} --un $fq_out 2> $metrics
            """
    }
}


/*
* STEP 5 - Align on human genome, output unaligned reads
*/

if (params.ctrl != "none"){
    process bowtie_align_to_human_genome {

        cpus = 18
        // cache false
        tag "$name"
        publishDir "${params.outdir}/human_removed", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".fastq") > 0)  "./$filename"
            }

        input:
            set val(name), file(reads) from fq_unaligned_ctrl_reads

        output:
            set val(name), file('*.fastq') into fq_unaligned_human_reads
            file("*.metrics") into human_aln_metrics

        script:
            fq_out = name+"_human_unal.fastq"
            metrics = name+".metrics"
            """
            bowtie2 -x ${params.hgindex} -U $reads --no-sq --threads ${task.cpus} --un $fq_out 2> $metrics
            """
    }
} else {
    process bowtie_align_to_human_genome {

        cpus = 18
        // cache false
        tag "$name"
        publishDir "${params.outdir}/human_removed", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".fastq") > 0)  "./$filename"
            }

        input:
            set val(name), file(reads) from trimmed_reads

        output:
            set val(name), file('*.fastq') into fq_unaligned_human_reads
            file("*.metrics") into human_aln_metrics

        script:
            fq_out = name+"_human_unal.fastq"
            metrics = name+".metrics"
            """
            bowtie2 -x ${params.hgindex} -U $reads --no-sq --threads ${task.cpus} --un $fq_out 2> $metrics
            """
    }
}







/*
* STEP 6 - Align on organellome database
*/


process bowtie_align_to_organellome_db {

    cpus = 18
    // cache false
    tag "$name"
    publishDir "${params.outdir}/alignments", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".sam") > 0)  "./$filename"
        }

    input:
        set val(name), file(reads) from fq_unaligned_human_reads

    output:
        set val(name), file('*.sam') into aligned_reads
        file("*.metrics") into organellome_aln_metrics

    script:
        sam_out = name+".sam"
        metrics = name+".metrics"
        """
        bowtie2 -x ${params.btindex} -U $reads --end-to-end --threads ${task.cpus} -S $sam_out -a 2> $metrics
        """
}

/*
* STEP 7 - Extract Mapped reads
*/

process extract_mapped_reads {

    tag "$name"
    // cache false
    publishDir "${params.outdir}/alignments", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".mapped.sam") > 0)  "./$filename"
        }

    input:
        set val(name), file(align) from aligned_reads

    output:
        set val(name), file('*.mapped.sam') into mapped_reads

    script:
        mapped_out = name+".mapped.sam"
        """
        samtools view -S -F4 $align > $mapped_out
        """
}

/*
* STEP 8 - Filter reads on quality and length
*/

process extract_best_reads {
    tag "$name"

    input:
        set val(name), file(sam) from mapped_reads

    output:
        set val(name), file("*.best.aligned.fa") into best_match

    script:
        """
        python $py_specie $sam
        """
}

/*
* STEP 9 - Align on NR database
*/

process diamond_align_to_nr {
    tag "$name"
    cpus = 18


    publishDir "${params.outdir}/nr_alignment", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".diamond.out") > 0)  "./$filename"
        }

    input:
        set val(name), file(best_fa) from best_match

    output:
        set val(name), file("*.diamond.out") into nr_aligned

    script:
        diamond_out = name+".diamond.out"
        """
        diamond blastx -d ${params.nrdb} -q $best_fa -o $diamond_out -f 6 -p ${task.cpus}
        """
}

/*
* STEP 10 - Assign LCA
*/

process lca_assignation {
    tag "$name"


    publishDir "${params.outdir}/taxonomy", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".basta.out") > 0)  "./$filename"
        }

    // beforeScript "set +u; source activate py27"
    // afterScript "set +u; source deactivate py27"

    input:
        set val(name), file(aligned_nr) from nr_aligned

    output:
        set val(name), file("*.basta.out") into lca_result

    script:
        basta_name = name+".basta.out"
        sorted_nr = name+"_diamond_nr.sorted"
        """
        sort -k3 -r -n $aligned_nr > $sorted_nr
        basta sequence $sorted_nr $basta_name prot -t ${params.bastamode} -m 1 -n ${params.bastanum} -i ${params.bastaid}
        """
}


/*
* STEP 11 - Generate Krona output
*/

process visual_results {
    tag "$name"


    publishDir "${params.outdir}/krona", mode: 'copy',
        saveAs: {filename ->  "./$filename"}

    // beforeScript "set +u; source activate py27"
    // afterScript "set +u; source deactivate py27"

    input:
        set val(name), file(basta_res) from lca_result

    output:
        set val(name), file("*.krona.html") into krona_res

    script:
        krona_out = name+".krona.html"
        """
        basta2krona.py $basta_res $krona_out
        """
}

/*
* STEP 12 - Generate run summary
*/

process multiqc {
    tag "$prefix"
    // cache false
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    // beforeScript "set +u; source activate py27"
    // afterScript "set +u; source deactivate py27"

    input:
        file (fastqc:'fastqc_before_trimming/*') from fastqc_results.collect()
        file ('adapter_removal/*') from adapter_removal_results.collect()
        file("fastqc_after_trimming/*") from fastqc_results_after_trim.collect()
        file('aligned_to_blank/*') from ctrl_aln_metrics.collect()
        file('aligned_to_human/*') from human_aln_metrics.collect()
        file('aligned_to_organellomeDB/*') from organellome_aln_metrics.collect()

    output:
        file '*multiqc_report.html' into multiqc_report
        file '*_data' into multiqc_data
        file '.command.err' into multiqc_stderr

    script:
        prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
        """
        multiqc -f -d fastqc_before_trimming adapter_removal fastqc_after_trimming aligned_to_blank aligned_to_human aligned_to_organellomeDB -c ${params.multiqc_conf}
        """
}
