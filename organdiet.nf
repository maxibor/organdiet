#!/usr/bin/env nextflow

// params.reads = "/home/maxime/Documents/ulug_depe/data/*_R{1,2}.fastq.gz"
params.reads = "$baseDir/data/*_R{1,2}.fastq.gz"
params.ctrl = "$baseDir/data/control/*_R{1,2}.fastq.gz"
params.outdir = "$baseDir/results"
params.btindex = "$baseDir/data/db/bowtie/organellome"
scriptdir = "$baseDir/bin/"
py_specie = scriptdir+"process_mapping.py"


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

    cpus = 18
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


    script:
    """
    AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --collapse --threads ${task.cpus}
    fastqc -q *.collapsed
    rename 's/(.collapsed)/\$1.fastq/' *
    """
}

process adapter_removal_ctrl {

    cpus = 18
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

/*
* STEP 3 - Build Bowtie DB of control
*/

process ctr_bowtie_db {
    cpus = 12



    input:
    file(read) from trimmed_ctrl

    output:
    file "ctrl_index*" into ctrl_index


    """
    sed '/^@/!d;s//>/;N' $read > ctrl.fa
    bowtie2-build --threads ${task.cpus} ctrl.fa ctrl_index
    """


}


/*
* STEP 4 - Align on control
*/

process bowtie_align_to_ctrl {

    cpus = 18
    tag "$name"
    publishDir "${params.outdir}/control_removed", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".fastq") > 0)  "./$filename"
        }

    input:
    set val(name), file(reads) from trimmed_reads
    file bt_index from ctrl_index

    output:
    set val(name), file('*.fastq') into fq_unaligned_ctrl_reads
    //something

    script:

    index_base = bt_index.toString().tokenize(' ')[0].tokenize('.')[0]
    sam_out = name+".sam"
    fq_out = name+"_unal.fastq"
    """
    bowtie2 -x $index_base -U $reads --no-sq --threads ${task.cpus} --un $fq_out
    """
}

/*
* STEP 5 - Align on organellome database
*/


process bowtie_align_to_organellome_db {

    cpus = 18
    tag "$name"
    publishDir "${params.outdir}/alignments", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".sam") > 0)  "./$filename"
        }

    input:
    set val(name), file(reads) from fq_unaligned_ctrl_reads

    output:
    set val(name), file('*.sam') into aligned_reads

    script:

    sam_out = name+".sam"
    """
    bowtie2 -x ${params.btindex} -U $reads --end-to-end --threads ${task.cpus} -S $sam_out
    """
}

/*
* STEP 6 - Extract Mapped reads
*/

process extract_mapped_reads {

    tag "$name"
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
* STEP 6 - Get specie composition
*/

process mapped_reads_to_species {
    tag "$name"


    publishDir "${params.outdir}/taxonomic_compositions", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".csv") > 0)  "./$filename"
        }
    input:
        set val(name), file(sam) from mapped_reads
    output:
        set val(name), file("*.csv") into taxo_compo



    script:
        """
        python3 $py_specie $sam
        """

}
