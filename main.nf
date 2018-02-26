#!/usr/bin/env nextflow


def helpMessage() {
    log.info"""
    =========================================
     OrganDiet version ${version}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run maxibor/organdiet --reads '*_R{1,2}.fastq.gz' --btindex db_basename --hgindex db_basename --nrdb db_basename
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --btindex                     Path to organellome database bowtie2 index
      --hgindex                     Path to human genome bowtie2 index

    Options:
      --aligner2                    Specifies the 2nd aligner to nt or nr db (respectively centrifuge or diamond). The proper db associated with aligner2 program must be specified. Defaults to ${params.aligner2}
      --adna                        Specifies if you have ancient dna (true) or modern dna (false). Defaults to ${params.adna}
      --ctrl_index                  Specifies control fastq sequencing data. Must be the same specified the same way as --reads. Defaults to ${params.ctrl}
      --bastamode                   Specifies the mode of LCA for BASTA. Only used if --aligner2 is set to diamond. Defaults to ${params.bastamode}
      --bastaid                     Specifies the identity lower threshold for BASTA LCA. Only used if --aligner2 is set to diamond. Defaults to ${params.bastaid}
      --bastanum                    Specifies the number of hits to retain for BASTA LCA. Only used if --aligner2 is set to diamond. Defaults to ${params.bastanum}
      --trimmingCPU                 Specifies the number of CPU used to trimming/cleaning by AdapterRemoval. Defaults to ${params.trimmingCPU}
      --bowtieCPU                   Specifies the number of CPU used by bowtie2 aligner. Defaults to ${params.bowtieCPU}
      --diamondCPU                  Specifies the number of CPU used by diamond aligner. Only used if --aligner2 is set to diamond. Defaults to ${params.diamondCPU}
      --centrifugeCPU               Specifies the number of CPU used by centrifuge aligner. Only used if --aligner2 is set to centrifuge. Default to ${params.centrifugeCPU}

    References:
      --nrdb                        Path to diamond nr db index. Must be specified if --aligner2 is set to diamond
      --centrifugedb                Path to centrifuge nt db index. Must be specified if --aligner2 is set to centrifuge

    Other options:
      --results                     Name of result directory. Defaults to ${params.results}

    """.stripIndent()
}

//Pipeline version
version = "0.2"

// File locations - default for testing
// params.reads = "$baseDir/data/*_R{1,2}.fastq.gz"
params.reads = "/home/maxime/Documents/data/anna_africa/*_{1,2}.fastq.gz"
// params.ctrl = "$baseDir/data/control/*_R{1,2}.fastq.gz"
params.ctrl = "none"


// Result directory
params.results = "$baseDir/results"

// Script and configurations
params.adna = true
params.multiqc_conf="$baseDir/.multiqc_config.yaml"
params.aligner2 = "diamond"
scriptdir = "$baseDir/bin/"
py_specie = scriptdir+"process_mapping.py"
recentrifuge = scriptdir+"recentrifuge/recentrifuge.py"
basta = scriptdir+"BASTA/bin/basta"
basta2krona = scriptdir+"BASTA/scripts/basta2krona.py"

// Databases locations
params.btindex = "/home/maxime/Documents/db/organellome/bowtie2index/organellome"
params.hgindex = "/home/maxime/Documents/db/hs_genome/Homo_sapiens_Ensembl_GRCh37/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"
params.nrdb = "/home/maxime/Documents/db/diamond/nr"
params.centrifugedb = "/home/maxime/Documents/db/centrifuge/nt/nt"
recentrifugeNodes = scriptdir+"recentrifuge/taxdump"

// BASTA (LCA) parameters
params.bastamode = "majority"
params.bastaid = "97"
params.bastanum = "5"

//CPU parameters
params.trimmingCPU = 12
params.bowtieCPU = 18
params.diamondCPU = 18
params.centrifugeCPU = 18

// Show help emssage
params.help = false
params.h = false
if (params.help || params.h){
    helpMessage()
    exit 0
}



// Header log info
log.info "========================================="
log.info " OrganDiet version ${version}"
log.info "========================================="
def summary = [:]
summary['Reads']        = params.reads
if (params.ctrl) summary['Control']    = params.ctrl
summary['DNA type']    = params.adna ? 'Ancient DNA' : 'Modern DNA'
summary['Organellome database']   = params.btindex
summary['Human genome']     = params.hgindex
summary['Aligner2'] = params.aligner2
if (params.aligner2 == "diamond") summary['Diamond DB'] = params.nrdb
if (params.aligner2 == "centrifuge") summary['Centrifuge DB']  = params.centrifugedb
if (params.aligner2 == "diamond"){
    summary["BASTA mode"] = params.bastamode
    summary["BASTA identity threshold"] = params.bastaid
    summary["BASTA hit threshold"] = params.bastanum
}
summary["CPU for Trimming"] = params.trimmingCPU
summary["CPU for Bowtie2"] = params.bowtieCPU
if (params.aligner2 == "diamond") summary["CPU for diamond"] = params.diamondCPU
if (params.alinger2 == "centrifuge") summary["CPU for centrifuge"] = params.centrifugeCPU
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


Channel
    .fromFilePairs( params.reads, size: 2 )
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}" }
	.into { raw_reads_fastqc; raw_reads_trimming }
if (params.ctrl != "none"){
    Channel
        .fromFilePairs(params.ctrl, size: 2)
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.ctrl}"}
        .into { raw_ctrl_fastqc; raw_ctrl_trimming }

}


/*
* STEP 1 - FastQC
*/
process fastqc {
    // cache false
    tag "$name"
    publishDir "${params.results}/fastqc", mode: 'copy',
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

if (params.ctrl != "none"){
    process fastqc_control {
        // cache false
        tag "$name"
        publishDir "${params.results}/fastqc", mode: 'copy',
           saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

        input:
            set val(name), file(reads) from raw_ctrl_fastqc

        script:
            """
            fastqc -q $reads
            """
    }
}


/*
 * STEP 2 - AdapterRemoval
 */


if (params.adna == true){
    process adapter_removal_ancient_dna {

        cpus = params.trimmingCPU
        // cache false
        tag "$name"
        publishDir "${params.results}/trimmed", mode: 'copy'

        input:
            set val(name), file(reads) from raw_reads_trimming

        output:
            set val(name), file('*.truncated.fastq') into truncated_reads
            set val(name), file('*.collapsed.fastq') into collapsed_reads
            set val(name), file("*.settings") into adapter_removal_results
            file '*_fastqc.{zip,html}' into fastqc_results_after_trim



        script:
            out1 = name+".pair1.truncated.fastq"
            out2 = name+".pair2.truncated.fastq"
            col_out = name+".collapsed.fastq"
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --collapse --output1 $out1 --output2 $out2 --outputcollapsed $col_out --threads ${task.cpus}
            fastqc -q *.collapsed*q
            """
    }

    if (params.ctrl != "none"){
        process adapter_removal_ctrl_ancient_dna {

            cpus = params.trimmingCPU
            // cache false
            tag "$name"
            publishDir "${params.results}/trimmed", mode: 'copy'

            input:
                set val(name), file(reads) from raw_ctrl_trimming

            output:
                set val(name), file('*.collapsed.fastq') into collapsed_reads_ctrl
                set val(name), file(out1), file(out2) into truncated_reads_ctrl


            script:
                out1 = name+".pair1.truncated.fastq"
                out2 = name+".pair2.truncated.fastq"
                col_out = name+".collapsed.fastq"
                """
                AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --collapse --output1 $out1 --output2 $out2 --outputcollapsed $col_out --threads ${task.cpus}
                """
        }
    }
} else {
    process adapter_removal_modern_dna {

        cpus = params.trimmingCPU
        // cache false
        tag "$name"
        publishDir "${params.results}/trimmed", mode: 'copy'

        input:
            set val(name), file(reads) from raw_reads_trimming

        output:
            // set val(name), file('*.truncated.fastq') into truncated_reads
            set val(name), file(out1), file(out2) into truncated_reads
            set val(name), file("*.settings") into adapter_removal_results
            file '*_fastqc.{zip,html}' into fastqc_results_after_trim



        script:
            out1 = name+".pair1.truncated.fastq"
            out2 = name+".pair2.truncated.fastq"
            """
            AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --output1 $out1 --output2 $out2 --threads ${task.cpus}
            fastqc -q *.truncated*
            """
    }

    if (params.ctrl != "none"){
        process adapter_removal_ctrl_modern_dna {

            cpus = params.trimmingCPU
            // cache false
            tag "$name"
            publishDir "${params.results}/trimmed", mode: 'copy'

            input:
                set val(name), file(reads) from raw_ctrl_trimming

            output:
                set val(name), file(out1), file(out2) into truncated_reads_ctrl

            script:
                out1 = name+".pair1.truncated.fastq"
                out2 = name+".pair2.truncated.fastq"
                """
                AdapterRemoval --basename $name --file1 ${reads[0]} --file2 ${reads[1]} --trimns --trimqualities --output1 $out1 --output2 $out2 --threads ${task.cpus}
                """
        }
    }
}





/*
* STEP 3 - Build Bowtie DB of control
*/
if (params.adna == true){
    if (params.ctrl != "none"){
        process ctr_bowtie_db_ancient_dna {
            cpus = params.bowtieCPU
            // cache false

            input:
                set val(name), file(col_read) from collapsed_reads_ctrl

            output:
                file "ctrl_index*" into ctrl_index

            script:
                """
                sed '/^@/!d;s//>/;N' $col_read > ctrl.fa
                bowtie2-build --threads ${task.cpus} ctrl.fa ctrl_index
                """
        }
    }
} else {
    if (params.ctrl != "none"){
        process ctr_bowtie_db_modern_dna {
            cpus = params.bowtieCPU
            // cache false

            input:
                set val(name), file(trun_read1), file(trun_read2) from truncated_reads_ctrl

            output:
                file "ctrl_index*" into ctrl_index

            script:
                merge_file = name+"_merged.fq"
                """
                cat ${trun_read1} > $merge_file
                cat ${trun_read2} >> $merge_file
                sed '/^@/!d;s//>/;N' $merge_file > ctrl.fa
                bowtie2-build --threads ${task.cpus} ctrl.fa ctrl_index
                """

        }
    }

}





/*
* STEP 4 - Align on control, output unaligned reads
*/

if (params.adna == true){
    if (params.ctrl != "none"){
        process bowtie_align_to_ctrl_ancient_dna {

            cpus = params.bowtieCPU
            // cache false
            tag "$name"
            publishDir "${params.results}/control_removed", mode: 'copy',
                saveAs: {filename ->
                    if (filename.indexOf(".fastq") > 0)  "./$filename"
                }

            input:
                set val(name), file(col_reads) from collapsed_reads
                set val(name), file(trun_read1), file(trun_read2) from truncated_reads
                file bt_index from ctrl_index.collect()

            output:
                set val(name), file('*.unal.fastq') into fq_unaligned_ctrl_reads
                file("*.metrics") into ctrl_aln_metrics

            script:
                sam_out = name+".sam"
                fq_out = name+".unal.fastq"
                metrics = name+".metrics"
                """
                bowtie2 -x ctrl_index -U $col_reads --no-sq --threads ${task.cpus} --un $fq_out 2> $metrics
                """

        }
    }
} else {
    if (params.ctrl != "none"){
        process bowtie_align_to_ctrl_modern_dna {

            cpus = params.bowtieCPU
            // cache false
            tag "$name"
            publishDir "${params.results}/control_removed", mode: 'copy',
                saveAs: {filename ->
                    if (filename.indexOf(".fastq") > 0)  "./$filename"
                }

            input:
                set val(name), file(trun_read1), file(trun_read2) from truncated_reads
                file bt_index from ctrl_index.collect()

            output:
                // set val(name), file('*.unal.fastq') into fq_unaligned_ctrl_reads
                file("*.metrics") into ctrl_aln_metrics
                set val(name), file(out1), file(out2) into fq_unaligned_ctrl_reads

            script:
                sam_out = name+".sam"
                fq_out = name+".unal.fastq"
                out1=name+".unal.1.fastq"
                out2=name+".unal.2.fastq"
                metrics = name+".metrics"
                """
                bowtie2 -x ctrl_index -1 $trun_read1 -2 $trun_read2 --no-sq --threads ${task.cpus} --un-conc $fq_out 2> $metrics
                """

        }
    }

}





/*
* STEP 5 - Align on human genome, output unaligned reads
*/

if (params.ctrl != "none"){
    if (params.adna == true){
        process bowtie_align_to_human_genome_from_ctrl_ancient_dna {

            cpus = params.bowtieCPU
            // cache false
            tag "$name"
            publishDir "${params.results}/human_removed", mode: 'copy',
                saveAs: {filename ->
                    if (filename.indexOf(".fastq") > 0)  "./$filename"
                }

            input:
                set val(name), file(reads) from fq_unaligned_ctrl_reads

            output:
                set val(name), file('*.human_unal.fastq') into fq_unaligned_human_reads
                file("*.metrics") into human_aln_metrics

            script:
                fq_out = name+".human_unal.fastq"
                metrics = name+".metrics"
                """
                bowtie2 -x ${params.hgindex} -U $reads --no-sq --threads ${task.cpus} --un $fq_out 2> $metrics
                """

        }
    } else {
        process bowtie_align_to_human_genome_from_ctrl_modern_dna {

            cpus = params.bowtieCPU
            // cache false
            tag "$name"
            publishDir "${params.results}/human_removed", mode: 'copy',
                saveAs: {filename ->
                    if (filename.indexOf(".fastq") > 0)  "./$filename"
                }

            input:
                // set val(name), file(reads) from fq_unaligned_ctrl_reads
                set val(name), file(trun_read1), file(trun_read2) from fq_unaligned_ctrl_reads

            output:
                // set val(name), file('*.human_unal.fastq') into fq_unaligned_human_reads
                set val(name), file(out1), file(out2) into fq_unaligned_human_reads
                file("*.metrics") into human_aln_metrics

            script:
                fq_out = name+".human_unal.fastq"
                out1 = name+".human_unal.1.fastq"
                out2 = name+".human_unal.2.fastq"
                metrics = name+".metrics"
                """
                bowtie2 -x ${params.hgindex} -1 $trun_read1 -2 $trun_read2 --no-sq --threads ${task.cpus} --un-conc $fq_out 2> $metrics
                """

        }
    }

} else {
    if (params.adna == true){
        process bowtie_align_to_human_genome_no_control_ancient_dna {

            cpus = params.bowtieCPU
            // cache false
            tag "$name"
            publishDir "${params.results}/human_removed", mode: 'copy',
                saveAs: {filename ->
                    if (filename.indexOf(".fastq") > 0)  "./$filename"
                }

            input:
            set val(name), file(col_reads) from collapsed_reads
            set val(name), file(trun_read1), file(trun_read2) from truncated_reads
            // set val(name), file(trun_reads) from truncated_reads.fromFilePairs("*.pair{1,2}.fastq")

            output:
                set val(name), file('*.human_unal.fastq') into fq_unaligned_human_reads
                file("*.metrics") into human_aln_metrics

            script:
                fq_out = name+".human_unal.fastq"
                metrics = name+".metrics"
                """
                bowtie2 -x ${params.hgindex} -U $col_reads --no-sq --threads ${task.cpus} --un $fq_out 2> $metrics
                """


        }
    } else {
        process bowtie_align_to_human_genome_no_control_modern_dna {

            cpus = params.bowtieCPU
            // cache false
            tag "$name"
            publishDir "${params.results}/human_removed", mode: 'copy',
                saveAs: {filename ->
                    if (filename.indexOf(".fastq") > 0)  "./$filename"
                }

            input:
            set val(name), file(trun_read1), file(trun_read2) from truncated_reads
            // set val(name), file(trun_reads) from truncated_reads.fromFilePairs("*.pair{1,2}.fastq")

            output:
                // set val(name), file('*.human_unal.fastq') into fq_unaligned_human_reads
                set val(name), file(out1), file(out2) into fq_unaligned_human_reads
                file("*.metrics") into human_aln_metrics

            script:
                fq_out = name+".human_unal.fastq"
                out1 = name+".human_unal.1.fastq"
                out2 = name+".human_unal.2.fastq"
                metrics = name+".metrics"
                """
                bowtie2 -x ${params.hgindex} -1 $trun_read1 -2 $trun_read2 --no-sq --threads ${task.cpus} --un-conc $fq_out 2> $metrics
                """
        }
    }

}







/*
* STEP 6 - Align on organellome database
*/

if (params.adna == true ){
    process bowtie_align_to_organellome_db_ancient_dna {

        cpus = params.bowtieCPU
        // cache false
        tag "$name"
        publishDir "${params.results}/alignments", mode: 'copy',
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
} else {
    process bowtie_align_to_organellome_db_modern_dna {

        cpus = params.bowtieCPU
        // cache false
        tag "$name"
        publishDir "${params.results}/alignments", mode: 'copy',
            saveAs: {filename ->
                if (filename.indexOf(".sam") > 0)  "./$filename"
            }

        input:
            // set val(name), file(reads) from fq_unaligned_human_reads
            set val(name), file(read1), file(read2) from fq_unaligned_human_reads

        output:
            set val(name), file('*.sam') into aligned_reads
            file("*.metrics") into organellome_aln_metrics

        script:
            sam_out = name+".sam"
            metrics = name+".metrics"
            """
            bowtie2 -x ${params.btindex} -1 $read1 -2 $read2 --end-to-end --threads ${task.cpus} -S $sam_out -a 2> $metrics
            """
    }
}


/*
* STEP 7 - Extract Mapped reads
*/

process extract_mapped_reads {

    tag "$name"
    // cache false
    publishDir "${params.results}/alignments", mode: 'copy',
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

if (params.aligner2 == "diamond"){
    process diamond_align_to_nr {
        tag "$name"
        cpus = params.diamondCPU


        publishDir "${params.results}/nr_alignment", mode: 'copy',
            saveAs: {filename ->  "./$filename"}

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
} else if (params.aligner2 == "centrifuge"){

    process centrifuge_align_to_nr{
        tag "$name"
        cpus = params.centrifugeCPU

        publishDir "${params.results}/nr_alignment", mode: 'copy',
            saveAs: {filename ->  "./$filename"}

        input:
            set val(name), file(best_fa) from best_match

        output:
            file("*.centrifuge.out") into nr_aligned

        script:
            centrifuge_out = name+".centrifuge.out"
            centrifuge_report = name+"_centrifuge_report.tsv"
            """
            centrifuge -x ${params.centrifugedb} -r $best_fa -p ${task.cpus} -f --report-file $centrifuge_report -S $centrifuge_out
            """
    }

    process recentrifuge {
        tag "${centrifuge_aligned[0].baseName}"

        beforeScript "set +u; source activate py36"
        afterScript "set +u; source deactivate py36"

         publishDir "${params.results}/krona", mode: 'copy',
            saveAs: {filename ->  "./$filename"}

        input:
            file centrifuge_aligned from nr_aligned.toList()

        output:
            file("recentrifuge_result.html") into recentrifuge_result

        script:
            allfiles = centrifuge_aligned.join(" -f ")
            """
            $recentrifuge -f $allfiles -n $recentrifugeNodes -o recentrifuge_result.html
            """
    }
}


if (params.aligner2 == "diamond"){

    /*
    * STEP 10 - Assign LCA
    */

    process lca_assignation {
        tag "$name"


        publishDir "${params.results}/taxonomy", mode: 'copy',
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
            $basta sequence $sorted_nr $basta_name prot -t ${params.bastamode} -m 1 -n ${params.bastanum} -i ${params.bastaid}
            """
    }


    /*
    * STEP 11 - Generate Krona output
    */

    process visual_results {
        tag "$name"


        publishDir "${params.results}/krona", mode: 'copy',
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
            $basta2krona $basta_res $krona_out
            """
    }
}


/*
* STEP 12 - Generate run summary
*/
if (params.adna == true){
    process multiqc_no_control {
        tag "$prefix"
        // cache false
        publishDir "${params.results}/MultiQC", mode: 'copy'

        // beforeScript "set +u; source activate py27"
        // afterScript "set +u; source deactivate py27"

        input:
            file (fastqc:'fastqc_before_trimming/*') from fastqc_results.collect()
            file ('adapter_removal/*') from adapter_removal_results.collect()
            file("fastqc_after_trimming/*") from fastqc_results_after_trim.collect()
            file('aligned_to_human/*') from human_aln_metrics.collect()
            file('aligned_to_organellomeDB/*') from organellome_aln_metrics.collect()

        output:
            file '*multiqc_report.html' into multiqc_report
            file '*_data' into multiqc_data
            file '.command.err' into multiqc_stderr

        script:
            prefix = fastqc[0].toString() - '_fastqc.html' - 'fastqc/'
            """
            multiqc -f -d fastqc_before_trimming adapter_removal fastqc_after_trimming aligned_to_human aligned_to_organellomeDB -c ${params.multiqc_conf}
            """

    }
} else {
    process multiqc_with_control {
        tag "$prefix"
        // cache false
        publishDir "${params.results}/MultiQC", mode: 'copy'

        // beforeScript "set +u; source activate py27"
        // afterScript "set +u; source deactivate py27"

        input:
            file (fastqc:'fastqc_before_trimming/*') from fastqc_results.collect()
            file ('adapter_removal/*') from adapter_removal_results.collect()
            file("fastqc_after_trimming/*") from fastqc_results_after_trim.collect()
            file('aligned_to_human/*') from human_aln_metrics.collect()
            file('aligned_to_organellomeDB/*') from organellome_aln_metrics.collect()
            file('aligned_to_blank/*') from ctrl_aln_metrics.collect()

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
}
