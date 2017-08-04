#!/usr/bin/env nextflow

/*
*params input 
*/

params.bam = "$baseDir/data/*.bam"
params.output_dir = "data" 
params.annotation  ="$baseDir/data/*.gtf"
params.type_data = "Illumina"
params.compare = null
params.help = false


//print usage
if (params.help) {
    log.info ''
    log.info 'featureCounts '
    log.info '-----------------------'
    log.info '.'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow feauturecount.nf --compare PATH/file.bam --type_data Illumina'
    log.info '            --output_dir toto --annotation file.gtf --bam  "PATH/*.bam"   '
    log.info 'Or '
    log.info '    nextflow feauturecount.nf --type_data Illumina   --output_dir toto    '
    log.info '             --annotation PATH/file.gtf   --bam "PATH/*.bam"                    '
    log.info ''
    log.info 'Options:'
    log.info '    --help                              Show this message and exit.'
    log.info '    --type_data                         Solid, Illumina or Proton data [default : Illumina]'
    log.info '    --annotation                        File of genome annotation in GTF .'
    log.info '    --output_dir                        Name Directorie of output result.'
    log.info '    --bam                               Input file BAM [ex :"././*.bam"].'
    log.info '    --compare                           Input file BAM . Option for compare BAM and return plot [need R]  .'
        
    exit 1
}

/*
*Create channel by params. input 
*/
file_annotation = file(params.annotation)

file_mapping =  Channel.fromPath(params.bam).map { file -> tuple(file.baseName, file) }


/*
* Parameter of featureCount analyse in fonction of data
*/
//parametre -f et R supprimer par rapport au script precedant 

if(params.type_data == "Solid" ){
    parameter = "-p -M -O --largestOverlap -s 1 "
}
if(params.type_data == "Illumina" ){
    parameter = "-p -M -O --largestOverlap -s 2 "
}
if(params.type_data == "Proton" ){
    parameter = "-M -O --largestOverlap -s 1 "
}



/*
*This program executes a refined paired-end with two files SAM  
*
*OPTION :
* 
*-T     Number of the threads. 1 by default.
*-p     Count fragments (read pairs) instead of individual reads
*-M     Multi-mapping reads will also be counted. For a multi-mapping read,
*       all its reported alignments will be counted
*-largestOverlap    Assign reads to a meta-feature/feature that has the largest number of overlapping bases.
*-t     Specify feature type in GTF annotation. `exon' by default.
*-g     Specify attribute type in GTF annotation. `gene_id' by default.
*-a     Name of an annotation file. GTF format by default.
*-o     Name of the output file including read counts.
*-O     Assign reads to all their overlapping meta-features
*-s     Perform strand-specific read counting. Possible values:  
*           0 (unstranded), 1 (stranded) and 2 (reversely stranded).
*/

process featureCounts {
    cpus 2
    tag{id}
    publishDir "result/featureCounts/$params.output_dir/", mode: "copy"

    input: 
    set id , file (mapping) from file_mapping
    each anno from file_annotation

    output: 

    file "${id}" into count 
    file "*.summary" into summary 

    """
    featureCounts -T ${task.cpus} \
    ${parameter} \
    -t exon -g gene_id \
    -a ${anno} \
    -o ${id} ${mapping} \
    """
}

/*
*Execute this option if you will compare BAM 
* 
*/

if(params.compare != null){


    bam_compare =  Channel.fromPath(params.compare).map { file -> tuple(file.baseName, file) }

    process featureCounts_compare {
        cpus 2
        tag{id}
        publishDir "result/featureCounts/$params.output_dir/", mode: "copy"

            input: 
            set id , file (mapping) from bam_compare
            each anno from file_annotation

            output: 

            file "${id}" into count_compare 
            file "*.summary" into summary_compare 

            """
            featureCounts -T ${task.cpus} \
            ${parameter} \
            -t exon -g gene_id \
            -a ${anno} \
            -o ${id} ${mapping} \
            """

        }

    /*
    *Process call script R for plot count gene in different alignement    
    */    

    process plot { 
        publishDir "result/featureCounts/$params.output_dir/", mode: "copy"
        
        input: 
        file count_1 from count
        each count_2 from count_compare


        output: 
        file "*plot" into plot 

        """
        Rscript '$baseDir/plotcount.R' ${count_1} ${count_2}
        """
        }
    }