#!/usr/bin/env nextflow

/*
*params input 
*/

params.file_bam = "$baseDir/data/*.bam"
params.path_featurecounts  ="/bin/featureCqounts" 
params.name_dir  ="data" 
params.annotation  ="$baseDir/data/*.gtf"
params.type_data = "Illumina"
params.file_bam_compare = null
params.help = false


//print usage
if (params.help) {
    log.info ''
    log.info 'featureCounts '
    log.info '-----------------------'
    log.info '.'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow feauturecount.nf --file_bam_compare ../file.bam --type_data Illumina'
    log.info '            --name_dir toto --annotation file.gtf --path_featurecounts /bin/featureCqounts '
    log.info ''
    log.info 'Options:'
    log.info '    --help                              Show this message and exit.'
    log.info '    --type_data                         Solid or Illumina data [default : Illumina]'
    log.info '    --annotation                        File of genome annotation in GTF .'
    log.info '    --name_dir                          Name Dir.'
    log.info '    --path_featurecounts                Path of your tool featureCounts [default : /bin/featureCqounts] .'
    log.info '    --file_bam                          Input file BAM [ex :"././*.bam"].'
    log.info '    --file_bam_compare                  Input file BAM . Option for compare BAM and return plot [need R]  .'
        
    exit 1
}


file_annotation = file(params.annotation)

path_featurecounts= file(params.path_featurecounts)

file_mapping =  Channel.fromPath(params.file_bam).map { file -> tuple(file.baseName, file) }




if(params.type_data == "Solid" ){
    parameter = 1
}
if(params.type_data == "Illumina" ){
    parameter = 2
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

process featureCounts_Solid {
    cpus 2
    tag{id}
    publishDir "result/featureCounts/$params.name_dir/", mode: "copy"

    input: 
    set id , file (mapping) from file_mapping
    each anno from file_annotation

    output: 

    file "${id}" into count 
    file "*.summary" into summary 

    """
    featureCounts -T ${task.cpus} \
    -p -M -O --largestOverlap -s ${parameter} \
    -t exon -g gene_id \
    -a ${anno} \
    -o ${id} ${mapping} \
    """
}

/*
*
*Execute this option if you will compare two BAM 
*First step execute the programme featureCounts 
*
*/

if(params.file_bam_compare != null){


    bam_compare =  Channel.fromPath(params.file_bam_compare).map { file -> tuple(file.baseName, file) }

    process featureCounts_Solid {
        cpus 4
        tag{id}

            input: 
            set id , file (mapping) from bam_compare
            each anno from file_annotation

            output: 

            file "${id}" into count_compare 
            file "*.summary" into summary_compare 

            """
            featureCounts -T ${task.cpus} \
            -p -M -O --largestOverlap -s ${parameter} \
            -t exon -g gene_id \
            -a ${anno} \
            -o ${id} ${mapping} \
            """

        }

    process plot { 
        publishDir "result/featureCounts/$params.name_dir/", mode: "copy"
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