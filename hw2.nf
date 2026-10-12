nextflow.enable.dsl=2


/*
 * Command-line parameters
 */
params.vcf = null
params.index = null
params.bed = null



/*
 * Extract one genomic region from the input VCF.
 *
 * One instance of this process should run for each
 * region in the BED file.
 */
process EXTRACT_REGION {

    tag "${prefix}"

    publishDir 'results/regions', mode: 'copy'


    input:

    /*
     * TODO 1:
     *
     * Declare ONE tuple input containing:
     *
     *   vcf      = VCF file
     *   index    = VCF index file
     *   region   = ordinary value
     *   prefix   = ordinary value
     *
     * Remember:
     *   path(...) is used for files
     *   val(...)  is used for ordinary values
     */

    // TODO 1: add the tuple input declaration here

    tuple path(vcf), path(index), val(region), val(prefix)

    output:

    /*
     * TODO 2:
     *
     * Output ONE tuple containing:
     *
     *   prefix
     *   "${prefix}.vcf"
     *
     * This tuple will become the input to MAKE_PLINK.
     */

    // TODO 2: add the tuple output declaration here
    tuple path("${prefix}.vcf"), val(prefix)

    script:

    """
    bcftools view \
        -r ${region} \
        ${vcf} \
        -Ov \
        -o ${prefix}.vcf
    """
}



/*
 * Convert one regional VCF into a binary PLINK dataset.
 *
 * There should be one MAKE_PLINK task corresponding
 * to every regional VCF produced by EXTRACT_REGION.
 */
process MAKE_PLINK {

    tag "${prefix}"

    publishDir 'results/regions', mode: 'copy'


    input:

    /*
     * TODO 3:
     *
     * Declare the tuple received from EXTRACT_REGION.
     *
     * It contains:
     *
     *   prefix = ordinary value
     *   vcf    = VCF file
     */

    // TODO 3: add the tuple input declaration here
    tuple val(prefix), path(vcf)

    output:

    /*
     * TODO 4:
     *
     * Output ONE tuple containing:
     *
     *   prefix
     *   "${prefix}.bed"
     *   "${prefix}.bim"
     *   "${prefix}.fam"
     */

    // TODO 4: add the tuple output declaration here
    tuple val(prefix), path("${prefix}.bed"), path("${prefix}.bim"), path("${prefix}.fam")

    script:

    """
    plink \
        --vcf ${vcf} \
        --make-bed \
        --out ${prefix}
    """
}



/*
 * Merge all of the regional PLINK datasets.
 *
 * Unlike the two processes above, this process should
 * run only ONCE, after all regional PLINK datasets exist.
 */
process MERGE_PLINK {

    publishDir 'results', mode: 'copy'


    input:

    /*
     * TODO 5:
     *
     * Declare one path input called:
     *
     *   plink_files
     *
     * This will receive the complete collection of
     * .bed, .bim, and .fam files.
     */

    // TODO 5: add the path input declaration here
    path(plink_files)

    output:

    tuple path('combined.bed'),
          path('combined.bim'),
          path('combined.fam')


    script:

    """
    ls *.bed | sed 's/\\.bed\$//' > merge_list.txt

    plink \
        --merge-list merge_list.txt \
        --make-bed \
        --out combined
    """
}



/*
 * Define how the processes are connected.
 */
workflow {


    /*
     * Command-line parameters initially enter Nextflow
     * as strings.
     *
     * file() tells Nextflow that these strings refer
     * to files that Nextflow should track.
     */
    vcf_file = file(params.vcf)
    index_file = file(params.index)



    /*
     * Read the BED file one line at a time.
     *
     * You are NOT responsible for the Groovy/string-
     * manipulation syntax in this section.
     *
     * splitText() causes each BED row to become a
     * separate item flowing through the channel.
     */
    regions_ch = Channel
        .fromPath(params.bed)
        .splitText()
        .filter { line -> line.trim() }
        .map { line ->

            /*
             * Split the BED row into its three fields.
             */
            def fields = line.trim().split('\t')

            def chr = fields[0]
            def left = fields[1] as int
            def right = fields[2] as int


            /*
             * Remove the "chr" prefix.
             *
             * Example:
             * chr12 -> 12
             */
            def chr_num = chr.replaceFirst('^chr', '')


            /*
             * BED uses a 0-based left coordinate.
             * VCF uses a 1-based coordinate.
             */
            def left_adj = left + 1


            /*
             * Construct the VCF-style region.
             *
             * Example:
             * 12:150001-500000
             */
            def region = "${chr_num}:${left_adj}-${right}"


            /*
             * Construct an informative prefix for
             * files produced for this region.
             *
             * Example:
             * chr12_150001_500000
             */
            def prefix = "${chr}_${left_adj}_${right}"


            /*
             * Each BED row produces ONE tuple containing:
             *
             *   VCF file
             *   VCF index
             *   region
             *   output prefix
             *
             * Since the BED file has multiple rows,
             * regions_ch will contain multiple tuples.
             */
            tuple(
                vcf_file,
                index_file,
                region,
                prefix
            )
        }



    /*
     * TODO 6:
     *
     * Run EXTRACT_REGION using regions_ch.
     *
     * Save the output channel as:
     *
     *   regional_vcfs
     */

    // TODO 6: call EXTRACT_REGION here
    regional_vcfs = EXTRACT_REGION(regions_ch)


    /*
     * TODO 7:
     *
     * Send regional_vcfs into MAKE_PLINK.
     *
     * Save the output channel as:
     *
     *   plink_sets
     */

    // TODO 7: call MAKE_PLINK here
    plink_sets = MAKE_PLINK(regional_vcfs)


    /*
     * plink_sets contains tuples like:
     *
     *   [prefix, bed, bim, fam]
     *
     * The code below removes the prefix and keeps
     * the three PLINK files.
     *
     * You are NOT responsible for the Groovy syntax
     * in the map block.
     */
    merge_files = plink_sets
        .map { prefix, bed, bim, fam ->
            [bed, bim, fam]
        }
        .flatten()


        /*
         * TODO 8:
         *
         * Add the Nextflow operator that:
         *
         *   - waits for all upstream items
         *   - gathers them into ONE collection
         *
         * Without this step, the files would continue
         * downstream one at a time.
         *
         * With this step, MERGE_PLINK can receive all
         * regional PLINK files together and run once.
         */

        // TODO 8: add the operator here
	.collect()


    /*
     * TODO 9:
     *
     * Run MERGE_PLINK using merge_files.
     */

    // TODO 9: call MERGE_PLINK here
    MERGE_PLINK(merge_files)	


    /*
     * TODO 10:
     *
     * Add brief comments below answering:
     *
     * 1. Why can several EXTRACT_REGION tasks run
     *    independently or in parallel?
     *		
     *		several EXTRACT_REGION tasks can run at the same time since they are just looking at one region of the genome from the VCF file, and since it should run for each region in the bed file running multiple tasks at once is good. 
     *
     * 2. Why must MAKE_PLINK for a particular region
     *    wait until its corresponding EXTRACT_REGION
     *    task has finished?
     *	
     *		MAKE_PLINK uses the output from the EXTRACT_REGION process so it needs the corresponding EXTRACT_REGION to finish so it has its input. 
     *
     *
     * 3. Why must MERGE_PLINK wait until ALL
     *    MAKE_PLINK tasks have finished?
     *
     *		The MERGE_PLINK process uses the plink files output from the MAKE_PLINK process so would be unable to start until the MAKE_PLINK has finished giving it its input. 
     *
     */


}

