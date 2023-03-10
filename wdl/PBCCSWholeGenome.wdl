version 1.0

######################################################################################
## A workflow that performs single sample variant calling on PacBio HiFi reads from
## one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to variant calling.
######################################################################################

import "tasks/PBUtils.wdl" as PB
import "tasks/Utils.wdl" as Utils
import "tasks/CallVariantsPBCCS.wdl" as VAR
import "tasks/Finalize.wdl" as FF

import "tasks/SampleLevelAlignedMetrics.wdl" as COV

workflow PBCCSWholeGenome {
    input {
        Array[File] aligned_bams
        Array[File] aligned_bais

        File? bed_to_compute_coverage

        File ref_map_file

        String participant_name

        String gcs_out_root_dir

        Boolean call_svs = true
        Boolean? pbsv_call_per_chr = true

        Boolean call_small_variants = true
        Boolean? call_small_vars_on_mitochondria = false
        File? sites_vcf
        File? sites_vcf_tbi

        Boolean? run_dv_pepper_analysis = true
        Int? dvp_threads = 32
        Int? dvp_memory = 128
        File? ref_scatter_interval_list_locator
        File? ref_scatter_interval_list_ids
    }

    parameter_meta {
        aligned_bams:       "GCS path to aligned BAM files"
        aligned_bais:       "GCS path to aligned BAM file indices"
        participant_name:   "name of the participant from whom these samples were obtained"

        ref_map_file:       "table indicating reference sequence and auxillary file locations"
        gcs_out_root_dir:   "GCS bucket to store the reads, variants, and metrics files"

        call_svs:               "whether to call SVs"
        pbsv_call_per_chr:      "when using PBSV, make calls per chromosome then merge (trade lower sensitivity for faster speed)"

        call_small_variants: "whether to call small variants"
        call_small_vars_on_mitochondria: "if false, will not attempt to call variants on mitochondria; if true, some samples might fail (caller feature) due to lack of signal"
        sites_vcf:     "for use with Clair"
        sites_vcf_tbi: "for use with Clair"

        run_dv_pepper_analysis:  "to turn on DV-Pepper analysis or not (non-trivial increase in cost and runtime)"
        ref_scatter_interval_list_locator: "A file holding paths to interval_list files; needed only when running DV-Pepper"
        ref_scatter_interval_list_ids:     "A file that gives short IDs to the interval_list files; needed only when running DV-Pepper"
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/PBCCSWholeGenome/~{participant_name}"

    # gather across (potential multiple) input CCS BAMs
    if (length(aligned_bams) > 1) {
        scatter (pair in zip(aligned_bams, aligned_bais)) {
            call Utils.InferSampleName {input: bam = pair.left, bai = pair.right}
        }
        call Utils.CheckOnSamplenames {input: sample_names = InferSampleName.sample_name}

        call Utils.MergeBams as MergeAllReads { input: bams = aligned_bams, prefix = participant_name }
    }

    File bam = select_first([MergeAllReads.merged_bam, aligned_bams[0]])
    File bai = select_first([MergeAllReads.merged_bai, aligned_bais[0]])

    call PB.PBIndex as IndexCCSUnalignedReads { input: bam = bam }
    File pbi = IndexCCSUnalignedReads.pbi

    call COV.SampleLevelAlignedMetrics as coverage {
        input:
            aligned_bam = bam,
            aligned_bai = bai,
            ref_fasta   = ref_map['fasta'],
            bed_to_compute_coverage = bed_to_compute_coverage
    }

    String dir = outdir + "/alignments"

    call FF.FinalizeToFile as FinalizeBam { input: outdir = dir, file = bam, name = "~{participant_name}.bam" }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = dir, file = bai, name = "~{participant_name}.bam.bai" }
    call FF.FinalizeToFile as FinalizePbi { input: outdir = dir, file = pbi, name = "~{participant_name}.bam.pbi" }

    if (defined(bed_to_compute_coverage)) { call FF.FinalizeToFile as FinalizeRegionalCoverage { input: outdir = dir, file = select_first([coverage.bed_cov_summary]) } }

    ####################################################################################################
    if (call_svs || call_small_variants) {

        # verify arguments are provided
        if (call_svs) {
            if (! defined(pbsv_call_per_chr)) {call Utils.StopWorkflow as pbsv_call_per_chr_not_provided {input: reason = "Calling SVs without specifying arg pbsv_call_per_chr"}}
        }
        if (call_small_variants) {
            if (! defined(call_small_vars_on_mitochondria)) {call Utils.StopWorkflow as call_small_vars_on_mitochondria_not_provided {input: reason = "Unprovided arg call_small_vars_on_mitochondria"}}
            if (! defined(run_dv_pepper_analysis)) {call Utils.StopWorkflow as run_dv_pepper_analysis_not_provided {input: reason = "Unprovided arg run_dv_pepper_analysis"}}
            if (! defined(dvp_threads)) {call Utils.StopWorkflow as dvp_threads_not_provided {input: reason = "Unprovided arg dvp_threads"}}
            if (! defined(ref_scatter_interval_list_locator)) {call Utils.StopWorkflow as ref_scatter_interval_list_locator_not_provided {input: reason = "Unprovided arg ref_scatter_interval_list_locator"}}
            if (! defined(ref_scatter_interval_list_ids)) {call Utils.StopWorkflow as ref_scatter_interval_list_ids_not_provided {input: reason = "Unprovided arg ref_scatter_interval_list_ids"}}
        }

        call VAR.CallVariants {
            input:
                bam               = bam,
                bai               = bai,
                sample_id         = participant_name,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                tandem_repeat_bed = ref_map['tandem_repeat_bed'],

                prefix = participant_name,

                call_svs = call_svs,
                pbsv_call_per_chr = select_first([pbsv_call_per_chr]),

                call_small_variants = call_small_variants,
                call_small_vars_on_mitochondria = select_first([call_small_vars_on_mitochondria]),
                sites_vcf = sites_vcf,
                sites_vcf_tbi = sites_vcf_tbi,

                run_dv_pepper_analysis = select_first([run_dv_pepper_analysis]),
                dvp_threads = select_first([dvp_threads]),
                dvp_memory = select_first([dvp_memory]),
                ref_scatter_interval_list_locator = select_first([ref_scatter_interval_list_locator]),
                ref_scatter_interval_list_ids = select_first([ref_scatter_interval_list_ids])
        }

        String svdir = outdir + "/variants/sv"
        String smalldir = outdir + "/variants/small"

        if (call_svs) {
            call FF.FinalizeToFile as FinalizePBSV { input: outdir = svdir, file = select_first([CallVariants.pbsv_vcf]) }
            call FF.FinalizeToFile as FinalizePBSVtbi { input: outdir = svdir, file = select_first([CallVariants.pbsv_tbi]) }

            call FF.FinalizeToFile as FinalizeSniffles { input: outdir = svdir, file = select_first([CallVariants.sniffles_vcf]) }
            call FF.FinalizeToFile as FinalizeSnifflesTbi { input: outdir = svdir, file = select_first([CallVariants.sniffles_tbi]) }
            call FF.FinalizeToFile as FinalizeSnifflesSnf { input: outdir = svdir, file = select_first([CallVariants.sniffles_snf]) }
        }

        if (call_small_variants) {
            call FF.FinalizeToFile as FinalizeClairVcf { input: outdir = smalldir, file = select_first([CallVariants.clair_vcf])}
            call FF.FinalizeToFile as FinalizeClairTbi { input: outdir = smalldir, file = select_first([CallVariants.clair_tbi])}

            call FF.FinalizeToFile as FinalizeClairGVcf { input: outdir = smalldir, file = select_first([CallVariants.clair_gvcf])}
            call FF.FinalizeToFile as FinalizeClairGTbi { input: outdir = smalldir, file = select_first([CallVariants.clair_gtbi])}

            if (select_first([run_dv_pepper_analysis])) {
                call FF.FinalizeToFile as FinalizeDVPepperVcf { input: outdir = smalldir, file = select_first([CallVariants.dvp_vcf])}
                call FF.FinalizeToFile as FinalizeDVPepperTbi { input: outdir = smalldir, file = select_first([CallVariants.dvp_tbi])}
                call FF.FinalizeToFile as FinalizeDVPepperGVcf { input: outdir = smalldir, file = select_first([CallVariants.dvp_g_vcf])}
                call FF.FinalizeToFile as FinalizeDVPepperGTbi { input: outdir = smalldir, file = select_first([CallVariants.dvp_g_tbi])}
                call FF.FinalizeToFile as FinalizeDVPEPPERPhasedVcf { input: outdir = smalldir, file = select_first([CallVariants.dvp_phased_vcf]), name = "~{participant_name}.deepvariant_pepper.phased.vcf.gz" }
                call FF.FinalizeToFile as FinalizeDVPEPPERPhasedTbi { input: outdir = smalldir, file = select_first([CallVariants.dvp_phased_tbi]), name = "~{participant_name}.deepvariant_pepper.phased.vcf.gz.tbi" }
            }
        }
    }

    output {
        File aligned_bam = FinalizeBam.gcs_path
        File aligned_bai = FinalizeBai.gcs_path
        File aligned_pbi = FinalizePbi.gcs_path

        Map[String, Float] aln_summary = {
            "aligned_num_reads" :coverage.aligned_num_reads,
            "aligned_num_bases" :coverage.aligned_num_bases,
            "aligned_frac_bases" :coverage.aligned_frac_bases,
            "aligned_est_fold_cov" :coverage.aligned_est_fold_cov,

            "aligned_read_length_mean" :coverage.aligned_read_length_mean,
            "aligned_read_length_median" :coverage.aligned_read_length_median,
            "aligned_read_length_stdev" :coverage.aligned_read_length_stdev,
            "aligned_read_length_N50" :coverage.aligned_read_length_N50,

            "average_identity" :coverage.average_identity,
            "median_identity" :coverage.median_identity
        }

        File? bed_cov_summary = FinalizeRegionalCoverage.gcs_path
        ########################################
        File? pbsv_vcf = FinalizePBSV.gcs_path
        File? pbsv_tbi = FinalizePBSVtbi.gcs_path

        File? sniffles_vcf = FinalizeSniffles.gcs_path
        File? sniffles_tbi = FinalizeSnifflesTbi.gcs_path
        File? sniffles_snf = FinalizeSnifflesSnf.gcs_path

        File? clair_vcf = FinalizeClairVcf.gcs_path
        File? clair_tbi = FinalizeClairTbi.gcs_path

        File? clair_gvcf = FinalizeClairGVcf.gcs_path
        File? clair_gtbi = FinalizeClairGTbi.gcs_path

        File? dvp_vcf = FinalizeDVPepperVcf.gcs_path
        File? dvp_tbi = FinalizeDVPepperTbi.gcs_path
        File? dvp_g_vcf = FinalizeDVPepperGVcf.gcs_path
        File? dvp_g_tbi = FinalizeDVPepperGTbi.gcs_path
        File? dvp_phased_vcf = FinalizeDVPEPPERPhasedVcf.gcs_path
        File? dvp_phased_tbi = FinalizeDVPEPPERPhasedTbi.gcs_path
    }
}
