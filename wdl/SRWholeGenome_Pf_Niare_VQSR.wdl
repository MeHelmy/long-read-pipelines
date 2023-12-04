version 1.0

######################################################################################
## A workflow that performs single sample variant calling on Illumina reads from
## one or more flow cells. The workflow merges multiple samples into a single BAM
## prior to variant calling.
######################################################################################

import "tasks/Utils.wdl" as Utils
import "tasks/SRUtils.wdl" as SRUTIL
import "tasks/Finalize.wdl" as FF
import "tasks/VariantUtils.wdl" as VARUTIL
import "tasks/Pf_Niare_HaplotypeCaller.wdl" as Niare_HC

workflow SRWholeGenome_Pf_Niare_VQSR {
    input {
        Array[File] aligned_bams
        Array[File] aligned_bais

        File ref_map_file

        String participant_name

        File vcf_calling_interval_list
        File genotype_gvcfs_intervals

        String gcs_out_root_dir

        File vqsr_sites_vcf
        File vqsr_sites_vcf_index

        Boolean call_vars_on_mitochondria = false
        String mito_contig = "chrM"
        Array[String] contigs_names_to_ignore = ["RANDOM_PLACEHOLDER_VALUE"]  ## Required for ignoring any filtering - this is kind of a hack - TODO: fix the task!
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/SRWholeGenome_Pf_Niare_VQSR/~{participant_name}"

    String bam_dir = outdir + "/alignments"
    String metrics_dir = outdir + "/metrics"
    String smalldir = outdir + "/variants/small"
    String recalibration_dir = outdir + "/variants/recalibration_files"

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

    ####################################################################################################
    # HC Call Variants:

    # Now we handle HaplotypeCaller data:
    call Niare_HC.CallVariantsWithHaplotypeCaller {
        input:
            bam               = bam,
            bai               = bai,
            sample_id         = participant_name,
            ref_fasta         = ref_map['fasta'],
            ref_fasta_fai     = ref_map['fai'],
            ref_dict          = ref_map['dict'],

            vcf_calling_interval_list = vcf_calling_interval_list,
            genotype_gvcfs_intervals = genotype_gvcfs_intervals,

            prefix = participant_name + ".haplotype_caller",

            mito_contig = ref_map['mt_chr_name'],
            contigs_names_to_ignore = contigs_names_to_ignore,
    }

    # Make sure our sample name is correct:
    call VARUTIL.RenameSingleSampleVcf as RenameRawHcVcf {
        input:
            vcf = CallVariantsWithHaplotypeCaller.output_vcf,
            vcf_index = CallVariantsWithHaplotypeCaller.output_vcf_index,
            prefix = participant_name + ".haplotype_caller.renamed",
            new_sample_name = participant_name
    }
    call VARUTIL.RenameSingleSampleVcf as RenameRawHcGvcf {
        input:
            vcf = CallVariantsWithHaplotypeCaller.output_gvcf,
            vcf_index = CallVariantsWithHaplotypeCaller.output_gvcf_index,
            prefix = participant_name + ".haplotype_caller.renamed",
            is_gvcf = true,
            new_sample_name = participant_name
    }

    ################################################################################################
    # VQSR:

    # Scatter by chromosome:
    Array[String] use_filter = if (call_vars_on_mitochondria) then contigs_names_to_ignore else flatten([[mito_contig], contigs_names_to_ignore])
    call Utils.MakeChrIntervalList as SmallVariantsScatterPrep {
        input:
            ref_dict = ref_map['dict'],
            filter = use_filter
    }

    # Call over the scattered intervals:
    scatter (c in SmallVariantsScatterPrep.chrs) {
        String contig_for_small_var = c[0]

        call VARUTIL.SubsetVCF as GetHcCallsForContig {
            input:
                vcf_gz = RenameRawHcVcf.new_sample_name_vcf,
                vcf_tbi = RenameRawHcVcf.new_sample_name_vcf_index,
                locus = contig_for_small_var,
                prefix = participant_name + "." + contig_for_small_var,
        }

        call Niare_HC.NormalizeVcfSplittingMultiallelics as NormalizeVcfPreVqsr {
            input:
                input_vcf = GetHcCallsForContig.subset_vcf,
                input_vcf_index = GetHcCallsForContig.subset_tbi,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                prefix = participant_name + "." + contig_for_small_var + ".norm"
        }

        call Niare_HC.VariantRecalibratorIndel as VariantRecalibratorIndel {
            input:
                input_vcf = NormalizeVcfPreVqsr.output_vcf,
                input_vcf_index = NormalizeVcfPreVqsr.output_vcf_index,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                sites_only_vcf = vqsr_sites_vcf,
                sites_only_vcf_index = vqsr_sites_vcf_index,
                prefix = participant_name + "." + contig_for_small_var + ".norm",
        }

        call Niare_HC.ApplyVqsrIndel as ApplyVqsrIndel {
            input:
                input_vcf = NormalizeVcfPreVqsr.output_vcf,
                input_vcf_index = NormalizeVcfPreVqsr.output_vcf_index,
                recal_file = VariantRecalibratorIndel.recalibration,
                recal_file_index = VariantRecalibratorIndel.recalibration_index,
                recal_tranches = VariantRecalibratorIndel.tranches,
                prefix = participant_name + "." + contig_for_small_var + ".norm",
        }

        call Niare_HC.VariantRecalibratorSnp as VariantRecalibratorSnp {
            input:
                input_vcf = ApplyVqsrIndel.output_vcf,
                input_vcf_index = ApplyVqsrIndel.output_vcf_index,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                sites_only_vcf = vqsr_sites_vcf,
                sites_only_vcf_index = vqsr_sites_vcf_index,
                prefix = participant_name + "." + contig_for_small_var + ".norm",
        }

        call Niare_HC.ApplyVqsrSnp as ApplyVqsrSnp {
            input:
                input_vcf = ApplyVqsrIndel.output_vcf,
                input_vcf_index = ApplyVqsrIndel.output_vcf_index,
                recal_file = VariantRecalibratorSnp.recalibration,
                recal_file_index = VariantRecalibratorSnp.recalibration_index,
                recal_tranches = VariantRecalibratorSnp.tranches,
                prefix = participant_name + "." + contig_for_small_var + ".norm.indel_recal",
        }

        call Niare_HC.MergeMultiAllelicSitesPostRecalibration as MergeMultiAllelicSitesPostRecalibration {
            input:
                input_vcf = ApplyVqsrSnp.output_vcf,
                input_vcf_index = ApplyVqsrSnp.output_vcf_index,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                prefix = participant_name + "." + contig_for_small_var,
        }
    }

    call SRUTIL.MergeVCFs as MergeVCFs {
        input:
            input_vcfs = MergeMultiAllelicSitesPostRecalibration.output_vcf,
            input_vcfs_indexes = MergeMultiAllelicSitesPostRecalibration.output_vcf_index,
            prefix = participant_name + ".recalibrated"
    }

    ################################################################################################

    # Create a Keyfile for finalization:
    File keyfile = MergeVCFs.output_vcf_index

    # Finalize the raw Joint Calls:
    call FF.FinalizeToFile as FinalizeRawHCVcf { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcVcf.new_sample_name_vcf }
    call FF.FinalizeToFile as FinalizeRawHCTbi { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcVcf.new_sample_name_vcf_index }
    call FF.FinalizeToFile as FinalizeHCGVcf { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcGvcf.new_sample_name_vcf }
    call FF.FinalizeToFile as FinalizeHCGTbi { input: outdir = smalldir, keyfile = keyfile, file = RenameRawHcGvcf.new_sample_name_vcf_index }
    call FF.FinalizeToFile as FinalizeHCBamOut { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.bamout }
    call FF.FinalizeToFile as FinalizeHCBaiOut { input: outdir = smalldir, keyfile = keyfile, file = CallVariantsWithHaplotypeCaller.bamout_index }
    call FF.FinalizeToFile as FinalizeRecalibratedVcf { input: outdir = smalldir, keyfile = keyfile, file = MergeVCFs.output_vcf }
    call FF.FinalizeToFile as FinalizeRecalibratedVcfIndex { input: outdir = smalldir, keyfile = keyfile, file = MergeVCFs.output_vcf_index }

    ################################
    # Finalize the VETS files:
    ############

    output {
        Boolean successfully_processed = true

        ########################################

        File? hc_g_vcf    = FinalizeHCGVcf.gcs_path
        File? hc_g_tbi    = FinalizeHCGTbi.gcs_path
        File? hc_bamout   = FinalizeHCBamOut.gcs_path
        File? hc_baiout   = FinalizeHCBaiOut.gcs_path
        File? hc_raw_vcf  = FinalizeRawHCVcf.gcs_path
        File? hc_raw_tbi  = FinalizeRawHCTbi.gcs_path
        File? hc_rescored_vcf = FinalizeRecalibratedVcf.gcs_path
        File? hc_rescored_tbi = FinalizeRecalibratedVcfIndex.gcs_path
    }
}
