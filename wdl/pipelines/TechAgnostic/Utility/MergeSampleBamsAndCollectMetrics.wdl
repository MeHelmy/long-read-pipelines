version 1.0

import "../../../tasks/Utility/Utils.wdl" as Utils
import "../../../tasks/Utility/PBUtils.wdl" as PB
import "../../../tasks/Utility/BAMutils.wdl" as BU

import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../../tasks/QC/SampleLevelAlignedMetrics.wdl" as COV

import "../../TechAgnostic/Utility/VerifyBamFingerprint.wdl" as QC0
import "../../TechAgnostic/Utility/LongReadsContaminationEstimation.wdl" as QC1
import "../../TechAgnostic/Utility/SexCheckNaive.wdl" as QC2
import "../../TechAgnostic/Utility/CountTheBeans.wdl" as QC3

workflow Work {
    meta {
        description: "Merge a sample's (potential) multiple SMRT-/flow-cells data, and collect alignment metrics."
    }
    parameter_meta {
        tech: "LR technology used for generating the data; accepted value: [ONT, Sequel, Revio]"
        bams_suspected_to_contain_dup_record: "Some ONT output files from basecall dirs have a strange duplicate issue."
        bed_to_compute_coverage: "BED file holding regions-of-interest for computing coverage over."
        bed_descriptor: "Description of the BED file, will be used in the file name so be careful naming things"

        fingerprint_store:  "Path to the GCS folder holding fingerprint VCFs"
        sample_id_at_store: "Sampld ID at the fingerprint store for uniquely locating the fingerprint VCF (assumes some naming convention)"
        vbid2_config_json: "Path to Json file specifying config to be passed to VBID2 workflow (cross-individual contamination); currently only supports GRCh38."
        expected_sex_type: "Expected sex type of the sample that generated the BAM; accepted value: [M, F, NA]; if provided, performs sex concordance check"
        check_postaln_methyl_tags: "if true, gather statistics and reads without MM/ML tags from the aligned BAM"
    }
    input {
        String gcs_out_dir

        # sample specific
        String sample_name
        Array[File] aligned_bams
        Array[File] aligned_bais

        String tech
        Boolean bams_suspected_to_contain_dup_record

        File? bed_to_compute_coverage
        String? bed_descriptor

        File ref_map_file

        # args for optional QC subworkflows
        String? fingerprint_store
        String? sample_id_at_store
        File? vbid2_config_json
        String? expected_sex_type
        Boolean check_postaln_methyl_tags = true
    }
    output {
        File aligned_bam = FinalizeBam.gcs_path
        File aligned_bai = FinalizeBai.gcs_path
        File? aligned_pbi = FinalizePbi.gcs_path

        Float coverage = AlignmentMetrics.coverage
        File coverage_per_chr = FinalizeMosdepSumm.gcs_path
        String nanoplots = FinalizeNanoplots.gcs_dir
        File? bed_cov_summary = FinalizeRegionalCoverage.gcs_path

        Map[String, Float] alignment_metrics = AlignmentMetrics.reads_stats

        Map[String, Float] sam_flag_stats = AlignmentMetrics.sam_flag_stats

        # contamination QC
        Float? contamination_est = VBID2.contamination_est

        # sex concordance QC
        Map[String, String]? inferred_sex_info = SexConcordance.inferred_sex_info

        # post-alignment methylation
        Map[String, String]? methyl_tag_simple_stats = NoMissingBeans.methyl_tag_simple_stats

        # fingerprint check
        Map[String, String]? fingerprint_check = VerifyBamFingerprint.fingerprint_check
    }

    if (defined(bed_to_compute_coverage)) {
        if (!defined(bed_descriptor)) {
            call Utils.StopWorkflow { input: reason = "Must provied descriptive name of the BED file if the file is provided."}
        }
    }
    if (defined(fingerprint_store) != defined(sample_id_at_store)) {
        call Utils.StopWorkflow as CannotFingerprint { input: reason = "fingerprint_store and sample_id_at_store must be specified together or omitted together" }
    }

    String outdir = sub(gcs_out_dir, "/$", "") + "/alignments"

    ###########################################################
    # some input validation
    scatter (pair in zip(aligned_bams, aligned_bais)) {
        call Utils.InferSampleName {input: bam = pair.left, bai = pair.right}
    }
    call Utils.CheckOnSamplenames {input: sample_names = InferSampleName.sample_name}
    if (InferSampleName.sample_name[0] != sample_name) {
        call Utils.StopWorkflow as SM_mismatch { input: reason = "Provided sample name and those encoded in the BAM(s) don't match."}
    }

    ###########################################################
    # gather across (potential multiple) input BAMs
    if (length(aligned_bams) > 1) {
        call BU.MergeBamsWithSamtools as MergeAllReads { input: bams = aligned_bams, out_prefix = sample_name }
    }

    File bam = select_first([MergeAllReads.merged_bam, aligned_bams[0]])
    File bai = select_first([MergeAllReads.merged_bai, aligned_bais[0]])

    ###########################################################
    # ont specific: sometimes there are duplicate reads
    if (('ONT'==tech) && bams_suspected_to_contain_dup_record) {
        call Utils.DeduplicateBam as RemoveONTDuplicates {
            input: aligned_bam = bam, aligned_bai = bai
        }
    }

    ###########################################################
    # save bam and index
    File use_this_bam = select_first([RemoveONTDuplicates.corrected_bam, bam])
    File use_this_bai = select_first([RemoveONTDuplicates.corrected_bai, bai])

    call FF.FinalizeToFile as FinalizeBam { input: outdir = outdir, file = use_this_bam, name = "~{sample_name}.bam" }
    call FF.FinalizeToFile as FinalizeBai { input: outdir = outdir, file = use_this_bai, name = "~{sample_name}.bam.bai" }

    ###########################################################
    # pacbio specific index
    if ('ONT'!=tech) {
        call PB.PBIndex as PBIndexSampleReads { input: bam = use_this_bam }
        call FF.FinalizeToFile as FinalizePbi { input: outdir = outdir, file = PBIndexSampleReads.pbi, name = "~{sample_name}.bam.pbi" }
    }
    ###########################################################
    call COV.SampleLevelAlignedMetrics as AlignmentMetrics {
        input:
            aligned_bam = use_this_bam,
            aligned_bai = use_this_bai,
            bed_to_compute_coverage = bed_to_compute_coverage,
            bed_descriptor = bed_descriptor
    }
    # save alignment metrics
    String aln_metrics_out = sub(outdir, "/$", "") + "/metrics/~{sample_name}"
    call FF.FinalizeToFile as FinalizeMosdepSumm { input: outdir = aln_metrics_out, file = AlignmentMetrics.coverage_per_chr }
    call FF.FinalizeToDir as FinalizeNanoplots { input: outdir = aln_metrics_out, files = AlignmentMetrics.nano_plots }

    if (defined(bed_to_compute_coverage)) {
        call FF.FinalizeToFile as FinalizeRegionalCoverage { input:
            outdir = aln_metrics_out, file = select_first([AlignmentMetrics.bed_cov_summary]),
            name = "~{sample_name}.mosdepth_coverage.coverage_over_bed.~{bed_descriptor}.summary.txt"
        }
    }
    ###################################################################################
    # more QCs and metrics
    # (optional) fingerprint verification
    if (defined(fingerprint_store)) {
        Map[String, String] ref_map = read_map(ref_map_file)
        call QC0.VerifyBamFingerprint { input:
            aligned_bam=use_this_bam,
            aligned_bai=use_this_bai,
            fp_store=select_first([fingerprint_store]),
            sample_id_at_store=select_first([sample_id_at_store]),
            ref_specific_haplotype_map=ref_map['haplotype_map']
        }
    }

    # (optional) contamination
    if (defined(vbid2_config_json)) {
        Map[String, String] vbid2_config = read_json(select_first([vbid2_config_json]))
        call QC1.LongReadsContaminationEstimation as VBID2 { input:
            bam=use_this_bam,
            bai=use_this_bai,
            ref_map_file=ref_map_file,
            tech = tech,
            gt_sites_bed = vbid2_config['genotyping_sites_bed'],
            is_hgdp_sites = vbid2_config['is_HGDP_sites'],
            is_100k_sites = vbid2_config['is_100K_sites'],
            disable_baq = vbid2_config['disable_BAQ'],
            disk_type = "SSD",
        }
    }

    # (optional) sex concordance
    if (defined(expected_sex_type)) {
        call QC2.SexCheckNaive as SexConcordance { input:
            bam=use_this_bam,
            bai=use_this_bai,
            expected_sex_type=select_first([expected_sex_type]),
            mosdepth_summary_txt=AlignmentMetrics.coverage_per_chr
        }
    }

    # (optional) verify methylation tags aren't missing
    if (check_postaln_methyl_tags) {
        call QC3.CountTheBeans as NoMissingBeans { input:
            bam=use_this_bam,
            bai=use_this_bai,
            bam_descriptor="SAMPLE_OUTPUT",
            gcs_out_root_dir=gcs_out_dir,
            use_local_ssd=false
        }
    }
}
