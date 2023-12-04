version 1.0

import "../../../tasks/Utility/GeneralUtils.wdl" as GU
import "../../../tasks/Utility/Utils.wdl"

import "../../../tasks/Alignment/AlignAndCheckFingerprintCCS.wdl" as major

import "../../TechAgnostic/Utility/LongReadsContaminationEstimation.wdl" as QC1
import "../../TechAgnostic/Utility/SexCheckNaive.wdl" as QC2
import "../../TechAgnostic/Utility/CountTheBeans.wdl" as QC3

import "../../../tasks/Utility/Finalize.wdl" as FF

workflow ProcessOnInstrumentDemuxedChunk {

    meta {
        desciption: "!!! WARN: THIS IS PROJECT-CENTER SPECIFIC !!! Given an on-instrument demultiplexed hifi_reads.bam, perform alignment and QC check."
    }

    parameter_meta {
        readgroup_id: "ID of a readgroup; used for storing outputs; no whitespaces allowed"
        bam_SM_field: "value to place in the SM field of the resulting BAM header's @RG lines"

        platform: "PacBio platform used for generating the data; accepted value: [Sequel, Revio]"

        fingerprint_store:  "Path to the GCS folder holding fingerprint VCFs"
        sample_id_at_store: "Sampld ID at the fingerprint store for uniquely locating the fingerprint VCF (assumes some naming convention)"
        vbid2_config_json: "Path to Json file specifying config to be passed to VBID2 workflow (cross-individual contamination); currently only supports GRCh38."
        expected_sex_type: "Expected sex type of the sample that generated the BAM; accepted value: [M, F, NA]; if provided, performs sex concordance check"
        check_postaln_methyl_tags: "if true, gather statistics and reads without MM/ML tags from the aligned BAM"
    }

    input {
        File  uBAM
        File? uPBI

        String readgroup_id

        String bam_SM_field

        String platform

        File ref_map_file

        String gcs_out_root_dir

        String disk_type = "SSD"

        # args for optional QC subworkflows
        String? fingerprint_store
        String? sample_id_at_store
        File? vbid2_config_json
        String? expected_sex_type
        Boolean check_postaln_methyl_tags = true
    }

    ###################################################################################
    # prep work

    # where to store final results
    String workflow_name = "ProcessOnInstrumentDemuxedChunk"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name

    if (defined(fingerprint_store) != defined(sample_id_at_store)) {
        call Utils.StopWorkflow { input: reason = "fingerprint_store and sample_id_at_store must be specified together or omitted together" }
    }

    ###################################################################################
    # major work
    String fp_store = select_first([fingerprint_store,  'None'])
    String fp_smid  = select_first([sample_id_at_store, 'None'])

    call major.AlignAndCheckFingerprintCCS {
        input:
            uBAM = uBAM,
            uPBI = uPBI,
            bam_sample_name = bam_SM_field,

            turn_off_fingperprint_check = !(defined(fingerprint_store)),
            fp_store = fp_store,
            sample_id_at_store = fp_smid,
            ref_map_file = ref_map_file
    }
    File aBAM = AlignAndCheckFingerprintCCS.aligned_bam
    File aBAI = AlignAndCheckFingerprintCCS.aligned_bai
    File aPBI = AlignAndCheckFingerprintCCS.aligned_pbi

    ###################################################################################
    # more QCs and metrics

    # (optional) contamination
    if (defined(vbid2_config_json)) {
        Map[String, String] vbid2_config = read_json(select_first([vbid2_config_json]))
        call QC1.LongReadsContaminationEstimation as VBID2 { input:
            bam=aBAM,
            bai=aBAI,
            ref_map_file=ref_map_file,
            tech = platform,
            gt_sites_bed = vbid2_config['genotyping_sites_bed'],
            is_hgdp_sites = vbid2_config['is_HGDP_sites'],
            is_100k_sites = vbid2_config['is_100K_sites'],
            disable_baq = vbid2_config['disable_BAQ'],
            disk_type = disk_type,
        }
    }

    # (optional) sex concordance
    if (defined(expected_sex_type)) {
        call QC2.SexCheckNaive as SexConcordance { input:
            bam=aBAM,
            bai=aBAI,
            expected_sex_type=select_first([expected_sex_type]),
            mosdepth_summary_txt=AlignAndCheckFingerprintCCS.coverage_per_chr
        }
    }

    # (optional) verify methylation tags aren't missing
    if (check_postaln_methyl_tags) {
        call QC3.CountTheBeans as NoMissingBeans { input:
            bam=aBAM,
            bai=aBAI,
            bam_descriptor="POST_ALN",
            gcs_out_root_dir=gcs_out_root_dir,
            use_local_ssd=disk_type=='LOCAL'
        }
    }

    ###################################################################################
    # finalize
    String bc_specific_aln_out    = outdir + '/alignments/' + readgroup_id
    String bc_specific_metric_out = outdir + "/metrics/"    + readgroup_id

    call FF.FinalizeToFile as FinalizeAlignedBam { input: outdir = bc_specific_aln_out, file = aBAM, name = readgroup_id + '.bam' }
    call FF.FinalizeToFile as FinalizeAlignedBai { input: outdir = bc_specific_aln_out, file = aBAI, name = readgroup_id + '.bai' }
    call FF.FinalizeToFile as FinalizeAlignedPbi { input: outdir = bc_specific_aln_out, file = aPBI, name = readgroup_id + '.pbi' }

    call FF.FinalizeToFile as FinalizePerChrCov  { input: outdir = bc_specific_metric_out, file = AlignAndCheckFingerprintCCS.coverage_per_chr }
    call FF.FinalizeToFile as FinalizeAlnMetrics { input: outdir = bc_specific_metric_out, file = AlignAndCheckFingerprintCCS.alignment_metrics_tar_gz }

    if (defined(fingerprint_store)) {
        call FF.FinalizeToFile as FinalizeFPDetails  { input: outdir = bc_specific_metric_out, file = select_first([AlignAndCheckFingerprintCCS.fingerprint_detail_tar_gz]) }
    }

    ###################################################################################

    call GU.GetTodayDate as today {}

    output {
        File aligned_bam = FinalizeAlignedBam.gcs_path
        File aligned_bai = FinalizeAlignedBai.gcs_path
        File aligned_pbi = FinalizeAlignedPbi.gcs_path
        Float wgs_cov = AlignAndCheckFingerprintCCS.wgs_cov

        File coverage_per_chr = FinalizePerChrCov.gcs_path

        File alignment_metrics_tar_gz = FinalizeAlnMetrics.gcs_path
        Map[String, Float] alignment_metrics = AlignAndCheckFingerprintCCS.alignment_metrics
        Map[String, Float] sam_flag_stats = AlignAndCheckFingerprintCCS.sam_flag_stats

        String movie = AlignAndCheckFingerprintCCS.movie

        # the following QC/metrics aren't always available
        # fingerprint QC
        String? fingerprint_check_result = AlignAndCheckFingerprintCCS.fp_status
        Float? fingerprint_check_LOD = AlignAndCheckFingerprintCCS.fp_lod_expected_sample
        File? fingerprint_check_tar_gz = FinalizeFPDetails.gcs_path

        # contamination QC
        Float? contamination_est = VBID2.contamination_est

        # sex concordance QC
        Map[String, String]? inferred_sex_info = SexConcordance.inferred_sex_info

        # post-alignment methylation
        Map[String, String]? methyl_tag_simple_stats = NoMissingBeans.methyl_tag_simple_stats

        String last_processing_date = today.yyyy_mm_dd
    }
}
