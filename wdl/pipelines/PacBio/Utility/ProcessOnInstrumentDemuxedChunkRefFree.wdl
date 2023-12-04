version 1.0

import "../../../tasks/Visualization/NanoPlot.wdl"

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/BAMutils.wdl" as BU
import "../../../tasks/Utility/GeneralUtils.wdl" as GU

import "../../../tasks/Utility/Finalize.wdl" as FF

import "../../TechAgnostic/Utility/CountTheBeans.wdl" as QC1
import "../../TechAgnostic/Utility/DystPeaker.wdl" as QC2
import "../../TechAgnostic/Utility/FASTQstats.wdl" as QC3

workflow ProcessOnInstrumentDemuxedChunkRefFree {

    meta {
        desciption: "!!! WARN: THIS IS PROJECT-CENTER SPECIFIC !!! Given an on-instrument demultiplexed hifi_reads.bam, perform ref-independent prep work."
    }

    parameter_meta {
        bam_descriptor: "a one-word description of the input BAM (used for saving the reads that don't have any MM/ML tags)"
        short_reads_threshold: "a length threshold below which reads are classified as short"
    }

    input {
        File uBAM

        String readgroup_id

        String bam_descriptor
        Int short_reads_threshold

        String gcs_out_root_dir
    }

    ###################################################################################
    # prep work

    # where to store final results
    String workflow_name = "ProcessOnInstrumentDemuxedChunkRefFree"
    String outdir = sub(gcs_out_root_dir, "/$", "") + "/" + workflow_name
    String outdir_ref_free = outdir + '/RefFree'

    ###################################################################################
    # stats
    call NanoPlot.NanoPlotFromUBam {input: uBAM = uBAM}

    call BU.BamToFastq {input: bam = uBAM, prefix = "does_not_matter"}

    ###################################################################################
    # finalize
    call BU.GetReadGroupInfo as RG {input: bam = uBAM, keys = ['PU']}
    String movie_name = RG.read_group_info['PU']
    String bc_specific_fastq_out = outdir_ref_free + '/' + readgroup_id

    call FF.FinalizeToFile as FinalizeFQ { input: outdir = bc_specific_fastq_out, file = BamToFastq.reads_fq, name =  readgroup_id + '.hifi.fq.gz' }

    ###################################################################################
    # more QCs and metrics

    call QC1.CountTheBeans { input: bam=uBAM, bam_descriptor=bam_descriptor, gcs_out_root_dir=gcs_out_root_dir }
    call QC2.DystPeaker { input: input_file=uBAM, input_is_bam=true, id=readgroup_id, short_reads_threshold=short_reads_threshold, gcs_out_root_dir=gcs_out_root_dir }
    call QC3.FASTQstats { input: reads=BamToFastq.reads_fq, file_type='FASTQ' }

    ###################################################################################

    call GU.GetTodayDate as today {}

    output {
        String movie = movie_name

        File hifi_fq = FinalizeFQ.gcs_path
        # todo merge these two together
        Map[String, Float] hifi_stats = NanoPlotFromUBam.stats_map
        Map[String, Float] hifi_fq_stats = FASTQstats.stats

        # read length metrics
        File read_len_hist = DystPeaker.read_len_hist
        Array[Int] read_len_peaks = DystPeaker.read_len_peaks
        Array[Int] read_len_deciles = DystPeaker.read_len_deciles
        Map[String, String] read_len_summaries = DystPeaker.read_len_summaries

        # methylation call rate stats
        Map[String, String] methyl_tag_simple_stats = CountTheBeans.methyl_tag_simple_stats

        String last_processing_date = today.yyyy_mm_dd
    }
}
