version 1.0

#############################################################################################################
## A workflow that performs joint calling on single-sample gVCFs from GATK4 HaplotypeCaller using GenomicsDB.
#############################################################################################################

import "tasks/SRJointGenotyping.wdl" as SRJOINT
import "tasks/VariantUtils.wdl" as VARUTIL
import "tasks/Utils.wdl" as UTILS
import "tasks/Finalize.wdl" as FF
import "tasks/Pf_Niare_HaplotypeCaller.wdl" as Niare_HC

workflow SRJointCallGVCFsWithGenomicsDB_Pf_Niare_VQSR {
    input {
        Array[File] gvcfs
        Array[File] gvcf_indices

        File ref_map_file

        File vqsr_sites_vcf
        File vqsr_sites_vcf_index

        String prefix

        String gcs_out_root_dir
    }

    parameter_meta {
        gvcfs:            "GCS paths to gVCF files"
        gvcf_indices:     "GCS paths to gVCF tbi files"
        ref_map_file:     "table indicating reference sequence and auxillary file locations"
        prefix:           "prefix for output joint-called gVCF and tabix index"
        gcs_out_root_dir: "GCS bucket to store the reads, variants, and metrics files"
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/SRJointCallGVCFsWithGenomicsDB_Pf_Niare_VQSR/~{prefix}"

    Map[String, String] ref_map = read_map(ref_map_file)

    # Create interval list over which to shard the processing:
    call UTILS.MakeChrIntervalList as MakeChrIntervalList {
        input:
            ref_dict = ref_map['dict'],
    }

    # Create sample-name map:
    call SRJOINT.CreateSampleNameMap as CreateSampleNameMap {
        input:
            gvcfs = gvcfs,
            prefix = prefix
    }

    # Shard by contig for speed:
    scatter (idx_1 in range(length(MakeChrIntervalList.contig_interval_list_files))) {

        String contig = MakeChrIntervalList.chrs[idx_1][0]
        File contig_interval_list = MakeChrIntervalList.contig_interval_list_files[idx_1]

        call Niare_HC.GenomicsDbImport as GenomicsDbImport {
            input:
                sample_name_map = CreateSampleNameMap.sample_name_map,
                interval_list = contig_interval_list,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                prefix = prefix + "." + contig,
        }

        # Shard again by contig chunk:
        call UTILS.SplitContigToIntervals as SplitContigToIntervals {
            input:
                ref_dict = ref_map['dict'],
                contig = contig
        }

        scatter (idx_2 in range(length(SplitContigToIntervals.individual_bed_files))) {

            File genotype_gvcfs_intervals = SplitContigToIntervals.individual_bed_files[idx_2]

            call Niare_HC.GenotypeGVCFs as GenotypeGVCFs {
                input:
                    input_gvcf_data = GenomicsDbImport.output_genomicsdb,
                    interval_list = genotype_gvcfs_intervals,
                    ref_fasta         = ref_map['fasta'],
                    ref_fasta_fai     = ref_map['fai'],
                    ref_dict          = ref_map['dict'],
                    prefix = prefix + "." + contig + ".raw",
            }
        }

        # Merge all raw VCFs:
        call VARUTIL.GatherVcfs as GatherVcfs {
            input:
                input_vcfs = GenotypeGVCFs.output_vcf,
                input_vcf_indices = GenotypeGVCFs.output_vcf_index,
                prefix = prefix + "." + contig + ".raw.merged",
        }

        # Normalize variants here:
        call Niare_HC.NormalizeVcfSplittingMultiallelics as NormalizeVcfSplittingMultiallelics {
            input:
                input_vcf = GatherVcfs.output_vcf,
                input_vcf_index = GatherVcfs.output_vcf_index,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                prefix = prefix + "." + contig + ".raw.merged.normalized",
        }

        # Run variant recalibrator and VQSR:
        call Niare_HC.VariantRecalibratorIndel as VariantRecalibratorIndel {
            input:
                input_vcf = NormalizeVcfSplittingMultiallelics.output_vcf,
                input_vcf_index = NormalizeVcfSplittingMultiallelics.output_vcf_index,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                sites_only_vcf = vqsr_sites_vcf,
                sites_only_vcf_index = vqsr_sites_vcf_index,
                prefix = prefix + "." + contig,
        }

        call Niare_HC.ApplyVqsrIndel as ApplyVqsrIndel {
            input:
                input_vcf = NormalizeVcfSplittingMultiallelics.output_vcf,
                input_vcf_index = NormalizeVcfSplittingMultiallelics.output_vcf_index,
                recal_file = VariantRecalibratorIndel.recalibration,
                recal_file_index = VariantRecalibratorIndel.recalibration_index,
                recal_tranches = VariantRecalibratorIndel.tranches,
                prefix = prefix + "." + contig + "raw",
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
                prefix = prefix + "." + contig,
        }

        call Niare_HC.ApplyVqsrIndel as ApplyVqsrSnp {
            input:
                input_vcf = ApplyVqsrIndel.output_vcf,
                input_vcf_index = ApplyVqsrIndel.output_vcf_index,
                recal_file = VariantRecalibratorSnp.recalibration,
                recal_file_index = VariantRecalibratorSnp.recalibration_index,
                recal_tranches = VariantRecalibratorSnp.tranches,
                prefix = prefix + "." + contig + "raw",
        }

        # Merge multi-allelic sites after recalibration:
        call Niare_HC.MergeMultiAllelicSitesPostRecalibration as MergeMultiAllelicSitesPostRecalibration {
            input:
                input_vcf = ApplyVqsrSnp.output_vcf,
                input_vcf_index = ApplyVqsrSnp.output_vcf_index,
                ref_fasta         = ref_map['fasta'],
                ref_fasta_fai     = ref_map['fai'],
                ref_dict          = ref_map['dict'],
                prefix = prefix + "." + contig + ".recalibrated",
        }
    }

    # Consolidate files:
    call VARUTIL.GatherVcfs as GatherRescoredVcfs {
        input:
            input_vcfs = MergeMultiAllelicSitesPostRecalibration.output_vcf,
            input_vcf_indices = MergeMultiAllelicSitesPostRecalibration.output_vcf_index,
            prefix = prefix + ".rescored.combined"
    }

    ################################
    # Finalize the regular output files:
    ############
    File keyfile = GatherRescoredVcfs.output_vcf_index

    call FF.FinalizeToFile as FinalizeVETSVCF { input: outdir = outdir, keyfile = keyfile, file = GatherRescoredVcfs.output_vcf }
    call FF.FinalizeToFile as FinalizeVETSTBI { input: outdir = outdir, keyfile = keyfile, file = GatherRescoredVcfs.output_vcf_index }


    output {
        File joint_recalibrated_vcf     = FinalizeVETSVCF.gcs_path
        File joint_recalibrated_vcf_tbi = FinalizeVETSTBI.gcs_path
    }
}

