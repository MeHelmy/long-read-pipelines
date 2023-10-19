version 1.0

import "../../../tasks/Utility/Utils.wdl"
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Phasing/StatisticalPhasing.wdl" as StatPhase
import "../../../tasks/Phasing/Hiphase.wdl"
import "../../../tasks/Phasing/SplitJointCallbySample.wdl"


workflow HybridPhase {
    meta{
        description : "..."
    }
    parameter_meta {
        one_chr_bams_from_all_samples:  "GCS path to subset BAM files"
        one_chr_bais_from_all_samples:  "GCS path to subset BAI file indices"
        one_chr_joint_vcf:  "path to subset joint vcf per chromosome"
        one_chr_joint_vcf_tbi:  "path to subset joint vcf index per chromosome"
        reference: "path to reference genome fasta file"
        reference_index: "path to reference genome index fai file"
        genetic_mapping_tsv_for_shapeit4: "path to the tsv file for the genetic mapping file address per chromosome"
        chromosome: "string for chromosome to be processed"
        prefix: "output file prefix, usually using chromosome"
        num_t: "integer for threads"
    }

    input {
        Array[File] one_chr_bams_from_all_samples
        Array[File] one_chr_bais_from_all_samples
        File one_chr_joint_vcf
        File one_chr_joint_vcf_tbi
        File one_chr_sv_vcf
        File one_chr_sv_vcf_tbi
        File reference
        File reference_index
        File genetic_mapping_tsv_for_shapeit4
        String chromosome
        String prefix
        Int num_t
    }
    
    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)

    scatter (bam_bai in zip(one_chr_bams_from_all_samples, one_chr_bais_from_all_samples)) {
        File bam = bam_bai.left
        File bai = bam_bai.right
        
        call Utils.InferSampleName { input: 
            bam = bam, 
            bai = bai}
        String sample_id = InferSampleName.sample_name

        call SplitJointCallbySample.SplitVCFbySample as Split { input:
            joint_vcf = one_chr_joint_vcf,
            region = chromosome,
            samplename = sample_id
        }

        call Hiphase.Hiphase as hiphase { input:
            bam = bam,
            bai = bai,
            unphased_snp_vcf = Split.single_sample_vcf,
            unphased_snp_tbi = Split.single_sample_vcf_tbi,
            unphased_sv_vcf = one_chr_sv_vcf,
            unphased_sv_tbi = one_chr_sv_vcf_tbi,
            ref_fasta = reference,
            ref_fasta_fai = reference_index,
            prefix = prefix
        }
    }
        
    call VU.MergePerChrVcfWithBcftools as MergeAcrossSamplesSNPs { input:
        vcf_input = hiphase.phased_snp_vcf,
        tbi_input = hiphase.phased_snp_vcf_tbi,
        pref = prefix
    }
    call VU.MergePerChrVcfWithBcftools as MergeAcrossSamplesSVs { input:
        vcf_input = hiphase.phased_sv_vcf,
        tbi_input = hiphase.phased_sv_vcf_tbi,
        pref = prefix
    }
    call StatPhase.Shapeit4 as scaffold { input:
        vcf_input = MergeAcrossSamplesSNPs.merged_vcf,
        vcf_index = MergeAcrossSamplesSNPs.merged_tbi,
        mappingfile = genetic_mapping_dict[chromosome],
        region = chromosome,
        num_threads = num_t
    }
    call StatPhase.Shapeit4_phaseSVs as SVphase { input:
        vcf_input = MergeAcrossSamplesSVs.merged_vcf,
        vcf_index = MergeAcrossSamplesSVs.merged_tbi,
        scaffold_vcf = scaffold.scaffold_vcf,
        mappingfile = genetic_mapping_dict[chromosome],
        region = chromosome,
        num_threads = num_t
    }

    output{
        File phased_scaffold = SVphase.final_phased_vcf
    }
}
