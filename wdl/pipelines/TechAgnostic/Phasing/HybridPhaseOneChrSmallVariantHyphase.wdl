version 1.0

import "../../../tasks/Utility/Utils.wdl" as U
import "../../../tasks/Utility/VariantUtils.wdl" as VU
import "../../../tasks/Phasing/StatisticalPhasing.wdl" as StatPhase
import "../../../tasks/Phasing/Hiphase.wdl"
import "../../../tasks/Phasing/SplitJointCallbySample.wdl"
import "../../../tasks/Utility/Finalize.wdl" as FF


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
        Array[File] wholegenome_bams_from_all_samples
        Array[File] wholegenome_bais_from_all_samples
        Array[File] wholegenome_sv_vcf
        Array[File] wholegenome_sv_vcf_tbi

        File one_chr_joint_vcf
        File one_chr_joint_vcf_tbi
        File reference
        File reference_index
        File genetic_mapping_tsv_for_shapeit4
        String chromosome
        String prefix
        String gcs_out_root_dir
        Int num_t
    }



    Map[String, String] genetic_mapping_dict = read_map(genetic_mapping_tsv_for_shapeit4)
    Int data_length = length(wholegenome_bams_from_all_samples)
    Array[Int] indexes= range(data_length)

    scatter (idx in indexes)  {
        File all_chr_bam = wholegenome_bams_from_all_samples[idx]
        File all_chr_bai = wholegenome_bais_from_all_samples[idx]
        File pbsv_vcf = wholegenome_sv_vcf[idx]
        File pbsv_vcf_tbi = wholegenome_sv_vcf_tbi[idx]

        ####### Subset Bam and SVs####
        call U.SubsetBam as SubsetBam { input:
            bam = all_chr_bam,
            bai = all_chr_bai,
            locus = chromosome
        }
        call VU.SubsetVCF as SubsetSVs { input:
            vcf_gz = pbsv_vcf,
            vcf_tbi = pbsv_vcf_tbi,
            locus = chromosome
        }
        ####### DONE ######


        call U.InferSampleName { input: 
            bam = all_chr_bam, 
            bai = all_chr_bai
        }
        String sample_id = InferSampleName.sample_name

        call SplitJointCallbySample.SplitVCFbySample as SP { input:
            joint_vcf = one_chr_joint_vcf,
            region = chromosome,
            samplename = sample_id
        }

        call Hiphase.HiphaseSVs as HP { input:
            bam = SubsetBam.subset_bam,
            bai = SubsetBam.subset_bai,
            unphased_snp_vcf = SP.single_sample_vcf,
            unphased_snp_tbi = SP.single_sample_vcf_tbi,
            unphased_sv_vcf = SubsetSVs.subset_vcf,
            unphased_sv_tbi = SubsetSVs.subset_tbi,
            ref_fasta = reference,
            ref_fasta_fai = reference_index,
            samplename = sample_id

        }
    }
    
    ### phase small variants as scaffold  
    call VU.MergePerChrVcfWithBcftools as MergeAcrossSamples { input:
        vcf_input = HP.phased_snp_vcf,
        tbi_input = HP.phased_snp_vcf_tbi,
        pref = prefix
    }
    
    call StatPhase.Shapeit4 as Shapeit4scaffold { input:
        vcf_input = MergeAcrossSamples.merged_vcf,
        vcf_index = MergeAcrossSamples.merged_tbi,
        mappingfile = genetic_mapping_dict[chromosome],
        region = chromosome,
        num_threads = num_t
    }
    call FF.FinalizeToFile as Finalizescaffold {
        input: outdir = gcs_out_root_dir, file = Shapeit4scaffold.scaffold_vcf
    }
    ##### phase structural variants
    call VU.MergePerChrVcfWithBcftools as MergeAcrossSamplesSVs { input:
        vcf_input = HP.phased_sv_vcf,
        tbi_input = HP.phased_sv_vcf_tbi,
        pref = prefix
    }
    call StatPhase.Shapeit4_phaseSVs as Shapeit4SVphase { input:
        vcf_input = MergeAcrossSamplesSVs.merged_vcf,
        vcf_index = MergeAcrossSamplesSVs.merged_tbi,
        scaffold_vcf = Shapeit4scaffold.scaffold_vcf,
        mappingfile = genetic_mapping_dict[chromosome],
        region = chromosome,
        num_threads = num_t
    }
    call FF.FinalizeToFile as FinalizeSVs {
        input: outdir = gcs_out_root_dir, file = Shapeit4SVphase.final_phased_vcf
    }
    output{
        File shapeit4scaffold = Shapeit4scaffold.scaffold_vcf
        File resource = Shapeit4scaffold.resouce_monitor_log
        File finalscaffold = Shapeit4SVphase.final_phased_vcf
        

    }
}
