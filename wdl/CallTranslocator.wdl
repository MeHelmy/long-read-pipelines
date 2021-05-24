version 1.0

##########################################################################################
# A workflow that runs Translocator
##########################################################################################

import "tasks/Translocator.wdl" as Tr

workflow Translocator {
    input {
        File aligned_bam
        File ref_map_file

        String prefix
        Float? min_het_af
    }

    Map[String, String] ref_map = read_map(ref_map_file)

    call Tr.Translocator as Translocator {
        input:
            aligned_bam = aligned_bam,
            ref_fasta = ref_map['fasta'],
            prefix = prefix,
            min_het_af = min_het_af
    }

    output {
        translocator_vcf = Translocator.vcf
    }
}
