version 1.0

import "Structs.wdl"

task BamToFq {
    input {
        File bam
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools sort -n ~{bam} | samtools bam2fq \
            -n \
            -s /dev/null \
            -1 ~{prefix}.end1.fq.gz \
            -2 ~{prefix}.end2.fq.gz \
            -0 ~{prefix}.unpaired.fq.gz
    >>>

    output {
        File fq_end1 = "~{prefix}.end1.fq.gz"
        File fq_end2 = "~{prefix}.end1.fq.gz"
        File fq_unpaired = "~{prefix}.unpaired.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.8"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task Bam2FqPicard {
    input {
        File bam
        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        java -Xms8192m -Xmx30768m -jar /usr/picard/picard.jar \
            SamToFastq \
            INPUT=~{bam} \
            FASTQ=~{prefix}.fastq \
            INTERLEAVE=true \
            NON_PF=true
    >>>

    output {
        File fastq = "~{prefix}.fastq"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }

}

task BwaMem2 {
    input {
        File name_sorted_fastq

        File ref_fasta
        File ref_fasta_index
        File ref_dict
        File ref_0123
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(name_sorted_fastq, "GB"))
                      + 4*ceil(size(ref_fasta, "GB"))
                      + 4*ceil(size(ref_fasta_index, "GB"))
                      + 4*ceil(size(ref_dict, "GB"))
                      + 4*ceil(size(ref_amb, "GB"))
                      + 4*ceil(size(ref_ann, "GB"))
                      + 4*ceil(size(ref_bwt, "GB"))
                      + 4*ceil(size(ref_pac, "GB"))
                      + 4*ceil(size(ref_0123, "GB"))

    command <<<
        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        let np=${np}-1

        # Breakdown of the arguments:
        # -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []
        # -p            Smart  pairing.  If  two  adjacent reads have the same
        #               name, they are considered to form a  read  pair.  This
        #               way, paired-end and single-end reads can be mixed in a
        #               single FASTA/Q stream.
        # -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]
        # -t INT        number of threads [1]
        # -Y            use soft clipping for supplementary alignments

        bwa-mem2 mem -K 100000000 -p -v 3 -t ${np} -Y ~{ref_fasta} ~{name_sorted_fastq} | \
            samtools view -1 - > ~{prefix}.bam

        samtools index -@${np} ~{prefix}.bam
    >>>

    output {
        File bam = "~{prefix}.bam"
        File bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/sr-utils:0.2.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MergeBamAlignment {
    input {
        File aligned_bam
        File unaligned_bam

        File ref_fasta
        File ref_fasta_index
        File ref_dict

        String prefix = "out"

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(aligned_bam, "GB"))
                      + 4*ceil(size(unaligned_bam, "GB"))
                      + 4*ceil(size(ref_fasta, "GB"))
                      + 4*ceil(size(ref_fasta_index, "GB"))
                      + 4*ceil(size(ref_dict, "GB"))

    command <<<
        set -euxo pipefail

        # Make sure we use all our proocesors:
        np=$(cat /proc/cpuinfo | grep ^processor | tail -n1 | awk '{print $NF+1}')
        let np=${np}-1

        java -Dsamjdk.compression_level=2 -Xms8192m -Xmx30768m -jar /usr/picard/picard.jar \
            MergeBamAlignment \
            VALIDATION_STRINGENCY=SILENT \
            EXPECTED_ORIENTATIONS=FR \
            ATTRIBUTES_TO_RETAIN=X0 \
            ATTRIBUTES_TO_REMOVE=NM \
            ATTRIBUTES_TO_REMOVE=MD \
            ALIGNED_BAM=~{aligned_bam} \
            UNMAPPED_BAM=~{unaligned_bam} \
            OUTPUT=~{prefix}.bam \
            REFERENCE_SEQUENCE=~{ref_fasta} \
            SORT_ORDER="unsorted" \
            IS_BISULFITE_SEQUENCE=false \
            ALIGNED_READS_ONLY=false \
            CLIP_ADAPTERS=false \
            MAX_RECORDS_IN_RAM=2000000 \
            ADD_MATE_CIGAR=true \
            MAX_INSERTIONS_OR_DELETIONS=-1 \
            PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
            PROGRAM_RECORD_ID="bwamem" \
            PROGRAM_GROUP_VERSION="${BWA_VERSION}" \
            PROGRAM_GROUP_COMMAND_LINE="bwa-mem2 mem -K 100000000 -p -v 3 -t 15 -Y" \
            PROGRAM_GROUP_NAME="bwa-mem2" \
            UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
            ALIGNER_PROPER_PAIR_FLAGS=true \
            UNMAP_CONTAMINANT_READS=true \
            ADD_PG_TAG_TO_READS=false
    >>>

    output {
        File bam = "~{prefix}.bam"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task MarkDuplicates {
    input {
        File input_bam

        String prefix

        # The program default for READ_NAME_REGEX is appropriate in nearly every case.
        # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
        # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing
        String? read_name_regex

        Float? sorting_collection_size_ratio

        RuntimeAttr? runtime_attr_override
    }

    Int compression_level = 2
    Int java_memory_size_mb = 30768

    Int disk_size = 1 + 4*ceil(size(input_bam, "GB"))

    # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
    # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
    # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"

    command {
        java -Dsamjdk.compression_level=~{compression_level} -Xms~{java_memory_size_mb}m -jar /usr/picard/picard.jar \
        MarkDuplicates \
        INPUT=~{input_bam} \
        OUTPUT=~{prefix}.bam \
        METRICS_FILE=~{prefix}.metrics.txt \
        VALIDATION_STRINGENCY=SILENT \
        ~{"READ_NAME_REGEX=" + read_name_regex} \
        ~{"SORTING_COLLECTION_SIZE_RATIO=" + sorting_collection_size_ratio} \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        ASSUME_SORT_ORDER="queryname" \
        CLEAR_DT="false" \
        ADD_PG_TAG_TO_READS=false
    }

    output {
        File bam = "~{prefix}.bam"
        File metrics = "~{prefix}.metrics.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}


# Generate Base Quality Score Recalibration (BQSR) model
task BaseRecalibrator {
    input {
        File input_bam
        File input_bam_index

        File ref_dict
        File ref_fasta
        File ref_fasta_index

        File known_sites_vcf
        File known_sites_index

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(input_bam, "GB"))
                      + 4*ceil(size(input_bam_index, "GB"))
                      + 2*ceil(size(ref_dict, "GB"))
                      + 2*ceil(size(ref_fasta, "GB"))
                      + 2*ceil(size(ref_fasta_index, "GB"))
                      + 2*ceil(size(known_sites_vcf, "GB"))
                      + 2*ceil(size(known_sites_index, "GB"))

    parameter_meta {
        input_bam: {
            localization_optional: true
        }
    }

    command {
        gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
            -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
            -Xloggc:gc_log.log -Xms5000m -Xmx5500m" \
            BaseRecalibrator \
            -R ~{ref_fasta} \
            -I ~{input_bam} \
            --use-original-qualities \
            -O ~{prefix}.txt \
            --known-sites ~{known_sites_vcf}
    }
    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          16,
        mem_gb:             32,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        1,
        docker:             "us.gcr.io/broad-gotc-prod/picard-cloud:2.26.10"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
    output {
        File recalibration_report = "~{prefix}.txt"
    }
}

