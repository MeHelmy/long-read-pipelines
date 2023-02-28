version 1.0

#######################################################
# This pipeline calls small variants using DeepVariant.
#######################################################

import "Structs.wdl"


workflow DeepVariant {

    meta {
        description: "Workflow for getting VCF and gVCF from DeepVariant. Note VCF is un-phased."
    }

    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        Int pepper_threads
        Int pepper_memory

        Int dv_threads
        Int dv_memory

        String zones = "us-central1-b us-central1-c"
    }

    parameter_meta {
        # when running large scale workflows, we sometimes see errors like the following
        #   A resource limit has delayed the operation: generic::resource_exhausted: allocating: selecting resources: selecting region and zone:
        #   no available zones: 2763 LOCAL_SSD_TOTAL_GB (738/30000 available) usage too high
        zones: "select which zone (GCP) to run this task"
    }

    call DV as deep_variant {
        input:
            bam = bam,
            bai = bai,
            ref_fasta = ref_fasta,
            ref_fasta_fai = ref_fasta_fai,
            threads = dv_threads,
            memory = dv_memory,
            zones = zones
    }

    output {
        File VCF        = deep_variant.VCF
        File VCF_tbi    = deep_variant.VCF_tbi

        File gVCF       = deep_variant.gVCF
        File gVCF_tbi   = deep_variant.gVCF_tbi
    }
}

task DV {

    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai

        Int threads
        Int memory
        String zones

        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bam, ".bam") + ".deepvariant"
    String output_root = "/cromwell_root/dv_output"

    Int bam_sz = ceil(size(bam, "GB"))
    Boolean is_big_bam = bam_sz > 100
    Int inflation_factor = if (is_big_bam) then 10 else 5
    Int minimal_disk = 1000
	Int disk_size = if inflation_factor * bam_sz > minimal_disk then inflation_factor * bam_sz else minimal_disk

    command <<<
        set -euxo pipefail

        export MONITOR_MOUNT_POINT="/cromwell_root"
        wget https://raw.githubusercontent.com/broadinstitute/long-read-pipelines/jts_kvg_sp_malaria/scripts/monitor/legacy/vm_local_monitoring_script.sh -O monitoring_script.sh
        chmod +x monitoring_script.sh
        ./monitoring_script.sh &> resources.log &
        monitoring_pid=$!

        num_core=$(cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l)

        mkdir -p "~{output_root}"

        /opt/deepvariant/bin/run_deepvariant \
            --model_type=WGS \
            --ref=~{ref_fasta} \
            --reads=~{bam} \
            --output_vcf="~{output_root}/~{prefix}.vcf.gz" \
            --output_gvcf="~{output_root}/~{prefix}.g.vcf.gz" \
            --num_shards="${num_core}"

        find "~{output_root}/" -print | sed -e 's;[^/]*/;|____;g;s;____|; |;g' \
            > "~{output_root}/dir_structure.txt"

        kill $monitoring_pid
    >>>

    output {
        File output_dir_structure = "~{output_root}/dir_structure.txt"

        File VCF        = "~{output_root}/~{prefix}.vcf.gz"
        File VCF_tbi    = "~{output_root}/~{prefix}.vcf.gz.tbi"

        File gVCF       = "~{output_root}/~{prefix}.g.vcf.gz"
        File gVCF_tbi   = "~{output_root}/~{prefix}.g.vcf.gz.tbi"

        File visual_report_html = "~{output_root}/~{prefix}.visual_report.html"

        File monitoring_log = "resources.log"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          threads,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             "google/deepvariant:1.4.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
