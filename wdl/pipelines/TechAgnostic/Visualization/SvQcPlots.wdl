version 1.0

import "../../structs/Structs.wdl"

workflow SV_QC_Metrics{

    input{
        Array[File] vcf_samples
    }
    scatter(input_vcf in vcf_samples){
        call bcfQuerySV{
            input:
                pbsv_vcf = input_vcf + ".pbsv.vcf.gz",
                sniffles_vcf = input_vcf + ".sniffles.vcf.gz",
                pav_vcf = input_vcf + ".pav.vcf.gz"
        }
    }


}


task bcfQuerySV{
    input{
        File pbsv_vcf
        File sniffles_vcf
        File pav_vcf
        RuntimeAttr? runtime_attr_override
    }

    String sample_basename = basename(pbsv_vcf, ".pbsv.vcf.gz") # assuming all vcf files have the same basename (sample name)
    String pbsv_stat_out = sample_basename+ ".pbsv.svlen"
    String sniffles_stat_out = sample_basename+ ".sniffles.svlen"
    String pav_stat_out = sample_basename+ ".pav.svlen"

    Int disk_size = size(pbsv_vcf) + size(sniffles_vcf) + size(pav_vcf) + 1000000000 # 1GB buffer
    Int minimal_disk_size = (ceil(size(pbsv_vcf, "GB") + size(sniffles_vcf, "GB") + size(pav_vcf, "GB") ) + 100 ) # 100GB buffer
    Int disk_size = if minimal_disk_size > 100 then minimal_disk_size else 100


    command{
        cat ~{pbsv_vcf} | bcftools query -i '(INFO/SVLEN>49 || INFO/SVLEN<-49) && FILTER=="PASS"' --format "%SVTYPE\t%SVLEN\n" > ~{pbsv_stat_out}
        cat ~{sniffles_vcf} | bcftools query -i '(INFO/SVLEN>49 || INFO/SVLEN<-49) && FILTER=="PASS"' --format "%SVTYPE\t%SVLEN\n" > ~{sniffles_stat_out}
        cat ~{pav_vcf} | bcftools query -i '(INFO/SVLEN>49 || INFO/SVLEN<-49) && FILTER=="PASS"' --format "%SVTYPE\t%SVLEN\n" > ~{pav_stat_out}
    }
    output{
        File pbsv_stat_out = pbsv_stat_out
        File sniffles_stat_out = sniffles_stat_out
        File pav_stat_out = pav_stat_out
    }
        #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             24,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-basic:latest"
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

