version 1.0

workflow MosdepthCoverageStats {

    meta {
        description: "Calculate coverage statistics using mosdepth."
    }
    parameter_meta {
        aligned_bam: "Aligned BAM file."
        aligned_bai: "Aligned BAM index file."
        bed_file: "BED file containing regions of interest."
        bin_length: "Length of bins to use for coverage calculation."
        preemptible_tries: "Number of times to retry a preempted task."
    }

    input {
        File aligned_bam
        File aligned_bai
        File? bed_file

        Int bin_length = 1000

        # Runtime parameters
        Int? preemptible_tries = 3
    }

    call MosDepthOverBed {
        input:
            bam = aligned_bam,
            bai = aligned_bai,
            bed = bed_file,
            bin_length = bin_length,
            preemptible_tries = preemptible_tries,
    }

    String basename = basename(aligned_bam, ".bam")

    call CoverageStats {
        input:
            mosdepth_regions = MosDepthOverBed.regions,
            basename_input = basename,
            preemptible_tries = preemptible_tries,
    }

    output {
        File cov_stat_summary_file      = CoverageStats.cov_stat_summary_file
        Map[String, String] cov_stat_summary = CoverageStats.cov_stat_summary
    }

}


task MosDepthOverBed {

    meta {
        description: "Calculate coverage using mosdepth."
    }
    parameter_meta {
        bam: "Aligned BAM file."
        bai: "Aligned BAM index file."
        bed: "BED file containing regions of interest."
        bin_length: "Length of bins to use for coverage calculation."
        preemptible_tries: "Number of times to retry a preempted task."
    }

    input {
        File bam
        File bai
        File? bed
        Int bin_length

        # Runtime parameters
        Int? preemptible_tries = 3
    }

    String mosdepth_by_region = select_first([bed, bin_length])

    Int disk_size = 2 * ceil(size(bam, "GB") + size(bai, "GB"))
    String basename = basename(bam, ".bam")
    String prefix = "~{basename}.coverage_over_bed"

    command {
        set -euxo pipefail

        # Create symbolic links for bam and bai in the current working directory
        ln -s ~{bam} ./~{basename}.bam
        ln -s ~{bai} ./~{basename}.bai

        mosdepth -t 4 -b ~{mosdepth_by_region} -n -x -Q 1 ~{prefix} ./~{basename}.bam
    }

    output {
        File global_dist      = "~{prefix}.mosdepth.global.dist.txt"
        File region_dist      = "~{prefix}.mosdepth.region.dist.txt"
        File regions          = "~{prefix}.regions.bed.gz"
        File regions_csi      = "~{prefix}.regions.bed.gz.csi"
    }

    runtime {
        cpu:                    4
        memory:                 8 + " GiB"
        disks: "local-disk " +  100 + " HDD"
        preemptible:            preemptible_tries
        maxRetries:             1
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-mosdepth:latest"
    }
}

task CoverageStats {

    meta {
        description: "Calculate coverage statistics from mosdepth output."
    }
    parameter_meta {
        mosdepth_regions: "Mosdepth output file containing coverage values."
        cov_col: "Column holding the coverage values."
        basename_input: "Basename to use for output files."
        preemptible_tries: "Number of times to retry a preempted task."
    }

    input {
        File mosdepth_regions
        Int cov_col = 4 # column holding the coverage values
        String? basename_input

        # Runtime parameters
        Int? preemptible_tries = 3
    }

    Int disk_size = 2*ceil(size(mosdepth_regions, "GB"))
    String basename = select_first([basename_input, basename(mosdepth_regions)])
    String prefix = "~{basename}.coverage_over_bed"

    command <<<
        set -euxo pipefail

        # Use Datamash to calculate summary statistics
        zcat ~{mosdepth_regions} | datamash -H \
        mean ~{cov_col} \
        q1 ~{cov_col} \
        median ~{cov_col} \
        q3 ~{cov_col} \
        iqr ~{cov_col} \
        sstdev ~{cov_col} \
        mad ~{cov_col} > ~{prefix}.cov_stat_summary.txt

        cat ~{prefix}.cov_stat_summary.txt

        # Replace floating-point numbers suffix ing field names with '_coverage
        sed -i '1s/([0-9]*\.[0-9]*)/_coverage/g' ~{prefix}.cov_stat_summary.txt

        # Calculate covrage percentage with greater than 4x coverage
        total_bases=$(zcat ~{mosdepth_regions} | wc -l)
        bases_above_4x=$(zcat ~{mosdepth_regions} | awk -v cov_col=~{cov_col} '$cov_col > 4' | wc -l)
        percent_above_4x=$(python3 -c "print($bases_above_4x/$total_bases)")

        # Append the percentage to the summary tsv file
        sed -i "1s/$/\tpercent_above_4x_coverage/" ~{prefix}.cov_stat_summary.txt
        sed -i "2s/$/\t$percent_above_4x/" ~{prefix}.cov_stat_summary.txt

        # Extract field names from the header
        header=$(head -n 1 ~{prefix}.cov_stat_summary.txt)
        IFS=$'\t' read -ra field_names <<< "$header"

        # Extract data values
        data=$(tail -n 1 ~{prefix}.cov_stat_summary.txt)
        data_values=($(echo "$data" | awk '{ for(i=1; i<=NF; i++) print $i }'))

        # Create the JSON object
        json="{"
        for (( i=0; i<${#field_names[@]}; i++ )); do
            json+="\"${field_names[i]}\":${data_values[i]}"
            if [[ $i -lt $(( ${#field_names[@]} - 1 )) ]]; then
                json+=", "
            fi
        done
        json+="}"

        # Write the JSON object to the output file
        echo "$json" > ~{prefix}.cov_stat_summary.json
    >>>

    output {
        File cov_stat_summary_file      = "~{prefix}.cov_stat_summary.txt"
        Map[String, Float] cov_stat_summary = read_json("~{prefix}.cov_stat_summary.json")
    }

    runtime {
        cpu:                    2
        memory:                 4 + " GiB"
        disks: "local-disk " +  100 + " HDD"
        preemptible:            preemptible_tries
        maxRetries:             1
        docker:                 "us.gcr.io/broad-dsp-lrma/lr-metrics:0.1.11"
    }
}
