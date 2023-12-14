version 1.0

workflow Bifrost{
    meta{
        description: "a workflow that construct bifrost colored de bruijn graph"
    }
    input{
        Array[File] input_fastas
        File reference_fasta
        String outputprefix
        Int k
    }
    call construct{input: fas = input_fastas, ref = reference_fasta, outputpref = outputprefix, kmersize = k}
    output{
        Array[File] output_files = construct.graph_files
        # File graph = construct.graph
        # File graph_index = construct.graph_index
        # File colors = construct.color_file
        # File inputfasta = construct.fasta

    }
}

task construct{
    input{
        Array[File] fas
        File ref
        String outputpref
        Int num_threads = 16
        Int kmersize
    }
    command <<<
    set -x pipefail
    Bifrost build -t ~{num_threads} -k ~{kmersize} -i -d -c -s ~{sep=" -s " fas} -r ~{ref} -o ~{outputpref}_Bfrost_graph
    cat ~{sep=" " fas} > all.fasta
    >>>

    output{
        Array[File] graph_files = glob("*")
        # File color_file="~{outputpref}_Bfrost_graph.bfg_colors"
        # File graph = "~{outputpref}_Bfrost_graph.gfa"
        # File graph_index = "~{outputpref}_Bfrost_graph.bfi"
        # File fasta = "all.fasta"
    }

    Int disk_size = 100 

    runtime {
        cpu: num_threads
        memory: "64 GiB"
        disks: "local-disk " + disk_size + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 10
        preemptible: 2
        maxRetries: 1
        docker: "us-central1-docker.pkg.dev/broad-dsp-lrma/fusilli/fusilli:devel" # "hangsuunc/bifrost:v1"
    }
}