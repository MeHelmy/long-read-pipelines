version 1.0

import "tasks/Quast.wdl" as Quast

workflow QuastEval {
    input {
        File reference
        File? reference_gff
        Array[File] assemblies
        Array[String]? labels

        Boolean icarus = false
    }

    call Quast.QuastBenchmark as Benchmark {
        input:
            ref = reference,
            ref_gff = reference_gff,
            icarus = icarus,

            assemblies = assemblies,
            labels = labels
    }

    output {
        File report_html = Benchmark.report_html
        File report_txt = Benchmark.report_txt
        File report_pdf = Benchmark.report_pdf
        File metrics_tsv = Benchmark.metrics_tsv

        Array[File] plots = Benchmark.plots
        Array[File] icarus_main = Benchmark.icarus_main
        Array[File] icarus_viewers = Benchmark.icarus_viewers

    }
}
