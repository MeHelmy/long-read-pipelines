version 1.0

import "tasks/Fusilli.wdl" as Fusilli
import "tasks/TrimGalore.wdl" as Trim

workflow FusilliBuildAndCleanSampleGraph {
    input {
        String sample_id

        File reads_fq1
        File? reads_fq2
    }

    call Trim.TrimGalore as TrimGalore {
        input:
            reads_fq1 = reads_fq1,
            reads_fq2 = reads_fq2
    }

    call Fusilli.BuildSampleGraph as BuildSampleGraph {
        input:
            sample_id = sample_id,
            reads_fq1 = TrimGalore.trimmed_fq1,
            reads_fq2 = TrimGalore.trimmed_fq2
    }

    Map[String, Float] clean_meta = read_json(BuildSampleGraph.clean_meta)

    output {
        File trimmed_fq1 = TrimGalore.trimmed_fq1
        File? trimmed_fq2 = TrimGalore.trimmed_fq2
        File sample_graph = BuildSampleGraph.sample_graph
        File sample_graph_colors = BuildSampleGraph.sample_graph_colors
        File kmer_counts = BuildSampleGraph.kmer_counts
        File kmer_spectrum = BuildSampleGraph.kmer_spectrum
        File cleaned_graph = BuildSampleGraph.cleaned_graph
        File cleaned_graph_colors = BuildSampleGraph.cleaned_graph_colors
        Float est_coverage = clean_meta['est_coverage']
        Int prune_threshold = round(clean_meta['threshold'])
        Int num_nodes_before = round(clean_meta['num_nodes_before'])
        Int num_nodes_after = round(clean_meta['num_nodes_after'])
    }
}
