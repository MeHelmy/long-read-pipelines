version 1.0

import "tasks/Structs.wdl"
import "tasks/FunctionalAnnotation.wdl" as FUNK
import "tasks/Finalize.wdl" as FF

workflow ExpandedDrugResistanceMarkerExtraction {
    input {
        String sample_name

        File vcf
        File snpeff_db
        File protein_drug_resistance_list
        File gene_drug_resistance_list

        String dir_prefix
        String gcs_out_root_dir

        Boolean do_functional_annotation = true
    }

    String outdir = sub(gcs_out_root_dir, "/$", "") + "/ExpandedDrugResistanceMarkerExtraction/~{dir_prefix}"

    if (do_functional_annotation) {
        call FUNK.FunctionallyAnnotateVariants { input: vcf = vcf, snpeff_db = snpeff_db }
    }

    call CallDrugResistanceMutations {
        input:
            vcf = select_first([FunctionallyAnnotateVariants.annotated_vcf, vcf]),
            protein_drug_resistance_list = protein_drug_resistance_list,
            gene_drug_resistance_list = gene_drug_resistance_list,
            prefix = sample_name
    }

    # Finalize data
    String dir = outdir + "/reports"

    call FF.FinalizeToFile as FinalizeDRReportAllMarkers { input: outdir = dir, file = CallDrugResistanceMutations.all_markers }
    call FF.FinalizeToFile as FinalizeDRReportProteinMarkers { input: outdir = dir, file = CallDrugResistanceMutations.protein_coding_markers }

    if (do_functional_annotation) {
        call FF.FinalizeToFile as FinalizeAnnotatedVCF { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.annotated_vcf]) }
        call FF.FinalizeToFile as FinalizeAnnotatedVCFIndex { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.annotated_vcf_index]) }
        call FF.FinalizeToFile as FinalizeSnpEffSummary { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.snpEff_summary]) }
        call FF.FinalizeToFile as FinalizeSnpEffGenes { input: outdir = dir, file = select_first([FunctionallyAnnotateVariants.snpEff_genes]) }
    }

    output {
        File drug_res_report_all = FinalizeDRReportAllMarkers.gcs_path
        File drug_res_report_prot_only = FinalizeDRReportProteinMarkers.gcs_path

        File? annotated_vcf = FinalizeAnnotatedVCF.gcs_path
        File? annotated_vcf_index = FinalizeAnnotatedVCFIndex.gcs_path
        File? snpEff_summary = FinalizeSnpEffSummary.gcs_path
        File? snpEff_genes = FinalizeSnpEffGenes.gcs_path
    }
}

task CallDrugResistanceMutations {
    input {
        File vcf
        File protein_drug_resistance_list
        File gene_drug_resistance_list

        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 10*ceil(size([vcf, protein_drug_resistance_list, gene_drug_resistance_list], "GB"))

    command <<<
        set -euxo pipefail

        python3 <<CODE

        import gzip

        gene_info = []
        with open("~{gene_drug_resistance_list}", 'r') as f:
            for line in f:
                gene_name, gene_id = line.strip().split()
                gene_info.append((gene_name, gene_id))

        print(f"Gene info:")
        for gene_name, gene_id in gene_info:
            print(f"\t{gene_name}\t{gene_id}")
        print()

        p_change_marker_info = []
        with open("~{protein_drug_resistance_list}", 'r') as f:
            for line in f:
                gene_name, gene_id, prot_change = line.strip().split()
                p_change_marker_info.append((gene_name, gene_id, prot_change))

        print(f"Protein Change Drug Res Info:")
        for gene_name, gene_id, prot_change in p_change_marker_info:
            print(f"\t{gene_name}\t{gene_id}\t{prot_change}")
        print()

        gene_drug_report_all = "~{prefix}.expanded_drug_report.ALL.tsv"
        gene_drug_report_prot = "~{prefix}.expanded_drug_report.PROTEIN_CHANGES_ONLY.tsv"

        def make_ann_dict(ann_info, annotation_fields):
            ann_dicts = []

            vcf_annotations = ann_info.split(",")
            for v in vcf_annotations:
                ann_dict = dict()
                vcf_ann_fields = v.split("|")

                if len(vcf_ann_fields) != len(annotation_fields):
                    print(vcf_ann_fields)
                    raise RuntimeError(f"Error: Field length of annotation is not the length of the number of fields: {len(vcf_ann_fields)} != expected({len(annotation_fields)})")
                for i in range(len(vcf_ann_fields)):
                    ann_dict[annotation_fields[i]] = vcf_ann_fields[i]

                ann_dicts.append(ann_dict)

            return ann_dicts

        def ann_impact_sort_key(x):
            s = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
            return s.index(x[0])

        annotations = []
        with gzip.open("~{vcf}", 'rt') as f:
            for line in f:
                if line.startswith("##INFO=<ID=ANN"):
                    needle = 'Functional annotations: '
                    i1 = line.find(needle)
                    i2 = line.find(' ">')
                    annotation_fields = line[i1+len(needle):i2].replace("'", "").split(" | ")
                    continue
                elif line.startswith("#"):
                    continue
                chrom, pos, idd, ref, alt, qual, flt, info, gt_f, gt_d = line.strip().split("\t")
                infos = {}
                for i in info.split(";"):
                    if "=" in i:
                        k, v = i.split("=")
                        infos[k] = v
                ann_dicts = make_ann_dict(infos["ANN"], annotation_fields)

                # Now check if our Gene level drug markers:
                gene_annotations = []
                for gene_name, gene_id in gene_info:
                    for ann_dict in ann_dicts:
                        if gene_id in ann_dict["Gene_Name"] or gene_id in ann_dict["Gene_ID"]:
                            # is it protein coding?
                            if len(ann_dict["HGVS.p"]) > 0:
                                gene_annotations.append((ann_dict["Annotation_Impact"], gene_name, gene_id, ann_dict['HGVS.p'], ann_dict["Annotation"]))
                            else:
                                gene_annotations.append((ann_dict["Annotation_Impact"], gene_name, gene_id, ann_dict['HGVS.c'], ann_dict["Annotation"]))

                if len(gene_annotations) > 0:
                    # Sort gene annotations to make the highest / worst effect at the front:
                    gene_annotations = sorted(gene_annotations, key=ann_impact_sort_key)
                    # Only output the greatest effect annotation
                    # only output the fields from gene_name onward - we dont need the impact
                    annotations.append(gene_annotations[0][1:])

                # Now check for our protein change string drug markers:
                for gene_name, gene_id, prot_change in p_change_marker_info:
                    for ann_dict in ann_dicts:
                        if gene_id in ann_dict["Gene_Name"] or gene_id in ann_dict["Gene_ID"]:
                            if len(ann_dict["HGVS.p"]) > 0 and ann_dict["HGVS.p"] == prot_change:
                                annotations.append((gene_name, gene_id, ann_dict['HGVS.p'], ann_dict["Annotation"]))

        with open(gene_drug_report_all, 'w') as f:
            f.write(f"Gene_Name\tGene_ID\tChange_String\tAnnotation_Type\n")
            for a in annotations:
                if a[3] != "synonymous_variant":
                    f.write("\t".join(a))
                    f.write("\n")

        with open(gene_drug_report_prot, 'w') as f:
            f.write(f"Gene_Name\tGene_ID\tChange_String\tAnnotation_Type\n")
            for a in annotations:
                if a[2].startswith("p.") and a[3] != "synonymous_variant":
                    f.write("\t".join(a))
                    f.write("\n")
        CODE

    >>>

    output {
        File all_markers = "~{prefix}.expanded_drug_report.ALL.tsv"
        File protein_coding_markers = "~{prefix}.expanded_drug_report.PROTEIN_CHANGES_ONLY.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/snpeff:5.1d--hdfd78af_0"
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
