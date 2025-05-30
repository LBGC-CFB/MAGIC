rule assembly:
    input: f"{config['dirs']['outdir']}/alignment/{{sample}}.minimap2.bam",
    output: f"{config['dirs']['outdir']}/assembly/{{construction}}/{{sample}}_assembly.gtf"
    params: 
        construction_gtf = lambda wildcards: f"{config['dirs']['constructiondir']}/{config["samples"][wildcards.sample]}/{config["samples"][wildcards.sample]}.gtf",
        alignment_mode = '-L' if config['options']['longread'] else ''
    threads:
        config['threads']
    container:
        "docker://aucam/magic:latest"
    shell:
        """
        stringtie \
            -p {threads} \
            -G {params.construction_gtf} \
            -o {output} \
            {input}
        """
        #{params.alignment_mode} \



rule merge:
    input: lambda wildcards: [f"{config['dirs']['outdir']}/assembly/{{construction}}/{sample}_assembly.gtf" for sample in constructions_dic[wildcards.construction]]
    output: 
        f"{config['dirs']['outdir']}/assembly/{{construction}}/{{construction}}_merged.gtf"
    params:
        merge_name = f"M_{{construction}}",
        construction_gtf = f"{config['dirs']['constructiondir']}/{{construction}}/{{construction}}.gtf",
        alignment_mode = '-L' if config['options']['longread'] else ''
    threads:
        config['threads']
    container:
        "docker://aucam/magic:latest"
    shell:
        """
        stringtie \
            -p {threads} \
            --merge -l {params.merge_name} \
            -G {params.construction_gtf} \
            -o {output} \
            {input}
        """
        #{params.alignment_mode} \



rule expression:
    input: 
        bam_file = f"{config['dirs']['outdir']}/alignment/{{sample}}.minimap2.bam",
        transcript_guide = lambda wildcards: f"{config['dirs']['outdir']}/assembly/{config["samples"][wildcards.sample]}/{config["samples"][wildcards.sample]}_merged.gtf"
    output: 
        output_exp = f"{config['dirs']['outdir']}/expression/{{sample}}_expression.gtf",
    params:
        alignment_mode = '-L' if config['options']['longread'] else ''
    threads:
        config['threads']
    container:
        "docker://aucam/magic:latest"
    shell:
        """
        stringtie \
            -p {threads} \
            -e \
            -G {input.transcript_guide}  \
            -o {output} \
            {input.bam_file}
        """
        #{params.alignment_mode} \