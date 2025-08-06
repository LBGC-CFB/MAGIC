rule assembly:
    input: 
        bam = os.path.join(config['dirs']['outdir'], "alignment", "{sample}.minimap2.bam")
    output: 
        gtf = os.path.join(config['dirs']['outdir'], "assembly", "{construction}", "{sample}_assembly.gtf")
    params: 
        construction_gtf = lambda wildcards: os.path.join(config['dirs']['constructiondir'], config["samples"][wildcards.sample], f"{config['samples'][wildcards.sample]}.gtf"),
        alignment_mode = '-L' if config['options']['longread'] else ''
    threads:
        config['threads']
    container:
        "docker://aucam/magic:latest"
    shell:
        r"""
        stringtie \
            -p {threads} \
            -G {params.construction_gtf} \
            -o {output.gtf} \
            {input.bam}
        """
        #{params.alignment_mode} \



rule merge:
    input: lambda wildcards: [os.path.join(config['dirs']['outdir'], "assembly", wildcards.construction, f"{sample}_assembly.gtf") for sample in constructions_dic[wildcards.construction]]
    output: 
        gtf_merged=os.path.join(config['dirs']['outdir'], "assembly", "{construction}", "{construction}_merged.gtf")
    params:
        merge_name = lambda wildcards: f"M_{wildcards.construction}",
        construction_gtf = lambda wildcards: os.path.join(config['dirs']['constructiondir'], wildcards.construction, f"{wildcards.construction}.gtf"),
        alignment_mode = '-L' if config['options']['longread'] else ''
    threads:
        config['threads']
    container:
        "docker://aucam/magic:latest"
    shell:
        r"""
        stringtie \
            -p {threads} \
            --merge -l {params.merge_name} \
            -G {params.construction_gtf} \
            -o {output.gtf_merged} \
            {input}
        """
        #{params.alignment_mode} \



rule expression:
    input: 
        bam_file = os.path.join(config['dirs']['outdir'], "alignment", "{sample}.minimap2.bam"),
        transcript_guide = lambda wildcards: os.path.join(config['dirs']['outdir'], "assembly", config["samples"][wildcards.sample], f"{config['samples'][wildcards.sample]}_merged.gtf")
    output: 
        expfile = os.path.join(config['dirs']['outdir'], "expression", "{sample}_expression.gtf")
    params:
        alignment_mode = '-L' if config['options']['longread'] else ''
    threads:
        config['threads']
    container:
        "docker://aucam/magic:latest"
    shell:
        r"""
        stringtie \
            -p {threads} \
            -e \
            -G {input.transcript_guide}  \
            -o {output.expfile} \
            {input.bam_file}
        """
        #{params.alignment_mode} \