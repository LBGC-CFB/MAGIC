rule minimap2_align:
    input: 
        fastq = f"{config['dirs']['indir']}/{{sample}}.fastq" if config['options']['longread'] else [f"{config['dirs']['indir']}/{{sample}}_L001_R1_001.fastq.gz", f"{config['dirs']['indir']}/{{sample}}_L001_R2_001.fastq.gz"],
        ref = lambda wildcards: f"{config['dirs']['constructiondir']}/{config["samples"][wildcards.sample]}/{config["samples"][wildcards.sample]}.fasta",
        bed = lambda wildcards: f"{config['dirs']['constructiondir']}/{config["samples"][wildcards.sample]}/{config["samples"][wildcards.sample]}.bed"
    output: 
        sam = temp(f"{config['dirs']['outdir']}/alignment/{{sample}}.minimap2.sam")
    params: 
        alignment_mode = 'splice' if config['options']['longread'] else 'splice:sr'
    threads: 
        config['threads']
    container:
        "docker://aucam/magic:latest"
    shell:
        """
        minimap2 \
            -ax {params.alignment_mode} \
            -t {threads} \
            {input.ref} \
            -o {output.sam} \
            {input.fastq}
        """
        #--junc-bed {input.bed} --junc-bonus 9 \



rule samtools_sort_index_bam:
    input: f"{config['dirs']['outdir']}/alignment/{{sample}}.minimap2.sam"
    output: 
        bam = f"{config['dirs']['outdir']}/alignment/{{sample}}.minimap2.bam",
        bai = f"{config['dirs']['outdir']}/alignment/{{sample}}.minimap2.bam.bai"
    threads:
        config['threads']
    container: 
        "docker://aucam/magic:latest"
    shell:
        """
        samtools sort \
            -@ {threads} \
            -o {output.bam} \
            {input} && \
                samtools index \
                -@ {threads} \
                {output.bam}
        """