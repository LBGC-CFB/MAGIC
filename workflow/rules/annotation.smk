rule SOSTAR:
	input:
		exp_file = expand(f"{config['dirs']['outdir']}/expression/{{sample}}_expression.gtf", sample=config['samples']),
		transcripts_ref = f"{config['dirs']['outdir']}/constructions_all.gtf"
	output: 
		output_file = f"{config['dirs']['outdir']}/MAGIC_SOSTAR_annotation_table_results.xlsx",
	params:
		input_folder = f"{config['dirs']['outdir']}/expression",
		output_folder = config['dirs']['outdir']
	shell:
		"""
		python3 scripts/MAGIC_SOSTAR.py \
			-I {params.input_folder} \
			-G {input.transcripts_ref} \
			-O {params.output_folder}
		"""