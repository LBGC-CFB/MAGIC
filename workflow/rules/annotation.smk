rule SOSTAR:
	input:
		exp_file = expand(os.path.join(config['dirs']['outdir'], "expression", "{sample}_expression.gtf"), sample=config['samples']),
		transcripts_ref=os.path.join(config['dirs']['outdir'], "constructions_all.gtf")
	output: 
		output_file=os.path.join(config['dirs']['outdir'], "MAGIC_SOSTAR_annotation_table_results.xlsx")
	params:
		input_folder=os.path.join(config['dirs']['outdir'], "expression"),
		output_folder = config['dirs']['outdir']
	shell:
		r"""
		python3 scripts/MAGIC_SOSTAR.py \
			-I {params.input_folder} \
			-G {input.transcripts_ref} \
			-O {params.output_folder}
		"""