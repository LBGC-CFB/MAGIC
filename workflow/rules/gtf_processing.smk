rule gtf_to_bed:
  input:
     gtf = os.path.join(config['dirs']['constructiondir'], "{construction}", "{construction}.gtf")
  output:
     bed = os.path.join(config['dirs']['constructiondir'], "{construction}", "{construction}.bed")
  run:
    chr = None
    with open(input.gtf, "r") as filin:
      with open(output.bed, "w") as filout:
        filout.write("# MAGIC: conversion of the gtf construction file into a bed file\n")
        for line in filin:
          if not line.startswith("#"):
            line = line.strip().split()
            #print(line)
            if line[2] == "transcript":
              if chr:
                filout.write(f"{chr}\t{startt}\t{endt}\t{id}\t723\t{strand}\t{startt}\t{endt}\t255,0,0\t{str(exon_count)}\t{str(','.join(exon_length))}\t{str(','.join(exon_start))}\n")
              exon_count = 0
              exon_length = []; exon_start = []
              chr = line[0]; startt = str(int(line[3])-1); endt = line[4]; id = line[11].split('"')[1]; strand = line[6]
            elif line[2] == "exon":
              exon_count += 1
              exon_length.append(str(int(line[4]) - int(line[3])+1))
              exon_start.append(str(int(line[3]) - int(startt)-1))



rule concat_gtf:
  input:
     gtf_all = lambda wildcards: expand(os.path.join(config['dirs']['constructiondir'], "{construction}", "{construction}.gtf"), construction=list(constructions_dic.keys()))
  output:
     gtf_merged = os.path.join(config['dirs']['outdir'], "constructions_all.gtf")
  run:
    with open(output.gtf_merged, "w") as filout:
      for file in input.gtf_all:
        with open(file, "r") as filin:
          content = filin.read()
          if not content.endswith("\n"):
            content += "\n"
          filout.write(content)