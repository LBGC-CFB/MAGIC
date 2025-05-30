import os
import time
import glob
import argparse
import pandas as pd
import numpy as np



def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-I', '--input', dest='input', type=str,
                        help='/path/to/gtf files')
    parser.add_argument('-O', '--output', dest='output', type=str,
                        help='/path/to/output directory')
    parser.add_argument('-G', '--gtf_coord', dest='gtf_coord', type=str,
                        help='gtf annotation file of construction')
    return parser.parse_args() 



def get_ex_nb(annot):
    """
    Get exon number of a splicing event according to our defined annotation
    """
    annot = str(annot)
    if "exo" in annot:
        return int(annot.split("exo")[0])
    elif "int" in annot:
        return int(annot.split("int")[0])
    elif annot.startswith("Δ"):
        if "," in annot:
            annot = annot.split(",")[0]
            get_ex_nb(annot)
        if "p" in annot or "q" in annot:
            i_start = 1; i_end = annot.index("p") if "p" in annot else annot.index("q")
            return int(annot[i_start:i_end])
        else:
            i_start = 1; i_end = annot.index("(") if "(" in annot else None
            return int(annot[i_start:i_end])
    elif annot.startswith("▼"):
        if "," in annot:
            annot = annot.split(",")[1]
            get_ex_nb(annot)
        if "q" in annot:
            i_start = 1; i_end = annot.index("q")
            return int(annot[i_start:i_end])
        elif "p" in annot:
            i_start = 1; i_end = annot.index("p")
            return int(annot[i_start:i_end])+1
        else:
            i_start = 1; i_end = annot.index("(") if "(" in annot else None
            return int(annot[i_start:i_end])
    elif annot.isdigit():
        return int(annot)


def check_intro(lst, dico_bed, gene):
    """
    Check if an intronization of exonic sequence is present in the intermediate list of events found and if so, return the correct annotation for that event
    """
    if lst:
        if len(lst)>1 and get_ex_nb(lst[0]) == get_ex_nb(lst[1]) and (lst[0][0] == lst[1][0] == "Δ"):
            if ("," in lst[0]) and ("," in lst[1]):
                return lst
            elif ("," in lst[0]) or("," in lst[1]):
                ex_nb = get_ex_nb(lst[0])
                ex_length = int(dico_bed[gene][ex_nb-1][1]) - int(dico_bed[gene][ex_nb-1][0]) +1
                for ele in lst[:2]:
                    if "," in ele:
                        ele=ele.split(",")
                        left=int(ele[0][ele[0].index("(")+1:ele[0].index(")")])
                        keep=ele[1]
                    else:
                        simple=int(ele[ele.index("(")+1:ele.index(")")])
                nb_p = ex_length - left
                nb_q = ex_length - simple
                nb_int = left - nb_q
                value = str(ex_nb) + "int(" + str(nb_int) +")[q"+str(nb_q)+",p"+str(nb_p)+"],"+keep
            else:
                ex_nb = get_ex_nb(lst[0])
                ix = [dico_bed[gene].index(ex) for ex in dico_bed[gene] if str(ex_nb) == ex[2]][0]
                ex_length = int(dico_bed[gene][ix][1]) - int(dico_bed[gene][ix][0]) +1
                nb_p = int(lst[0][lst[0].index("(")+1:lst[0].index(")")] if "p" in lst[0] else lst[1][lst[1].index("(")+1:lst[1].index(")")])
                nb_q = int(lst[0][lst[0].index("(")+1:lst[0].index(")")] if "q" in lst[0] else lst[1][lst[1].index("(")+1:lst[1].index(")")])
                nb_int = ex_length - ((ex_length-nb_p)+(ex_length-nb_q))
                value = str(ex_nb) + "int(" + str(nb_int) +")[q"+str(nb_p-nb_int)+",p"+str(nb_q-nb_int)+"]"
            if len(lst) > 2:
                return [value] + check_intro(lst[2:], dico_bed, gene)
            else:
                return [value]
        elif len(lst) == 1:
            return [lst[0]]
        else:
            return [lst[0]] + check_intro(lst[1:], dico_bed, gene)
    else:
        return lst


def check_next(lst, previous_values = []):
    """
    Check if the events in the intermediate list of found events are contigous, and if so goup these events in a new annotation name.
    """
    # If the following element is to be grouped with the current element
    if len(lst)>1 and get_ex_nb(lst[0])+1 == get_ex_nb(lst[1]) and (lst[0][0] == lst[1][0] == "Δ") and ("q" not in lst[1]) and ("p" not in lst[0]) and (("p" or "q" not in lst[1]) and ("(" not in lst[1])) and (("(" not in lst[0]) and ("p" or "q" not in lst[0])):
        # If "p" or "q" in the current or next element, add it to the list of elements to keep
        if len(lst) >1 and ("p" in lst[0] or "q" in lst[0]):
            return check_next(lst[1:], previous_values + [lst[0]])
        # Call the recursive function with the following element if previous values to keep
        elif(previous_values):
            return check_next(lst[1:], previous_values)
        # Else take the current element as the element to keep
        else:
            return check_next(lst[1:], [lst[0]])
    
    # If the list contains at least 2 elements and previous values to keep, call the recursive function with the next element and build the display of the current element
    elif len(lst)>1 and previous_values:
        # Take the value of the elements to keep, without the initial Δ and add a "-" then the value of the current element
        value = "Δ(" + ",".join(map(lambda x: x[1:], previous_values)) + "-" + str(lst[0][1:]) + ")"
        return [value] + check_next(lst[1:])
    
    # If no element to keep, display the current element
    elif len(lst)>1:
        return [lst[0]] + check_next(lst[1:])
    
    # If last element, build according to the previous conditions
    elif len(lst) == 1 and previous_values:
        value = "Δ(" + ",".join(map(lambda x: x[1:], previous_values)) + "-" + str(lst[0][1:]) + ")"
        return [value] 
    else:
        return [lst[0]]



def annotate_isoforms(indir, ref_file, outdir):
    """
    Get annotation per isoform in a final table containing: isoform ID, isoform coordinates, gene name, reference annotation, descriptive annotation, expression values per patient and number of occurrences in the patient cohort
    """
    #get dictionnary from bed transcript reference in form
    #dico = {gene: [[start_exon1, end_exon1, nb_exon1], [start_exon2, end_exon2, nb_exon2], [start_exon3, end_exon3, nb_exon3]], ...}
    dico_bed = {}
    dico_bed_UTR = {}
    with open(ref_file, "r") as gtfile:
        for line in gtfile:
            lines = line.split()
            construction = lines[0]
            if lines[2] == "transcript":
                dico_bed.setdefault(construction, [])
                dico_bed_UTR.setdefault(construction, [])
            if lines[2] == "exon":
                dic_attr = {attr.split()[0]: attr.split()[1].replace('"',"") for attr in line.split("\t")[-1].split(";") if attr.split()} 
                dico_bed[construction].append([lines[3], lines[4], dic_attr["exon_number"]])
            if lines[2] == "UTR":
                dico_bed_UTR[construction].append([lines[3], lines[4]])
    
    #get gtf files list
    os.chdir(indir)
    gtf_files = sorted(glob.glob("*.gtf"))
    
    #create final dataframe to save
    df = pd.DataFrame.from_dict({"transcript_id": [], "construction": [], "gene": [], "annot_ref": [], "annot_find": []})
    
    #get expression values from StringTie gtf files
    samples = []
    for gtf in gtf_files:
        filename = gtf.split("_exp")[0]
        samples.append(filename)
        gene = None; list_exon_found = None; exon_count = None
        with open(gtf, "r") as filin:
            for line in filin:
                lines = line.split()
                #initialise transcript and get gene information
                if lines[2] == "transcript":
                    if gene:
                        list_exon_found = sorted(list_exon_found, key=get_ex_nb)
                        list_exon_found = check_intro(list_exon_found, dico_bed, gene)
                        list_exon_nb_found = list(set([get_ex_nb(exon) for exon in list_exon_found if not ("exo" in exon)]))
                        list_annot_keep = [ex for ex in list_exon_found if not ex.isdigit() is True]
                        list_exon_nb_ref = [int(ex[-1]) for ex in dico_bed[gene]]
                        
                        list_del_exon = sorted([i for i in list_exon_nb_ref if i not in list_exon_nb_found])
                        list_annot_keep.extend("Δ" + str(ex) for ex in list_del_exon)
                        list_annot_sorted = sorted(list_annot_keep, key = get_ex_nb)

                        if sorted(list_exon_found, key=get_ex_nb) == list(map(str, list_exon_nb_ref)):
                            list_annot_final = ["-".join(str(x) for x in list_exon_nb_ref)]
                        else:
                            list_annot_final = check_next(list_annot_sorted)
                        list_annot_final = "-".join(list_annot_final)

                        if (str(list_annot_final) in list(df["annot_find"])) and (gene == df.loc[list(df["annot_find"]).index(list_annot_final)]["construction"]):
                            if filename in list(df.columns):
                                df.loc[(df["annot_find"] == list_annot_final) & (df["construction"] == gene), "transcript_id"] += (" & " + info[4]) if info[4] not in (sum([ele.split(" & ") for ele in list(df["transcript_id"])], [])) else ""
                                if np.isnan(df.loc[(df["annot_find"] == list_annot_final) & (df["construction"] == gene), filename]).bool():
                                    df.loc[(df["annot_find"] == list_annot_final) & (df["construction"] == gene), filename] = round(float(cov),1)
                                else:
                                    df.loc[(df["annot_find"] == list_annot_final) & (df["construction"] == gene), filename] += round(float(cov),1)
                            else:
                                df.loc[(df["annot_find"] == list_annot_final) & (df["construction"] == gene), filename] = round(float(cov),1)
                        
                        elif info[4] not in list(df["transcript_id"]):
                            df_tmp = pd.DataFrame.from_dict({"transcript_id": [info[4]], "construction": [info[0]], "gene": [gene.split("_")[1]], "annot_ref": ["-".join(str(x) for x in list_exon_nb_ref)], "annot_find": [list_annot_final], filename: [round(float(cov),1)]})
                            df = pd.concat([df, df_tmp], ignore_index=True)
                        elif info[4] in list(df["transcript_id"]):
                            df.loc[list(df["transcript_id"]).index(info[4]), filename] = round(float(cov),1)

                    dic_attr = {attr.split()[0]: attr.split()[1].replace('"',"") for attr in line.split("\t")[-1].split(";") if attr.split()}
                    if lines[0] in dico_bed.keys():
                        gene = lines[0]; transcript_ID = dic_attr["transcript_id"]; cov = dic_attr["cov"]; exon_count_ref = dico_bed[gene][-1][-1]#exon_count_ref = len(dico_bed[gene]); cov = dic_attr["cov"]
                        info = [lines[0], lines[3], lines[4], lines[6], transcript_ID]
                       
                       
                        list_exon_found = []; exon_count = 0
                    else:
                        gene = None; continue
                

                if lines[2] == "exon":
                    dic_attr = {attr.split()[0]: attr.split()[1].replace('"',"") for attr in line.split("\t")[-1].split(";") if attr.split()}
                    if dic_attr["transcript_id"] == transcript_ID:
                        exon_count += 1
                        start_exon_stg = int("+"+lines[3]); end_exon_stg = int("+"+lines[4])
                        count = 0
                        for exon_ref in dico_bed[gene]:
                            count += 1
                            start_exon_ref = int("+"+exon_ref[0]); end_exon_ref = int("+"+exon_ref[1]); exon_nb = int(exon_ref[2]); sign = "+"
                            #retrieve exon numbering if in exon
                            if start_exon_stg in range(start_exon_ref, end_exon_ref+1):
                                #case exon == exon
                                if (start_exon_stg == start_exon_ref) and (end_exon_stg == end_exon_ref):
                                    list_exon_found.append(f"{exon_nb}"); break
                                #case del start exon
                                if (start_exon_stg > start_exon_ref) and (end_exon_stg == end_exon_ref):
                                    list_exon_found.append(f"Δ{exon_nb}p({start_exon_stg-start_exon_ref})"); break
                                #case del end exon
                                if (start_exon_stg == start_exon_ref) and (end_exon_stg < end_exon_ref):
                                    list_exon_found.append(f"Δ{exon_nb}q({end_exon_ref-end_exon_stg})"); break
                                #case del start and ins end exon
                                if (start_exon_stg > start_exon_ref) and (end_exon_stg > end_exon_ref):
                                    list_exon_found.append(f"Δ{exon_nb}p({start_exon_stg-start_exon_ref}),▼{exon_nb}q({end_exon_stg-end_exon_ref})"); break
                                #case del start and end exon
                                if (start_exon_stg > start_exon_ref) and (end_exon_stg < end_exon_ref):
                                    list_exon_found.append(f"Δ{exon_nb}p({start_exon_stg-start_exon_ref}),Δ{exon_nb}q({end_exon_ref-end_exon_stg})"); break
                                #case ins end exon or full exonisation of intron
                                if (start_exon_stg == start_exon_ref) and (end_exon_stg > end_exon_ref):
                                    #case full exonisation of intron
                                    if exon_nb != exon_count_ref:
                                        for val in dico_bed[gene]:
                                            val_sign = [int(sign+v) for v in val]
                                            if end_exon_stg in val_sign:
                                                for ins in list(range(exon_nb, int(val[2]))):
                                                    list_exon_found.append(f"▼{ins}")
                                                list_exon_found.append(f"{val[2]}"); break
                                    #case ins end exon 
                                        else:
                                            list_exon_found.append(f"▼{exon_nb}q({end_exon_stg-end_exon_ref})"); break
                                        break
                                    elif exon_nb == exon_count_ref:
                                        list_exon_found.append(f"▼{exon_nb}({end_exon_stg-end_exon_ref})"); break
                            elif end_exon_stg in range(start_exon_ref, end_exon_ref+1):
                                #case ins start exon
                                if (start_exon_stg < start_exon_ref) and (end_exon_stg == end_exon_ref):
                                    list_exon_found.append(f"▼{exon_nb-1}p({start_exon_ref-start_exon_stg})"); break
                                #case ins start and del end exon
                                if (start_exon_stg < start_exon_ref) and (end_exon_stg < end_exon_ref):
                                    list_exon_found.append(f"▼{exon_nb-1}p({start_exon_ref-start_exon_stg}),Δ{exon_nb}q({end_exon_ref-end_exon_stg})"); break
                            #retrieve intron numbering if in intron
                            elif (len(dico_bed[gene])>1) and (start_exon_stg in range(int(dico_bed[gene][0][1]), int(dico_bed[gene][-1][1])+1)):
                            #elif (str(exon_nb) != str(exon_count_ref)) and (str(exon_nb) != str(dico_bed[gene][0][-1])) and (len(dico_bed[gene])>1):
                                if lines[6] == "+":
                                    ix = [dico_bed[gene].index(ex) for ex in dico_bed[gene] if str(exon_nb) == ex[2]][0]
                                    start_intron_ref = end_exon_ref+1; end_intron_ref = int(dico_bed[gene][ix+1][0])-1; intron_nb = int(exon_ref[2])
                                elif lines[6] == "-":
                                    ix = [dico_bed[gene].index(ex) for ex in dico_bed[gene] if str(exon_nb) == ex[2]][0]
                                    start_intron_ref = end_exon_ref+1; end_intron_ref = int("-"+dico_bed[gene][ix+1][1])-1; intron_nb = int(exon_ref[2])
                                #case exonisation part of intron
                                if start_exon_stg in range(start_intron_ref, end_intron_ref+1):
                                    if start_exon_stg > start_intron_ref and end_exon_stg < end_intron_ref:
                                        list_exon_found.append(f"{intron_nb}exo({end_exon_stg-start_exon_stg+1})[p{start_exon_stg-start_intron_ref},q{end_intron_ref-end_exon_stg}]"); break
                                #case ins start and end exon 
                                if end_exon_stg in range(start_intron_ref, end_intron_ref+1):
                                    if (start_exon_stg < start_exon_ref) and (end_exon_stg > end_exon_ref):
                                        if exon_nb == exon_count_ref:
                                            list_exon_found.append(f"▼{intron_nb-1}p({start_exon_ref-start_exon_stg}),▼{exon_nb}({end_exon_stg-end_exon_ref})"); break
                                        elif exon_nb == 1:
                                            list_exon_found.append(f"▼{exon_nb}({start_exon_ref-start_exon_stg}),▼{exon_nb}q({end_exon_stg-end_exon_ref})"); break
                                        elif exon_nb != 1:
                                            list_exon_found.append(f"▼{intron_nb-1}p({start_exon_ref-start_exon_stg}),▼{exon_nb}q({end_exon_stg-end_exon_ref})"); break
                            #exonisation of intron before first or after last exon
                            elif(str(exon_nb) == str(dico_bed[gene][0][-1])) or (str(exon_nb) == str(dico_bed[gene][-1][-1])):
                                if [str(start_exon_stg), str(end_exon_stg)] not in dico_bed_UTR[gene]:
                                    #if (start_exon_stg < int(dico_bed[gene][0][0])) and (end_exon_stg < int(dico_bed[gene][0][1])):
                                    if (start_exon_stg > int(dico_bed_UTR[gene][0][1])) and (end_exon_stg < int(dico_bed[gene][0][0])):
                                        list_exon_found.append(f"{str(int(dico_bed[gene][0][-1])-1)}exo({end_exon_stg-start_exon_stg+1})"); break
                                    elif (start_exon_stg > int(dico_bed[gene][-1][0])) and (end_exon_stg < int(dico_bed_UTR[gene][-1][1])):
                                        list_exon_found.append(f"{dico_bed[gene][-1][-1]}exo({end_exon_stg-start_exon_stg+1})"); break
                            else:
                                continue
                    else:
                        continue

            if gene:
                list_exon_found = sorted(list_exon_found, key=get_ex_nb)
                list_exon_found = check_intro(list_exon_found, dico_bed, gene)
                list_exon_nb_found = list(set([get_ex_nb(exon) for exon in list_exon_found if not ("exo" in exon)]))
                list_annot_keep = [ex for ex in list_exon_found if not ex.isdigit() is True]
                list_exon_nb_ref = [int(ex[-1]) for ex in dico_bed[gene]]
                
                list_del_exon = sorted([i for i in list_exon_nb_ref if i not in list_exon_nb_found])
                list_annot_keep.extend("Δ" + str(ex) for ex in list_del_exon)
                list_annot_sorted = sorted(list_annot_keep, key = get_ex_nb)
                if sorted(list_exon_found, key=get_ex_nb) == list(map(str, list_exon_nb_ref)):
                    list_annot_final = ["-".join(str(x) for x in list_exon_nb_ref)]
                else:
                    list_annot_final = check_next(list_annot_sorted)
                list_annot_final = "-".join(list_annot_final)
                
                if (str(list_annot_final) in list(df["annot_find"])) and (gene == df.loc[list(df["annot_find"]).index(list_annot_final)]["construction"]):
                    if filename in list(df.columns):
                        df.loc[(df["annot_find"] == list_annot_final) & (df["construction"] == gene), "transcript_id"] += (" & " + info[4]) if info[4] not in (sum([ele.split(" & ") for ele in list(df["transcript_id"])], [])) else ""
                        if np.isnan(df.loc[(df["annot_find"] == list_annot_final) & (df["construction"] == gene), filename]).bool():
                            df.loc[(df["annot_find"] == list_annot_final) & (df["construction"] == gene), filename] = round(float(cov),1)
                        else:
                            df.loc[(df["annot_find"] == list_annot_final) & (df["construction"] == gene), filename] += round(float(cov),1)
                    else:
                        df.loc[(df["annot_find"] == list_annot_final) & (df["construction"] == gene), filename] = round(float(cov),1)
                elif info[4] not in list(df["transcript_id"]):
                    df_tmp = pd.DataFrame.from_dict({"transcript_id": [info[4]], "construction": [info[0]], "gene": [gene.split("_")[1]], "annot_ref": ["-".join(str(x) for x in list_exon_nb_ref)], "annot_find": [list_annot_final], filename: [round(float(cov),1)]})
                    df = pd.concat([df, df_tmp], ignore_index=True)
                elif info[4] in list(df["transcript_id"]):
                    df.loc[list(df["transcript_id"]).index(info[4]), filename] = round(float(cov),1)

    df["occurence"] = (df.iloc[:,5:].fillna(0).astype(bool).sum(axis=1)).tolist()


    for sample in samples:
        dico_calc = {}; sum_calc = 0
        for index, row in df.iterrows():
            if not np.isnan(row[sample]):
                dico_calc[row["annot_find"]] = str(row[sample])+"_"+str(index)
                sum_calc += row[sample]
        for ann in dico_calc.keys():
            df.loc[int(dico_calc[ann].split("_")[-1]), "P_"+sample] = round((float(dico_calc[ann].split("_")[0])/sum_calc)*100,2)

    with pd.ExcelWriter(outdir + "/MAGIC_SOSTAR_annotation_table_results.xlsx") as writer:
        df.to_excel(writer, index=False)



def main():
    start_time = time.time()
    
    # get arguments
    args = get_arguments()
    indir = os.path.realpath(args.input)
    outdir = os.path.realpath(args.output)
    ref_gtf = os.path.realpath(args.gtf_coord)

    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # annotate isoforms
    annotate_isoforms(indir=indir, ref_file=ref_gtf, outdir=outdir)

    step = f"--- Program {os.path.basename(__file__)} executed in {round((float(time.time() - start_time) / 60), 2)} minutes ---"; print(step)



if __name__ == "__main__":
    main()