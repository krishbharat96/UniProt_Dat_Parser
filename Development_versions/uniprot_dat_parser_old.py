import re

prot_dict = dict()
def create_de(de_dict):
    de_string = ""
    if not de_dict is None:
        de_string = de_dict["Actual_prot_name"][0]
        if len(de_dict["Actual_prot_name"]) > 1:
            for n in range(len(de_dict["Actual_prot_name"])):
                if not (n == 0):
                    de_string = de_string + " (" + de_dict["Actual_prot_name"][n] + ")"
        if len(de_dict["Alternate_names"]) > 0:
            for n in range(len(de_dict["Alternate_names"])):
                de_string = de_string + " (" + de_dict["Alternate_names"][n] + ")"
        if len(de_dict["Contains"]) > 0:
            for n in range(len(de_dict["Contains"])):
                de_string = de_string + " [Cleaved into: " + de_dict["Contains"][n] + "]"
    return de_string

def get_sv(sv_string):
    sv_int = ""
    sv_split = sv_string.strip().split(",")
    for item in sv_split:
        if "version" in item:
            replace_arr = [".", "sequence version"]
            sv_int = item
            for r in replace_arr:
                sv_int = sv_int.replace(r, "")
    return sv_int.strip()

def elim_brackets(string):
    string = re.sub(r'\[.*\]', '', string).strip()
    string = re.sub(r'\([^()]*\)', '', string).strip()
    string = re.sub(r'\{[^()]*\}', '', string).strip()
    return string

def create_iso(isostring):
    iso_dict = dict()
    if not (isostring.strip() == ""):
        name_arr = []
        isoid_arr = []
        isoseq_arr = []
        for line in isostring.splitlines():
            l_name_arr = process_dat_line(line)
            if "Name=" in line:
                name_init = l_name_arr[1]
                name = name_init.split(";")[0]
                replace_arr = [";", "Name="]
                for r in replace_arr:
                    name = name.replace(r, "")
                name_arr.append(elim_brackets(name))
            if "IsoId=" in line:
                # Get rid of last character
                isoid = ""
                info_arr = l_name_arr[1].split(";")
                for info in info_arr:
                    if "IsoId=" in info:
                        isoid = info.replace("IsoId=", "").strip()
                        isoid_arr.append(isoid)
            if "Sequence=" in line:
                seq_type = ""
                inf_arr = l_name_arr[1].split(";")
                for inf in inf_arr:
                    if "Sequence=" in inf:
                        seq_type = inf.replace("Sequence=", "").strip()
                        isoseq_arr.append(seq_type)
        for n in range(len(name_arr)):
            iso_dict.update({name_arr[n]:{"isoid":isoid_arr[n], "isoseq":isoseq_arr[n]}})

    # print iso_dict
    return iso_dict

prot_dict = dict()
def process_de(arr):
    output_str = ""
    dict_de = dict()
    dict_de.update({"Actual_prot_name": [], "Alternate_names":[], "Contains":[]})
    marked_contains = False
    prev_rec = False
    prev_alt = False
    for item in arr:
        if (marked_contains == False):
            if "RecName:" in item:
                edited_item = item
                replace_arr = ["RecName: ", "Full=", ";"]
                for r in replace_arr:
                    edited_item = edited_item.replace(r, "")
                dict_de["Actual_prot_name"].append(edited_item.strip())
                prev_alt = False
                prev_rec = True
            if "AltName:" in item:
                edited_item = item
                replace_arr = ["AltName: ", "Full=", ";"]
                for r in replace_arr:
                    edited_item = edited_item.replace(r, "")
                dict_de["Alternate_names"].append(edited_item.strip())
                prev_alt = True
                prev_rec = False
            if "Contains:" in item:
                replace_arr = ["RecName: ", "Full=", ";", "Contains:"]
                edited_item = item
                for r in replace_arr:
                    edited_item = edited_item.replace(r, "")
                if not (edited_item.strip() == ""):
                    dict_de["Contains"].append(edited_item.strip())
                marked_contains = True
            if not ":" in item:
                if (prev_alt == True) and (prev_rec == False):
                    act_term = item.split("=")[1].strip()
                    edited_item = act_term
                    replace_arr = ["RecName: ", "Full=", ";"]
                    for r in replace_arr:
                        edited_item = edited_item.replace(r, "")
                    dict_de["Alternate_names"].append(edited_item)
                else:
                    act_term = item.split("=")[1].strip()
                    edited_item = act_term
                    replace_arr = ["RecName: ", "Full=", ";"]
                    for r in replace_arr:
                        edited_item = edited_item.replace(r, "")
                    dict_de["Actual_prot_name"].append(edited_item)
        else:
            try:
                c_item1 = item.split(":")[1].strip()
                c_item2 = c_item1.split("=")[1].strip()
                c_item3 = c_item2.replace(";", "")
                dict_de["Contains"].append(c_item3.strip())
            except:
                continue
            marked_contains = False
    return create_de(dict_de)
        
def process_dat_line(lin):
    init_arr = lin.split("  ")
    new_arr = []
    for item in init_arr:
        if not (item.strip() == ""):
            new_arr.append(item.strip())
    return new_arr

def process_protein(prot_string):
    id_prot = ""
    de_arr = []
    on_isoform = False
    isostring = ""
    sv_int = ""
    gene_arr = []
    pe_val = ""
    prot_name = ""
    db_name = ""
    chembl_id = ""
    os = ""
    for line in prot_string.splitlines():
        t_split = process_dat_line(line)
        if (t_split[0] == "ID") and (len(t_split) > 1):
            prot_name = t_split[1]
            if "Reviewed" in t_split[2]:
                db_name = "SwissProt"
            else:
                db_name = "TrEMBL"
        if (t_split[0] == "AC") and (len(t_split) > 1):
            semi_split = t_split[1].split(";")
            prot_id = semi_split[0].strip()
            if (id_prot == ""):
                id_prot = prot_id
        if (t_split[0] == "DE") and (len(t_split) > 1):
            de_arr.append(t_split[1])
        if (t_split[0] == "CC") and ("-!- ALTERNATIVE PRODUCTS:" in t_split[1]):
            on_isoform = True
            isostring = line
            continue
        if (on_isoform == True):
            if not "-!-" in line:
                isostring = isostring + "\n" + line
            else:
                on_isoform = False
        if (t_split[0] == "DT") and (len(t_split) > 1):
            if "sequence version" in t_split[1]:
                sv_int = get_sv(t_split[1])
        if (t_split[0] == "GN") and (len(t_split) > 1):
            gene_str = t_split[1].split(";")
            for item in gene_str:
                if "Name" in item:
                    gene_act = elim_brackets(item).replace("Name=", "").strip()
                    gene_arr.append(gene_act)
        if (t_split[0] == "PE") and (len(t_split) > 1):
            pe_split = t_split[1].split(":")
            pe_val = pe_split[0].strip()
        if (t_split[0] == "DR") and (len(t_split) > 1):
            if "CHEMBL" in t_split[1] and "ChEMBL" in t_split[1]:
                c_split = t_split[1].split(";")
                for c_item in c_split:
                    if "CHEMBL" in c_item:
                        chembl_id = c_item.strip()
        if (t_split[0] == "OS"):
            os = t_split[1].replace(".", "").strip()
    prot_dict.update({id_prot:{"Name":prot_name, "Description":process_de(de_arr), "Isoform":create_iso(isostring), "SV":sv_int, "PE":pe_val, "Gene":gene_arr, "Database":db_name, "Chembl_id":chembl_id, "Organism": os}})

def parse(filename):
    file = open(filename, "r")
    str1 = ""
    for line in file:
        if line:
            t_split = process_dat_line(line)
            if (t_split[0] == "ID"):
                if not (str1.strip() == ""):
                    process_protein(str1)
                str1 = line
                # print t_split[1], t_split[2]
            else:
                str1 = str1 + line
    process_protein(str1)
    # print prot_dict
    return prot_dict
    

