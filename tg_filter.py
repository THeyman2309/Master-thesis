# create hashtable to store proband data.
proband_table = dict()

# open output genome-wide analysis patient 1 by tandem-genotypes

#with open("./proband/proband_tg_rmsk_filt.txt", "r") as fr:
with open("straglr_out_1dif.txt", "r") as fr:
    for line in fr:
        # ignore header
        if not line.startswith("#"):
            # split into relevant units
            chrom, start, end, unit, gene, pos, lenght_1, lenght_2, r_1, r_2 = line.split()
            # convert some units to integers to allow computations
            start, end, lenght_1, lenght_2 = int(start), int(end), int(lenght_1), int(lenght_2)
            # exclude intergenic regions
            if pos != "intergenic":
                # if this TR hasn't occured on this chromosome, make new entry
                if f"{chrom}: {unit}" not in proband_table:
                    proband_table[f"{chrom}: {unit}"] = [{"position": [start, end, pos, gene], "proband": [lenght_1, lenght_2]}]
                # if this TR has occured on this chromosome, add new position to existing entry
                else:
                    # store information per chromosome-motif combination
                    proband_table[f"{chrom}: {unit}"].append({"position": [start, end, pos, gene], "proband": [lenght_1, lenght_2]})

fr.close()


# open output genome-wide analysis father of patient 1 by tandem-genotypes

#with open("./father/father_tg_rmsk_filt.txt", "r") as fr:
with open("straglr_out.txt", "r") as fr:
    for line in fr:
        # ignore header
        if not line.startswith("#"):
            # split into relevant units
            chrom, start, end, unit, gene, pos, lenght_1, lenght_2, r_1, r_2 = line.split()
            # convert some units to integers to allow computations
            start, end, lenght_1, lenght_2 = int(start), int(end), int(lenght_1), int(lenght_2)
            # this repeat is only relevant if it was also detected in the proband
            if f"{chrom}: {unit}" in proband_table:
                # identify positon in the proband
                for index, occurence in enumerate(proband_table[f"{chrom}: {unit}"]):
                    # give 100 bp room for erronous rapportation of begin and end positions
                    if abs(occurence["position"][0] - start) < 100 or abs(occurence["position"][1] - end) < 100:
                        # add a key (father) to the dictionary associated with this specific TR, to store info
                        proband_table[f"{chrom}: {unit}"][index]["father"] = [lenght_1, lenght_2]


# perform exactly the same process for the output of the genome-wide analysis of the mother of patient 1.

#with open("./mother/mother_tg_rmsk_filt.txt", "r") as fr:
with open("straglr_out.txt", "r") as fr:
    for line in fr:
        if not line.startswith("#"):
            # split into relevant units
            chrom, start, end, unit, gene, pos, lenght_1, lenght_2, r_1, r_2 = line.split()
            start, end, lenght_1, lenght_2 = int(start), int(end), int(lenght_1), int(lenght_2)

            # this repeat is only relevant if it was also detected in the proband
            if f"{chrom}: {unit}" in proband_table:
                # identify positon in the proband
                for index, occurence in enumerate(proband_table[f"{chrom}: {unit}"]):
                    # give 100 bp room for erronous rapportation
                    if abs(occurence["position"][0] - start) < 100 or abs(occurence["position"][1] - end) < 100:
                        # add a key (mother) to the dictionary associated with this specific TR, to store info
                        proband_table[f"{chrom}: {unit}"][index]["mother"] = [lenght_1, lenght_2]

# After storing the information on all TRs that were detected in the patient, the filtering can start.
# Make list to store TRs that meet the filter criteria
TRs_of_interest = []

# For every stored chromosome-motif combination
for key in proband_table.keys():
    # For every TR with this chromosome-motif combination
    for entry in proband_table[key]:
        # Retrieve chromosome and unit
        chrom, unit = key.split(": ")
        # Only include TRs where the largest expansion in the patient compared to the reference is -
        # at least 20 repeated units.
        if max(entry["proband"]) > 20:
            # Include all TRs that only occured in proband (entry only contains the keys "position" and "proband").
            if len(entry) == 2:
                TRs_of_interest.append([chrom, entry["position"][0], entry["position"][1], unit,
                                        entry["position"][2], entry["position"][3],
                                        entry["proband"][0], entry["proband"][1],
                                        "-", "-",
                                        "-", "-",
                                        "de novo"])

            else:
                # If the TR was detected in the patient (always the case) and both parents.
                # Determine the maximal and minum TR length detected in each individual
                if "father" in entry and "mother" in entry:
                    min_fath = min(entry["father"])
                    min_mother = min(entry["mother"])
                    max_parent = max(entry["father"] + entry["mother"])
                    min_proband = min(entry["proband"])
                    max_proband = max(entry["proband"])

                    # If the smallest allele detected in the patient is bigger than the smallest allele of both parents.
                    # This section allows for the detection of TRs that could have a recessive inheritance pattern.
                    if min_fath + 20 < min_proband and min_mother + 20 < min_proband:
                        TRs_of_interest.append([chrom, entry["position"][0], entry["position"][1], unit,
                                                entry["position"][2], entry["position"][3],
                                                entry["proband"][0], entry["proband"][1],
                                                entry["mother"][0], entry["mother"][1],
                                                entry["father"][0], entry["father"][1],
                                                "bigger than smallest allele of both parents"])
                    # If the largest allele is at least 20 expansions bigger than the largest parental allele,
                    # the TR is preserved.
                    elif max_proband - 20 > max_parent:
                        TRs_of_interest.append([chrom, entry["position"][0], entry["position"][1], unit,
                                                entry["position"][2], entry["position"][3],
                                                entry["proband"][0], entry["proband"][1],
                                                entry["mother"][0], entry["mother"][1],
                                                entry["father"][0], entry["father"][1],
                                                "bigger than largest parental allele"])

                # If the tandem-repeat was not detected in both parents
                else:
                    # If one parent is missing, the third key in the entry of the TR is the parent where the TR
                    # was detected. We store the name of the parent where the TR was detected.
                    parent = [i for i in entry.keys()][2]
                    max_parent = max(entry[parent])
                    min_parent = min(entry[parent])
                    min_proband = min(entry["proband"])
                    max_proband = max(entry["proband"])

                    # If the smallest allele of the patient is bigger than the smallest allele of the parent
                    if min_parent + 20 < min_proband:
                        # If the TR was detected in the mother
                        if parent == "mother":
                            TRs_of_interest.append([chrom, entry["position"][0], entry["position"][1], unit,
                                                    entry["position"][2], entry["position"][3],
                                                    entry["proband"][0], entry["proband"][1],
                                                    entry["mother"][0], entry["mother"][1],
                                                    "-", "-",
                                                    "bigger than smallest maternal allele"])
                        # If the TR was detected in the father
                        else:
                            TRs_of_interest.append([chrom, entry["position"][0], entry["position"][1], unit,
                                                     entry["position"][2], entry["position"][3],
                                                     entry["proband"][0], entry["proband"][1],
                                                     "-", "-",
                                                     entry["father"][0], entry["father"][1],
                                                     "bigger than smallest parental allele"])

                    # If the TR was detected in only one parent and the smallest allele of the patient was not larger
                    # than the smallest allele of the parent, but the largest allele of the patient was larger than
                    # the largest allele of the parent.
                    elif max_proband - 20 > max_parent:
                        if parent == "mother":
                            TRs_of_interest.append([chrom,  entry["position"][0], entry["position"][1], unit,
                                                    entry["position"][2], entry["position"][3],
                                                    entry["proband"][0], entry["proband"][1],
                                                    entry["mother"][0], entry["mother"][1],
                                                    "-", "-",
                                                    "bigger than largest maternal allele"])

                        # If the TR was detected in the father
                        else:
                            TRs_of_interest.append([chrom, entry["position"][0], entry["position"][1], unit,
                                                    entry["position"][2], entry["position"][3],
                                                    entry["proband"][0], entry["proband"][1],
                                                    "-", "-",
                                                    entry["father"][0], entry["father"][1],
                                                    "bigger than largest paternal allele"])


# Write outputfile

with open("filter_output.txt", "w") as fw:
    # Write header for the output file.
    fw.write("#CHROM\tSTART\tEND\tMOTIF\tLOCATION\tGENE\tPROB_1\tPROB_2\tMOTH_1\tMOTH_2\tFATH_1\tFATH_2\tCOMMENT")
    for TR in TRs_of_interest:
        fw.write("\n")
        # Join only works on a list of strings.
        fw.write('\t'.join([str(i) for i in TR]))
fw.close()
