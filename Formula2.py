n_homozygous, sum_hom, n_heterozygous, sum_het, sum_dif, n_different, tg = 0, 0, 0, 0, 0, 0, 0
with open("hg_t2t_calls.txt", "r") as file:
    # for every tandem-repeat
    for line in file:
        if not line.startswith("#"):
            line = line.split()
            for person in range(6):
                # if both tools called it as heterozygous
                if "/" in line[person] and "/" in line[person + 6]:
                    n_heterozygous += 2
                    hg38 = [float(x) for x in line[person].split("/")]
                    T2T = [float(x) for x in line[person + 6].split("/")]
                    sum_het += abs(hg38[0] - T2T[0]) + abs(hg38[1]  - T2T[1])
                    if person >= 3:
                        tg += 1

                # if both called it as homozygous
                elif line[person] != "n.d." and line[person + 6] != "n.d" and "/" not in line[person] and "/" not in line[person + 6]:
                    n_homozygous += 2
                    hg38 = float(line[person])
                    T2T = float(line[person + 6])
                    sum_hom += 2 * abs(hg38 - T2T)

                # if only T2T called it as heterozygous
                elif line[person] != "n.d." and "/" in line[person + 6]:
                    n_different += 2
                    hg38 = float(line[person])
                    T2T = [float(x) for x in line[person + 6].split("/")]
                    sum_dif += abs(hg38 - T2T[0]) + abs(hg38 - T2T[1])

                # if only hg38 called it as heterozygous
                elif line[person + 6] != "n.d" and "/" in line[person]:
                    n_different += 2
                    hg38 = [int(x) for x in line[person].split("/")]
                    T2T = float(line[person + 6])
                    sum_dif += abs(hg38[0] - T2T) + abs(hg38[1] - T2T)




    print("On average:",sum_hom/n_homozygous, f"from {n_homozygous} observations")
    print("On average:",sum_het/n_heterozygous, f"from {n_heterozygous} observations")
    print("contribution tandem-genotypes : ", tg)
    print("On average when they differed:", sum_dif / n_different, f"from {n_different} observations")
    print("Combined on average:", (sum_het + sum_hom +sum_dif) / (n_heterozygous + n_homozygous + n_different), f"from {n_heterozygous + n_homozygous + n_different} observations")
