n_homozygous, sum_hom, n_heterozygous, sum_het = 0, 0, 0, 0
with open("vals.txt", "r") as file:
    for line in file:
        line = line.split()
        for person in range(3):
            if "/" in line[person]:
                n_heterozygous += 2
                strag = [float(x) for x in line[person].split("/")]
                tandem_gen = [float(x) for x in line[person + 3].split("/")]
                sum_het += abs(strag[0] - tandem_gen[0]) + abs(strag[1]  - tandem_gen[1])
            elif line[person] != "n.d." and line[person] != "n.d":
                n_homozygous += 2
                strag = float(line[person])
                tandem_gen = [int(x) for x in line[person + 3].split("/")]
                sum_hom += abs(strag - tandem_gen[0]) + abs(strag  - tandem_gen[1])

    print("On average:",sum_hom/n_homozygous, f"from {n_homozygous} observations")
    print("On average:",sum_het/n_heterozygous, f"from {n_heterozygous} observations")
    print("Combined on average:", (sum_het + sum_hom) / (n_heterozygous + n_homozygous), f"from {n_heterozygous + n_homozygous} observations")
