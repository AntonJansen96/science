from typing import List


def genMissingAtomRep(
    fname: str,
    resRemark=465,
    atomRemark=470,
    chains: List[str] = ["A", "B", "C", "D", "E"],
    residues: List[str] = ["all"],
) -> None:
    """Generates .vmd representation file for visualzing the missing atoms in a .pdb file."""

    def list2str(array: list) -> str:
        xstr = ""
        for val in array:
            xstr += f"{val} "
        return xstr[:-1]

    out = open(f"{fname[:-4]}_missing.vmd", "w+")
    count = 1

    for line in open(fname).read().splitlines():
        if f"REMARK {resRemark} " in line and "SSSEQI" not in line:
            for chain in chains:
                if f" {chain} " in line:

                    line = line.split()
                    if not (line[2] in residues or residues[0] == "all"):
                        break

                    out.write(
                        f'atomselect macro l{count} "chain {line[3]} and resid {int(line[4])}"\n'
                    )
                    count += 1

                    break

        elif f"REMARK {atomRemark} " in line:
            for chain in chains:
                if f" {chain} " in line:

                    line = line.split()
                    if not (line[2] in residues or residues[0] == "all"):
                        break

                    out.write(
                        f'atomselect macro l{count} "chain {line[3]} and resid {int(line[4])} and name '
                    )
                    for idx in range(5, len(line)):
                        out.write(f"{line[idx]} ")
                    out.write('"\n')
                    count += 1

                    break

    xstr = ""
    for x in range(1, count):
        xstr += f"l{x} or "
    xstr = xstr[: len(xstr) - 4]
    out.write(f'atomselect macro missing "{xstr}"\n')

    rawstring = f"""
mol modselect     0 0 chain {list2str(chains)}
mol modstyle      0 0 NewCartoon 0.3 20 4.1 0
mol modmaterial   0 0 AOChalky
mol modcolor      0 0 ColorID 10

mol addrep 0
mol modselect     1 0 (chain {list2str(chains)}) and missing
mol modstyle      1 0 NewCartoon 0.3 20 4.1 0
mol modmaterial   1 0 AOChalky
mol modcolor      1 0 ColorID 1

mol addrep 0
mol modselect     2 0 (chain {list2str(chains)}) and missing
mol modstyle      2 0 Licorice 0.3 20 4.1 0
mol modmaterial   2 0 AOChalky
mol modcolor      2 0 ColorID 1
"""

    out.write(rawstring)
    out.close()
