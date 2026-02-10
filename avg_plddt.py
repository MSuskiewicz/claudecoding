import sys, os, glob
from openpyxl import Workbook

AA = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q',
      'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
      'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

folder = sys.argv[1]
wb = Workbook()
ws = wb.active
ws.append(["PDB", "Residues", "pLDDT_1", "pLDDT_2", "pLDDT_3", "pLDDT_4", "Average"])

for pdb in sorted(glob.glob(os.path.join(folder, "*.pdb"))):
    # Extract CA atoms from chain A: list of (one-letter code, pLDDT)
    res = [(AA.get(l[17:20].strip(), "X"), float(l[60:66]))
           for l in open(pdb) if l[:4] == "ATOM" and l[21] == "A" and l[12:16].strip() == "CA"]
    # Sliding window of 4 with highest average pLDDT
    i = max(range(len(res) - 3), key=lambda i: sum(r[1] for r in res[i:i+4]))
    win = res[i:i+4]
    seq = "".join(r[0] for r in win)
    scores = [r[1] for r in win]
    ws.append([os.path.basename(pdb), seq] + scores + [sum(scores) / 4])

wb.save(os.path.join(folder, "plddt_results.xlsx"))
print(f"Saved {os.path.join(folder, 'plddt_results.xlsx')}")
