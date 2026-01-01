# Please cite the following article if you use this code:
# Gheibzadeh MS, Manyumwa CV, Tastan Bishop Ö, Shahbani Zahiri H, Parkkila S, Zolfaghari Emameh R.
# Genome Study of α-, β-, and γ-Carbonic Anhydrases from the Thermophilic Microbiome of Marine Hydrothermal Vent Ecosystems.
# Biology (Basel). 2023 May 25;12(6):770. doi: 10.3390/biology12060770.


import re, json
from collections import Counter
from dataclasses import dataclass, field, asdict
from typing import List, Tuple, Optional

# -----------------------------
# Pattern definitions
# -----------------------------
AMINO = r"[ACDEFGHIKLMNPQRSTVWY]"
GAP_MIN, GAP_MAX = 0, 300

ALPHA_PATTERN = re.compile(rf"(H{AMINO}H{AMINO}{{16}}H)")
BETA_PATTERN  = re.compile(rf"(C{AMINO}D{AMINO}R){AMINO}{{{GAP_MIN},{GAP_MAX}}}(H{AMINO}{{2}}C)")
GAMMA_PATTERN = re.compile(rf"(N{AMINO}Q{AMINO}{{5}}H){AMINO}{{{GAP_MIN},{GAP_MAX}}}(H{AMINO}{{4}}H)")

# -----------------------------
# Data structures
# -----------------------------
@dataclass
class MotifMatch:
    family: str
    span: Tuple[int,int]
    segments: List[str]
    gap_len: Optional[int] = None

@dataclass
class ClassificationResult:
    id: str
    sequence: str
    family: str
    matches: List[MotifMatch] = field(default_factory=list)

# -----------------------------
# Motif search functions
# -----------------------------
def find_alpha(seq: str) -> List[MotifMatch]:
    return [MotifMatch("Alpha", m.span(1), [m.group(1)]) for m in ALPHA_PATTERN.finditer(seq)]

def find_beta(seq: str) -> List[MotifMatch]:
    matches=[]
    for m in BETA_PATTERN.finditer(seq):
        gap_len = m.span(2)[0]-m.span(1)[1]
        matches.append(MotifMatch("Beta",(m.span(1)[0],m.span(2)[1]),[m.group(1),m.group(2)],gap_len))
    return matches

def find_gamma(seq: str) -> List[MotifMatch]:
    matches=[]
    for m in GAMMA_PATTERN.finditer(seq):
        gap_len = m.span(2)[0]-m.span(1)[1]
        matches.append(MotifMatch("Gamma",(m.span(1)[0],m.span(2)[1]),[m.group(1),m.group(2)],gap_len))
    return matches

def classify_sequence(seq_id:str, seq:str)->ClassificationResult:
    clean = re.sub(r"[^A-Z]","",seq.upper())
    alpha,beta,gamma = find_alpha(clean),find_beta(clean),find_gamma(clean)
    families=[]
    if alpha: families.append("Alpha")
    if beta:  families.append("Beta")
    if gamma: families.append("Gamma")
    family="None" if not families else families[0] if len(families)==1 else "Ambiguous"
    return ClassificationResult(seq_id,clean,family,alpha+beta+gamma)

# -----------------------------
# Reading FASTA
# -----------------------------
def read_fasta_file(filepath):
    records,seq_id,buf=[],None,[]
    with open(filepath) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if seq_id: records.append((seq_id,"".join(buf)))
                seq_id=line[1:].strip(); buf=[]
            else: buf.append(line)
    if seq_id: records.append((seq_id,"".join(buf)))
    return records

# -----------------------------
# Writing FASTA outputs
# -----------------------------
def write_fasta(filename, records, line_length=70):
    with open(filename, "w") as f:
        for sid, seq in records:
            f.write(f">{sid}\n")
            # Split sequence into fixed-length lines
            for i in range(0, len(seq), line_length):
                f.write(seq[i:i+line_length] + "\n")

# -----------------------------
# Main execution
# -----------------------------
def main():
    input_path="CA.fasta"   # Input file
    records=read_fasta_file(input_path)
    results=[classify_sequence(sid,seq) for sid,seq in records]

    # Classification for FASTA files
    alpha=[(r.id,r.sequence) for r in results if r.family=="Alpha"]
    beta=[(r.id,r.sequence) for r in results if r.family=="Beta"]
    gamma=[(r.id,r.sequence) for r in results if r.family=="Gamma"]
    none=[(r.id,r.sequence) for r in results if r.family=="None"]

    write_fasta("alpha.txt",alpha)
    write_fasta("beta.txt",beta)
    write_fasta("gamma.txt",gamma)
    write_fasta("none.txt",none)

    # Overall statistics
    counts=Counter([r.family for r in results])
    total=len(results)
    print("Statistics:")
    for fam,count in counts.items():
        print(f"{fam}: {count} ({count/total:.1%})")

    # Save statistics to text file
    with open("statistics.txt", "w") as f:
        f.write("Statistics:\n")
        for fam, count in counts.items():
            f.write(f"{fam}: {count} ({count/total:.1%})\n")

    # Save detailed results to JSON
    with open("detailed_results.json","w") as f:
        json.dump([asdict(r) for r in results], f, ensure_ascii=False, indent=2)

if __name__=="__main__":
    main()
