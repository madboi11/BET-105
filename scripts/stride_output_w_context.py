import sys
import os


HELIX_TYPES = {"AlphaHelix", "310Helix", "PiHelix"}


def parse_stride(stride_file):
    residues = []
    with open(stride_file) as f:
        for line in f:
            if not line.startswith("ASG"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            resname = parts[1]
            chain = parts[2]
            try:
                resnum = int(parts[3])
            except ValueError:
                continue
            ss = parts[6]
            residues.append((chain, resnum, resname, ss))
    return residues


def main():
    stride_file = sys.argv[1]
    output_file = sys.argv[2]
    target = sys.argv[3]

    if os.path.getsize(stride_file) == 0:
        open(output_file, "w").close()
        return

    residues = parse_stride(stride_file)

    with open(output_file, "w") as out:
        out.write("chain\tprev_resnum\tprev_resname\tprev_ss\tcenter_resnum\tcenter_resname\tcenter_ss\tnext_resnum\tnext_resname\tnext_ss\n")
        for i in range(1, len(residues) - 1):
            prev = residues[i - 1]
            center = residues[i]
            nxt = residues[i + 1]

            if center[2] != target:
                continue
            if prev[0] != center[0] or center[0] != nxt[0]:
                continue
            if center[3] not in HELIX_TYPES:
                continue
            if prev[3] not in HELIX_TYPES:
                continue
            if nxt[3] not in HELIX_TYPES:
                continue

            out.write(
                f"{center[0]}\t{prev[1]}\t{prev[2]}\t{prev[3]}\t"
                f"{center[1]}\t{center[2]}\t{center[3]}\t"
                f"{nxt[1]}\t{nxt[2]}\t{nxt[3]}\n"
            )


if __name__ == "__main__":
    main()
