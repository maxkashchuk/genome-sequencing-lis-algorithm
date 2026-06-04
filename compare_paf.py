import sys

def parse_paf(filename):
    alignments = {}
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 11:
                continue

            strand = parts[4]
            if strand == '-':
                continue
            
            q_name = parts[0]
            q_start, q_end = int(parts[2]), int(parts[3])
            t_name = parts[5]
            t_start, t_end = int(parts[7]), int(parts[8])
            
            length = q_end - q_start
            if q_name not in alignments or length > (alignments[q_name]['q_end'] - alignments[q_name]['q_start']):
                alignments[q_name] = {
                    't_name': t_name,
                    'q_start': q_start, 'q_end': q_end,
                    't_start': t_start, 't_end': t_end
                }
    return alignments

def main():
    my_paf = "./paf_results_new/alignments.paf"
    mini_paf = "./paf_results_new/minimap2_align.paf"

    print(f"Loading {my_paf}...")
    my_aln = parse_paf(my_paf)
    
    print(f"Loading {mini_paf}...")
    mini_aln = parse_paf(mini_paf)

    common_reads = set(my_aln.keys()).intersection(set(mini_aln.keys()))
    
    print(f"\n--- Comparison Results ---")
    print(f"Reads aligned by our program : {len(my_aln)}")
    print(f"Reads aligned by minimap2    : {len(mini_aln)}")
    print(f"Common aligned reads         : {len(common_reads)}")

    good_matches = 0
    tolerance = 100

    for read in common_reads:
        m1 = my_aln[read]
        m2 = mini_aln[read]
        
        if (m1['t_name'] == m2['t_name'] and 
            abs(m1['t_start'] - m2['t_start']) <= tolerance and 
            abs(m1['t_end'] - m2['t_end']) <= tolerance):
            good_matches += 1

    print(f"Matches with close coordinates ({tolerance}bp tolerance): {good_matches}")
    if len(common_reads) > 0:
        print(f"Accuracy relative to minimap2: {(good_matches / len(common_reads)) * 100:.2f}%")

if __name__ == "__main__":
    main()