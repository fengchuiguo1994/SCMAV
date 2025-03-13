import sys
import json

with open("{0}/summary.json".format(sys.argv[1]), 'r') as f:
    data = json.load(f)

print("ATAC_Estimated_cells\t{0}".format(format(data["cells_detected"], ',')))
print("ATAC_doublet_rate\t{0:.4f}".format(data["estimated_gelbead_doublet_rate"]))
print("ATAC_fragment_in_peaks_rate\t{0:.4f}".format(data["frac_cut_fragments_in_peaks"]))
print("ATAC_fragment_overlapping_peaks_rate\t{0:.4f}".format(data["frac_fragments_overlapping_peaks"]))
print("ATAC_fragment_overlapping_TSS_rate\t{0:.4f}".format(data["frac_fragments_overlapping_tss"]))
print("ATAC_median_fragment_count_percell\t{0}".format(data["median_fragments_per_cell"]))
