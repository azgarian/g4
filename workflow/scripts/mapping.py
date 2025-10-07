import os
from helper_func import make_windows, map_xr_ds
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, required=True)
parser.add_argument("-r", "--region_name", type=str, required=True)
parser.add_argument("-t", "--tss", type=bool, default=False)
parser.add_argument("--g4-bed", dest="g4_bed", type=str, default=None, help="Override G4 BED path for windowing")
args = parser.parse_args()

file_path = args.input
region_name = args.region_name
tss = args.tss

if os.path.exists(f"results/master/{region_name}/mapped_all.csv"):
    print(f"results/master/{region_name}/mapped_all.csv already exists!", flush=True)
else:
    print(f"Processing {region_name}...", flush=True)
    info = region_name.split("_")
    if tss:
        if len(info) >= 3:
            region = make_windows(
                "resources/ref_genomes/hg38/hg38_proteinCoding_genes.bed",
                quant=False, quant_num=4,
                region_size=(int(info[-2]) - 1) * int(info[-1]),
                interval_length=int(info[-1]), tss=True,
            )
            map_xr_ds(file_path, region, "results/master", region_name, int(info[-1]))
    else:
        if len(info) >= 4:
            bed_path = args.g4_bed or f"results/g4_miner/oqs_{"_".join(info[1:-2])}.bed"
            region = make_windows(
                bed_path,
                quant=False, quant_num=4,
                region_size=(int(info[-2]) - 1) * int(info[-1]),
                interval_length=int(info[-1]),
            )
            map_xr_ds(file_path, region, "results/master", region_name, int(info[-1]))
