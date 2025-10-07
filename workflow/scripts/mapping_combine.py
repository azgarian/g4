import pandas as pd
import argparse
import yaml
import gc
from helper_func import make_windows

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--region_name", type=str, required=True)
parser.add_argument("-p", "--mapped_paths", nargs="+", type=str, required=True)
parser.add_argument("-s", "--sim_paths", nargs="+", type=str, required=True)
parser.add_argument("-c", "--config_file", type=str, required=True)
args = parser.parse_args()
region_name = args.region_name
mapped_paths = args.mapped_paths
sim_paths = args.sim_paths

with open(args.config_file, "r") as f:
    cfg = yaml.safe_load(f)
samples = cfg.get("samples", {})

print(f"Combining {region_name}...", flush=True)

info = region_name.split("_")

if "oqs" in region_name:
    base_df = make_windows(
        f"results/g4_miner/oqs_{"_".join(info[1:-2])}.bed",
        quant=False,
        quant_num=4,
        region_size=(int(info[-2]) - 1) * int(info[-1]),
        interval_length=int(info[-1]),
    )
elif "tss" in region_name:
    base_df = make_windows(
        f"resources/ref_genomes/hg38/hg38_proteinCoding_genes.bed",
        quant=False,
        quant_num=4,
        region_size=(int(info[-2]) - 1) * int(info[-1]),
        interval_length=int(info[-1]),
        tss=True,
    )
else:
    raise ValueError(f"Invalid region name: {region_name}")

base_plus = base_df.copy()
base_plus["dam_strand"] = "+"
base_minus = base_df.copy()
base_minus["dam_strand"] = "-"
base_df = pd.concat([base_plus, base_minus], ignore_index=True)
merge_keys = ["chrom", "start", "end", "name", "score", "strand", "dam_strand", "window"]
base_df = base_df.sort_values(merge_keys).reset_index(drop=True)

expected = len(samples.keys())
if len(mapped_paths) == expected and len(sim_paths) == expected:
    # Read files one by one and concatenate incrementally to save memory
    all_values = None
    
    # Process mapped files
    for file in mapped_paths:
        df = pd.read_csv(file, index_col=False, dtype='float32')  # Use float32 to save memory
        df.reset_index(drop=True, inplace=True)
        if all_values is None:
            all_values = df
        else:
            all_values = pd.concat([all_values, df], axis=1)
        del df  # Free memory immediately
        gc.collect()  # Force garbage collection
    
    # Process sim files
    for file in sim_paths:
        df = pd.read_csv(file, index_col=False, dtype='float32')  # Use float32 to save memory
        df.reset_index(drop=True, inplace=True)
        all_values = pd.concat([all_values, df], axis=1)
        del df  # Free memory immediately
        gc.collect()  # Force garbage collection

    # Concatenate with base_df only once
    mapped_all = pd.concat([base_df, all_values], axis=1)

    rename_pairs = {i: samples[i]["name"] for i in samples.keys()}
    rename_pairs.update({f"{i}_sim": f"{samples[i]['name']}_sim" for i in samples.keys()})
    mapped_all.rename(columns=rename_pairs, inplace=True)
    
    # Use chunked writing for large files to save memory
    chunk_size = 10000
    if len(mapped_all) > chunk_size:
        mapped_all.to_csv(f"results/master/{region_name}/mapped_all.csv", index=False, chunksize=chunk_size)
    else:
        mapped_all.to_csv(f"results/master/{region_name}/mapped_all.csv", index=False)
    
    # Clean up memory
    del mapped_all, all_values
    gc.collect()
    print(f"results/master/{region_name}/mapped_all.csv is created!", flush=True)
else:
    print(f"not all files created properly! expected {expected} mapped and {expected} sim, got {len(mapped_paths)} and {len(sim_paths)}", flush=True)