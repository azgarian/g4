"""Helper utilities for window generation, BED binning, and coverage metrics."""

from typing import Any, Callable, Optional
import pandas as pd
import numpy as np
from datetime import datetime
from zoneinfo import ZoneInfo
from functools import wraps
import inspect
import os
import subprocess
import bioframe as bf


def timing_decorator(timezone: str = 'Europe/Istanbul') -> Callable[[Callable[..., Any]], Callable[..., Any]]:
    """Decorator that logs start and finish timestamps for function execution.

    If the decorated function has a parameter named `file`, its value is included
    in the log messages for additional context.

    Parameters:
    - timezone: IANA timezone identifier to format timestamps.
    """
    def decorator(func: Callable[..., Any]) -> Callable[..., Any]:
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            # Get the start time
            start_time = datetime.now(ZoneInfo(timezone))

            # Inspect the function's signature to find the 'file' argument
            try:
                sig = inspect.signature(func)
                bound_args = sig.bind(*args, **kwargs)
                bound_args.apply_defaults()
                file_value = bound_args.arguments.get('file', None)
            except Exception:
                file_value = None

            # Extract the value of the 'file' argument if it exists
            # Log the start message with the function name and file value
            if file_value:
                print(start_time.strftime('%Y-%m-%d %H:%M:%S'), f"{func.__name__} ({file_value}) started.", flush=True)
            else:
                print(start_time.strftime('%Y-%m-%d %H:%M:%S'), f"{func.__name__} started.", flush=True)

            # Execute the function
            result = func(*args, **kwargs)

            # Get the end time
            end_time = datetime.now(ZoneInfo(timezone))
            duration = end_time - start_time

            # Log the end message with the function name and file value
            if file_value:
                print(end_time.strftime('%Y-%m-%d %H:%M:%S'), f"{func.__name__} ({file_value}) finished. Duration:", duration, flush=True)
            else:
                print(end_time.strftime('%Y-%m-%d %H:%M:%S'), f"{func.__name__} finished. Duration:", duration, flush=True)

            return result
        return wrapper
    return decorator


@timing_decorator()
def split_bed_into_bins(bed_file: str, interval_length: int = 10, output_file: Optional[str] = None) -> Optional[pd.DataFrame]:
    """
    Split a BED file into bins with defined interval length.

    Parameters:
    -----------
    bed_file : str
        Path to the input BED file
    interval_length : int
        Length of each bin in base pairs (default: 10)
    output_file : str, optional
        Path to save the output BED file. If None, returns DataFrame

    Returns:
    --------
    pandas.DataFrame or None
        DataFrame with binned regions if output_file is None, otherwise saves to file
    """

    if interval_length <= 0:
        raise ValueError("interval_length must be a positive integer")

    # Read the BED file
    regions = pd.read_table(bed_file, header=None, sep="\t")

    # Set column names based on BED format
    if len(regions.columns) >= 6:
        regions.columns = ["chrom", "start", "end", "name", "score", "strand"]
    elif len(regions.columns) >= 4:
        regions.columns = ["chrom", "start", "end", "name"] + [f"col_{i}" for i in range(4, len(regions.columns))]
    else:
        regions.columns = ["chrom", "start", "end"] + [f"col_{i}" for i in range(3, len(regions.columns))]

    # Ensure integer coordinates
    regions["start"] = regions["start"].astype(int)
    regions["end"] = regions["end"].astype(int)

    # Create binned regions
    binned_regions = []

    for idx, row in regions.iterrows():
        chrom = row['chrom']
        start = int(row['start'])
        end = int(row['end'])

        if end <= start:
            continue

        # Create bins for this region
        bin_starts = np.arange(start, end, interval_length)
        bin_ends = np.minimum(bin_starts + interval_length, end)

        # Create DataFrame for this region's bins
        region_bins = pd.DataFrame({
            'chrom': chrom,
            'start': bin_starts,
            'end': bin_ends,
            'original_name': row.get('name', f'region_{idx}'),
            'bin_id': range(len(bin_starts)),
            'original_start': start,
            'original_end': end
        })

        # Add additional columns if they exist
        for col in ['score', 'strand']:
            if col in row:
                region_bins[col] = row[col]

        binned_regions.append(region_bins)

    # Combine all binned regions
    if binned_regions:
        result = pd.concat(binned_regions, ignore_index=True)

        # Create unique bin names
        result['name'] = result['original_name'] + '_' + result['bin_id'].astype(str)

        # Reorder columns
        cols = ['chrom', 'start', 'end', 'name']
        if 'score' in result.columns:
            cols.append('score')
        if 'strand' in result.columns:
            cols.append('strand')
        # cols.extend(['original_name', 'bin_id', 'original_start', 'original_end'])

        result = result[cols]

        # Save to file if specified
        if output_file:
            result.to_csv(output_file, sep='\t', index=False, header=False)
            print(f"Binned BED file saved to: {output_file}")
            return None
        else:
            return result
    else:
        print("No regions found in the BED file")
        return pd.DataFrame()

# Example usage:
# result = split_bed_into_bins("input.bed", interval_length=50, output_file="output_binned.bed")


@timing_decorator()
def make_windows(file: str, quant: bool = False, quant_num: int = 4, region_size: int = 1000, interval_length: int = 10, interval_num: Optional[int] = None, tss: bool = False) -> pd.DataFrame:
    """Generate fixed windows centered on regions or around TSS/TES.

    Expects an input file with 6 columns: chrom, start, end, name, signal, strand.

    Parameters:
    - file: Path to the BED-like file.
    - quant: If True, bin `signal` into `quant_num` quantiles.
    - quant_num: Number of quantiles for binning when `quant=True`.
    - region_size: Total span covered by the windowed region (bp).
    - interval_length: Size of each bin (bp). Overridden by `interval_num` if provided.
    - interval_num: If provided, set interval_length = region_size / interval_num.
    - tss: If True, generate windows around TSS and TES; otherwise center on region midpoints.
    """
    if interval_num:
        interval_length = region_size / interval_num
    region_size = int(region_size)
    interval_length = int(interval_length)
    if region_size <= 0 or interval_length <= 0:
        raise ValueError("region_size and interval_length must be positive integers")

    half_bins = region_size // (2 * interval_length)
    num_bins = 2 * half_bins + 1  # always odd, centered at 0

    regions = pd.read_table(file, header=None)
    regions.columns = ["chrom", "start", "end", "name", "signal", "strand"]

    if quant:
        regions['quantile'] = pd.qcut(regions["signal"], q=quant_num, labels=False)
    else:
        regions['quantile'] = regions["signal"]

    regions_org = regions[["chrom", "start", "end", "quantile", "strand"]].copy().reset_index(drop=True)

    if not tss:
        # Centered mode (original)
        regions_org["name"] = "region_" + regions_org.index.astype(str)
        regions_org.columns = ["chrom", "start", "end", "quantile", "strand", "name"]
        center = (regions_org["start"] + regions_org["end"]) // 2
        starts = np.tile(center.values - half_bins * interval_length, (num_bins, 1)).T + np.arange(num_bins) * interval_length
        ends = starts + interval_length
        window = np.arange(-half_bins, half_bins + 1)
        df = pd.DataFrame({
            "chrom": np.repeat(regions_org["chrom"].values, num_bins),
            "start": starts.flatten(),
            "end": ends.flatten(),
            "name": np.repeat(regions_org["name"].values, num_bins),
            "score": np.repeat(regions_org["quantile"].values, num_bins),
            "strand": np.repeat(regions_org["strand"].values, num_bins),
            "window": np.tile(window, len(regions_org))
        })
        return df

    # TSS/TES mode
    plus = regions_org["strand"] == "+"
    minus = regions_org["strand"] == "-"

    # TSS
    tss_coord = np.where(plus, regions_org["start"], regions_org["end"])  # kept for readability
    tss_name = ["tss_{}".format(i) for i in range(len(regions_org))]
    tss_starts = np.empty((len(regions_org), num_bins), dtype=int)
    tss_ends = np.empty((len(regions_org), num_bins), dtype=int)
    for i, w in enumerate(np.arange(-half_bins, half_bins + 1)):
        # The 0th window is centered on the TSS/TES coordinate
        tss_starts[plus, i] = regions_org.loc[plus, "start"] + w * interval_length - interval_length // 2
        tss_ends[plus, i]   = regions_org.loc[plus, "start"] + w * interval_length + interval_length // 2
        tss_starts[minus, i] = regions_org.loc[minus, "end"] - w * interval_length - interval_length // 2
        tss_ends[minus, i]   = regions_org.loc[minus, "end"] - w * interval_length + interval_length // 2
    tss_df = pd.DataFrame({
        "chrom": np.repeat(regions_org["chrom"].values, num_bins),
        "start": tss_starts.flatten(),
        "end": tss_ends.flatten(),
        "name": np.repeat(tss_name, num_bins),
        "score": np.repeat(regions_org["quantile"].values, num_bins),
        "strand": np.repeat(regions_org["strand"].values, num_bins),
        "window": np.tile(np.arange(-half_bins, half_bins + 1), len(regions_org))
    })

    # TES
    tes_coord = np.where(plus, regions_org["end"], regions_org["start"])  # kept for readability
    tes_name = ["tes_{}".format(i) for i in range(len(regions_org))]
    tes_starts = np.empty((len(regions_org), num_bins), dtype=int)
    tes_ends = np.empty((len(regions_org), num_bins), dtype=int)
    for i, w in enumerate(np.arange(-half_bins, half_bins + 1)):
        tes_starts[plus, i] = regions_org.loc[plus, "end"] + w * interval_length - interval_length // 2
        tes_ends[plus, i]   = regions_org.loc[plus, "end"] + w * interval_length + interval_length // 2
        tes_starts[minus, i] = regions_org.loc[minus, "start"] - w * interval_length - interval_length // 2
        tes_ends[minus, i]   = regions_org.loc[minus, "start"] - w * interval_length + interval_length // 2
    tes_df = pd.DataFrame({
        "chrom": np.repeat(regions_org["chrom"].values, num_bins),
        "start": tes_starts.flatten(),
        "end": tes_ends.flatten(),
        "name": np.repeat(tes_name, num_bins),
        "score": np.repeat(regions_org["quantile"].values, num_bins),
        "strand": np.repeat(regions_org["strand"].values, num_bins),
        "window": np.tile(np.arange(-half_bins, half_bins + 1), len(regions_org))
    })

    return pd.concat([tss_df, tes_df], ignore_index=True)


def calculate_rpkm(region_df: pd.DataFrame, read_counts_column: str, total_reads: int) -> pd.Series:
    """Calculate RPKM values for the given region DataFrame.

    Parameters:
        region_df (DataFrame): Must contain integer 'start' and 'end' columns defining region lengths.
        read_counts_column (str): Column containing read counts for each region.
        total_reads (int): Total number of reads in the dataset (line count of the file).

    Returns:
        Series: A pandas Series with RPKM values.
    """
    if total_reads <= 0:
        raise ValueError("total_reads must be a positive integer")

    # Region lengths in kilobases
    region_lengths_kb = (region_df["end"] - region_df["start"]) / 1000  # Length in kilobases
    if (region_lengths_kb <= 0).any():
        raise ValueError("All regions must have positive length")

    # Total reads in millions
    total_reads_million = total_reads / 1_000_000  # Convert to millions

    # RPKM calculation
    rpkm_values = region_df[read_counts_column] / (region_lengths_kb * total_reads_million)
    return rpkm_values

def exact_damage_site(df: pd.DataFrame, method: str, strand: str) -> pd.DataFrame:
    """Return a copy of df with coordinates adjusted to exact damage sites.

    Parameters:
        method: 'xr' or 'ds'.
        strand: '+' or '-'.
    """
    m = method.lower()
    if m not in {"xr", "ds"}:
        raise ValueError("method must be either 'xr' or 'ds'")
    if strand not in {"+", "-"}:
        raise ValueError("strand must be '+' or '-'")

    out = df.copy()
    if m == "xr":
        if strand == "+":
            out["start"] = out["end"] - 8
            out["end"] = out["end"] - 6
        else:  # '-'
            out["end"] = out["start"] + 8
            out["start"] = out["start"] + 6
    else:  # ds
        if strand == "+":
            out["start"] = out["end"] - 6
            out["end"] = out["end"] - 4
        else:  # '-'
            out["end"] = out["start"] + 6
            out["start"] = out["start"] + 4
    return out

def process_xr_ds(curr_file, expanded_regions, name, method, strand):

        df = bf.read_table(curr_file, schema="bed")
        df = exact_damage_site(df, method, strand)
        overlapped = bf.count_overlaps(expanded_regions, df)

        overlapped.rename(columns={"count": name}, inplace=True)

        result = subprocess.run(['wc', '-l', curr_file], capture_output=True, text=True)
        total_reads = int(result.stdout.split()[0])

        # with open(curr_file, 'r') as f:
        #     total_reads = sum(1 for _ in f)

        if "_plus.bed" in curr_file:
            opposite_strand = curr_file.replace("_plus.bed", "_minus.bed")
        else:
            opposite_strand = curr_file.replace("_minus.bed", "_plus.bed")

        result = subprocess.run(['wc', '-l', opposite_strand], capture_output=True, text=True)
        total_reads += int(result.stdout.split()[0])

        # with open(opposite_strand, 'r') as f_opp:
        #     total_reads += sum(1 for _ in f_opp)

        overlapped[name] = calculate_rpkm(
            region_df=overlapped,
            read_counts_column=name,
            total_reads=total_reads
        )

        return overlapped

@timing_decorator()
def map_xr_ds(file, expanded_regions, output="results/master", region_name="", interval=10):

    os.makedirs(f"{output}/{region_name}", exist_ok=True)
    name = os.path.basename(file).replace(".bed", "")

    if "XR" in name:
        method = "xr"
    elif "DS" in name:
        method = "ds"
    else:
        raise ValueError(f"Invalid method: {name}")

    for ext in ["", "_sim"]:
        if os.path.exists(f"{output}/{region_name}/{name}{ext}_mapped.csv"):
            print(f"{output}/{region_name}/{name}{ext}_mapped.csv already exists!", flush=True)
            continue
        else:
            curr_name = name + ext
            curr_file = file + ext + "_plus.bed"
            overlapped_plus = process_xr_ds(curr_file, expanded_regions, curr_name, method, "+")
            overlapped_plus.loc[:, "window"] *= interval
            overlapped_plus.loc[:, "dam_strand"] = "+"

            curr_file = file + ext + "_minus.bed"
            overlapped_minus = process_xr_ds(curr_file, expanded_regions, curr_name, method, "-")
            overlapped_minus.loc[:, "window"] *= interval * -1
            overlapped_minus.loc[:, "dam_strand"] = "-"

            overlapped = pd.concat([overlapped_plus, overlapped_minus], ignore_index=True)
            merge_keys = ["chrom", "start", "end", "name", "score", "strand", "dam_strand", "window"]
            overlapped = overlapped.sort_values(merge_keys).reset_index(drop=True)

            # if os.path.exists(f"{output}/{region_name}/{region_name}.csv") == False:
            #     overlapped[["chrom", "start", "end", "name", "score", "strand", "dam_strand", "window"]].to_csv(f"{output}/{region_name}/{region_name}.csv")

            # mapped_agg = overlapped[["strand", "dam_strand", "window", curr_name]].groupby(["strand", "dam_strand", "window"]).agg("mean").reset_index()
            overlapped[curr_name].to_csv(f"{output}/{region_name}/{name}{ext}_mapped.csv", index=False)
