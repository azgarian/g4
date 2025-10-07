import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def g4_to_uv(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["strand_info"] = np.where(
        ((df["strand"] == "+") & (df["dam_strand"] == "+"))
        | ((df["strand"] == "-") & (df["dam_strand"] == "-")),
        "on_G4",
        "opposite_G4",
    )
    return df


def organize_agg(df):
    df_new = df[["HXC_12_rep1","HDC_0_rep1", "HXC_12_rep1_sim", "HDC_0_rep1_sim","HX64_12_rep1","HD64_0_rep1","HX64_12_rep1_sim","HD64_0_rep1_sim","window","strand_info"]].groupby(["window","strand_info"]).mean().reset_index()
    df_new["C_rr"] = df_new["HXC_12_rep1"] / df_new["HDC_0_rep1"]
    df_new["C_sim_rr"] = df_new["HXC_12_rep1_sim"] / df_new["HDC_0_rep1_sim"]
    df_new["C_ds_rs"] = df_new["HDC_0_rep1"] / df_new["HDC_0_rep1_sim"]
    df_new["C_xr_rs"] = df_new["HXC_12_rep1"] / df_new["HXC_12_rep1_sim"]
    df_new["C_rr_rs"] = df_new["C_rr"] / df_new["C_sim_rr"]
    df_new["64_rr"] = df_new["HX64_12_rep1"] / df_new["HD64_0_rep1"]
    df_new["64_sim_rr"] = df_new["HX64_12_rep1_sim"] / df_new["HD64_0_rep1_sim"]
    df_new["64_ds_rs"] = df_new["HD64_0_rep1"] / df_new["HD64_0_rep1_sim"]
    df_new["64_xr_rs"] = df_new["HX64_12_rep1"] / df_new["HX64_12_rep1_sim"]
    df_new["64_rr_rs"] = df_new["64_rr"] / df_new["64_sim_rr"]

    return df_new

def g4_plot(df, column_list, dam_name, p_ticks, p_ticklabs, outdir):

    fig, axes = plt.subplots(3, 3, figsize=(25, 10), sharex=True)
    strand_colors = {"on_G4": "#E63946", "opposite_G4": "#457B9D"}

    g=sns.lineplot( ax=axes[0,0],
        data=df,
        x="window", y=column_list[0],
        hue="strand_info",
        palette=strand_colors,
    )
    g.set_ylabel("Damage RPKM", fontsize=14)
    g.set_xlabel("")
    g.set_xticks(p_ticks)
    g.set_xticklabels(p_ticklabs)
    g.set_ylim(0.09, 0.41)

    g=sns.lineplot( ax=axes[0,1],
        data=df,
        x="window", y=column_list[1],
        hue="strand_info",
        palette=strand_colors,
        legend=False,
    )
    g.set_title(f'{dam_name}', fontsize=14)
    g.set_ylabel("Repair RPKM", fontsize=14)
    g.set_xlabel("")
    g.set_xticks(p_ticks)
    g.set_xticklabels(p_ticklabs)
    g.set_ylim(0.09, 0.41)
    
    g=sns.lineplot( ax=axes[0,2],
        data=df,
        x="window", y=column_list[2],
        hue="strand_info",
        palette=strand_colors,
        legend=False,
    )
    g.set_ylabel("Repair Rate", fontsize=14)
    g.set_xlabel("")
    g.set_xticks(p_ticks)
    g.set_xticklabels(p_ticklabs)
    g.set_ylim(0.5, 2)
    
    g=sns.lineplot( ax=axes[1,0],
        data=df,
        x="window", y=column_list[3],
        hue="strand_info",
        palette=strand_colors,
        legend=False,
    )
    g.set_ylabel("Simulated Damage RPKM", fontsize=14)
    g.set_xlabel("")
    g.set_xticks(p_ticks)
    g.set_xticklabels(p_ticklabs)
    g.set_ylim(0.09, 0.41)
    
    g=sns.lineplot( ax=axes[1,1],
        data=df,
        x="window", y=column_list[4],
        hue="strand_info",
        palette=strand_colors,
        legend=False,
    )
    g.set_ylabel("Simulated Repair RPKM", fontsize=14)
    g.set_xlabel("")
    g.set_xticks(p_ticks)
    g.set_xticklabels(p_ticklabs)
    g.set_ylim(0.09, 0.41)

    g=sns.lineplot( ax=axes[1,2],
        data=df,
        x="window", y=column_list[5],
        hue="strand_info",
        palette=strand_colors,
        legend=False,
    )
    g.set_ylabel("Simulated Repair Rate", fontsize=14)
    g.set_xlabel("")
    g.set_xticks(p_ticks)
    g.set_xticklabels(p_ticklabs)
    g.set_ylim(0.5, 2)
    
    g=sns.lineplot( ax=axes[2,0],
        data=df,
        x="window", y=column_list[6],
        hue="strand_info",
        palette=strand_colors,
        legend=False,
    )
    g.set_ylabel("Damage / Simulation", fontsize=14)
    g.set_xlabel("G4 center", fontsize=14)
    g.set_xticks(p_ticks)
    g.set_xticklabels(p_ticklabs)
    g.set_ylim(0.5, 2)
    
    g=sns.lineplot( ax=axes[2,1],
        data=df,
        x="window", y=column_list[7],
        hue="strand_info",
        palette=strand_colors,
        legend=False,
    )
    g.set_ylabel("Repair / Simulation", fontsize=14)
    g.set_xlabel("G4 center", fontsize=14)
    g.set_xticks(p_ticks)
    g.set_xticklabels(p_ticklabs)
    g.set_ylim(0.5, 2)
    
    g=sns.lineplot( ax=axes[2,2],
        data=df,
        x="window", y=column_list[8],
        hue="strand_info",
        palette=strand_colors,
        legend=False,
    )
    g.set_ylabel("Normalized Repair Rate", fontsize=14)
    g.set_xlabel("G4 center", fontsize=14)
    g.set_xticks(p_ticks)
    g.set_xticklabels(p_ticklabs)
    g.set_ylim(0.5, 2)

    dam = dam_name.lower().replace(" ","_")
    plt.savefig(f"{outdir}/{dam}.png", dpi=200, bbox_inches='tight')
    plt.savefig(f"{outdir}/{dam}.svg", bbox_inches='tight')
    plt.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True, help="Region set name (directory under results/master)")
    ap.add_argument(
        "--input", required=False, help="Optional CSV path; overrides results/master/{name}/mapped_all.csv"
    )
    ap.add_argument("--outdir", required=True, help="Output directory for aggregated TSV")
    args = ap.parse_args()

    in_csv = args.input or f"results/master/{args.name}/mapped_all.csv"
    df = pd.read_csv(in_csv)

    df = g4_to_uv(df)
    df_repair_agg = organize_agg(df)

    # ================= Plotting (ported from notebook) =================
    # Cosmetics: ticks and labels based on window count and interval length
    def g4_plot_cosmetics(df_in: pd.DataFrame):
        win_num = int(df_in["window"].max())
        interval_length = int(df_in.iloc[0]["end"] - df_in.iloc[0]["start"]) if {"start", "end"}.issubset(df_in.columns) else 1
        boundary = win_num * interval_length / 1000
        p_ticks = [
            int(df_in["window"].min()),
            int(df_in["window"].min() / 2),
            int(df_in["window"].median()),
            int(df_in["window"].max() / 2),
            int(df_in["window"].max()),
        ]
        p_ticklabs = [
            f"-{boundary}kb",
            f"-{boundary/2:.1f}kb",
            0,
            f"-{boundary/2:.1f}kb",
            f"+{boundary}kb",
        ]
        return p_ticks, p_ticklabs

    p_ticks, p_ticklabs = g4_plot_cosmetics(df)

    list_64 = ["HD64_0_rep1", "HX64_12_rep1", "64_rr", "HD64_0_rep1_sim", "HX64_12_rep1_sim", "64_sim_rr", "64_ds_rs", "64_xr_rs", "64_rr_rs"]
    list_CPD = ["HDC_0_rep1", "HXC_12_rep1", "C_rr", "HDC_0_rep1_sim", "HXC_12_rep1_sim", "C_sim_rr", "C_ds_rs", "C_xr_rs", "C_rr_rs"]

    g4_plot(df_repair_agg, list_CPD, "All CPD", p_ticks, p_ticklabs, args.outdir)
    g4_plot(df_repair_agg, list_64, "All 64", p_ticks, p_ticklabs, args.outdir)

    # Renaming and colors from notebook
    metric_rename = {
        'HD64_0_rep3_norm': '0m 64',
        'HD64_15_rep1_norm': '15m 64',
        'HD64_30_rep2_norm': '30m 64',
        'HD64_60_rep2_norm': '60m 64',
        'HD64_240_rep1_norm': '240m 64',
        'HD64_480_rep1_norm': '480m 64',
        'HDC_0_rep3_norm': '0m CPD',
        'HDC_15_rep1_norm': '15m CPD',
        'HDC_30_rep2_norm': '30m CPD',
        'HDC_60_rep2_norm': '60m CPD',
        'HDC_240_rep1_norm': '240m CPD',
        'HDC_480_rep1_norm': '480m CPD',
    }

    metric_color = {
        '0m 64': '#ef476f',
        '15m 64': '#ffd166',
        '30m 64': '#06d6a0',
        '60m 64': '#118ab2',
        '240m 64': '#073b4c',
        '480m 64': 'purple',
        '0m CPD': '#ef476f',
        '15m CPD': '#ffd166',
        '30m CPD': '#06d6a0',
        '60m CPD': '#118ab2',
        '240m CPD': '#073b4c',
        '480m CPD': 'purple',
    }

    legend_order_2025_64 = ['0m 64', '15m 64', '30m 64', '60m 64', '240m 64', '480m 64']
    legend_order_2025_CPD = ['0m CPD', '15m CPD', '30m CPD', '60m CPD', '240m CPD', '480m CPD']

    # Build aggregated frames for plotting following the notebook flow
    target_list = ['strand_info', 'window'] + [c for c in df.columns if c.startswith('H')]
    df_agg_plot = df[target_list].groupby(["strand_info", "window"]).mean(numeric_only=True).reset_index()
    df_agg_real = df_agg_plot.copy()
    df_agg_sim = df_agg_plot.copy()

    for column in [c for c in df.columns if c.startswith('H') and ('sim' not in c) and (f"{c}_sim" in df.columns)]:
        df_agg_plot[f"{column}_norm"] = df_agg_plot[column] / df_agg_plot[f"{column}_sim"].replace(0, np.nan)
        df_agg_real[f"{column}_norm"] = df_agg_real[column]
        df_agg_sim[f"{column}_norm"] = df_agg_sim[f"{column}_sim"]

        for frame in (df_agg_plot, df_agg_real, df_agg_sim):
            for drop_col in (column, f"{column}_sim"):
                if drop_col in frame.columns:
                    frame.drop(columns=[drop_col], inplace=True)

    def melt_and_rename(df_in: pd.DataFrame) -> pd.DataFrame:
        df_melted = pd.melt(df_in, id_vars=['strand_info', 'window'], var_name='metric', value_name='value')
        df_melted['metric'] = df_melted['metric'].map(metric_rename)
        df_melted = df_melted[~df_melted["metric"].isna()]
        return df_melted

    df_agg_melted = melt_and_rename(df_agg_plot)
    df_agg_real_melted = melt_and_rename(df_agg_real)
    df_agg_sim_melted = melt_and_rename(df_agg_sim)

    def plot_damage(df_melted: pd.DataFrame, damage: str, signal: str, tag: str):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 4), sharey=True)
        legend_order = legend_order_2025_64 if damage == '64' else legend_order_2025_CPD

        sns.lineplot(
            data=df_melted[(df_melted['metric'].str.contains(damage)) & (df_melted['strand_info'] == "on_G4")],
            x='window', y='value', hue='metric', hue_order=legend_order, palette=metric_color, ax=ax1
        )
        ax1.set_xlabel('G4 center (bp)')
        ax1.set_ylabel(signal)
        if ax1.legend_:
            ax1.legend_.remove()
        ax1.set_xticks(p_ticks)
        ax1.set_xticklabels(p_ticklabs)
        ax1.set_title("G4 Strand")

        sns.lineplot(
            data=df_melted[(df_melted['metric'].str.contains(damage)) & (df_melted['strand_info'] == "opposite_G4")],
            x='window', y='value', hue='metric', hue_order=legend_order, palette=metric_color, ax=ax2
        )
        ax2.set_xlabel('G4 center (bp)')
        ax2.set_ylabel(signal)
        if ax2.legend_:
            ax2.legend_.remove()
        ax2.set_xticks(p_ticks)
        ax2.set_xticklabels(p_ticklabs)
        ax2.set_title("Opposite G4 Strand")

        handles, labels = ax1.get_legend_handles_labels()
        fig.legend(handles, labels, loc='center right')
        plt.tight_layout()
        plt.subplots_adjust(right=0.85)

        base = Path(args.outdir) / f"{args.name}_{damage}_{tag}"
        fig.savefig(f"{base}.png", dpi=200, bbox_inches='tight')
        fig.savefig(f"{base}.svg", bbox_inches='tight')
        plt.close(fig)

    # Generate and save figures (PNG and SVG), mirroring the notebook calls
    plot_damage(df_agg_melted, '64', 'Normalized Damage signal', 'normalized')
    plot_damage(df_agg_melted, 'CPD', 'Normalized Damage signal', 'normalized')
    plot_damage(df_agg_real_melted, '64', 'Damage signal', 'real')
    plot_damage(df_agg_real_melted, 'CPD', 'Damage signal', 'real')
    plot_damage(df_agg_sim_melted, '64', 'Simulated Damage signal', 'sim')
    plot_damage(df_agg_sim_melted, 'CPD', 'Simulated Damage signal', 'sim')


if __name__ == "__main__":
    main()


