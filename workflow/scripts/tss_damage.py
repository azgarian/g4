
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

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

    # 'HD64_0_rep1_norm': '0m 64 (2022 rep1)',
    # 'HD64_0_rep2_norm': '0m 64 (2022 rep2)',
    # 'HDC_0_rep1_norm': '0m CPD (2022 rep1)',
    # 'HDC_0_rep2_norm': '0m CPD (2022 rep2)',
    # 'HX64_12_rep1_norm': '12m 64 XR (2022 rep1)',
    # 'HX64_12_rep2_norm': '12m 64 XR (2022 rep2)',
    # 'HXC_12_rep1_norm': '12m CPD XR (2022 rep1)',
    # 'HXC_12_rep2_norm': '12m CPD XR (2022 rep2)',
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

df = pd.read_csv("results/master/tss_201_100/mapped_all.csv", index_col=None)

df['ori'] = (df['strand'] == df['dam_strand']).map({True: 'NTS', False: 'TS'})

def plot_cosmetics(df):
    win_num = df["window"].max()
    interval_length = df.iloc[0]["end"] - df.iloc[0]["start"]
    boundary = win_num * interval_length / 1000
    p_ticks = [df["window"].min(),df["window"].min()/2,df["window"].median(),df["window"].max()/2,df["window"].max()]
    p_ticklabs = [f"-{boundary}kb",f"-{boundary/2:.1f}kb",0,f"-{boundary/2:.1f}kb",f"+{boundary}kb"]

    return p_ticks, p_ticklabs

p_ticks, p_ticklabs = plot_cosmetics(df)


df_tss = df[(df['name'].str.contains('tss'))]
df_tes = df[(df['name'].str.contains('tes'))]

target_list = ['ori', 'window'] + list(df.columns[df.columns.str.startswith('H')])
df_tss_agg = df_tss[target_list].groupby(["ori", "window"]).mean(numeric_only=True).reset_index()
df_tes_agg = df_tes[target_list].groupby(["ori", "window"]).mean(numeric_only=True).reset_index()

df_tss_agg_real = df_tss_agg.copy()
df_tes_agg_real = df_tes_agg.copy()
df_tss_agg_sim = df_tss_agg.copy()
df_tes_agg_sim = df_tes_agg.copy()

for column in list(df_tss.columns[df_tss.columns.str.startswith('H') & ~df_tss.columns.str.contains('sim')]):
    df_tss_agg[f"{column}_norm"] = df_tss_agg[column] / df_tss_agg[f"{column}_sim"]
    df_tss_agg_real[f"{column}_norm"] = df_tss_agg_real[column]
    df_tss_agg_sim[f"{column}_norm"] = df_tss_agg_sim[f"{column}_sim"]

    df_tss_agg_real.drop(columns=[column, f"{column}_sim"], inplace=True)
    df_tss_agg_sim.drop(columns=[column, f"{column}_sim"], inplace=True)
    df_tss_agg.drop(columns=[column, f"{column}_sim"], inplace=True)

for column in list(df_tes.columns[df_tes.columns.str.startswith('H') & ~df_tes.columns.str.contains('sim')]):
    df_tes_agg[f"{column}_norm"] = df_tes_agg[column] / df_tes_agg[f"{column}_sim"]
    df_tes_agg_real[f"{column}_norm"] = df_tes_agg_real[column]
    df_tes_agg_sim[f"{column}_norm"] = df_tes_agg_sim[f"{column}_sim"]

    df_tes_agg_real.drop(columns=[column, f"{column}_sim"], inplace=True)
    df_tes_agg_sim.drop(columns=[column, f"{column}_sim"], inplace=True)
    df_tes_agg.drop(columns=[column, f"{column}_sim"], inplace=True)


def melt_and_rename(df):
    df_melted = pd.melt(df, id_vars=['ori', 'window'], var_name='metric', value_name='value')
    df_melted['metric'] = df_melted['metric'].map(metric_rename)
    df_melted = df_melted[~df_melted["metric"].isna()]
    return df_melted

df_tss_melted = melt_and_rename(df_tss_agg)
df_tes_melted = melt_and_rename(df_tes_agg)
df_tss_real_melted = melt_and_rename(df_tss_agg_real)
df_tes_real_melted = melt_and_rename(df_tes_agg_real)
df_tss_sim_melted = melt_and_rename(df_tss_agg_sim)
df_tes_sim_melted = melt_and_rename(df_tes_agg_sim)


def plot_tss_damage(df_tss_melted, df_tes_melted, damage, ori, signal):

    # Create figure with two subplots side by side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 4), sharey=True)

    if ori == 'TS':
        title = 'Transcribed Strand'
    else:
        title = 'Non-transcribed Strand'

    if damage == '64':
        legend_order = legend_order_2025_64
    else:
        legend_order = legend_order_2025_CPD


    fig.suptitle(title, fontsize=14, fontweight='bold', x=0.45)

    # Plot TSS
    sns.lineplot(data=df_tss_melted[(df_tss_melted['metric'].str.contains(damage)) & (df_tss_melted['ori'] == ori)], 
                x='window', y='value', hue='metric', hue_order=legend_order, palette=metric_color, ax=ax1)
    ax1.set_xlabel('Distance from TSS (bp)')
    ax1.set_ylabel(signal)
    ax1.legend_.remove()  # Remove individual legend
    ax1.set_xticks(p_ticks)
    ax1.set_xticklabels(p_ticklabs)


    # Plot TES
    sns.lineplot(data=df_tes_melted[(df_tes_melted['metric'].str.contains(damage)) & (df_tes_melted['ori'] == ori)], 
                x='window', y='value', hue='metric', hue_order=legend_order, palette=metric_color, ax=ax2)
    ax2.set_xlabel('Distance from TES (bp)')
    ax2.set_ylabel(signal)
    # ax2.legend_.remove()  # Remove individual legend
    ax2.set_xticks(p_ticks)
    ax2.set_xticklabels(p_ticklabs)

    plt.tight_layout()
    plt.subplots_adjust(right=0.85)  # Make room for the legend

    plt.show()

plot_tss_damage(df_tss_melted, df_tes_melted, '64', 'TS', 'Normalized Damage signal')
plot_tss_damage(df_tss_melted, df_tes_melted, 'CPD', 'TS', 'Normalized Damage signal')
plot_tss_damage(df_tss_melted, df_tes_melted, '64', 'NTS', 'Normalized Damage signal')
plot_tss_damage(df_tss_melted, df_tes_melted, 'CPD', 'NTS', 'Normalized Damage signal')

plot_tss_damage(df_tss_real_melted, df_tes_real_melted, '64', 'TS', 'Damage signal')
plot_tss_damage(df_tss_real_melted, df_tes_real_melted, 'CPD', 'TS', 'Damage signal')
plot_tss_damage(df_tss_real_melted, df_tes_real_melted, '64', 'NTS', 'Damage signal')
plot_tss_damage(df_tss_real_melted, df_tes_real_melted, 'CPD', 'NTS', 'Damage signal')

plot_tss_damage(df_tss_sim_melted, df_tes_sim_melted, '64', 'TS', 'Simulated Damage signal')
plot_tss_damage(df_tss_sim_melted, df_tes_sim_melted, 'CPD', 'TS', 'Simulated Damage signal')
plot_tss_damage(df_tss_sim_melted, df_tes_sim_melted, '64', 'NTS', 'Simulated Damage signal')
plot_tss_damage(df_tss_sim_melted, df_tes_sim_melted, 'CPD', 'NTS', 'Simulated Damage signal')