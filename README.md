## Minimal Snakemake Template for HPC/SLURM

This template helps our lab members structure and run data pipelines locally and on SLURM clusters using Snakemake. It includes:
- A ready-to-run example pipeline
- A working SLURM profile using the cluster-generic executor plugin
- Clear hooks for adding your own rules, environments, and configs

### Repository layout
```
Snakefile                       # Main workflow
config/
  ├─ config.yaml                # Project config (samples, params)
  └─ slurm/
     └─ config.yaml            # SLURM profile for cluster runs
workflow/
  ├─ rules/
  │  ├─ common.smk             # Shared helpers (e.g., is_odd)
  │  └─ example_rules.smk      # Example rules (optional include)
  └─ envs/
     └─ primary_env.yaml       # Example conda environment
logs/
  └─ cluster/                  # SLURM stderr by job (from profile)
results/                        # Outputs (created by running the workflow)
README.md
```

### Quick start
- **Create and activate the env** (includes Snakemake and the cluster plugin):
```bash
cd /cta/users/cazgari/pipelines/snakemake_template
conda env create -f workflow/envs/primary_env.yaml
conda activate snakemake_example_env
```
- **Dry-run the workflow** (no files created):
```bash
snakemake -n -r -p
```
- **Run locally** (1 core):
```bash
snakemake --cores 1
```
- **Run on SLURM** (using the bundled profile):
```bash
snakemake --profile config/slurm
```
Add `--use-conda` if you start using per-rule conda environments
(already included in `config/slurm/config.yaml`).

### What this example pipeline does
The example workflow generates two outputs per sample from `config/config.yaml`:
- `results/{sample}_number_of_lines.txt` from `results/{sample}_raw.txt`
- `results/{sample}_length.txt` computed from `results/{sample}_raw.txt`

Rules in `Snakefile`:
- `example_rule`: creates `results/{sample}_raw.txt` with a small message
- `example_rule2`: writes the line count of the raw file to `..._number_of_lines.txt`
- `example_rule3`: computes whether the raw file line count is odd and writes it to `..._length.txt`
- `rule all`: requests both outputs for every sample defined in `config/config.yaml`

You can also inspect `workflow/rules/example_rules.smk` for a modular version of these rules.

### Configure your samples and parameters
Edit `config/config.yaml`:
```yaml
samples:
  - sample1
  - sample2

params:
  threads: 4
  memory: "8G"
```
- **samples**: list the sample identifiers you want to process
- **params**: optional project-wide defaults you can reference from rules

### SLURM profile (cluster execution)
The profile in `config/slurm/config.yaml` is preconfigured for a working HPC setup using the cluster-generic executor plugin. Key fields:
- **default-resources**: sets `account` and `partition` used in the `sbatch` command
- **cluster-generic-submit-cmd**: how jobs are submitted; uses `{threads}` and `{resources.*}`
- **logs**: SLURM stderr is written to `logs/cluster/{rule}_{jobid}.err`

Update to match your cluster:
```yaml
default-resources: [account='YOUR_ACCOUNT', partition='YOUR_PARTITION']
```
Per-rule overrides (inside a rule):
```python
rule example_rule:
    output: "results/{sample}_raw.txt"
    threads: 8                          # Overrides threads
    resources:                          # Overrides resources (account, partition, etc)
        account="your_account"
        partition="long"
    shell: "echo hello > {output}"

```


Run on cluster:
```bash
snakemake --profile config/slurm
```

### Using conda environments
Add a per-rule `conda:` directive to pin software reproducibly. Example:
```python
rule example_rule:
    conda: "workflow/envs/primary_env.yaml" # environment in which the rule will run
    output: "results/{sample}_raw.txt"
    shell: "echo hello > {output}"
```
Then run with `--use-conda`:
```bash
snakemake --cores 1 --use-conda
# or
snakemake --profile config/slurm # already provided inside profile: software-deployment-method: "conda"
```
Create additional env files in `workflow/envs/` as needed.

### Adapting this template to your pipeline
1. **List your samples** in `config/config.yaml`.
2. **Write rules** for each processing step. Use clear, non-overlapping outputs.
3. **Chain rules** via `input:`/`output:` so Snakemake infers dependencies.
4. **Update `rule all`** to expand the final targets you care about, e.g.:
```python
rule all:
    input:
        expand("results/{sample}.final.bam", sample=config["samples"]) 
```
5. **Test** with a dry-run: `snakemake -n -r -p`.
6. **Run locally**, then move to **SLURM** with `--profile config/slurm`.

Tips for robust rules:
- Prefer deterministic, unique output names per rule to avoid ambiguity
- Use `wildcard_constraints` (see `Snakefile`) to restrict wildcard matches
- Use `log:` and `benchmark:` to capture execution details
- Pass small computed parameters via `params:` (see `example_rule3` using `is_odd` from `workflow/rules/common.smk`)

### Common operations
- **Dry-run with reasons and commands**:
```bash
snakemake -n -r -p
```
- **Force a single rule for one sample**:
```bash
snakemake results/sample1_number_of_lines.txt -R example_rule2 -p
```
- **Unlock a stuck working directory**:
```bash
snakemake --unlock
```
- **Rerun incomplete tasks**:
```bash
snakemake --rerun-incomplete
```
- **Clean outputs created by rules** (preview, then execute):
```bash
snakemake --delete-all-output -n
snakemake --delete-all-output
```

### Troubleshooting
- **AmbiguousRuleException**: Ensure two rules never write the same file. Use distinct output patterns and/or stricter `wildcard_constraints`. Keep `rule all` targets aligned with actual final outputs.
- **Missing input files**: Check rule `input:` paths and the dependency chain. Use `-n -r` to see which rule should create a missing file.
- **Profile submit errors**: Verify `account` and `partition` in `config/slurm/config.yaml`. Confirm that `sbatch` is available and you can submit a minimal job.
- **Conda not activated on cluster**: Always pass `--use-conda` when using `conda:` directives and ensure the compute nodes can access your conda installation.

### Where things go
- Rule logs: `logs/{rule}_{sample}.log`
- Benchmarks: `logs/{rule}_{sample}.benchmark.log`
- SLURM stderr: `logs/cluster/{rule}_{jobid}.err`
- Outputs: `results/`

This template is meant to be simple but extensible. Start small, keep outputs unique, and let Snakemake manage the dependency graph for you. If you get stuck, run a dry-run with `-n -r -p` and read the log paths shown for the failing rule.
