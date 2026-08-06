# Installing q2-gglasso

## Conda Environment

Create a dedicated conda environment for q2-gglasso.

```{note}
Two upstream renames landed in QIIME 2 **2026.4** and change every install
command you may have seen in older tutorials:

* the **`amplicon` distribution is now called `qiime2`**, so the channel and file
  paths contain `/qiime2/` rather than `/amplicon/`;
* the environment files are named **`rachis-*`** rather than `qiime2-*`, because
  the framework package was rebranded from `qiime2` to `rachis`. A compatibility
  shim keeps `import qiime2` working, so existing scripts do not need changes.
```

Pick the file matching your platform. `linux-64` and `osx-64` are available;
there is no `osx-arm64` build of this distribution, so Apple Silicon users should
run the `osx-64` build under Rosetta or use Docker.

```bash
# Create the QIIME 2 2026.7 environment (linux-64 shown)
conda env create \
  --name qiime2-2026.7 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.7/qiime2/released/rachis-qiime2-linux-64-conda.yml

# Activate the environment
conda activate qiime2-2026.7

# Clone and install q2-gglasso
git clone https://github.com/Vlasovets/q2-gglasso.git
cd q2-gglasso
python -m pip install --no-cache-dir -r requirements.txt
pip install -e .

# Refresh QIIME 2 cache
qiime dev refresh-cache
```

```{note}
If `conda env create` fails with

    package deblur-1.1.1 requires sortmerna 2.0, but none of the providers
    can be installed

you have hit a known defect in the upstream 2026.7 `linux-64` file: it pins
`zlib=1.3.2` while every `sortmerna` 2.0 build requires `zlib <1.3`. `deblur` is
not used anywhere in this tutorial. Download the environment file, delete the
`deblur`, `q2-deblur` and `sortmerna` lines, and create the environment from your
edited copy.
```

```{note}
`python setup.py install` used to appear in these instructions. It is deprecated
and redundant with `pip install -e .` — use the latter only.
```

## Docker Installation

Docker image of q2-gglasso is available through Docker Hub:

```bash
# Pull the Docker image
docker pull ovlasovets/q2-gglasso:latest

# Run q2-gglasso container with volume mapping for data
docker run -it -v $(pwd):/data ovlasovets/q2-gglasso:latest

# Alternative: Run specific analysis with data directory
docker run -it -v /path/to/your/data:/data ovlasovets/q2-gglasso:latest qiime gglasso --help
```

```{note}
The published `:latest` image is built on the retired `amplicon` base image and
has not yet been rebuilt for 2026.7. Prefer the conda instructions above until it
has been.
```

### Verification

To verify that q2-gglasso is correctly installed:

```bash
# Check that gglasso is available
qiime gglasso --help
```

You should see all six actions listed: `build-groups`, `calculate-covariance`,
`pca`, `solve-problem`, `summarize` and `transform-features`.

The installation is now complete! You can proceed to explore the plugin's functionality.
