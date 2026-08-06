# Installing q2-classo

## Conda Environment

q2-classo installs into the same QIIME 2 2026.7 environment as q2-gglasso — if
you already created it while following the previous page, activate it and skip
straight to the clone step.

```{note}
Two upstream renames landed in QIIME 2 **2026.4**: the `amplicon` distribution is
now called `qiime2` (so paths contain `/qiime2/`, not `/amplicon/`), and the
environment files are named `rachis-*` because the framework package was
rebranded from `qiime2` to `rachis`. A compatibility shim keeps `import qiime2`
working.
```

```bash
# Create the QIIME 2 2026.7 environment if you do not already have it
# (linux-64 shown; osx-64 also available, osx-arm64 is not)
conda env create \
  --name qiime2-2026.7 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.7/qiime2/released/rachis-qiime2-linux-64-conda.yml

# Activate the environment
conda activate qiime2-2026.7

# Clone and install q2-classo
git clone https://github.com/Vlasovets/q2-classo.git
cd q2-classo
python -m pip install --no-cache-dir -r requirements.txt
pip install -e .

# Refresh QIIME 2 cache
qiime dev refresh-cache
```

```{note}
Older instructions pointed at a repository called `q2-classo-latest` under the
`bio-datascience` organisation. The canonical repository is now
[`Vlasovets/q2-classo`](https://github.com/Vlasovets/q2-classo). Likewise,
`python setup.py install` is deprecated and redundant with `pip install -e .`.
```

## Docker Installation

Docker image of q2-classo is available through Docker Hub:

```bash
# Pull the Docker image
docker pull ovlasovets/q2-classo:latest

# Run q2-classo container with volume mapping for data
docker run -it -v $(pwd):/data ovlasovets/q2-classo:latest

# Alternative: Run specific analysis with data directory
docker run -it -v /path/to/your/data:/data ovlasovets/q2-classo:latest qiime classo --help
```

```{warning}
The published `:latest` image is built on the retired `amplicon` base image and
has not yet been rebuilt for 2026.7. Prefer the conda instructions above until it
has been.
```

### Verification

To verify that q2-classo is correctly installed:

```bash
# Check that classo is available
qiime classo --help
```

You should see all eight actions: `add-covariates`, `add-taxa`, `classify`,
`generate-data`, `predict`, `regress`, `summarize` and `transform-features`.

```{note}
In releases before this one, `classify` was mistakenly registered under the name
`regress`, so `qiime classo --help` showed two actions with the same name. If you
see that, your q2-classo predates the fix.
```
