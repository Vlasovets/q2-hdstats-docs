#!/usr/bin/env python
"""Render static figures from the committed .qzv visualizations.

Three of the figures in the shotgun section are not drawn by this repository at all --
they are the visualizations q2-gglasso and q2-classo emit, screenshotted from the very
`.qzv` files the chapter links for download. That is the point: a reader who opens
`motus-sgl-summary.qzv` at https://view.qiime2.org sees the same image the page shows,
not a redrawing of it that could drift from the plugin's actual output.

The alternative -- reimplementing each plot in matplotlib from an exported TSV -- is what
the other two figures in that section do, and it is the right choice *only* where the
plugin emits no such view (the eBIC path and the mu/rank map are computed by this
analysis, not by q2-gglasso). Where the plugin does draw the thing, screenshot the
plugin.

Inputs are the .qzv files under docs/_static/qzv/, which are committed. No cluster and no
QIIME 2 are needed -- only a headless browser, because a .qzv is a web page:

    pip install playwright && playwright install chromium
    python analysis/scripts/render_qzv_figures.py

Notes on the two settings that matter:

  * The bokeh toolbar (zoom, reset, save) and the plotly modebar are interactive chrome.
    They carry no information in a static page and are hidden before capture.
  * `device_scale_factor=2` renders at 2x. The precision-matrix heatmap is 100x100 cells
    with per-feature tick labels; at 1x the mOTU identifiers are unreadable in print.
"""
import argparse
import pathlib
import sys
import tempfile
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
QZV_DIR = ROOT / "docs" / "_static" / "qzv"
DEFAULT_OUT = ROOT / "docs" / "images" / "png" / "generated"

# Interactive chrome, removed before capture.
HIDE_CSS = ".bk-toolbar, .bk-logo, .modebar { display: none !important; }"

# (qzv, bokeh tab label or None for a plotly sub-page, sub-page, output name)
FIGURES = [
    ("motus-sgl-summary.qzv", "Estimated inverse covariance", "index.html",
     "motus-qzv-precision.png"),
    ("motus-latent-pca.qzv", "Single plot", "index.html", "motus-qzv-pca.png"),
    ("motus-regress-summary.qzv", None, "cv-graph.html", "motus-qzv-cv.png"),
]


def render(pw, qzv, tab, page, out):
    tmp = tempfile.mkdtemp()
    with zipfile.ZipFile(qzv) as z:
        z.extractall(tmp)
    # A .qzv is <uuid>/data/<the actual web page>.
    data = next(p for p in pathlib.Path(tmp).rglob("data") if p.is_dir())

    browser = pw.chromium.launch(args=["--no-sandbox", "--disable-dev-shm-usage"])
    pg = browser.new_page(viewport={"width": 1400, "height": 1000}, device_scale_factor=2)
    pg.goto(f"file://{data / page}", wait_until="networkidle", timeout=90_000)
    pg.wait_for_timeout(4000)  # bokeh lays out after networkidle

    if tab is not None:
        for el in pg.query_selector_all(".bk-tab"):
            if el.inner_text().strip() == tab:
                el.click()
                pg.wait_for_timeout(3500)
                break
        else:
            sys.exit(f"{qzv.name}: no tab named {tab!r}")

    pg.add_style_tag(content=HIDE_CSS)
    pg.wait_for_timeout(1200)

    # Capture the plot, not the padded page. Inactive bokeh panels stay in the DOM but
    # hidden, and screenshotting a hidden element hangs waiting for it to be stable.
    best = None
    for sel in (".js-plotly-plot", ".bk-Figure"):
        for el in pg.query_selector_all(sel):
            if not el.is_visible():
                continue
            box = el.bounding_box()
            if box and box["width"] > 300 and (best is None or
                                               box["width"] * box["height"] > best[1]):
                best = (el, box["width"] * box["height"])
        if best:
            break
    if best is None:
        sys.exit(f"{qzv.name}: found no visible plot element to capture")
    best[0].screenshot(path=str(out))
    browser.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=pathlib.Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        sys.exit("playwright is required: pip install playwright && "
                 "playwright install chromium")

    missing = [n for n, _, _, _ in FIGURES if not (QZV_DIR / n).is_file()]
    if missing:
        sys.exit(f"missing visualizations in {QZV_DIR}: {missing}")

    with sync_playwright() as pw:
        for name, tab, page, out_name in FIGURES:
            out = args.outdir / out_name
            render(pw, QZV_DIR / name, tab, page, out)
            print(f"  {out.name:30} {out.stat().st_size / 1024:6.0f} KB  <- {name}")
    print(f"  -> {args.outdir}")


if __name__ == "__main__":
    main()
