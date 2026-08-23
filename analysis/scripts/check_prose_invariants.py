#!/usr/bin/env python
"""Gate a prose-only rewrite: everything that is not a sentence must be unchanged.

A style pass rewrites sentences. It must not move a number, edit a command, rename a
heading, drop a citation or repoint a link. Those are exactly the changes that are
invisible in review -- a reworded paragraph is obvious, a 0.30 that became 0.3 is not --
and this book's whole claim is that its numbers come from committed tables.

Compares each file against a git revision (default HEAD) and reports every difference in:

  headings          anchors are generated from them, so other pages link to them
  code fences       commands readers run verbatim; outputs are transcripts
  directive heads   ```{figure} path, ```{csv-table} title -- the type and its argument
  directive options :name:, :file:, :width:, :header-rows:, :delim:
  inline code       --p-gamma, mclr.qza, FeatureTable[Frequency]
  numbers           every numeral outside a code fence, in order
  citations         {cite}`key` -- resolved against references.bib at build time
  link targets      [text](target) and <target>
  math              $...$ and $$...$$
  admonition types  {note} must not become {warning} for emphasis

Prose inside directives IS rewritable, so ```{note} bodies are compared only for the
items above, not for wording.

    python analysis/scripts/check_prose_invariants.py                  # all chapters
    python analysis/scripts/check_prose_invariants.py --rev main docs/chapters/x.md
    python analysis/scripts/check_prose_invariants.py --verbose        # show matches too

Exits non-zero if any invariant moved. A finding is not automatically a defect -- fixing
a genuinely wrong number is legitimate -- but it must be a decision, not a side effect.
"""
import argparse
import collections
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]

FENCE = re.compile(r"^(\s*)(`{3,}|~{3,})(.*)$")
HEADING = re.compile(r"^\s{0,3}#{1,6}\s+\S")
OPTION = re.compile(r"^\s*:[A-Za-z0-9_-]+:")
# Numbers: integers, decimals, thousands separators, exponents, signed.
NUMBER = re.compile(r"[-+]?\d[\d,]*(?:\.\d+)?(?:[eE][-+]?\d+)?")
INLINE_CODE = re.compile(r"`([^`\n]+)`")
CITE = re.compile(r"\{cite[a-z:]*\}`([^`]+)`")
MDLINK = re.compile(r"\[[^\]]*\]\(([^)]+)\)")
AUTOLINK = re.compile(r"<(https?://[^>]+)>")
ROLE = re.compile(r"\{(ref|numref|doc|eq|term)\}`([^`]+)`")
MATH = re.compile(r"\$\$?[^$]+\$\$?")


def parse(text):
    """Split a MyST document into the parts a style pass may and may not touch."""
    out = {k: [] for k in ("headings", "code", "directive_heads", "options",
                           "inline_code", "numbers", "citations", "links", "roles",
                           "math", "admonitions")}
    lines = text.split("\n")
    prose_lines = []
    stack = []  # open fences: (indent, marker, is_directive)

    for line in lines:
        m = FENCE.match(line)
        if m:
            indent, marker, info = m.group(1), m.group(2), m.group(3).strip()
            if stack and marker[0] == stack[-1][1][0] and len(marker) >= len(stack[-1][1]) \
                    and not info:
                stack.pop()
                continue
            if not stack or stack[-1][2]:
                # Opening a new fence. Nested fences inside a directive are common
                # (a code block inside a {note}), so only track directive nesting.
                is_directive = info.startswith("{")
                if is_directive:
                    out["directive_heads"].append(info)
                    dtype = info[1:].split("}")[0]
                    out["admonitions"].append(dtype)
                else:
                    out["code"].append(f"<<<fence {info}>>>")
                stack.append((indent, marker, is_directive))
                continue
            stack.append((indent, marker, info.startswith("{")))
            continue

        in_code = bool(stack) and not all(s[2] for s in stack)
        if in_code:
            out["code"].append(line)
            continue

        if OPTION.match(line) and stack:
            out["options"].append(line.strip())
            continue

        if HEADING.match(line):
            out["headings"].append(line.strip())
            continue

        prose_lines.append(line)

    prose = "\n".join(prose_lines)
    # Roles and citations before inline code, so their backticked payload is not also
    # counted as an inline code span.
    out["citations"] = CITE.findall(prose)
    out["roles"] = [f"{{{a}}}`{b}`" for a, b in ROLE.findall(prose)]
    stripped = ROLE.sub(" ", CITE.sub(" ", prose))
    out["inline_code"] = INLINE_CODE.findall(stripped)
    out["links"] = MDLINK.findall(prose) + AUTOLINK.findall(prose)
    out["math"] = [m.strip() for m in MATH.findall(prose)]
    # Numbers from prose with inline code and math removed: a flag like --p-n-lambda1
    # or $\lambda_1$ is not a quantitative claim.
    numeric_src = MATH.sub(" ", INLINE_CODE.sub(" ", stripped))
    numeric_src = MDLINK.sub(" ", numeric_src)
    out["numbers"] = NUMBER.findall(numeric_src)
    return out


def git_show(rev, relpath):
    try:
        return subprocess.run(["git", "show", f"{rev}:{relpath}"], cwd=ROOT,
                              capture_output=True, text=True, check=True).stdout
    except subprocess.CalledProcessError:
        return None


def diff_section(name, before, after, ordered):
    if ordered:
        if before == after:
            return []
        msgs = []
        cb, ca = collections.Counter(before), collections.Counter(after)
        for item in (cb - ca).elements():
            msgs.append(f"    - removed: {item[:150]}")
        for item in (ca - cb).elements():
            msgs.append(f"    + added:   {item[:150]}")
        if not msgs:
            msgs.append(f"    ~ reordered ({len(before)} items, same multiset)")
        return msgs
    cb, ca = collections.Counter(before), collections.Counter(after)
    if cb == ca:
        return []
    return ([f"    - removed: {i[:150]}" for i in (cb - ca).elements()] +
            [f"    + added:   {i[:150]}" for i in (ca - cb).elements()])


ORDERED = {"headings", "code", "directive_heads", "options", "numbers", "admonitions"}

# Headings are allowed to change -- the chattiest prose in this book is in its headings,
# and only 8 links in 443 target an anchor. But a renamed heading silently breaks any
# link that did target it, so heading changes are REVIEW (reported, not fatal) and are
# backed by anchor_report(), which is fatal.
REVIEW = {"headings"}


def slug(heading):
    """Approximate the anchor MyST generates for a heading."""
    text = re.sub(r"^\s*#+\s*", "", heading)
    text = CITE.sub("", text)
    text = re.sub(r"[`*_]", "", text)
    text = re.sub(r"\$[^$]*\$", lambda m: m.group(0), text)
    text = text.lower()
    text = re.sub(r"[^a-z0-9\s-]", "", text)
    return re.sub(r"\s+", "-", text.strip())


def anchor_report(files):
    """Every [text](target#anchor) must resolve to a heading that exists.

    Link targets are resolved against the file on disk, not against the checked set --
    checking one file must not report every link it makes into the rest of the book as
    broken.
    """
    heads = {}

    def headings_of(path):
        if path not in heads:
            heads[path] = ({slug(h) for h in parse(path.read_text())["headings"]}
                           if path.is_file() else None)
        return heads[path]

    broken = []
    for f in files:
        for target in MDLINK.findall(f.read_text()):
            if "#" not in target or target.startswith(("http://", "https://", "mailto:")):
                continue
            path_part, _, anchor = target.partition("#")
            dest = (f.parent / path_part).resolve() if path_part else f.resolve()
            found = headings_of(dest)
            if found is None:
                broken.append((f, target, "target file does not exist"))
            elif anchor not in found:
                broken.append((f, target, "no heading generates this anchor"))
    return broken


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="*", help="files to check (default: docs/chapters/**.md)")
    ap.add_argument("--rev", default="HEAD", help="git revision to compare against")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if args.paths:
        files = [pathlib.Path(p).resolve() for p in args.paths]
    else:
        files = sorted((ROOT / "docs" / "chapters").rglob("*.md"))

    total_findings = 0
    checked = 0
    for f in files:
        rel = f.relative_to(ROOT).as_posix()
        old = git_show(args.rev, rel)
        if old is None:
            print(f"  NEW  {rel} (not in {args.rev}, skipped)")
            continue
        new = f.read_text()
        checked += 1
        if old == new:
            if args.verbose:
                print(f"  ==   {rel} (unchanged)")
            continue
        b, a = parse(old), parse(new)
        hard, review = [], []
        for key in b:
            msgs = diff_section(key, b[key], a[key], key in ORDERED)
            if not msgs:
                continue
            (review if key in REVIEW else hard).append(f"  [{key}]")
            (review if key in REVIEW else hard).extend(msgs)
        if hard:
            total_findings += 1
            print(f"\nFAIL {rel}")
            print("\n".join(hard))
            if review:
                print("\n".join(review))
        elif review:
            print(f"\nREVIEW {rel} (invariants held; headings changed)")
            if args.verbose:
                print("\n".join(review))
        else:
            print(f"  OK   {rel} (prose changed, invariants held)")

    broken = anchor_report(files)
    if broken:
        print(f"\nFAIL broken anchor links ({len(broken)}):")
        for f, target, why in broken:
            print(f"    {f.relative_to(ROOT).as_posix()} -> {target}  ({why})")

    print(f"\n{checked} file(s) compared against {args.rev}; "
          f"{total_findings} with moved invariants; {len(broken)} broken anchors.")
    return 1 if (total_findings or broken) else 0


if __name__ == "__main__":
    sys.exit(main())
