# Prose style for this book

The reference page is
[`docs/chapters/04_highdim_atacama/07_motus_shotgun.md`](docs/chapters/04_highdim_atacama/07_motus_shotgun.md).
When this document and that page disagree, the page wins.

The goal is formal technical documentation: precise, unhurried, and free of the tics that
make writing read as generated. Formality is not impersonality. The book gives the reader
instructions, and it should keep giving them in the imperative and the second person. An
earlier revision of the reference page removed every "you" from it and the result read
*more* synthetic, not less, because the operational advice turned into agentless
obligation ("$\gamma$ should be reported") where the neighbouring chapters say "Report
$\gamma$ with any network you derive from it."

## Part 1 — Hard invariants

**A style pass changes sentences. It changes nothing else.** This book's value is that its
numbers come from committed tables; a rewrite that quietly moves one is worse than no
rewrite. The following must survive byte-identical unless a factual error is found:

| Invariant | Why |
|---|---|
| Every numeral, unit and percentage | They are generated from `analysis/results/tables/` and asserted by the figure scripts |
| Everything inside fenced code blocks | Commands are run verbatim by readers; outputs are transcripts |
| Inline code spans naming flags, files, types, columns | `--p-gamma`, `mclr.qza`, `FeatureTable[Frequency]` |
| ``{cite}`key` `` keys | Resolved against `docs/references.bib` at build time |
| Link targets, `{ref}`, `{numref}`, `{doc}`, and `:name:` anchors | Cross-page links break silently in the rendered HTML |
| Directive types and options | `:file:`, `:width:`, `:header-rows:`, `:delim:` |
| Table structure — same rows, same columns, same order | |

**Headings are the exception, and they are where much of the chattiest prose lives**
("What each chapter adds", "What you should have now", "Before you start"). Rewrite them.
Only 8 of the book's 443 links target a heading anchor, so the cost is low — but a renamed
heading silently breaks any link that did target it. If you rename a heading, grep the
book for its anchor (`#what-you-cannot-conclude`) and update every inbound link in the
same change.

`analysis/scripts/check_prose_invariants.py` enforces all of the above against `git HEAD`,
and resolves every anchor link in the book. Run it before you consider a file done.

Admonition **types** (`{note}`, `{tip}`, `{important}`, `{warning}`) stay as they are. Do
not promote a note to a warning for emphasis, and do not demote one to reduce drama. The
prose inside them is rewritten like any other prose.

If you believe a number or a claim is *wrong*, do not fix it silently and do not rewrite
around it. Leave it and report it.

## Part 2 — What to remove

Each rule below is followed by a real before/after from the reference page's revision.

### 2.1 The document must not narrate itself

The reader can see the page. Sentences about what the page is doing, why it exists, or
what it will do later are the single loudest generated-text signal.

> ~~That difference turns out to decide which parameter values work, and it is the reason
> this page exists rather than a sentence saying "shotgun tables also work".~~
> **The lower depth determines which parameter values yield a usable model.**

> ~~That result is relevant to the classification analysis below and is discussed under
> [Limitations](#limitations-and-interpretation).~~
> **The bearing of this conclusion on the classification analysis is given under
> [Limitations](#limitations-and-interpretation).**

> ~~Profiling requires 13 GiB of reads, and is therefore summarised rather than reproduced
> here.~~
> **Profiling requires 13 GiB of sequencing reads and the 2.9 GiB mOTUs reference
> database. The corresponding pipeline stages are `analysis/slurm/32`–`35`; the steps are
> given here in summary.**

Cross-references to *other* pages are not self-narration — keep them. What goes is
narration about *this* page's own structure and editorial choices.

### 2.2 No bare claims of importance

Asserting that something matters is weaker than showing it. If the next sentence
demonstrates the point, delete the announcement.

> ~~The consequence is substantive.~~ Under prevalence filtering, the extended BIC selects
> the empty graph…
> **The two filters select different models.** Under prevalence filtering, …

Also drop: "Importantly," "It is worth noting that," "Crucially," "It is important to
understand."

### 2.3 No X-not-Y antithesis as decoration

The frame is fine when the contrast is the actual content ("Rank features by abundance,
not by prevalence" — the whole point is the comparison). It is padding when Y is a
strawman nobody proposed.

### 2.4 No metaphor where a number is available

> ~~Three latent dimensions absorb 55% of the network.~~
> **Three latent dimensions remove 264 of the 481 edges (55%).**

> ~~the latent block explains away very little~~
> **the latent block removes few edges**

### 2.5 No counted openers that create ordinal back-references

> ~~Two things differ: a feature is a marker-gene-defined species…, and the number of reads
> is far smaller. The second difference determines which parameter values yield a usable
> model.~~
> **The differences are the feature definition and the depth per feature: a feature here
> is a marker-gene-defined species…, and each count rests on far fewer reads. The lower
> depth determines which parameter values yield a usable model.**

A count is fine when it is stable and the list follows immediately — "Three artifacts are
distributed with this book:" is good prose, not a tic. The problem is only the count that
forces the reader to hold an ordinal across paragraphs.

### 2.6 Bold is for defined terms and list labels, not for stress

Two bolded phrases within four lines for emphasis is a tic. `**Repeated measures.**` as a
bullet label is correct. `**36 runs from 18 infants**` mid-sentence is not.

### 2.7 Do not name the book's own tiers as if they were conditions

> ~~differs from the $\gamma = 0.3$ used elsewhere in this tier~~
> **differs from the $\gamma = 0.3$ used for the Atacama data**

Name the dataset, the chapter, or the model. "This tier" means nothing to a reader who
arrived from a search engine.

## Part 3 — What to keep and what to restore

### 3.1 Imperative voice for anything operational

This is the rule most likely to be over-applied in the wrong direction. Agentless deontic
passive is the failure mode, not the target.

> ~~$\gamma$ should be reported alongside any network derived from it.~~
> **Report $\gamma$ with any network you derive from it.**

> ~~this should be stated in any methods description that relies on it~~
> **State this in any methods description that relies on it.**

> ~~If `--p-lambda2-*` is left unset, a five-point default path is substituted and the
> single fit becomes a model-selection run.~~
> **Leave `--p-lambda2-*` unset and the solver substitutes a five-point default path,
> which turns the single fit into a model-selection run.**

Keep an actor in the sentence: the solver substitutes, the action requires, mOTUs emits,
`add-taxa` reconstructs. Passive voice is correct where the actor is genuinely irrelevant
or unknown — "192 runs, of which half are amplicon" needs no agent.

### 3.2 Second person, for the reader's actions

Keep "you" for what the reader does: what they run, choose, need installed, will see.
Do not use it for what the *data* does.

### 3.3 State the requirement rather than withholding it

> ~~This requires a property that not every taxonomy possesses.~~
> **This requires distinct leaf labels, which not all taxonomies provide.**

### 3.4 Structural conventions of the reference page

- Every command block is followed by its verbatim output block (`Saved … to: …`).
- A command with many flags gets an `**Explanation:**` bullet list covering the
  non-obvious ones.
- Figures carry a `:name:` and a caption that states what the figure shows and the numbers
  in it.
- Quantitative claims name their source when it is not the adjacent table.

### 3.5 Punctuation

Do not use `---` as a horizontal rule between sections. Headings already separate
sections, and a transition immediately after a heading is a MyST error
("Document or section may not begin with a transition") that fails the build under
`--warningiserror`.

Em-dashes are standard formal punctuation — keep them where they set off a genuine aside,
and do not manufacture more. Prefer a full stop to a semicolon joining two independent
clauses in an elliptical antithesis:

> ~~The covariance estimate is unaffected; the cross-validated regression is not.~~
> **The covariance estimate is unaffected. The cross-validated regression is affected.**

## Part 4 — Before you finish a file

1. Read every heading aloud against the original. Identical?
2. Diff the code blocks. Byte-identical?
3. Search your output for: "it is worth", "importantly", "crucially", "this section",
   "this page", "this chapter will", "we will see", "as we", "let's", "simply", "just",
   "note that", "in order to", "leverage", "delve", "robust" (unless statistical),
   "seamless", "powerful", "comprehensive".
4. Count the numerals. Same set, same order?
5. Does any sentence tell the reader what the document is doing rather than what the
   software or the data does?
