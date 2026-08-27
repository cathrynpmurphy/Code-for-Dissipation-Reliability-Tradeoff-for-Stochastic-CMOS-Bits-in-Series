# FigureSources

TikZ/pgfplots sources for every figure in the CMOS chain paper, plus the thin
Python "prep" scripts that turn raw exported data into space-free `.dat`
tables. `make` builds all five. The intent is that this directory plus the raw
data reproduces every figure, so it can be deposited on Zenodo as the
figure-generation record.

    raw data (../PostReviewFigures/Figure N Data/*.txt)
        -> prep_figN.py     deterministic, no styling decisions
        -> figN_*.dat       plain columns, no spaces in filenames
        -> figN.tex         standalone + \input{figstyle}
        -> figN.pdf         pure vector, exact width, tight MediaBox

## The style contract

`figstyle.tex` is the single source of truth. Every `figN.tex` does
`\input{figstyle}` and uses its styles. Do not restyle locally; change
`figstyle.tex` and rebuild all five.

**Scale must be 1.** This is the whole point. Author each figure at exactly the
width it occupies on the page, so `\includegraphics` never rescales it and a
9pt label in one figure is a 9pt label in every other. The five figures
currently ship at five different scale factors -- 0.770, 0.797, 0.807, 0.948,
0.843 -- which is why nothing matches today. A 9pt label lands at 6.9pt in
Fig. 1 and 8.5pt in Fig. 4.

Measured page geometry, `revtex4-2 [reprint,aps,pre]`:

| quantity | pt | in |
|---|---|---|
| `\columnwidth` | 246.0 | 3.404 |
| `\textwidth` (figure*) | 510.0 | 7.057 |

Class type sizes: body 10, `\small` 9, `\footnotesize` 8, `\scriptsize` 7.

**Type and weight**, from `figstyle.tex`: axis labels 9pt (`\figaxislabel`),
tick labels 8pt (`\figticklabel`), in-figure annotation 8pt (`\figannot`),
panel letters 9pt bold (`figpanellabel`). Axis/frame 0.6pt, data curves 1.0pt,
emphasised curves 1.4pt. Ticks inward on all four sides, minor ticks on, no
grid. Colormap viridis. Palette and node colours are named in `figstyle.tex`.

**Fonts: Computer Modern**, i.e. the LaTeX default. The paper is `revtex4-2`
with no font package, so it typesets in CM. Never load `newtx`, `times`, or
`mathptmx` in a figure.

**Pinning the width.** Put `\figcanvas{<w>}{<h>}` as the first line inside the
`tikzpicture`; it pins the bounding box so the emitted PDF is exactly that
size regardless of content.

**Verifying the width.** `pdfinfo` reports big points (72/in), not LaTeX pt
(72.27/in). A correct column figure reads **245.08 bp**; a correct full-width
figure reads **508.09 bp**. Check those numbers, not 246 and 510.

**Pure vector.** `pdfimages -list figN.pdf` must list no image objects.

## Target sizes per figure

| fig | env | width | canvas | include with |
|---|---|---|---|---|
| 1 | `figure*` | 510pt | 510 x 128pt | `[width=\textwidth]` |
| 2 | `figure` | see note | see note | `[width=...]` |
| 3 | `figure` | 246pt | 246 x ~250pt | `[width=\columnwidth]` |
| 4 | `figure` | 246pt | 246 x ~190pt | `[width=\columnwidth]` |
| 5 | `figure*` | 510pt | 510 x ~200pt | `[width=\textwidth]` |

Note on Fig. 2: its current aspect ratio is 204 x 427pt. Blown up to the full
246pt column it would stand **515pt tall**, two thirds of a page. Do not do
that. Either keep it at its present 162.6pt width (height 340pt) and include it
with `[width=162.6pt]` -- scale is still 1, which is all the contract requires
-- or re-proportion the layout while rebuilding so it is less severely
vertical. Re-proportioning is preferred; see `PROMPT_fig2.md`.

A narrower-than-column figure is fine. The contract is scale = 1 plus the
shared type spec, not one universal width.
