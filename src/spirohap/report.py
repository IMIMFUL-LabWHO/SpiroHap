from __future__ import annotations

import html
import math
from pathlib import Path
from typing import Any


def _table(
    rows: list[dict[str, Any]], columns: list[tuple[str, str]], limit: int | None = None
) -> str:
    shown = rows if limit is None else rows[:limit]
    head = "".join(f"<th>{html.escape(label)}</th>" for _, label in columns)
    body = []
    for row in shown:
        cells = "".join(
            f"<td>{html.escape(str(row.get(key, '')))}</td>" for key, _ in columns
        )
        body.append(f"<tr>{cells}</tr>")
    if not body:
        body.append(f'<tr><td colspan="{len(columns)}" class="muted">No rows</td></tr>')
    suffix = ""
    if limit is not None and len(rows) > limit:
        suffix = f'<p class="muted">Showing {limit} of {len(rows)} rows; see the TSV file for all rows.</p>'
    return f'<div class="table-wrap"><table><thead><tr>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>{suffix}'


def _network_svg(catalog: list[dict[str, Any]], edges: list[dict[str, Any]]) -> str:
    if not catalog:
        return '<p class="muted">No haplotypes passed filters.</p>'
    width, height = 760, 460
    center_x, center_y = width / 2, height / 2
    radius = min(width, height) * 0.34
    positions: dict[str, tuple[float, float]] = {}
    for index, row in enumerate(catalog):
        angle = -math.pi / 2 + 2 * math.pi * index / max(1, len(catalog))
        if len(catalog) == 1:
            positions[row["haplotype_id"]] = (center_x, center_y)
        else:
            positions[row["haplotype_id"]] = (
                center_x + radius * math.cos(angle),
                center_y + radius * math.sin(angle),
            )
    parts = [
        f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Minimum spanning haplotype network">'
    ]
    for edge in edges:
        x1, y1 = positions[edge["source"]]
        x2, y2 = positions[edge["target"]]
        parts.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" class="edge"/>'
        )
        mx, my = (x1 + x2) / 2, (y1 + y2) / 2
        parts.append(
            f'<text x="{mx:.1f}" y="{my - 5:.1f}" class="edge-label">{edge["distance"]}</text>'
        )
    for row in catalog:
        node = row["haplotype_id"]
        x, y = positions[node]
        human = int(row["human_samples"])
        tick = int(row["tick_samples"])
        control = int(row["negative_control_samples"])
        if control:
            color = "#d97706"
        elif human and tick:
            color = "#7c3aed"
        elif human:
            color = "#dc2626"
        elif tick:
            color = "#16803c"
        else:
            color = "#64748b"
        node_radius = 12 + min(18, 3 * math.sqrt(max(1, int(row["sample_count"]))))
        parts.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{node_radius:.1f}" fill="{color}" class="node"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{y + 4:.1f}" class="node-label">{html.escape(node)}</text>'
        )
    parts.append("</svg>")
    parts.append(
        '<p class="legend"><span class="human">human only</span><span class="tick">tick only</span>'
        '<span class="both">human + tick</span><span class="control">negative control</span></p>'
    )
    return "".join(parts)


def write_report(
    path: Path,
    summary: dict[str, Any],
    qc_rows: list[dict[str, Any]],
    catalog_rows: list[dict[str, Any]],
    sample_haplotype_rows: list[dict[str, Any]],
    association_rows: list[dict[str, Any]],
    network_edges: list[dict[str, Any]],
) -> None:
    cards = "".join(
        f'<div class="card"><div class="number">{html.escape(str(value))}</div><div>{html.escape(label)}</div></div>'
        for label, value in (
            ("samples", summary["sample_count"]),
            ("passing samples", summary["passing_samples"]),
            ("unique 23S haplotypes", summary["haplotype_count"]),
            ("retained reads", summary["retained_reads"]),
        )
    )
    document = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>SpiroHap report</title>
<style>
:root {{ color-scheme: light; --ink:#172033; --muted:#64748b; --line:#dbe3ee; --bg:#f5f7fb; }}
* {{ box-sizing:border-box; }}
body {{ margin:0; font:14px/1.5 system-ui,-apple-system,Segoe UI,sans-serif; color:var(--ink); background:var(--bg); }}
main {{ max-width:1180px; margin:0 auto; padding:32px 22px 70px; }}
h1 {{ margin:0 0 4px; font-size:30px; }} h2 {{ margin-top:34px; font-size:21px; }} h3 {{ margin-top:24px; }}
.subtitle,.muted {{ color:var(--muted); }}
.notice {{ background:#fff7ed; border:1px solid #fed7aa; border-left:5px solid #f97316; padding:14px 16px; border-radius:8px; margin:22px 0; }}
.cards {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr)); gap:12px; margin:22px 0; }}
.card {{ background:white; border:1px solid var(--line); border-radius:10px; padding:16px; box-shadow:0 1px 2px #0000000a; }}
.number {{ font-size:25px; font-weight:700; }}
.table-wrap {{ overflow:auto; background:white; border:1px solid var(--line); border-radius:9px; }}
table {{ border-collapse:collapse; width:100%; min-width:700px; }}
th,td {{ padding:9px 11px; border-bottom:1px solid var(--line); text-align:left; vertical-align:top; }}
th {{ background:#edf2f7; position:sticky; top:0; }} tr:last-child td {{ border-bottom:0; }}
svg {{ width:100%; height:auto; background:white; border:1px solid var(--line); border-radius:9px; }}
.edge {{ stroke:#94a3b8; stroke-width:2; }} .edge-label {{ fill:#475569; font-size:12px; text-anchor:middle; paint-order:stroke; stroke:white; stroke-width:4px; }}
.node {{ stroke:white; stroke-width:3; }} .node-label {{ fill:white; text-anchor:middle; font-size:10px; font-weight:700; }}
.legend span {{ margin-right:16px; }} .legend span:before {{ content:""; display:inline-block; width:10px; height:10px; border-radius:50%; margin-right:5px; }}
.human:before {{ background:#dc2626; }} .tick:before {{ background:#16803c; }} .both:before {{ background:#7c3aed; }} .control:before {{ background:#d97706; }}
code {{ background:#e8edf4; padding:2px 5px; border-radius:4px; }}
</style>
</head>
<body><main>
<h1>SpiroHap 23S amplicon report</h1>
<p class="subtitle">Run created {html.escape(summary["created_at"])} · reference: {html.escape(summary["reference_name"])}</p>
<div class="notice"><strong>Interpretation boundary.</strong> These are reconstructed alleles of one 23S amplicon, not whole-genome strains. A bacterium can carry more than one rRNA operon, and a one-locus association cannot establish a human-infecting lineage. Treat low-frequency, single-event, unreplicated, or control-associated calls as provisional.</div>
<div class="cards">{cards}</div>
<h2>Haplotype catalogue</h2>
{_table(catalog_rows, [("haplotype_id", "ID"), ("length", "bp"), ("sample_count", "samples"), ("human_samples", "human"), ("tick_samples", "tick"), ("negative_control_samples", "negative controls"), ("total_reads", "assigned reads"), ("nearest_haplotype", "nearest"), ("nearest_distance", "distance"), ("closest_reference", "closest reference"), ("reference_distance", "ref distance")])}
<h2>Minimum-spanning network</h2>
<p class="muted">Edges minimize total edit distance; labels are nucleotide edit distances. This is a descriptive haplotype network, not a phylogenetic tree.</p>
{_network_svg(catalog_rows, network_edges)}
<h2>Sample-level calls</h2>
{_table(sample_haplotype_rows, [("sample_id", "sample"), ("source", "source"), ("haplotype_id", "haplotype"), ("assigned_reads", "reads"), ("fraction", "fraction"), ("forward_reads", "+ strand"), ("reverse_reads", "− strand"), ("replicate_status", "replicate"), ("flags", "flags")], limit=250)}
<h2>Sample QC</h2>
{_table(qc_rows, [("sample_id", "sample"), ("source", "source"), ("status", "status"), ("raw_reads", "raw"), ("passed_mean_q", "Q pass"), ("primer_matched", "primer matched"), ("passed_length", "length pass"), ("retained_reads", "retained"), ("ambiguous_reads", "ambiguous"), ("warnings", "warnings")], limit=250)}
<h2>Exploratory human–tick comparison</h2>
<p class="muted">Fisher tests use only passing non-control samples labelled exactly <code>human</code> or <code>tick</code>. With three humans, estimates are extremely uncertain; q-values are descriptive screening statistics, not evidence of causality.</p>
{_table(association_rows, [("haplotype_id", "haplotype"), ("human_positive", "human +"), ("human_total", "human n"), ("tick_positive", "tick +"), ("tick_total", "tick n"), ("odds_ratio_corrected", "odds ratio"), ("fisher_p", "p"), ("bh_q", "BH q")])}
<h2>Files to inspect next</h2>
<p><code>haplotypes.fasta</code>, <code>sample_haplotypes.tsv</code>, <code>read_assignments.tsv</code>, <code>haplotype_variants.tsv</code>, <code>candidate_events.tsv</code>, <code>haplotype_network.tsv</code>, and <code>run_manifest.json</code>.</p>
</main></body></html>"""
    path.write_text(document, encoding="utf-8")
