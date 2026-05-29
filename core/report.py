"""
Self-contained HTML tuning report.

Summarises a TIR-Tuner session: the original vs final sequence (with a
colour-coded diff), the list of applied edits, the TIR before/after, the
designed mutagenesis primers, and an embedded ΔTIR heatmap image (PNG base64)
if one is supplied. No external assets — opens in any browser.
"""

import base64
import difflib
import html
from typing import List, Optional

_DIFF_COLORS = {
    'equal': None, 'replace': '#fde68a', 'insert': '#bbf7d0', 'delete': '#fecaca',
}


def _diff_html(a: str, b: str) -> str:
    sm = difflib.SequenceMatcher(None, a, b, autojunk=False)
    ao, bo, tags = [], [], []
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == 'equal':
            for k in range(i2 - i1):
                ao.append(a[i1 + k]); bo.append(b[j1 + k]); tags.append('equal')
        elif tag == 'replace':
            for k in range(max(i2 - i1, j2 - j1)):
                ao.append(a[i1 + k] if k < i2 - i1 else '-')
                bo.append(b[j1 + k] if k < j2 - j1 else '-')
                tags.append('replace')
        elif tag == 'delete':
            for k in range(i2 - i1):
                ao.append(a[i1 + k]); bo.append('-'); tags.append('delete')
        elif tag == 'insert':
            for k in range(j2 - j1):
                ao.append('-'); bo.append(b[j1 + k]); tags.append('insert')

    def cell(ch, tag):
        color = _DIFF_COLORS[tag]
        disp = ch if ch != '-' else '·'
        return (f'<span style="background:{color}">{disp}</span>'
                if color else disp)

    width = 80
    blocks = []
    for s in range(0, len(tags), width):
        e = min(s + width, len(tags))
        o = ''.join(cell(ao[k], tags[k]) for k in range(s, e))
        nw = ''.join(cell(bo[k], tags[k]) for k in range(s, e))
        blocks.append(f'<div>old {s + 1:>5}: {o}</div>'
                      f'<div>new {s + 1:>5}: {nw}</div><div>&nbsp;</div>')
    return ('<pre style="font-family:Consolas,monospace;font-size:12px;'
            'white-space:pre;line-height:1.4">' + ''.join(blocks) + '</pre>')


def build_report(*, title: str, original: str, final: str,
                 applied: List[dict], primers: Optional[List[dict]] = None,
                 baseline_tir: Optional[float] = None,
                 final_tir: Optional[float] = None,
                 heatmap_png: Optional[bytes] = None,
                 meta: Optional[dict] = None) -> str:
    esc = html.escape
    parts = [f"""<!doctype html><html><head><meta charset="utf-8">
<title>{esc(title)}</title>
<style>
 body{{font-family:'Segoe UI',system-ui,sans-serif;color:#1f2a37;margin:30px;background:#f7f9fc}}
 h1{{font-size:22px}} h2{{font-size:15px;color:#2563eb;margin-top:26px}}
 .card{{background:#fff;border:1px solid #dbe3ee;border-radius:10px;padding:16px 18px;margin:12px 0}}
 table{{border-collapse:collapse;width:100%;font-size:13px}}
 th,td{{border:1px solid #e5ebf3;padding:6px 9px;text-align:left}}
 th{{background:#f1f5fb;color:#475569}}
 .metric{{display:inline-block;margin-right:26px}}
 .metric b{{font-size:20px;color:#0f172a}}
 code{{background:#f1f5f9;padding:1px 5px;border-radius:4px}}
 .muted{{color:#64748b}}
</style></head><body>"""]
    parts.append(f'<h1>{esc(title)}</h1>')

    # Metrics
    mh = []
    if baseline_tir is not None and final_tir is not None:
        fold = (final_tir / baseline_tir) if baseline_tir else float('inf')
        foldtxt = '∞' if fold == float('inf') else f'{fold:.2f}×'
        mh.append(f'<span class="metric">Wild-type TIR<br><b>{baseline_tir:.1f}</b></span>')
        mh.append(f'<span class="metric">Final TIR<br><b>{final_tir:.1f}</b></span>')
        mh.append(f'<span class="metric">Fold change<br><b>{foldtxt}</b></span>')
    mh.append(f'<span class="metric">Length<br><b>{len(original)} → {len(final)} nt</b></span>')
    parts.append(f'<div class="card">{"".join(mh)}</div>')

    if meta:
        rows = ''.join(f'<tr><th>{esc(str(k))}</th><td>{esc(str(v))}</td></tr>'
                       for k, v in meta.items())
        parts.append(f'<div class="card"><h2>Run</h2><table>{rows}</table></div>')

    # Applied edits
    if applied:
        rows = ''.join(
            f'<tr><td>{i+1}</td><td>{esc(a.get("desc",""))}</td></tr>'
            for i, a in enumerate(applied))
        parts.append('<div class="card"><h2>Applied edits</h2>'
                     f'<table><tr><th>#</th><th>Edit</th></tr>{rows}</table></div>')

    # Diff
    parts.append('<div class="card"><h2>Sequence comparison</h2>'
                 '<p class="muted">amber = substitution · '
                 'green = insertion · red = deletion · gaps shown as ·</p>'
                 + _diff_html(original, final) + '</div>')

    # Primers
    if primers:
        cols = list(primers[0].keys())
        head = ''.join(f'<th>{esc(c)}</th>' for c in cols)

        def _cell(col, val):
            v = esc(str(val))
            if 'primer' in col.lower():
                return f'<td><code>{v}</code></td>'
            return f'<td>{v}</td>'

        body = ''.join(
            '<tr>' + ''.join(_cell(c, p.get(c, '')) for c in cols) + '</tr>'
            for p in primers)
        parts.append('<div class="card"><h2>Mutagenesis primers</h2>'
                     f'<table><tr>{head}</tr>{body}</table></div>')

    # Heatmap
    if heatmap_png:
        b64 = base64.b64encode(heatmap_png).decode('ascii')
        parts.append('<div class="card"><h2>ΔTIR heatmap</h2>'
                     f'<img style="max-width:100%" src="data:image/png;base64,{b64}"></div>')

    parts.append('<p class="muted">Generated by TIRex · scored with OSTIR</p>')
    parts.append('</body></html>')
    return ''.join(parts)
