"""Render the README images of show=True: python -I docs/images/render_show.py docs/images"""
import io, re, sys, contextlib
from html import escape
import tcrkit

COL = {'31': '#d73a49', '32': '#22863a', '33': '#b08800', '34': '#0366d6',
       '35': '#8f3fbf', '36': '#1b9aaa', '90': '#8a8f98'}
BOLD = ' font-weight="bold"'
FG, BG = '#24292e', '#ffffff'

def svg(text, path, cw=8.4, lh=18):
    lines = text.rstrip('\n').split('\n')
    plain = [re.sub(r'\033\[\d*m', '', l) for l in lines]
    W = int(max(map(len, plain)) * cw + 32); H = int(len(lines) * lh + 24)
    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
           f'<rect width="100%" height="100%" rx="6" fill="{BG}" stroke="#d0d7de"/>',
           '<g font-family="ui-monospace,SFMono-Regular,Menlo,Consolas,monospace" font-size="14" xml:space="preserve">']
    for i, line in enumerate(lines):
        y = 20 + (i + 1) * lh - 4
        spans, color = [], None
        for part in re.split(r'(\033\[\d*m)', line):
            m = re.fullmatch(r'\033\[(\d*)m', part)
            if m:
                color = COL.get(m.group(1)); continue
            if part:
                b = BOLD if color else ''
                spans.append(f'<tspan fill="{color or FG}"{b}>{escape(part)}</tspan>')
        out.append(f'<text x="16" y="{y}">{"".join(spans)}</text>')
    out.append('</g></svg>')
    open(path, 'w').write('\n'.join(out) + '\n')

def capture(f):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        f()
    return buf.getvalue()

d = sys.argv[1]
beta = tcrkit.stitch('TRBV19', 'TRBJ2-7', 'CASSIRSSYEQYF', chain='beta', squeeze=True)['sequence_aa']
svg(capture(lambda: tcrkit.anarci(beta, allow={'A', 'B'}, domains=True, show={'width': 80})), f'{d}/show_tcr.svg')
a2 = tcrkit.mhc_sequence('HLA-A*02:01', squeeze=True)['sequence']
svg(capture(lambda: tcrkit.annotate_mhc(a2, domains=True, show={'width': 80})), f'{d}/show_mhc.svg')
