#!/usr/bin/env python3
"""Reduce ITSME assembly graphs to nodes belonging to master-summary loci.

Python >=3.9, standard library only. Does not run assembly or BLAST.
Usage:
  python3 itsme_reduce_graph.py --run-dir itsme_SAMPLE
  python3 itsme_reduce_graph.py --run-dir itsme_SAMPLE --neighbors 1
  python3 itsme_reduce_graph.py --run-dir itsme_SAMPLE --graph graph.fastg \
      --contig-paths contigs.paths --output-dir reduced

Default output: RUN/final/summary_graphs (must be new or empty).
Produces separate reduced graphs for pooled and taxon-bin assemblies, plus
locus_nodes.tsv and report.json. Also writes graph_index.html and one
*.sequences.tsv per graph, linking graph files to sequence IDs and taxonomy.
Graph node IDs are local to each assembly. No extra flags are required.
Defaults to FASTG; --format gfa selects GFA1. Explicit --graph detects format
from content, irrespective of extension, and is only used for pooled loci.

Keeps whole nodes, both FASTG orientations, and all original links between
retained nodes (an induced subgraph), not only traversed path links. Neighbor
expansion is undirected. Does not trim node ends or prove biological phase.
Collapsed alternative paths are excluded: only the representative is retained.
Native contig paths may have semicolon-separated disconnected pieces; all
pieces are retained without inventing links. Output paths list source orientation,
which may be reversed relative to the final oriented FASTA.

Missing/ambiguous mappings abort before writing by default. --allow-missing
writes a clearly marked partial result and lists omitted loci. Empty summaries
produce an empty report, not an error. Existing nonempty output is never replaced.
"""
import argparse
import csv
import json
import html
import re
import sys
from collections import defaultdict
from pathlib import Path


def table(path):
    with path.open(newline='', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f, delimiter=',' if path.suffix == '.csv' else '\t')
        return list(reader), reader.fieldnames or []


def tokens(value):
    pieces = []
    for piece in value.strip().split(';'):
        nodes = []
        for token in piece.strip().split(','):
            token = token.strip()
            if not re.fullmatch(r'[^\s,;]+[+-]', token):
                raise ValueError('Invalid oriented graph path: ' + value)
            nodes.append(token)
        pieces.append(nodes)
    return pieces


def paths_file(path):
    result = {}; current = None
    with path.open() as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('NODE_'):
                current = line
                result[current] = []
            elif current:
                # SPAdes splits a discontinuous path over lines ending in ';'.
                for piece in line.rstrip(';').split(';'):
                    if piece.strip():
                        result[current].extend(tokens(piece))
            else:
                raise ValueError('Unrecognized SPAdes path header in ' + str(path))
    return result


def fastg_id(label):
    m = re.fullmatch(r'EDGE_(\d+)_length_\d+_cov_[^\s:\',;]+(\')?', label)
    if not m:
        raise ValueError('Unsupported FASTG node label: ' + label)
    return m[1]


class Graph:
    def __init__(self, path):
        self.path = path
        self.nodes = set(); self.adj = defaultdict(set); self.records = []
        self.lengths = {}
        with path.open() as f:
            first = next((s for s in f if s.strip()), '')
        self.format = 'fastg' if first.startswith('>') else 'gfa'
        if self.format == 'fastg':
            header = None; sequence = []
            def save():
                if header is None:
                    return
                parts = header.rstrip(';').split(':', 1)
                label = parts[0]; node = fastg_id(label)
                targets = parts[1].split(',') if len(parts) == 2 else []
                self.nodes.add(node)
                for target in targets:
                    other = fastg_id(target)
                    self.adj[node].add(other); self.adj[other].add(node)
                seq = ''.join(sequence)
                self.records.append((node, label, targets, seq))
                self.lengths[node] = len(seq)
            with path.open() as f:
                for line in f:
                    if line.startswith('>'):
                        save(); header = line[1:].strip(); sequence = []
                    else:
                        sequence.append(line.strip())
                save()
        else:
            with path.open() as f:
                for raw in f:
                    fields = raw.rstrip('\n').split('\t')
                    if not raw.strip():
                        continue
                    if fields[0] not in {'H', 'S', 'L', 'P', 'J'}:
                        raise ValueError('Unsupported GFA record type: ' + fields[0])
                    if fields[0] == 'H' and any(x.startswith('VN:Z:2') for x in fields):
                        raise ValueError('GFA2 is not supported')
                    self.records.append(fields)
                    if fields[0] == 'S':
                        self.nodes.add(fields[1])
                        seq = fields[2] if len(fields) > 2 else '*'
                        if seq != '*':
                            self.lengths[fields[1]] = len(seq)
                        else:
                            ln = next((x[5:] for x in fields[3:] if x.startswith('LN:i:')), '')
                            self.lengths[fields[1]] = int(ln) if ln.isdigit() else 0
                    elif fields[0] in {'L', 'J'}:
                        self.adj[fields[1]].add(fields[3]); self.adj[fields[3]].add(fields[1])
        if not self.nodes:
            raise ValueError('No graph nodes found in ' + str(path))
        unknown = set(self.adj) - self.nodes
        if unknown:
            raise ValueError('Graph links reference missing nodes in ' + str(path))

    def write(self, dest, keep):
        with dest.open('w') as f:
            if self.format == 'fastg':
                for node, label, targets, sequence in self.records:
                    if node not in keep:
                        continue
                    targets = [t for t in targets if fastg_id(t) in keep]
                    f.write('>' + label + (':' + ','.join(targets) if targets else '') + ';\n')
                    for i in range(0, len(sequence), 80):
                        f.write(sequence[i:i+80] + '\n')
            else:
                for row in self.records:
                    # Drop P records: they describe the unfiltered assembly, not
                    # necessarily the selected loci. Exact locus paths are in TSV.
                    if (row[0] == 'H' or
                        row[0] == 'S' and row[1] in keep or
                        row[0] in {'L', 'J'} and row[1] in keep and row[3] in keep):
                        f.write('\t'.join(row) + '\n')



def graph_view_data(graph, keep, graph_assignments):
    """Return only the information needed by the browser-side graph viewer."""
    node_loci = defaultdict(list)
    for row in graph_assignments:
        for node in row['nodes']:
            if node in keep:
                node_loci[node].append(row['locus'])

    nodes = []
    for node in sorted(keep, key=lambda x: (not str(x).isdigit(),
                                            int(x) if str(x).isdigit() else str(x))):
        nodes.append({
            'id': str(node),
            'length': graph.lengths.get(node, 0),
            'loci': sorted(set(node_loci.get(node, []))),
            'context_only': node not in node_loci,
        })

    edges = []
    seen = set()
    for a in keep:
        for b in graph.adj.get(a, ()):
            if b not in keep:
                continue
            key = tuple(sorted((str(a), str(b))))
            if key in seen:
                continue
            seen.add(key)
            edges.append([str(a), str(b)])
    return {'nodes': nodes, 'edges': edges}




SEQUENCE_FIELDS = ['locus', 'source_bin', 'graph_file', 'source', 'length_bp',
                   'locus_type', 'mean_depth', 'taxonomy_status',
                   'consensus_taxonomy', 'node_path']


def write_table(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def write_graph_index(out, reports, missing, graph_views):
    """Interactive graph/sequence index. Taxonomy is copied, never inferred."""
    esc = lambda value: html.escape(str(value if value is not None else ''))
    sections = []

    for idx, item in enumerate(reports):
        rows = item['sequences']
        body = ''.join('<tr>' + ''.join('<td>' + esc(r.get(k, '')) + '</td>'
                       for k in SEQUENCE_FIELDS if k != 'graph_file') + '</tr>' for r in rows)
        headings = ''.join('<th>' + esc(k.replace('_', ' ')) + '</th>'
                           for k in SEQUENCE_FIELDS if k != 'graph_file')

        display_bin = item['source_bin']
        if display_bin.startswith('bin_'):
            display_bin = 'Bin ' + display_bin[4:].replace('_', ' ')
        elif display_bin == 'pooled':
            display_bin = 'Pooled assembly'

        locus_chips = ''.join(
            '<button class="locus-chip" data-view="' + str(idx) + '" data-locus="' +
            esc(r.get('locus', '')) + '">' + esc(r.get('locus', '')) +
            (' · ' + esc(r.get('consensus_taxonomy', '')) if r.get('consensus_taxonomy') else '') +
            '</button>' for r in rows
        )

        sections.append(
            '<section class="graph-section"><h2>' + esc(display_bin) + '</h2>'
            '<p><strong>Graph:</strong> <a href="' + esc(item['file']) + '">' + esc(item['file']) +
            '</a> &nbsp;·&nbsp; <a href="' + esc(item['sequence_table']) + '">Sequence table</a>'
            ' &nbsp;·&nbsp; ' + str(len(rows)) + ' reconstructed sequence(s)'
            ' &nbsp;·&nbsp; ' + str(item['retained_nodes']) + ' retained graph nodes.</p>'
            '<div class="locus-chips">' + locus_chips + '</div>'
            '<div class="viewer-wrap">'
            '<div class="viewer-toolbar">'
            '<button type="button" onclick="resetView(' + str(idx) + ')">Reset view</button>'
            '<button type="button" onclick="fitView(' + str(idx) + ')">Fit</button>'
            '<span class="legend"><span class="swatch locus"></span>reconstructed locus node '
            '<span class="swatch shared"></span>shared by loci '
            '<span class="swatch context"></span>neighbor/context node</span>'
            '</div>'
            '<canvas id="graph-' + str(idx) + '" class="graph-canvas" width="1200" height="650"></canvas>'
            '<div id="node-info-' + str(idx) + '" class="node-info">'
            'Click a node to see its ID, length, and reconstructed sequence assignment.</div>'
            '</div>'
            '<details><summary>Sequence/taxonomy table</summary>'
            '<div class="scroll"><table><thead><tr>' + headings + '</tr></thead><tbody>' +
            body + '</tbody></table></div></details></section>'
        )

    omitted = ('<h2>Unmapped sequences</h2><ul>' + ''.join(
        '<li>' + esc(r['locus']) + ': ' + esc(r['reason']) + '</li>' for r in missing) +
        '</ul>') if missing else ''

    view_json = json.dumps([graph_views.get(item['graph'], {'nodes': [], 'edges': []})
                            for item in reports]).replace('</', '<\\/')

    page = r"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<title>ITSME interactive graph-to-sequence viewer</title>
<style>
:root{--ink:#182330;--muted:#64748b;--line:#ccd4dd;--panel:#f7f9fc;--accent:#2563eb;
--shared:#7c3aed;--context:#94a3b8}
*{box-sizing:border-box} body{font:15px system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
margin:30px;color:var(--ink);background:#fff} h1{font-size:27px;margin-bottom:8px}
h2{font-size:21px;margin:0 0 6px}.subtle{color:var(--muted)} section{margin:34px 0 44px}
a{color:#0759a4}.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;margin-top:12px}
th,td{padding:9px;border:1px solid var(--line);text-align:left;vertical-align:top}
th{background:#eaf0f6}td{overflow-wrap:anywhere;min-width:100px}
.viewer-wrap{border:1px solid var(--line);border-radius:10px;overflow:hidden;background:var(--panel)}
.viewer-toolbar{padding:9px 11px;border-bottom:1px solid var(--line);display:flex;gap:8px;
align-items:center;flex-wrap:wrap;background:white}.viewer-toolbar button,.locus-chip{border:1px solid #b8c3d1;
background:white;border-radius:6px;padding:6px 9px;cursor:pointer}.viewer-toolbar button:hover,.locus-chip:hover{
background:#edf4ff}.locus-chip.active{outline:2px solid var(--accent);background:#eef4ff}
.legend{margin-left:auto;color:var(--muted);font-size:13px}.swatch{display:inline-block;width:10px;height:10px;
border-radius:50%;margin:0 4px 0 10px;background:var(--accent)}.swatch.shared{background:var(--shared)}
.swatch.context{background:var(--context)}.graph-canvas{display:block;width:100%;height:650px;background:#fff;
cursor:grab}.graph-canvas.dragging{cursor:grabbing}.node-info{min-height:42px;border-top:1px solid var(--line);
padding:10px 12px;background:white}.locus-chips{display:flex;gap:6px;flex-wrap:wrap;margin:10px 0}
details{margin-top:12px}summary{cursor:pointer;font-weight:600}code{background:#eef2f7;padding:1px 4px;border-radius:4px}
@media(max-width:800px){body{margin:16px}.graph-canvas{height:500px}.legend{margin-left:0;width:100%}}
</style></head><body>
<h1>ITSME graph-to-sequence viewer</h1>
<p class="subtle">Each section is the reduced assembly graph for one source assembly. Nodes used by reconstructed
sequences are highlighted; additional neighbor nodes requested with <code>--neighbors</code> are shown as context.
Click a sequence label to isolate its path nodes, or click a node for details.</p>
<p>Sequence IDs match the master summary. Source bin numbers and locus numbers are independent.
Taxonomy/status are copied from the summary. Shared nodes can belong to several reconstructions, and node IDs are
local to each graph.</p>
""" + ''.join(sections) + omitted + r"""
<script>
const GRAPH_DATA = """ + view_json + r""";
const STATE = [];

function makeState(i){
  const canvas=document.getElementById('graph-'+i), ctx=canvas.getContext('2d');
  const data=GRAPH_DATA[i], byId=new Map(data.nodes.map(n=>[String(n.id),n]));
  const nodes=data.nodes.map((n,k)=>{
    const a=2*Math.PI*k/Math.max(1,data.nodes.length);
    const ring=180+Math.min(220,Math.sqrt(data.nodes.length)*12);
    return {...n,x:canvas.width/2+Math.cos(a)*ring,y:canvas.height/2+Math.sin(a)*ring,
            vx:0,vy:0};
  });
  const nodeMap=new Map(nodes.map(n=>[String(n.id),n]));
  const edges=data.edges.map(e=>[nodeMap.get(String(e[0])),nodeMap.get(String(e[1]))]).filter(e=>e[0]&&e[1]);
  const st={canvas,ctx,nodes,nodeMap,edges,scale:1,ox:0,oy:0,drag:null,pan:false,lastX:0,lastY:0,
            selectedNode:null,selectedLocus:null,frames:0};
  STATE[i]=st;
  bindEvents(i);
  simulate(i, Math.min(220, 70+nodes.length));
  fitView(i);
}

function nodeRadius(n){
  const len=Math.max(1,Number(n.length)||1);
  return Math.max(5,Math.min(15,4+Math.log10(len+1)*2.2));
}
function nodeColor(n){
  if(n.context_only) return '#94a3b8';
  if((n.loci||[]).length>1) return '#7c3aed';
  return '#2563eb';
}
function screenToWorld(st,x,y){return {x:(x-st.ox)/st.scale,y:(y-st.oy)/st.scale};}
function draw(i){
  const st=STATE[i], c=st.canvas, ctx=st.ctx, dpr=window.devicePixelRatio||1;
  const rect=c.getBoundingClientRect(), targetW=Math.round(rect.width*dpr), targetH=Math.round(rect.height*dpr);
  if(c.width!==targetW||c.height!==targetH){c.width=targetW;c.height=targetH;fitView(i);return;}
  ctx.clearRect(0,0,c.width,c.height); ctx.save(); ctx.translate(st.ox,st.oy); ctx.scale(st.scale,st.scale);
  ctx.lineWidth=1/st.scale; ctx.strokeStyle='#c7d0db';
  for(const [a,b] of st.edges){ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.stroke();}
  for(const n of st.nodes){
    const r=nodeRadius(n);
    const dim=st.selectedLocus && !(n.loci||[]).includes(st.selectedLocus);
    ctx.globalAlpha=dim?0.12:1; ctx.beginPath();ctx.arc(n.x,n.y,r,0,Math.PI*2);
    ctx.fillStyle=nodeColor(n);ctx.fill();
    if(st.selectedNode===n){ctx.lineWidth=3/st.scale;ctx.strokeStyle='#111827';ctx.stroke();}
    if(!dim && st.scale>0.65){
      ctx.globalAlpha=0.9;ctx.fillStyle='#111827';ctx.font=(11/st.scale)+'px system-ui';
      ctx.fillText(String(n.id),n.x+r+3/st.scale,n.y+3/st.scale);
    }
  }
  ctx.restore();ctx.globalAlpha=1;
}
function simulate(i,steps){
  const st=STATE[i], ns=st.nodes, es=st.edges;
  if(ns.length>450){draw(i);return;} // keep very large reduced graphs responsive
  for(let iter=0;iter<steps;iter++){
    for(let a=0;a<ns.length;a++){
      for(let b=a+1;b<ns.length;b++){
        const A=ns[a],B=ns[b],dx=B.x-A.x,dy=B.y-A.y,d2=dx*dx+dy*dy+25,d=Math.sqrt(d2);
        const f=1800/d2,fx=f*dx/d,fy=f*dy/d;A.vx-=fx;A.vy-=fy;B.vx+=fx;B.vy+=fy;
      }
    }
    for(const [A,B] of es){
      const dx=B.x-A.x,dy=B.y-A.y,d=Math.sqrt(dx*dx+dy*dy)||1, target=65;
      const f=(d-target)*0.018,fx=f*dx/d,fy=f*dy/d;A.vx+=fx;A.vy+=fy;B.vx-=fx;B.vy-=fy;
    }
    let cx=0,cy=0;for(const n of ns){cx+=n.x;cy+=n.y}cx/=Math.max(1,ns.length);cy/=Math.max(1,ns.length);
    for(const n of ns){n.vx+=(600-cx)*0.0008;n.vy+=(325-cy)*0.0008;n.vx*=0.82;n.vy*=0.82;n.x+=n.vx;n.y+=n.vy;}
  }
  draw(i);
}
function fitView(i){
  const st=STATE[i]; if(!st||!st.nodes.length)return;
  const xs=st.nodes.map(n=>n.x),ys=st.nodes.map(n=>n.y),minX=Math.min(...xs),maxX=Math.max(...xs),
        minY=Math.min(...ys),maxY=Math.max(...ys),w=Math.max(50,maxX-minX),h=Math.max(50,maxY-minY);
  st.scale=Math.min((st.canvas.width-80)/w,(st.canvas.height-80)/h,2.2);
  st.ox=st.canvas.width/2-(minX+maxX)/2*st.scale;st.oy=st.canvas.height/2-(minY+maxY)/2*st.scale;draw(i);
}
function resetView(i){const st=STATE[i];st.selectedNode=null;st.selectedLocus=null;
  document.querySelectorAll('.locus-chip[data-view="'+i+'"]').forEach(x=>x.classList.remove('active'));
  document.getElementById('node-info-'+i).textContent='Click a node to see its ID, length, and reconstructed sequence assignment.';
  fitView(i);}
function hitNode(st,x,y){
  const p=screenToWorld(st,x,y);let best=null,bd=Infinity;
  for(const n of st.nodes){const d=Math.hypot(n.x-p.x,n.y-p.y);if(d<nodeRadius(n)+5/st.scale&&d<bd){best=n;bd=d;}}
  return best;
}
function bindEvents(i){
  const st=STATE[i],c=st.canvas;
  c.addEventListener('mousedown',e=>{const r=c.getBoundingClientRect(),x=(e.clientX-r.left)*(c.width/r.width),
    y=(e.clientY-r.top)*(c.height/r.height),n=hitNode(st,x,y);st.lastX=x;st.lastY=y;
    if(n){st.drag=n;st.selectedNode=n;}else st.pan=true;c.classList.add('dragging');draw(i);});
  window.addEventListener('mousemove',e=>{if(!st.drag&&!st.pan)return;const r=c.getBoundingClientRect(),
    x=(e.clientX-r.left)*(c.width/r.width),y=(e.clientY-r.top)*(c.height/r.height);
    if(st.drag){const p=screenToWorld(st,x,y);st.drag.x=p.x;st.drag.y=p.y;}
    else{st.ox+=x-st.lastX;st.oy+=y-st.lastY;}st.lastX=x;st.lastY=y;draw(i);});
  window.addEventListener('mouseup',()=>{if(st.drag){showNode(i,st.drag)}st.drag=null;st.pan=false;c.classList.remove('dragging');});
  c.addEventListener('wheel',e=>{e.preventDefault();const r=c.getBoundingClientRect(),
    x=(e.clientX-r.left)*(c.width/r.width),y=(e.clientY-r.top)*(c.height/r.height),before=screenToWorld(st,x,y),
    factor=e.deltaY<0?1.12:0.89;st.scale=Math.max(0.08,Math.min(6,st.scale*factor));
    st.ox=x-before.x*st.scale;st.oy=y-before.y*st.scale;draw(i);},{passive:false});
  c.addEventListener('click',e=>{const r=c.getBoundingClientRect(),x=(e.clientX-r.left)*(c.width/r.width),
    y=(e.clientY-r.top)*(c.height/r.height),n=hitNode(st,x,y);if(n){st.selectedNode=n;showNode(i,n);draw(i);}});
}
function showNode(i,n){
  const loci=(n.loci||[]);document.getElementById('node-info-'+i).innerHTML=
    '<strong>Node '+escapeHtml(n.id)+'</strong> &nbsp;·&nbsp; '+Number(n.length||0).toLocaleString()+
    ' bp &nbsp;·&nbsp; '+(n.context_only?'context/neighbor node':
    'used by: '+loci.map(escapeHtml).join(', '));
}
function escapeHtml(s){return String(s).replace(/[&<>"']/g,m=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[m]));}
document.querySelectorAll('.locus-chip').forEach(btn=>btn.addEventListener('click',()=>{
  const i=Number(btn.dataset.view),locus=btn.dataset.locus,st=STATE[i];
  const next=st.selectedLocus===locus?null:locus;st.selectedLocus=next;
  document.querySelectorAll('.locus-chip[data-view="'+i+'"]').forEach(x=>x.classList.toggle('active',next&&x.dataset.locus===next));
  document.getElementById('node-info-'+i).innerHTML=next?
    '<strong>'+escapeHtml(next)+'</strong>: highlighted nodes are those contributing to this reconstructed sequence.':
    'Click a node to see its ID, length, and reconstructed sequence assignment.';
  draw(i);
}));
for(let i=0;i<GRAPH_DATA.length;i++)makeState(i);
window.addEventListener('resize',()=>STATE.forEach((_,i)=>draw(i)));
</script></body></html>"""
    (out / 'graph_index.html').write_text(page, encoding='utf-8')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--run-dir', required=True, type=Path)
    p.add_argument('--summary', type=Path, help='Override final/master_summary.tsv or root master_summary.csv')
    p.add_argument('--output-dir', type=Path)
    p.add_argument('--graph', type=Path, help='Explicit pooled graph; format is detected from content')
    p.add_argument('--contig-paths', type=Path, help='Explicit pooled SPAdes contigs.paths')
    p.add_argument('--format', choices=['fastg', 'gfa'], default='fastg')
    p.add_argument('--neighbors', type=int, default=0)
    p.add_argument('--allow-missing', action='store_true')
    a = p.parse_args()
    if a.neighbors < 0:
        p.error('--neighbors must be >= 0')
    root = a.run_dir.resolve(); final = root / 'final'
    summary = a.summary or final / 'master_summary.tsv'
    if not a.summary and not summary.exists():
        summary = root / 'master_summary.csv'
    rows, fields = table(summary)
    if 'contig' not in fields:
        raise ValueError('Summary must have a contig column')
    metadata = {}
    for row in rows:
        locus = row.get('contig', '').strip()
        if locus in metadata and metadata[locus] != row:
            raise ValueError('Conflicting summary rows for ' + locus)
        metadata[locus] = row
    loci = list(dict.fromkeys(r['contig'].strip() for r in rows if r.get('contig', '').strip()))
    out = (a.output_dir or final / 'summary_graphs').resolve()
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError('Output must be a new or empty directory: ' + str(out))
    mappings = defaultdict(list)
    source_map = final / 'locus_source_map.tsv'
    if source_map.exists():
        for r in table(source_map)[0]:
            if r.get('path_status') != 'COLLAPSED_SAME_TAXONOMY':
                mappings[r['locus']].append(r)
    bins = []
    for b in sorted((root / 'taxon_bins').glob('*')):
        if b.is_dir():
            prefix = 'ITSME_' + b.name.upper().replace('-', '_') + '_'
            bins.append((prefix, b))
    graphs = {}; native_paths = {}; assignments = []; missing = []
    for locus in loci:
        try:
            mapped = mappings.get(locus, [])
            if len(mapped) > 1:
                raise ValueError('Multiple representative source paths')
            source = mapped[0]['source_path'] if mapped else locus
            matches = [(prefix, b) for prefix, b in bins if source.startswith(prefix)]
            if len(matches) > 1:
                raise ValueError('Ambiguous taxon-bin prefix')
            label = 'pooled'; assembly = root / 'assembly' / 'spades'; original = source
            if matches:
                prefix, b = matches[0]; label = 'bin_' + b.name
                original = source[len(prefix):]; assembly = b / 'spades'
                metrics = b / 'bridge_rescue.tsv'
                if metrics.exists():
                    accepted = [int(r['round']) for r in table(metrics)[0]
                                if r.get('decision') == 'accepted_reassembled']
                    if accepted:
                        assembly = b / 'bridge_rescue' / ('round_%02d' % max(accepted)) / 'spades'
            elif source.startswith('ITSME_') and not source.startswith(('ITSME_PATH', 'ITSME_LOCUS')):
                raise ValueError('No matching taxon-bin directory for ' + source)
            graph_path = a.graph if label == 'pooled' and a.graph else None
            if graph_path is None:
                names = (['assembly_graph.fastg'] if a.format == 'fastg' else
                         ['assembly_graph_with_scaffolds.gfa', 'assembly_graph.gfa'])
                candidates = [assembly / name for name in names]
                if label == 'pooled':
                    candidates += [final / name for name in names]
                graph_path = next((q for q in candidates if q.is_file()), None)
                if graph_path is None:
                    raise ValueError('Missing source graph for ' + label)
            graph_path = graph_path.resolve()
            if label not in graphs:
                graphs[label] = Graph(graph_path)
            graph = graphs[label]
            if mapped:
                pieces = tokens(mapped[0]['node_path'])
            else:
                path_file = a.contig_paths if label == 'pooled' and a.contig_paths else assembly / 'contigs.paths'
                if path_file not in native_paths:
                    native_paths[path_file] = paths_file(path_file)
                pieces = native_paths[path_file].get(original)
                if not pieces:
                    raise ValueError('No exact contig path for ' + original)
            nodes = {token[:-1] for piece in pieces for token in piece}
            absent = nodes - graph.nodes
            if absent:
                raise ValueError('Path nodes absent from graph: ' + ','.join(sorted(absent)[:10]))
            details = {key: metadata[locus].get(key, '') for key in
                       ('length_bp', 'locus_type', 'mean_depth', 'taxonomy_status', 'consensus_taxonomy')}
            assignments.append(dict(locus=locus, source=source, graph=label,
                                    source_bin=matches[0][1].name if matches else 'pooled',
                                    graph_file=label + '.summary.' + graph.format, **details,
                                    source_graph=str(graph_path),
                                    node_path=';'.join(','.join(piece) for piece in pieces),
                                    nodes=nodes))
        except (ValueError, OSError, KeyError) as exc:
            missing.append(dict(locus=locus, reason=str(exc)))
    if missing and not a.allow_missing:
        for entry in missing:
            print('UNMAPPED ' + entry['locus'] + ': ' + entry['reason'], file=sys.stderr)
        raise ValueError('Incomplete mapping; no outputs written. Fix inputs or use --allow-missing.')
    out.mkdir(parents=True, exist_ok=True)
    reports = []
    graph_views = {}
    for label, graph in graphs.items():
        selected = set().union(*(r['nodes'] for r in assignments if r['graph'] == label))
        if not selected:
            continue
        keep = set(selected); frontier = set(selected)
        for _ in range(a.neighbors):
            frontier = set().union(*(graph.adj[n] for n in frontier)) - keep
            keep.update(frontier)
        filename = label + '.summary.' + graph.format
        graph.write(out / filename, keep)
        sequence_rows = [{k: row.get(k, '') for k in SEQUENCE_FIELDS}
                         for row in assignments if row['graph'] == label]
        sequence_table = label + '.sequences.tsv'
        write_table(out / sequence_table, sequence_rows, SEQUENCE_FIELDS)
        graph_assignments = [r for r in assignments if r['graph'] == label]
        graph_views[label] = graph_view_data(graph, keep, graph_assignments)
        reports.append(dict(graph=label, file=filename, source=str(graph.path),
                            source_bin=sequence_rows[0]['source_bin'],
                            sequence_table=sequence_table, sequences=sequence_rows,
                            original_nodes=len(graph.nodes), locus_nodes=len(selected),
                            retained_nodes=len(keep)))
    with (out / 'locus_nodes.tsv').open('w', newline='') as f:
        fields = ['locus', 'source', 'graph', 'source_graph', 'node_path'] + [
            k for k in SEQUENCE_FIELDS if k not in ('locus', 'source', 'node_path')]
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore'); w.writeheader(); w.writerows(assignments)
    write_graph_index(out, reports, missing, graph_views)
    report = dict(status='partial' if missing else 'complete', summary=str(summary),
                  requested_loci=len(loci), mapped_loci=len(assignments), neighbors=a.neighbors,
                  omitted_loci=missing, graphs=reports)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f"{report['status']}: mapped {len(assignments)}/{len(loci)} loci; output: {out}")
    for item in reports:
        print(f"  {item['file']}: {item['retained_nodes']}/{item['original_nodes']} nodes")
        for row in item['sequences']:
            print(f"    {row['locus']} | {row['length_bp']} bp | {row['consensus_taxonomy'] or 'taxonomy unavailable'}")
    print('Open graph_index.html to interactively view the reduced graphs and graph-to-sequence mapping.')
    if missing:
        print('WARNING: incomplete subgraphs; see report.json for omitted loci.', file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        sys.exit(1)
