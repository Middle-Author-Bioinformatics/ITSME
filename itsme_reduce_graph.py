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
                self.records.append((node, label, targets, ''.join(sequence)))
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



SEQUENCE_FIELDS = ['locus', 'source_bin', 'graph_file', 'source', 'length_bp',
                   'locus_type', 'mean_depth', 'taxonomy_status',
                   'consensus_taxonomy', 'node_path']


def write_table(path, rows, fields):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def write_graph_index(out, reports, missing):
    """Human-readable sequence identities; taxonomy is copied, never inferred."""
    esc = lambda value: html.escape(str(value if value is not None else ''))
    sections = []
    for item in reports:
        rows = item['sequences']
        body = ''.join('<tr>' + ''.join('<td>' + esc(r.get(k, '')) + '</td>'
                       for k in SEQUENCE_FIELDS if k != 'graph_file') + '</tr>' for r in rows)
        headings = ''.join('<th>' + esc(k.replace('_', ' ')) + '</th>'
                           for k in SEQUENCE_FIELDS if k != 'graph_file')
        sections.append('<section><h2>' + esc(item['source_bin']) + '</h2><p>Graph: <a href="' +
                        esc(item['file']) + '">' + esc(item['file']) + '</a> | <a href="' +
                        esc(item['sequence_table']) + '">Sequence table</a> | ' + str(len(rows)) +
                        ' reported sequence(s); ' + str(item['retained_nodes']) +
                        ' retained graph nodes.</p><div class="scroll"><table><thead><tr>' +
                        headings + '</tr></thead><tbody>' + body + '</tbody></table></div></section>')
    omitted = ('<h2>Unmapped sequences</h2><ul>' + ''.join(
        '<li>' + esc(r['locus']) + ': ' + esc(r['reason']) + '</li>' for r in missing) + '</ul>') if missing else ''
    page = """<!doctype html><html lang="en"><meta charset="utf-8">
<title>ITSME graph-to-sequence index</title><style>
body{font:15px system-ui,sans-serif;margin:30px;color:#182330}h1{font-size:26px}
section{margin:28px 0}.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%}
th,td{padding:10px;border:1px solid #ccd4dd;text-align:left;vertical-align:top}
th{background:#eaf0f6}td{overflow-wrap:anywhere;min-width:100px}a{color:#0759a4}
</style><h1>Which sequences does each graph represent?</h1>
<p>Sequence IDs match the master summary. Source bin numbers and locus numbers are independent.
Taxonomy and status below are copied from the summary; blank fields mean unavailable.</p>
<p>Graphs contain shared sequence pieces, not necessarily one separate component per reported sequence.
A shared node can belong to several reconstructions. Node IDs are local to each graph.
Paths use source orientation, which may differ from the final FASTA orientation.
Additional neighbor nodes, if requested, are context—not classified sequences.</p>
""" + ''.join(sections) + omitted + '</html>'
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
        reports.append(dict(graph=label, file=filename, source=str(graph.path),
                            source_bin=sequence_rows[0]['source_bin'],
                            sequence_table=sequence_table, sequences=sequence_rows,
                            original_nodes=len(graph.nodes), locus_nodes=len(selected),
                            retained_nodes=len(keep)))
    with (out / 'locus_nodes.tsv').open('w', newline='') as f:
        fields = ['locus', 'source', 'graph', 'source_graph', 'node_path'] + [
            k for k in SEQUENCE_FIELDS if k not in ('locus', 'source', 'node_path')]
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore'); w.writeheader(); w.writerows(assignments)
    write_graph_index(out, reports, missing)
    report = dict(status='partial' if missing else 'complete', summary=str(summary),
                  requested_loci=len(loci), mapped_loci=len(assignments), neighbors=a.neighbors,
                  omitted_loci=missing, graphs=reports)
    (out / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f"{report['status']}: mapped {len(assignments)}/{len(loci)} loci; output: {out}")
    for item in reports:
        print(f"  {item['file']}: {item['retained_nodes']}/{item['original_nodes']} nodes")
        for row in item['sequences']:
            print(f"    {row['locus']} | {row['length_bp']} bp | {row['consensus_taxonomy'] or 'taxonomy unavailable'}")
    print('Open graph_index.html to see the graph-to-sequence mapping.')
    if missing:
        print('WARNING: incomplete subgraphs; see report.json for omitted loci.', file=sys.stderr)


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError, csv.Error) as exc:
        print('ERROR: ' + str(exc), file=sys.stderr)
        sys.exit(1)
