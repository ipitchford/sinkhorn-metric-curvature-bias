# SPDX-License-Identifier: MIT
"""Build refs.bib (every \\bibitem of manuscript.tex) from the saved registry records in refs/records/:
Crossref REST records (api.crossref.org/works/<DOI>), the DataCite record (CSL-JSON) of the cited Zenodo release and the arXiv API
records, all retrieved 4 November 2026. Three web-only works (DLMF; two author-posted preprints) are entered from the
manuscript's own reference list with their URLs. Usage: python3 refs/make_refs.py (from the package root)."""
import json, re, sys, xml.etree.ElementTree as ET
from pathlib import Path
R = Path(__file__).resolve().parents[1]; REC = R / 'refs' / 'records'
tex = (R / 'manuscript.tex').read_text()
items = re.findall(r'\\bibitem\{([^}]+)\}(.*)', tex)
def esc(s): return re.sub(r'(?<!\\)([&%#_])', r'\\\1', ' '.join(str(s).split()))
def entry(typ, key, fields, src):
    return f'% Source: {src}\n@{typ}{{{key},\n' + ',\n'.join(f'  {k} = {{{esc(v) if k not in ("doi", "url", "eprint") else v}}}' for k, v in fields if v) + '\n}\n\n'
ns = {'a': 'http://www.w3.org/2005/Atom'}
arx = {}
for e in [x for f in ('arxiv_records.xml', 'arxiv_saidi_2610.36488.xml') for x in ET.parse(REC / f).getroot().findall('a:entry', ns)]:
    vid = e.find('a:id', ns).text.split('/abs/')[1]; arx[re.sub(r'v\d+$', '', vid)] = e
out = ['% refs.bib: every reference of manuscript.tex, built by refs/make_refs.py from refs/records/ (retrieved 4 November 2026).\n\n']
for key, rest in items:
    d = re.search(r'https://doi\.org/(\S+?)\}', rest); a = re.search(r'https://arxiv\.org/abs/(\S+?)\}', rest)
    if d and d.group(1).startswith('10.5281/zenodo.'):
        at = json.loads((REC / f'datacite_{d.group(1).split("/")[1]}.json').read_text())  # CSL-JSON from doi.org (DataCite)
        if at['DOI'].lower() != d.group(1).lower(): sys.exit('DOI mismatch ' + key)
        out.append(entry('misc', key, [('author', 'Anonymous'), ('title', at['title']), ('howpublished', 'Evidence Press research release; Zenodo'),
            ('version', at.get('version')), ('year', at['issued']['date-parts'][0][0]), ('doi', d.group(1))], f'DataCite (CSL-JSON), doi:{at["DOI"]}'))
    elif d:
        m = json.loads((REC / f'crossref_{d.group(1).replace("/", "_")}.json').read_text()); m = m.get('message', m)
        if m['DOI'].lower() != d.group(1).lower(): sys.exit('DOI mismatch ' + key)
        names = ' and '.join(f"{x.get('family', '')}, {x.get('given', '')}".strip(', ') for x in m.get('author') or []) or None
        yr = (m.get('published-print') or m.get('published') or m.get('issued'))['date-parts'][0][0]
        typ = {'journal-article': 'article', 'book-chapter': 'incollection', 'proceedings-article': 'inproceedings', 'book': 'book', 'monograph': 'book'}.get(m['type'], 'misc')
        f = [('author', names), ('title', re.sub(r'<[^>]+>', '', (m.get('title') or [''])[0]))]
        cont = (m.get('container-title') or [''])[0]
        if cont and typ != 'book': f.append(('booktitle' if typ in ('incollection', 'inproceedings') else 'journal', cont))
        for k, s in (('volume', 'volume'), ('number', 'issue'), ('pages', 'page'), ('publisher', 'publisher')): f.append((k, m.get(s)))
        if m.get('article-number') and not m.get('page'): f.append(('pages', m['article-number']))
        f += [('year', str(yr)), ('doi', m['DOI'])]
        if a: f.append(('eprint', a.group(1)))
        out.append(entry(typ, key, f, f'Crossref, doi:{m["DOI"]}'))
    elif a:
        e = arx[a.group(1)]; vid = e.find('a:id', ns).text.split('/abs/')[1]
        out.append(entry('misc', key, [('author', ' and '.join(x.find('a:name', ns).text for x in e.findall('a:author', ns))), ('title', e.find('a:title', ns).text),
            ('year', e.find('a:published', ns).text[:4]), ('eprint', re.sub(r'v\d+$', '', vid)), ('archiveprefix', 'arXiv'), ('url', 'https://arxiv.org/abs/' + re.sub(r'v\d+$', '', vid))], f'arXiv API record, {vid}'))
    else:
        urls = re.findall(r'\\url\{([^}]+)\}', rest); text = re.sub(r'\\url\{[^}]+\}', '', rest)
        out.append(entry('misc', key, [('note', re.sub(r'\\emph\{([^}]*)\}', r'\1', text).replace('\\S', 'Section ').strip(' .')), ('url', urls[0])], 'manuscript reference list (web-only work)'))
(R / 'refs.bib').write_text(''.join(out)); print(len(items), 'references written')
