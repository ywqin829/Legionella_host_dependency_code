"""PubMed E-utilities helper: search + verify candidate references.
Every PMID is esummary-verified (title/journal/year/authors). Iron rule: no unverified citation."""
import time, urllib.parse, urllib.request, json, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def get(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    for _ in range(3):
        try:
            return urllib.request.urlopen(req, timeout=30).read()
        except Exception:
            time.sleep(1.5)
    return b''

def esearch(query, n=5):
    u = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmax=%d&retmode=json&term=%s'
         % (n, urllib.parse.quote(query)))
    return json.loads(get(u))['esearchresult']['idlist']

def esummary(pmids):
    out = {}
    for i in range(0, len(pmids), 50):
        chunk = pmids[i:i+50]
        u = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=%s'
             % ','.join(chunk))
        d = json.loads(get(u))['result']
        for pid in chunk:
            r = d.get(pid, {})
            out[pid] = {
                'title': r.get('title', ''),
                'journal': r.get('fulljournalname', ''),
                'year': r.get('pubdate', '')[:4],
                'authors': ', '.join(a['name'] for a in r.get('authors', [])[:3]),
            }
        time.sleep(1)
    return out

QUERIES = {
    'SidM/DrrA AMPylation Rab1': 'DrrA SidM AMPylation Rab1 Legionella',
    'SidD deAMPylase': 'SidD deAMPylase Rab1 Legionella',
    'LepB Rab1 GAP': 'LepB Rab1 GAP Legionella',
    'AnkX phosphocholination Rab1': 'AnkX phosphocholination Rab1 Legionella',
    'Lem3 dephosphocholinase': 'Lem3 phosphocholine Rab Legionella',
    'SetA Rab1 glycosylation': 'SetA Legionella glycosylation Rab1',
    'RalF ARF1 GEF': 'RalF Legionella ARF1 GEF',
    'SidE family ubiquitination': 'SidE phosphoribosyl ubiquitination Legionella',
    'SidJ metaeffector': 'SidJ Legionella metaeffector',
    'Lgt translation inhibition': 'Lgt1 Legionella translation inhibition eEF1A',
    'MavC MvcA UBE2N': 'MavC MvcA UBE2N Legionella',
    'LCV biogenesis Rab1': 'Legionella containing vacuole biogenesis Rab1',
    'host dependency CRISPR screen': 'genome-wide CRISPR screen intracellular bacteria host factor',
}

for topic, q in QUERIES.items():
    ids = esearch(q, 4)
    if not ids:
        print(f"== {topic}: NO HITS")
        continue
    summ = esummary(ids)
    print(f"== {topic}: {q}")
    for pid in ids:
        r = summ.get(pid, {})
        print(f"   PMID {pid} | {r.get('year','')} | {r.get('journal','')[:30]} | {r.get('title','')[:70]} | {r.get('authors','')[:40]}")
    time.sleep(1)
