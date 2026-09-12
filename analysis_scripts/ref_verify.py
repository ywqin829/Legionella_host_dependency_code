"""Verify every PMID to be cited via PubMed esummary (title/journal/year/authors four-field)."""
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

def esummary(pmids):
    out = {}
    for i in range(0, len(pmids), 50):
        u = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=%s'
             % ','.join(pmids[i:i+50]))
        d = json.loads(get(u))['result']
        for pid in pmids[i:i+50]:
            r = d.get(pid, {})
            out[pid] = {'title': r.get('title', ''), 'journal': r.get('fulljournalname', ''),
                        'year': r.get('pubdate', '')[:4],
                        'authors': ', '.join(a['name'] for a in r.get('authors', [])[:4])}
        time.sleep(1)
    return out

# all PMIDs: existing manuscript refs + new mechanism candidates
PMIDS = [
 '35759189','31540829','15117975','24244168','28835546','28842705','19942850','30323948',  # existing
 '33469029','21734656','23821544','32440549','37076474','35204756','31123136','18281405',  # new mechanisms
 '31994202','16824952','27120161','20163564','31669129','29112293',  # + SidE, LCV, cross-pathogen candidates
]
s = esummary(PMIDIS := PMIDS)
for pid in PMIDS:
    r = s.get(pid, {})
    ok = 'VERIFIED' if r.get('title') else 'FAIL'
    print(f"{ok} PMID {pid} | {r.get('year','')} | {r.get('journal','')[:32]} | {r.get('title','')[:75]}")
    print(f"      authors: {r.get('authors','')[:60]}")
