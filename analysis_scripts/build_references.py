"""Build the verified reference table (all PMIDs esummary-verified: title/journal/year/authors)."""
import time, urllib.parse, urllib.request, json, sys, io
import pandas as pd
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def get(url):
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    return urllib.request.urlopen(req, timeout=30).read()

def esummary(pmids):
    u = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=%s' % ','.join(pmids)
    return json.loads(get(u))['result']

# final reference PMIDs (topic-tagged for manuscript placement)
REFS = [
    ('Wibawa2022_methods', '35759189'),
    ('Jeng2019_crispr', '31540829'),
    ('Kagan2004_rab1sec22b', '15117975'),
    ('Machner2006_rab', '16824952'),
    ('Suh2010_drrA', '19942850'),
    ('Du2021_drrA_ampyl', '33469029'),
    ('Tan2011_sidD', '21734656'),
    ('Mishra2013_lepB', '23821544'),
    ('Ernst2020_ankX', '32440549'),
    ('Kaspers2023_lem3', '37076474'),
    ('Wang2018_setA', '30323948'),
    ('Belyi2022_setA_review', '35204756'),
    ('FollyKlan2013_ralF', '24244168'),
    ('Qiu2016_sidE', '27049943'),
    ('Black2019_sidJ', '31123136'),
    ('Belyi2008_lgt', '18281405'),
    ('Yan2020_mavC', '31994202'),
    ('Steiner2017_atl3', '28835546'),
    ('Prevost2017_lamtor1', '28842705'),
    ('Arasaki2010_lcv', '20163564'),
    ('Rolando2024_redundancy', '39636264'),
    ('Napier2019_salmonella', '31594818'),
    ('Im2026_setA', '41894332'),
    ('Gan2019_mavC', '30420781'),
    ('Fu2024_ubiquitin', '39009586'),
    ('Kubori2024_rab10', '38771316'),
    ('Chadha2024_bap31', '39661521'),
    ('Mount2024_effector_interactome', '39526800'),
    ('Mallama2017_typeII_secretion', '28138020'),
    ('Joseph2021_metaeffectors', '33499048'),
]
rows = []
for i in range(0, len(REFS), 50):
    chunk = REFS[i:i+50]
    d = esummary([p for _, p in chunk])
    for tag, pid in chunk:
        r = d.get(pid, {})
        rows.append({'tag': tag, 'pmid': pid,
                     'authors': ', '.join(a['name'] for a in r.get('authors', [])[:6]),
                     'year': r.get('pubdate', '')[:4],
                     'title': r.get('title', ''),
                     'journal': r.get('fulljournalname', ''),
                     'verified': 'YES' if r.get('title') else 'FAIL'})
    time.sleep(1)

df = pd.DataFrame(rows)
df.to_csv('${LEGIONELLA_ROOT}/13_postfreeze_interpretation/references_verified.tsv', sep='\t', index=False)
# BibTeX-ish
with open('${LEGIONELLA_ROOT}/15_manuscript/references_verified.bib', 'w', encoding='utf-8') as f:
    for _, r in df.iterrows():
        f.write("@misc{%s,\n  author = {%s},\n  year = {%s},\n  title = {%s},\n  journal = {%s},\n  pmid = {%s},\n  verified = {%s}\n}\n\n"
                % (r['tag'], r['authors'], r['year'], r['title'], r['journal'], r['pmid'], r['verified']))
print(df[['pmid','year','journal','title','verified']].to_string(index=False))
print("\nwrote references_verified.tsv + references_verified.bib")
