from Bio import Entrez
from Bio.Entrez import efetch, esummary
from datetime import date, timedelta
import json
import requests
from habanero import Crossref
import habanero
from bs4 import BeautifulSoup

Entrez.email = 'your.email@example.com'

def search(query, week):
    handle = Entrez.esearch(db='pubmed',
                            sort='pub_date',
                            retmax='200',
                            retmode='xml',
                            mindate= str(date.today() - timedelta(weeks=week)),
                            maxdate= str(date.today()),
                            term=query)
    results = Entrez.read(handle)
    return results

def find_pmid(doi):
    handle = Entrez.esearch(db='pubmed',
                            sort='pub_date',
                            retmax='200',
                            retmode='xml',
                            term=doi)
    results = Entrez.read(handle)
    return results

def fetch_abstract(pmid):
    while True:
        try:
            handle = efetch(db='pubmed', id=pmid, retmode='xml')
            record = Entrez.read(handle)
            break
        except Exception:
            pass
    abstract = ''
    try:
        abstract = str(record['PubmedArticle'][0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0])
    except Exception:
        abstract = ''
    return abstract

def fetch_summary(pmid):
    while True:
        try:
            handle = esummary(db='pubmed', id=pmid, retmode='json')
            break
        except Exception:
            pass
    return handle.read()

def find_authors(doi):
    cr = habanero.Crossref()
    result = cr.works(ids = doi)
    if result and 'message' in result:
        doi_info = result['message']
        authors = []
        try:
            for i in range(len(doi_info.get('author', []))):
                authors += [doi_info['author'][i].get('given','') + " " + doi_info['author'][i].get('family','')]
            return ', '.join(authors)
        except Exception:
            return 'Authors Not Found'
    else:
        return 'Authors Not Found'

def find_matches(trainees, authors):
    matches = ''
    if trainees == None:
        return ''
    train_match = False
    for trainee in trainees:
        if trainee[0] in authors and trainee[1] in authors:
            train_match = True
        if train_match:
            if matches != '':
                matches += ', '
            matches += trainee[0] + ' ' + trainee[1]
    if(matches == ''):
        matches = 'No'
    else:
        matches = "Yes: " + matches
    return matches

def build_pub(doi, trainees, backup_authors, backup_journal, backup_pub_date):
    cr = Crossref()
    result = cr.works(ids=doi)
    if result and 'message' in result:
        doi_info = result['message']
        authors = ''
        try:
            for i in range(len(doi_info.get('author', []))):
                if(i == len(doi_info.get('author', [])) - 1 or len(doi_info.get('author', [])) == 1):
                    authors += doi_info['author'][i].get('given','') + " " + doi_info['author'][i].get('family','')
                else:
                    authors += doi_info['author'][i].get('given','') + " " + doi_info['author'][i].get('family','') + ", "
        except Exception:
            authors = backup_authors

        page = doi_info.get('page','')
        abstract = doi_info.get('abstract','')
        issue = doi_info.get('issue','')
        volume = doi_info.get('volume','')
        issn = doi_info.get('ISSN',[None])[0] if doi_info.get('ISSN') else ''
        journal = doi_info.get('container-title',[backup_journal])[0]
        try:
            pub_date = '/'.join(str(d) for d in doi_info['published']['date-parts'][0])
        except Exception:
            pub_date = backup_pub_date
        try:
            abstract = abstract.split('<jats:p>')[1].replace('</jats:p>','')
        except Exception:
            abstract = abstract

        return [doi_info.get('title',[''])[0], authors, '', find_matches(trainees, authors), str(doi_info.get('type','')).replace('-',' '), journal, issue,
                volume, page, doi, 'None', 'None', issn, pub_date, 'None', abstract, '']
    else:
        return None

def fetchData(results, funding, trainees, headers, RunLog):
    array = []
    for id in results['IdList']:
        response = json.loads(fetch_summary(id))
        doi = ''
        pmcid = ''
        for i in range(len(response['result'][id]['articleids'])):
            if(response['result'][id]['articleids'][i]['idtype'] == 'pmc'):
                pmcid = response['result'][id]['articleids'][i]['value']
            if(response['result'][id]['articleids'][i]['idtype'] == 'doi'):
                doi = response['result'][id]['articleids'][i]['value']

        if doi == '':
            try:
                doi = str(response['result'][id]['references'][0]['refsource']).split('doi: ')[1]
            except Exception:
                doi = ''

        authors_fullname = find_authors(doi) if doi else 'Authors Not Found'

        if(authors_fullname == 'Authors Not Found'):
            link = 'https://search.crossref.org/search/works?q=' + response['result'][id]['articleids'][3]['value'] +'&from_ui=yes'
            p = requests.get(link, headers = headers, timeout=30)
            soup = BeautifulSoup(p.text, features = "html.parser")
            author_text = soup.find(attrs={'class':'expand'})
            if author_text:
                authors = author_text.text.strip().replace('Authors: ','').split(' | ')
                authors_fullname = ', '.join(authors)

        try:
            journal = response['result'][id]['pubtype'][0]
        except Exception:
            journal = ''

        array.append([response['result'][id]['title'], authors_fullname, '', find_matches(trainees, authors_fullname), journal,
                      response['result'][id]['source'], response['result'][id]['volume'], response['result'][id]['issue'], response['result'][id]['pages'],
                      doi, id, pmcid, response['result'][id]['issn'], response['result'][id]['pubdate'], funding, str(fetch_abstract(id)),''])

    RunLog.append_row(['Number of Hits for: ' + funding, str(len(array)), str(date.today())])
    RunLog.columns_auto_resize(0,3)

    return array

def search_pubmed_by_title(title, trainees, headers):
    try:
        handle = Entrez.esearch(db='pubmed', term=title, retmax='5', retmode='xml')
        results = Entrez.read(handle)
    except Exception:
        return None

    def _norm(s):
        if not s:
            return ''
        return ''.join(ch.lower() for ch in s if ch.isalnum())

    for pmid in results.get('IdList', []):
        try:
            resp = json.loads(fetch_summary(pmid))
        except Exception:
            continue

        r = resp['result'].get(pmid, {})
        resp_title = r.get('title', '')
        if _norm(resp_title) == _norm(title) or _norm(title) in _norm(resp_title) or _norm(resp_title) in _norm(title):
            doi = ''
            pmcid = ''
            for aid in r.get('articleids', []):
                if aid.get('idtype') == 'doi':
                    doi = aid.get('value', '')
                if aid.get('idtype') == 'pmc':
                    pmcid = aid.get('value', '')

            authors_fullname = ''
            if doi:
                authors_fullname = find_authors(doi)

            if authors_fullname == 'Authors Not Found' or authors_fullname == '':
                try:
                    authors_list = r.get('authors', [])
                    authors_fullname = ', '.join([a.get('name') if isinstance(a, dict) else str(a) for a in authors_list])
                except Exception:
                    authors_fullname = 'Authors Not Found'

            journal = ''
            try:
                journal = r.get('pubtype', [''])[0]
            except Exception:
                journal = ''

            entry = [
                r.get('title', ''),
                authors_fullname,
                '',
                find_matches(trainees, authors_fullname),
                journal,
                r.get('source', ''),
                r.get('volume', ''),
                r.get('issue', ''),
                r.get('pages', ''),
                doi,
                pmid,
                pmcid,
                r.get('issn', ''),
                r.get('pubdate', ''),
                'None',
                str(fetch_abstract(pmid)),
                ''
            ]
            return entry

    return None

def checkDupe(array, master_array):
    for item in array:
        dup = False
        for i in range(len(master_array)):
            if(item[10] in master_array[i][10]):
                master_array[i][14] += ', ' + item[14]
                dup = True
        if not dup:
            master_array.append(item)
    return master_array


def resolve_google_entries(person_array, codes, headers):
    """For each publication in person_array, attempt to resolve missing PMID/PMCID
    and scan the DOI landing page for funding codes, updating the funding field.
    """
    from bs4 import BeautifulSoup
    for i in range(len(person_array)):
        # resolve PMID/PMCID if missing
        if person_array[i][10] == 'None':
            try:
                pubmed = find_pmid(person_array[i][9])
            except Exception:
                pubmed = {'IdList': []}

            for id in pubmed.get('IdList', []):
                try:
                    name = json.loads(fetch_summary(id))
                except Exception:
                    continue

                pmc = ''
                doi = ''
                for j in range(len(name['result'][id].get('articleids', []))):
                    aid = name['result'][id]['articleids'][j]
                    if aid.get('idtype') == 'doi':
                        doi = aid.get('value', '')
                    if aid.get('idtype') == 'pmc':
                        pmc = aid.get('value', '')

                if doi == '':
                    try:
                        doi = str(name['result'][id]['references'][0]['refsource']).split('doi: ')[1]
                    except Exception:
                        doi = ''

                if doi == person_array[i][9]:
                    person_array[i][10] = id
                    person_array[i][11] = pmc
                else:
                    try:
                        pubmed2 = find_pmid(person_array[i][0])
                    except Exception:
                        pubmed2 = {'IdList': []}
                    if str(name['result'][id].get('title','')) == person_array[i][0]:
                        person_array[i][10] = id
                        person_array[i][11] = pmc
                break

        # scan DOI landing page for funding codes
        try:
            if person_array[i][9]:
                p = requests.get('https://doi.org/' + person_array[i][9], headers = headers, timeout=30)
            else:
                p = None
        except Exception:
            p = None

        if p is None:
            continue

        soup = BeautifulSoup(p.text, features = "html.parser")
        text = soup.get_text()

        for item in codes:
            if item.replace("[Grants and Funding]",'') in text:
                if person_array[i][13] == 'None':
                    person_array[i][13] = item.replace("[Grants and Funding]",'')
                else:
                    person_array[i][13] += ', ' + item.replace("[Grants and Funding]",'')

    return person_array

def build_list(list_of_lists, google_scholar):
    person_array = []
    for entry in google_scholar:
        not_dupe = True
        for array in list_of_lists:
            for item in array:
                if item[9] in entry[9]:
                    not_dupe = False
        if not_dupe:
            person_array.append(entry)

    new_person_array = []
    for i in range(len(person_array)):
        duplicate = False
        for j in range(len(new_person_array)):
            if(person_array[i][9] in new_person_array[j][9]):
                duplicate = True
        if(not duplicate):
            new_person_array.append(person_array[i])

    return new_person_array

def compare_people(new_people,old_people):
    diff_removed = []
    diff_added = []

    def normalize_person(person):
        if isinstance(person, (list, tuple)):
            return ', '.join(str(part).strip() for part in person)
        return str(person).strip()

    new_people = [normalize_person(person) for person in new_people]
    old_people = [normalize_person(person) for person in old_people]

    for person in new_people:
        if person not in old_people:
            diff_added.append('"' + person + '"')
    str_diff1 = ' '.join(diff_added)
    for person in old_people:
        if person not in new_people:
            diff_removed.append('"' + person + '"')
    str_diff2 = ' '.join(diff_removed)
    return str_diff1, str_diff2


def search_crossref_by_title(title, trainees):
    """Search Crossref by title and return a publication entry if a fuzzy match is found.
    Returns None if no match found.
    """
    if not title or not title.strip():
        return None

    cr = Crossref()
    try:
        result = cr.works(query=title, limit=5)
    except Exception:
        return None

    def _norm(s):
        if not s:
            return ''
        return ''.join(ch.lower() for ch in s if ch.isalnum())

    if not result or 'message' not in result or 'items' not in result['message']:
        return None

    for item in result['message'].get('items', []):
        item_title = item.get('title', [''])[0] if isinstance(item.get('title'), list) else item.get('title', '')
        if not item_title:
            continue

        if _norm(item_title) == _norm(title) or _norm(title) in _norm(item_title) or _norm(item_title) in _norm(title):
            doi = item.get('DOI', '')
            if not doi:
                continue

            authors = ''
            try:
                authors_list = item.get('author', [])
                authors = ', '.join([f"{a.get('given', '')} {a.get('family', '')}".strip() for a in authors_list if a.get('family')])
            except Exception:
                authors = 'Authors Not Found'

            try:
                journal = item.get('container-title', [''])[0] if isinstance(item.get('container-title'), list) else item.get('container-title', '')
            except Exception:
                journal = ''

            try:
                pub_date_parts = item.get('published', {}).get('date-parts', [[]])[0]
                pub_date = '/'.join(str(d) for d in pub_date_parts) if pub_date_parts else ''
            except Exception:
                pub_date = ''

            try:
                abstract = item.get('abstract', '')
            except Exception:
                abstract = ''

            volume = item.get('volume', '')
            issue = item.get('issue', '')
            pages = item.get('page', '')
            issn = item.get('ISSN', [None])[0] if item.get('ISSN') else ''

            return [item_title, authors, '', find_matches(trainees, authors), item.get('type', '').replace('-', ' '),
                    journal, issue, volume, pages, doi, 'None', 'None', issn, pub_date, 'None', abstract, '']

    return None
