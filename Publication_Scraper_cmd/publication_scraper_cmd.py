# Simplified version of the publication scraper that scrapes publication data from PubMed using the PMID of one publication, and writes it to a CSV file.


####################################################################################### IMPORTS

####### Entrez (pubmed) api import
from Bio import Entrez # pip3 install biopython
from Bio.Entrez import efetch, read, esummary

###### misc imports
from datetime import date, timedelta
import requests
from fake_useragent import UserAgent # pip3 install fake_useragent
from bs4 import BeautifulSoup # pip3 install beautifulsoup4
import json
import sys
import csv

###### crossref api import
import habanero # pip3 install habanero

####################################################################################### DEFINED FUNCTIONS

def fetch_funding(pmid):
    while True:
        try:
            handle = efetch(db='pubmed', id=pmid, retmode='xml')
            record = Entrez.read(handle)
            break
        except: 
            pass
    try:
        funding = ''
        for i in range(len(record['PubmedArticle'][0]["MedlineCitation"]["Article"]["GrantList"])):
            funding += str(record['PubmedArticle'][0]["MedlineCitation"]["Article"]["GrantList"][i]["GrantID"])
            if(i != len(record['PubmedArticle'][0]["MedlineCitation"]["Article"]["GrantList"]) - 1):
                funding += ', '
    except:
        funding = ''
    return funding

# using the pmid of the publication, fetch an abstract (isn't used in the program)
def fetch_abstract(pmid):
    while True:
        try:
            handle = efetch(db='pubmed', id=pmid, retmode='xml')
            record = Entrez.read(handle)
            break
        except: 
            pass
    try:
        abstract = str(record['PubmedArticle'][0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0])
    except:
        abstract = ''
    return abstract

# New Array: [Title, Raw Authors, Display Authors, Trainees, Pubtype, Journal, Issue, Volume, Pages, DOI, PMID, PMCID, ISSN, Pubdate, Funding, Abstraact, Notes]

# using the pmid of the publication, fetch a summary of it
def fetch_summary(pmid):
    while True:
        try:
            handle = esummary(db='pubmed', id=pmid, retmode='json')
            break
        except: 
            pass
    return handle.read()

def find_authors(doi):
    cr = habanero.Crossref()
    result = cr.works(ids = doi)
    if result and 'message' in result:
        doi_info = result['message']

        authors = []

        for i in range(len(doi_info['author'])):
            authors += [doi_info['author'][i]['given'] + " " + doi_info['author'][i]['family']]

        return ', '.join(authors)
    else:
        return 'Authors Not Found'

# parsing the results of a search() function. 
# finding the title, authors, Journal, doi, pmid, pmcid, publication date, funding, and trainees of publication.
# returns publication data. 
def fetchData(pmid, headers):

    response = json.loads(fetch_summary(pmid))

    doi = ''
    pmcid = ''

    for i in range(len(response['result'][pmid]['articleids'])):

        if(response['result'][pmid]['articleids'][i]['idtype'] == 'pmc'):
            pmcid = response['result'][pmid]['articleids'][i]['value']

        if(response['result'][pmid]['articleids'][i]['idtype'] == 'doi'):
            doi = response['result'][pmid]['articleids'][i]['value']

    if doi == '':
        try:
            doi = str(response['result'][pmid]['references'][0]['refsource']).split('doi: ')[1]
        except:
            doi = 'Could Not find DOI'

    try:
        authors_fullname = find_authors(doi)
    except:
        authors_fullname = 'Authors Not Found'

    if(authors_fullname == 'Authors Not Found'):

        link = 'https://search.crossref.org/search/works?q=' + response['result'][id]['articleids'][3]['value'] +'&from_ui=yes'

        p = requests.get(link, headers = headers, timeout=10)
        soup = BeautifulSoup(p.text, features = "html.parser")

        author_text = soup.find(attrs={'class':'expand'})

        authors = author_text.text.strip().replace('Authors: ','').split(' | ')
        authors_fullname = ', '.join(authors)

    try:
        journal = response['result'][pmid]['pubtype'][0]
    except:
        journal = ''

    
    # Entry: [Title, Raw Authors, Display Authors, Trainees, Pubtype, Journal, Volume, Issue, Pages, DOI, PMID, PMCID, ISSN, Pubdate, Funding, Abstract, Notes]

    return [response['result'][pmid]['title'], authors_fullname, '', '', journal, 
                    response['result'][pmid]['source'], response['result'][pmid]['volume'], response['result'][pmid]['issue'], response['result'][pmid]['pages'], 
                    doi, pmid, pmcid, response['result'][pmid]['issn'], response['result'][pmid]['pubdate'], str(fetch_funding(pmid)), str(fetch_abstract(pmid)),'']

####################################################################################### Initializing Variables 

print("Total arguments:", len(sys.argv))
print("Script name:", sys.argv[0])
print("Arguments:", sys.argv[1:])

pmid = sys.argv[1]

####################################################################################################### PUBMED/FUNDING

print("Searching PUBMED...............................................................")

# setting useragent for requests
ua = UserAgent()
headers = {
        "User-Agent":
        str(ua.chrome)
        }

# If all codes were filtered out, exit early
if not pmid:
    print("No new PMIDs to process")
    exit()

publication_data = fetchData(pmid, headers)

####################################################################################################################### Write to CSV

header = ['Title', 'Raw Authors', 'Display Authors', 'Trainees', 'Journal', 'Pubtype', 'Issue', 'Volume', 'Pages', 'DOI', 'PMID', 'PMCID', 'ISSN', 'Pubdate', 'Funding', 'Abstract', 'Notes']


with open('pmid.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(header)
    writer.writerow(publication_data)

print("Done\nResults written to pmid.csv\n")
