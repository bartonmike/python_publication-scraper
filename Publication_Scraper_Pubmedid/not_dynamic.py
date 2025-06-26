# NAME: publication_scraper.py
#
# create_paper_scraper.py MUST BE RUN BEFORE THIS CODE (only once to initialize the google sheet)
#
# Python Version Tested With: 3.12.2 / 3.12.3 / 3.9.18
#
# Author: Carter Deal
#         Anderson Lab
#         Oregon State University
#
# Date: December 2024
# 
# Description: This program takes in google scholar profile codes, along with names, and pulls data from the 
#              publications within a certain amount of weeks. It utilizes the biopython and habanero to fetch the entrez
#              and Crossref libraries respectively, using the PMID and name of the article to fetch the data. It also uses the three superfund funding codes
#              to locate publications in the pubmed database. It writes the data it fetches to a google sheet thats name was initialized in the create_paper_scraper.py 
#              program and written to a txt file named sheet_name.txt. This program develops a list of unique publications across all entered people. It also has added trainee 
#              functionality, where one can provide a list of trainees and the program will automatically search for those names in the publication authors, putting any matches onto the sheet.
#
# Requirements: 
#               Python Modules: 
#                       biopython
#                       fake_useragent
#                       beautifulsoup4
#                       habanero
#                       selenium
#                       gspread
#                       oauth2client
#                       lxml_html_clean
#
#               files and other requirements:
#                       client_key.json file in same directory, containing info on google developer bot
#                       create_paper_scraper.py ran, which creates all the needed txt files
#                       bot whose details are containted in the client_key.json has edit access to given google sheet
#              
# Inputs: Google scholar profile codes and names, trainee names, weeks ago one wants to search.
#
#
# Outputs: All outputs are on google sheet
#           Table: Title, Authors, Journal, DOI, PMID, PMCID, Pubdate, Funding, Trainees, Notes
#           Changelog: removed and added for the inputs 
#           Run Log: entries found for each person/funding
#
#
# Note: Trainees must be populated or the program will error


# what if first and last name match?
# what if last name and first initial match? (with that being the only information given)
# what if two people have the same alias?

####################################################################################### IMPORTS

####### Entrez (pubmed) api import
from Bio import Entrez # pip install biopython
from Bio.Entrez import efetch, read, esummary

###### misc imports
from datetime import date, timedelta
import requests
from fake_useragent import UserAgent # pip3 install fake_useragent
from bs4 import BeautifulSoup # pip3 install beautifulsoup4
import json
from habanero import Crossref # pip3 install habanero
import time
import os

####### google scholar imports
from selenium import webdriver # pip3 install selenium
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

###### gspread imports
import gspread # pip3 install gspread
from oauth2client.service_account import ServiceAccountCredentials # pip3 install oauth2client

###### crossref api import
import habanero

###### progress bar
from tqdm import tqdm, trange

import re

def is_full_pubmed_date(date_str):
    # Matches "YYYY MonthName DD" (e.g., "2024 Jun 15")
    return bool(re.match(r"^\d{4} [A-Za-z]+ \d{1,2}$", date_str.strip()))


from datetime import datetime

def convert_history_date(history_date_str):
    # Example: "2024/06/25 12:34"
    try:
        dt = datetime.strptime(history_date_str.split()[0], "%Y/%m/%d")
        print(dt.strftime("%Y %b %d"))
        return dt.strftime("%Y %b %d")  # e.g., "2024 Jun 25"
    except Exception:
        return history_date_str  # fallback, just in case
    

def get_best_pubdate(pubdata):

    # Try pubdate
    pubdate = pubdata.get('pubdate', '').strip()
    epubdate = pubdata.get('epubdate', '').strip()
    history = pubdata.get('history', [])
    ret_value = pubdate
    if pubdate:
        if is_full_pubmed_date(pubdate):
            ret_value = pubdate
        else:
            # Try to convert if not in correct format
            converted = convert_history_date(pubdate)
            if is_full_pubmed_date(converted):
                ret_value = converted

    elif epubdate:
            if is_full_pubmed_date(epubdate):
                ret_value = epubdate
            else:
                converted = convert_history_date(epubdate)
                if is_full_pubmed_date(converted):
                    ret_value = converted

        # Try history[3]['date']
    elif history:
        history_date = None
        if isinstance(history, list) and len(history) > 3:
            history_date = history[3].get('date', '').strip()
        elif isinstance(history, dict) and 3 in history:
            history_date = history[3].get('date', '').strip()
        if history_date:
            converted = convert_history_date(history_date)
            if is_full_pubmed_date(converted):
                ret_value = converted

    else:
        # If no pubdate, epubdate, or history date, return empty string
        ret_value = ''
    # If all fail, return empty string
    return ret_value

####################################################################################### DEFINED FUNCTIONS

# function that searches through the pubmed entrez database based on a minimum date 
# that takes a number of weeks ago decided by the user, returns a list of publication 
# results
def search(query, week):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='pub_date',
                            retmax='200',
                            retmode='xml',
                            # datetype = 'pdat',
                            mindate= str(date.today() - timedelta(weeks=week)),
                            maxdate= str(date.today()),
                            term=query)
    results = Entrez.read(handle)
    return results

# function that searches the pubmed entrez database using the doi of a publication
# and returns the results from pubmed. If a publication with the doi exists, only one publication 
# is returned, if there isn't then there is either many publications returned or none at all. 
def find_pmid(doi):
    Entrez.email = 'your.email@example.com'
    handle = Entrez.esearch(db='pubmed',
                            sort='pub_date',
                            retmax='200',
                            retmode='xml',
                            term=doi)
    results = Entrez.read(handle)
    return results

def fetch_funding(pmid):
    #print(pmid)
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
                funding += ','
        #print (funding)
    except:
        funding = ''
    return funding

# using the pmid of the publication, fetch an abstract (isn't used in the program)
def fetch_abstract(pmid):
    #print(pmid)
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
    #print(pmid)
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

def find_matches(trainees, authors):
        matches = ''

        # cross reference the trainee list to see if there are any trainees in the authors
        for trainee in trainees:
            train_match = False
            if trainee[0] in authors and trainee[1] in authors:
                train_match = True  

            if train_match:
                if matches != '':
                    matches += ', '
                matches += trainee[0] + ' ' + trainee[1]
        # if there are trainees, list them
        if(matches == ''):
            matches = 'No'
        else:
            matches = "Yes: " + matches

        return matches

def build_pub(doi, trainees, backup_authors, backup_journal, backup_pub_date):

    # header = ['Title', 'Raw Authors', 'Display Authors', 'Trainees', 'Journal', 'Issue', 'Volume', 'Pages', 'DOI', 'PMID', 'PMCID', 'ISSN', 'Pubdate', 'Funding', 'Abstract', 'Notes']
    # Create a Crossref API client
    cr = habanero.Crossref()

    # Search for a DOI using a query string
    result = cr.works(ids=doi)

    # Access the DOI information
    if result and 'message' in result:
        doi_info = result['message']

        authors = ''

        try:
            for i in range(len(doi_info['author'])):
                if(i == len(doi_info['author']) - 1 or len(doi_info['author']) == 1):
                    authors += doi_info['author'][i]['given'] + " " + doi_info['author'][i]['family']
                else:
                    authors += doi_info['author'][i]['given'] + " " + doi_info['author'][i]['family'] + ", "
        except:
            authors = backup_authors

        try:
            page = doi_info['page']
        except:
            page = ''

        try:
            abstract = doi_info['abstract']
        except:
            abstract = ''
        
        try:
            issue = doi_info['issue']
        except:
            issue = ''

        try:
            volume = doi_info['volume']
        except:
            volume = ''
        
        try:
            issn = doi_info['ISSN'][0]
        except:
            issn = ''

        try:
            journal = doi_info['container-title'][0]
        except:
            journal = backup_journal

        try: 
            pub_date = '/'.join(str(date) for date in doi_info['published']['date-parts'][0])
        except: 
            pub_date = backup_pub_date

        try:
            abstract = abstract.split('<jats:p>')[1].replace('</jats:p>','')
        except:
            abstract = abstract
        
        return [doi_info['title'][0], authors, '', find_matches(trainees, authors), str(doi_info['type']).replace('-',' '), journal, issue, 
                volume, page, doi, 'None', 'None', issn, pub_date, 'None', abstract, '']
    else:
        print("DOI not found.")

# parsing the results of a search() function. 
# finding the title, authors, Journal, doi, pmid, pmcid, publication date, funding, and trainees in each publicaiton in
# results. 
# returns an array with each publications data. 
def fetchData(results, funding, trainees, headers):
    array = []
    # for each publication in the results of pubmed
    tqdm_results = tqdm(results, desc='Fetching PUBMED Data', unit='publication')
    for id in tqdm_results:

        id = id.replace('\r','')

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

        # Entry: [Title, Raw Authors, Display Authors, Trainees, Pubtype, Journal, Volume, Issue, Pages, DOI, PMID, PMCID, ISSN, Pubdate, Funding, Abstraact, Notes]
        # append to array an entry

        try:
            journal = response['result'][id]['pubtype'][0]
        except:
            journal = ''
        array.append([response['result'][id]['title'], authors_fullname, '', find_matches(trainees, authors_fullname), journal, 
                      response['result'][id]['source'], response['result'][id]['volume'], response['result'][id]['issue'], response['result'][id]['pages'], 
                      doi, id, pmcid, response['result'][id]['issn'], get_best_pubdate(response['result'][id]), str(fetch_funding(id)), str(fetch_abstract(id)),''])
        # array.append([response['result'][id]['title'], authors_fullname, '', '', journal, 
        #         response['result'][id]['source'], response['result'][id]['volume'], response['result'][id]['issue'], response['result'][id]['pages'], 
        #         doi, id, pmcid, response['result'][id]['issn'], response['result'][id]['pubdate'], '', str(fetch_abstract(id)),''])
    
    #RunLog.append_row(['Number of Hits for: ' + funding, str(len(array)), str(date.today())])
    RunLog.columns_auto_resize(0,3)

    return array

# this function takes a result array built from google scholar results and 
# cross references it with previous lists and itself, checking for any duplicates
# and emitting them from the final list that will be returned
def build_list(list_of_lists, google_scholar):

    person_array = []

    # eliminate same articles between people (DOI match)
    for entry in google_scholar:
        not_dupe = True
        for array in list_of_lists:
            for item in array:
                if item[9] in entry[9]:
                    not_dupe = False
        if not_dupe:
            person_array.append(entry)

    # eliminate same articles within the person (DOI match)
    new_person_array = []
    for i in range(len(person_array)):
        duplicate = False
        for j in range(len(new_person_array)):
            if(person_array[i][9] in new_person_array[j][9]):
                duplicate = True
        if(not duplicate):
            new_person_array.append(person_array[i])

    return new_person_array

# This function checks for pubmed duplicates within the pubmed results based on pmid, if there are
# then emit them
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

# function that compares two arrays of strings and finds the removed and added values of the new string in comparison
# to the old string
def compare_people(new_people,old_people):
    diff_removed = []
    diff_added = []
    # added values
    for person in new_people:
        not_found = 1
        for old_person in old_people:
            if person in old_person:
                not_found  = 0
            
        if not_found:
            diff_added.append('"' + person + '"')

    str_diff1 = ' '.join(diff_added)
    # removed values
    for person in old_people:
        not_found = 1
        for new_person in new_people:
            if person in new_person:
                not_found  = 0
            
        if not_found:
            diff_removed.append('"' + person + '"')

    str_diff2 = ' '.join(diff_removed)

    return str_diff1, str_diff2 # return both added and removed items each in a single string

def normalize_author_key(first, middle, last):
    middle = (middle or '').strip().lower()
    return f"{first.strip().lower()}_{middle}_{last.strip().lower()}"

####################################################################################### Initiazing Variables 
#Google sheets initialization
#Authorize the API
scope = [
    'https://www.googleapis.com/auth/spreadsheets',
    'https://www.googleapis.com/auth/drive'
    ]
file_name = 'client_key.json'
creds = ServiceAccountCredentials.from_json_keyfile_name(file_name,scope)
client = gspread.authorize(creds)

f = open('sheet_name.txt','r')
sheet_name = f.readlines()[0].strip()
print(sheet_name)
# opening all the different worksheets in the google sheet
sheet = client.open(sheet_name).sheet1 
prompts = client.open(sheet_name).worksheet('Prompts')
ChangeLog = client.open(sheet_name).worksheet('Change Log')
RunLog = client.open(sheet_name).worksheet('Run Log')

disp_authors_start = len(sheet.col_values(3)) + 2

# opening and stripping all the files storing the previous values
with open("people.txt", 'r', encoding = 'utf-8') as file:
    old_people = file.readlines()

with open("trainees.txt", 'r', encoding = 'utf-8') as file:
    old_trainees = file.readlines()

# strip the lines of the text
for i in range(len(old_people)):
    old_people[i] = old_people[i].strip()

# strip the lines of the text
for i in range(len(old_trainees)):
    old_trainees[i] = old_trainees[i].strip()

# fetching people list from the sheet
people = prompts.cell(2,1).value

if(people != None):
    people = people.split('\n')

    # strip the lines of the text
    for i in range(len(people)):
        people[i] = people[i].strip()

# fetching trainees from the sheet
trainees = prompts.cell(4,1).value

if(trainees != None):
    trainees = trainees.split('\n')

    # strip the lines of the text
    for i in range(len(trainees)):
        trainees[i] = trainees[i].strip()

# # fetching the week from the sheet
# week = int(prompts.cell(6,1).value)

##################################################################################################### Added/Removed
# finding added and removed values from the people/trainees from the sheet 
# with the ones stored in the files

# people_diff1, people_diff2 = compare_people(people, old_people)
train_diff1, train_diff2 = compare_people(trainees, old_trainees)

str_diff1 = ''
str_diff2 = ''

week_diff = False

if(train_diff1 != ''):
    str_diff1 += '\n\tTrainees: ' + train_diff1
if(train_diff2 != ''):
    str_diff2 += '\n\tTrainees: ' + train_diff2

# if trainees changed update the file
if(train_diff2.strip() != '' or train_diff1.strip() != ''):
    c = open("trainees.txt", "w")
    c.write('\n'.join(trainees))
    c.close()


####################################################################################################### PUBMED/FUNDING

print("Searching PUBMED...............................................................")

# strip first and last name from the trainees
for i in range(len(trainees)):
    trainees[i] = trainees[i].strip().split(' ')

# setting useragent for requests
ua = UserAgent()
headers = {
        "User-Agent":
        str(ua.chrome)
        }

codes = prompts.cell(2,2).value.split('\n')

# After the codes initialization and before the code_prog declaration:
# Get all PMIDs from column K (skipping header)
existing_pmids = sheet.col_values(11)[1:]  # 11 is column K (0-based index)
time.sleep(1)  # Rate limit handling

# Remove empty values and strip whitespace
existing_pmids = [pmid.strip() for pmid in existing_pmids if pmid.strip()]

print(existing_pmids)
print(codes)

# Filter codes array to only include values not in existing_pmids
codes = [code for code in codes if code not in existing_pmids]




# If all codes were filtered out, exit early
if not codes:
    print("No new PMIDs to process")
    exit()

code_prog = tqdm(codes)

# creating a pubmed array that stores all the values of the different funding
pubmed_array = []


code_array = fetchData(codes, None, trainees, headers)

pubmed_array = checkDupe(pubmed_array, code_array)

# creates a master list of all the different people/pubmed results
list_of_lists  = []
list_of_lists.append(pubmed_array)

print("Done\n\n")

####################################################################################################################### Authors Names Sorting

# Helper function to update a specific cell in the authors worksheet
# Retries until successful due to potential API rate limits
def update_author_general(sheet_authors, value, row, column):
        while True:
            try:
                sheet_authors.update_cell(row - 1, column, value)
                break
            except:
                pass

# Helper function to add a new author row to the worksheet
# Includes first name, middle initial, last name, aliases, and initial publication count
def new_author(sheet_authors, author_dict, name, aliases, publication):
    while True:
        try:
            sheet_authors.append_row([author_dict['Sheets Index'] - 2, name[0], author_dict['Middle Initial'], name[len(name) - 1], aliases, 1, publication])
            break
        except:
            pass

# Creates a dictionary to store author information
# middle initial, index in spreadsheet, number of publications, and aliases
def create_author(initial, index, pub_num, aliases, doi):
    author = {}
    author['Middle Initial'] = initial
    author['Sheets Index'] = index
    author['Pub_Number'] = pub_num
    author['Aliases'] = aliases
    author['DOIs'] = doi
    return author

f = open('sheet_name.txt','r')
sheet_name = f.readlines()[0].strip()

print("Fetching Authors..........................\n")

# Fetch spaced last names from Prompts!A2 (comma-separated)
spaced_last_names_raw = prompts.cell(2, 1).value
if spaced_last_names_raw:
    spaced_last_names = [name.strip() for name in spaced_last_names_raw.split(',')]
else:
    spaced_last_names = []

def parse_author_name_with_spaced_last(author_name, spaced_last_names):
    """
    Splits author_name into first, middle (if any), and last name,
    using spaced_last_names to detect special last names.
    Returns (first, middle, last, key_last).
    """
    for spaced_last in spaced_last_names:
        if spaced_last in author_name:
            idx = author_name.find(spaced_last)
            first = author_name.split(" ")[0]
            last = ' '.join(author_name.split(" ")[-2:])
            rest = ' '.join(author_name.split(" ")[1:-2])
            # If there is anything after the spaced last name, treat as middle/initial (rare)
            if rest:
                middle = rest
            else:
                middle = ''
            return first, middle, last, last
    # Default: split on first and last
    parts = author_name.split()
    if len(parts) == 1:
        return parts[0], '', '', parts[0]
    elif len(parts) == 2:
        return parts[0], '', parts[1], parts[1]
    else:
        return parts[0], parts[1][0].upper(), parts[-1], parts[-1]

# try:
# Get access to the Known Authors worksheet
sheet_authors = client.open(sheet_name).worksheet('Known Authors')

# Dictionary to store all known authors
known_authors = {}

# # Start index after header row
# index = 2

# Get all values from the authors worksheet
full_authors = sheet_authors.get_all_values()
time.sleep(1)  # Rate limit handling

# Remove header row
full_authors.pop(0)

# Track used sheet indices
used_sheet_indices = set()
max_sheet_index = 0

# Process existing authors from spreadsheet
for i in range(len(full_authors)):
    row_values = full_authors[i]
    if(row_values == []):
        break
    # Sheets Index is in the first column (should be int)
    try:
        sheets_index = int(row_values[0])
    except Exception:
        continue
    used_sheet_indices.add(sheets_index)
    if sheets_index > max_sheet_index:
        max_sheet_index = sheets_index

    # Create key from first, middle, last name
    first = row_values[1].strip()
    middle = row_values[2].strip()
    last = row_values[3].strip()
    key = normalize_author_key(first, middle, last)
    known_authors[key] = create_author(
        middle,
        sheets_index,
        int(row_values[5]),
        row_values[4].split(', '),
        row_values[6].split(', ')
    )

# Helper to get next available Sheets Index
def get_next_sheets_index():
    # Try to find the lowest unused index from 1 to max_sheet_index
    for idx in range(1, max_sheet_index + 1):
        if idx not in used_sheet_indices:
            used_sheet_indices.add(idx)
            return idx
    # Otherwise, use the next after max
    next_idx = max_sheet_index + 1
    used_sheet_indices.add(next_idx)
    return next_idx

# Get access to main publication worksheet
while True:
    try:
        sheet = client.open(sheet_name).worksheet('Publication Info')
        break
    except:
        pass

# Initialize display authors array
Display_authors = ['Display Authors', '']

for person in list_of_lists:
    tqdm_publications = tqdm(person)
    tqdm_publications.set_description("Processing Authors")

    for publication in tqdm_publications:
        authors_str = publication[1]
        doi = publication[9]
        
        if not authors_str:
            Display_authors.append('')
            continue

        for author in authors_str.split(', '):
            # Use the new parser for spaced last names
            first, middle, last, key_last = parse_author_name_with_spaced_last(author, spaced_last_names)
            norm_key = normalize_author_key(first, middle, key_last)

            

            # Check for incomplete names (initials)
            parts = author.split()
            is_incomplete = (
                len(parts) == 2 and
                (len(parts[0].replace('.', '')) <= 2 or len(parts[1].replace('.', '')) <= 2)
            )

            found = False

            # Try direct match in known_authors
            if norm_key in known_authors:
                author_entry = known_authors[norm_key]
                # If middle initials don't match and both are present, treat as different person
                if author_entry['Middle Initial'] and middle and author_entry['Middle Initial'] != middle:
                    pass  # skip, look for other matches
                else:
                    # Update middle initial if missing
                    if not author_entry['Middle Initial'] and middle:
                        author_entry['Middle Initial'] = middle
                    # Update publication count and DOI
                    author_entry['Pub_Number'] += 1
                    if doi not in author_entry['DOIs']:
                        author_entry['DOIs'].append(doi)
                    # Add alias if new
                    if author not in author_entry['Aliases']:
                        author_entry['Aliases'].append(author)
                    found = True

            # If not found and is incomplete, search aliases for match
            if not found and is_incomplete:
                for k, v in known_authors.items():
                    if any(author == alias for alias in v['Aliases']):
                        # Found by alias
                        v['Pub_Number'] += 1
                        if doi not in v['DOIs']:
                            v['DOIs'].append(doi)
                        if author not in v['Aliases']:
                            v['Aliases'].append(author)
                        found = True
                        break

            # If still not found, create new author entry
            if not found:
                new_index = get_next_sheets_index()
                known_authors[norm_key] = create_author(middle, new_index, 1, [author], [doi])


import pprint

with open('Authors_dict.txt', 'w', encoding = 'utf-8') as dict:
    pprint.pprint(known_authors, dict)

# After processing all authors and before writing Authors_dict.txt

# 1. Delete everything after the first row in the "Known Authors" sheet
sheet_authors.resize(rows=1)

# 2. Prepare the data to append
authors_rows = []
for idx, (norm_key, author) in enumerate(known_authors.items(), start=1):
    # norm_key is "first_last"
    split_key = norm_key.split('_')
    first_name = split_key[0].capitalize()
    last_name = split_key[-1].capitalize() if len(split_key) > 1 else ''
    authors_rows.append([
        idx,                           # Sequential Index
        first_name,                    # First Name
        author['Middle Initial'],      # Middle Initial
        last_name,                     # Last Name
        ', '.join(author['Aliases']),  # Aliases
        author['Pub_Number'],          # # of Publications
        ', '.join(author['DOIs']) if isinstance(author['DOIs'], list) else author['DOIs']  # DOIs
    ])

# 3. Append all rows at once
# Calculate how many rows and columns you need
needed_rows = len(authors_rows) + 1  # +1 for header row
needed_cols = 7  # Index, First, Middle, Last, Aliases, # of Publications, DOIs

# Resize the sheet to fit the data
sheet_authors.resize(rows=needed_rows, cols=needed_cols)

# Now append the data
sheet_authors.append_rows(authors_rows, value_input_option='USER_ENTERED')

####################################################################################################################### Write to Google Sheets
final_list = []

# format list of lists to be one huge list for printing purposes
for array in list_of_lists:
    final_list += array

for publication in final_list:
    publication[2] = str('=display_authors_one_line(CONCATENATE("B", TO_TEXT(ROW(B' + str(disp_authors_start - 1) + '))))')
    disp_authors_start += 1

# write results to the google sheet
header = ['Title', 'Raw Authors', 'Display Authors', 'Trainees', 'Journal', 'Pubtype', 'Issue', 'Volume', 'Pages', 'DOI', 'PMID', 'PMCID', 'ISSN', 'Pubdate', 'Funding', 'Abstract', 'Notes']

titles = sheet.row_values(1)
time.sleep(1)

if(str(titles) == "[]"):
    sheet.insert_row(header,1)
    time.sleep(1)


sheet.append_rows(final_list, value_input_option='USER_ENTERED')
time.sleep(1)
sheet.columns_auto_resize(0, 9)
