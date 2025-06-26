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


# using the pmid of the publication, fetch an abstract (isn't used in the program)
def fetch_abstract(pmid):
    print(pmid)
    while True:
        try:
            handle = efetch(db='pubmed', id=pmid, retmode='xml')
            record = Entrez.read(handle)
            break
        except: 
            pass
    abstract = str(record['PubmedArticle'][0]["MedlineCitation"]["Article"]["Abstract"]["AbstractText"][0])
    return abstract

# New Array: [Title, Raw Authors, Display Authors, Trainees, Pubtype, Journal, Issue, Volume, Pages, DOI, PMID, PMCID, ISSN, Pubdate, Funding, Abstraact, Notes]

# using the pmid of the publication, fetch a summary of it
def fetch_summary(pmid):
    print(pmid)
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
        train_match = False
        for trainee in trainees:
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
            doi = str(response['result'][id]['references'][0]['refsource']).split('doi: ')[1]

        authors_fullname = find_authors(doi)

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
                      doi, id, pmcid, response['result'][id]['issn'], response['result'][id]['pubdate'], funding, str(fetch_abstract(id)),''])
    
    RunLog.append_row(['Number of Hits for: ' + funding, str(len(array)), str(date.today())])
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

with open("date.txt", 'r', encoding = 'utf-8') as file:
    old_week = file.readlines()

# strip the lines of the text
for i in range(len(old_people)):
    old_people[i] = old_people[i].strip()

# strip the lines of the text
for i in range(len(old_trainees)):
    old_trainees[i] = old_trainees[i].strip()

# strip the lines of the text
for i in range(len(old_week)):
    old_week[i] = old_week[i].strip()

# fetching people list from the sheet
people = prompts.cell(2,1).value
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

# fetching the week from the sheet
week = int(prompts.cell(6,1).value)

##################################################################################################### Added/Removed
# finding added and removed values from the people/trainees from the sheet 
# with the ones stored in the files

people_diff1, people_diff2 = compare_people(people, old_people)
train_diff1, train_diff2 = compare_people(trainees, old_trainees)

str_diff1 = ''
str_diff2 = ''

week_diff = False

if(old_week == []):
    old_week.append(None)

if(old_week[0] != str(week)):
    week_diff = True

# formatting the output strings
if(people_diff1 != '' or train_diff1 != '' or week_diff):
    str_diff1 = 'Added:'
if(people_diff2 != '' or train_diff2 != '' or week_diff):
    str_diff2 = 'Removed:'

if(people_diff1 != ''):
    str_diff1 += '\n\tPrompts: ' + people_diff1
if(people_diff2 != ''):
    str_diff2 += '\n\tPrompts: ' + people_diff2

if(train_diff1 != ''):
    str_diff1 += '\n\tTrainees: ' + train_diff1
if(train_diff2 != ''):
    str_diff2 += '\n\tTrainees: ' + train_diff2

if(week_diff):
    str_diff1 += '\n\tWeeks: ' + str(week)
    str_diff2 += '\n\tWeeks: ' + str(old_week[0])

# if there are changes and, if so, what changes exactly
if(str_diff1.strip() != '' and str_diff2.strip() != ''):
    # appending to the document the changes
    ChangeLog.append_row([str_diff1 + '\n' + str_diff2, str(date.today())])
    ChangeLog.columns_auto_resize(0, 2)
elif(str_diff1.strip() != ''):
    # appending to the document the changes
    ChangeLog.append_row([str_diff1, str(date.today())])
    ChangeLog.columns_auto_resize(0, 2)
elif(str_diff2.strip() != ''):
    # appending to the document the changes
    ChangeLog.append_row([str_diff2, str(date.today())])
    ChangeLog.columns_auto_resize(0, 2)

# if people changed update the file
if(people_diff1.strip() != '' or people_diff2.strip() != ''):
    c = open("people.txt", "w")
    c.write('\n'.join(people))
    c.close()

# if trainees changed update the file
if(train_diff2.strip() != '' or train_diff1.strip() != ''):
    c = open("trainees.txt", "w")
    c.write('\n'.join(trainees))
    c.close()

# if week changed update the file
if(week_diff):
    c = open("date.txt", "w")
    c.write(str(week))
    c.close()

####################################################################################################### PUBMED/FUNDING

print("Searching PUBMED...............................................................")
# split the google scholar code from the name of the person
for i in range(len(people)):
    people[i] = people[i].split(', ')

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
code_prog = tqdm(codes)

# creating a pubmed array that stores all the values of the different funding
pubmed_array = []

# go through all the funding codes and find results on pubmed, check for duplicates and update accordingly
for code in code_prog:
    code_prog.set_description("Current Code: " + code)
    code_prog.refresh()

    results = search(code, week)

    code_array = fetchData(results, code.replace("[Grants and Funding]",''), trainees, headers)

    pubmed_array = checkDupe(pubmed_array, code_array)

# creates a master list of all the different people/pubmed results
list_of_lists  = []
list_of_lists.append(pubmed_array)

print("Done\n\n")

############################################################################################################# Google Scholar

print("Going through Google Scholar Profiles.........................................")

options = webdriver.ChromeOptions()
options.add_experimental_option("excludeSwitches", ['enable-logging'])
options.add_argument('--headless=new')
options.add_argument('--log-level=3')
#options.add_argument('--no-sandbox')
options.set_capability("browserVersion", "117")
options.add_argument("--disable-dev-shm-usage") #overcome limited resource problems
driver = webdriver.Chrome(options = options)

tqdm_people = tqdm(people)

# for each person listed in the sheet
for term in tqdm_people:
    tqdm_people.set_description("Current Person: " + term[1])
    tqdm_people.refresh()
    # initializing chrome tab and searching


    # building url 
    url = 'https://scholar.google.com/citations?view_op=list_works&hl=en&hl=en&user=' + term[0] + '&sortby=pubdate'

    # driver.get(url)

    # wait = WebDriverWait(driver, 10)

    driver.get(url)

    get_url = driver.current_url
    #wait.until(EC.url_to_be(url))

    if get_url == url:
        page_source = driver.page_source

    # finding each tile with an article on the person's google scholar page
    elements = driver.find_elements(By.CLASS_NAME,"gsc_a_at")

    google_scholar = []

    i = 0
    # for each publication
    for element in elements:
        # locate the link
        link = element.get_attribute('href')

        # open the link to the google scholar summary page for the publication
        p = requests.get(link, headers = headers, timeout=10)
        soup = BeautifulSoup(p.text, features = "html.parser")

        # fetch the body paragraph with all the data
        parts = soup.find_all(attrs= {'class':'gs_scl'})
        # locate the title of the article with link to source
        parent = soup.find(attrs={'class':'gsc_oci_title_link'})

        # get title and source link
        print(link)
        try:
            title = parent.get_text()
            url = parent['href']
        except:
            parent = soup.find(attrs={'id':'gsc_oci_title'})
            title = parent.get_text().replace('<div id="gsc_oci_title">', '').replace('</div>', '')

        
        out_of_date = False

        # for each section in the summary. Find author, pub date, journal
        for part in parts:
            text = part.find(attrs={'class':'gsc_oci_field'})
            values = part.find(attrs={'class': 'gsc_oci_value'})

            if 'Authors' in text:
                Authors = values.get_text()
            if 'Publication date' in text:
                Publication_date = values.get_text()
            if 'Journal' in text:
                Journal = values.get_text()

        # this section determines whether or not the publication date of the article
        # is within the bounds of the dates given, it is tricky since the pubdate isn't always
        # the exact date, it could only be the year or year/month, therefore we need to test for each
        # scenario.
        date_today = str(date.today() - timedelta(weeks=week)).replace('-','/')

        pubcount = Publication_date.count('/')

        withindate = False
        if(Publication_date == ''): # if publication date doesn't exist
            withindate = False

        elif(pubcount == 0): # if pubdate is only the year
            if(date_today >= Publication_date):
                withinddate = True

        elif (pubcount == 1): # if pubdate is year/month
            if(len(Publication_date) < 7):
                Publication_date.replace('/','/0')

            if(date_today.replace('/','') >= Publication_date.replace('/','')):
                withindate = True

        elif (pubcount == 2): # if pubdate is year/month/day
            Split_pub = Publication_date.split('/')

            if(len(Split_pub[1]) < 2):
                Split_pub[1] = '0' + Split_pub[1]
            
            if(len(Split_pub[2]) < 2):
                Split_pub[2] = '0' + Split_pub[2]

            Split_pub = ''.join(Split_pub)
            date_today = date_today.replace('/','')

            if(Split_pub >= date_today):
                withindate = True

        # if it is within the date, add to array, else break from loop as google scholar is sorted by date published
        if (withindate):
            # cr = Crossref()
            # result = cr.works(query = title)
            # doi = result['message']['items'][0]['DOI'].removesuffix('.s001')

            link = 'https://search.crossref.org/search/works?q=' + title +'&from_ui=yes'

            p = requests.get(link, headers = headers)
            soup = BeautifulSoup(p.text, features = "html.parser")

            #doi_text = soup.find(attrs={'class':'icon-external-link'})
            doi_text = soup.find(attrs={'class':'item-links'}).text

            doi_text = doi_text.strip().split('\n')[0].replace('https://doi.org/','')
            print("\n")
            print(doi_text)
            print(link)
            print("\n")            
            google_scholar.append(build_pub(doi_text, trainees, Authors, Journal, Publication_date))
        else: 
            out_of_date = True
            break

        i += 1

    # add hit number of person to runlog
    RunLog.append_row(['Number of Hits for: ' + term[1], str(len(google_scholar)), str(date.today())])
    RunLog.columns_auto_resize(0,3)

    # build an official unique list from the person, checking for duplicates from other lists
    # and the list itself
    person_array = build_list(list_of_lists, google_scholar)

############################################################################################################################### Google Scholar Authors/PMID/PMCID/
    # for each unique entry
    for i in range(len(person_array)):
        # if there is no PMID, attempt to find it and the PMCID
        if person_array[i][10] == 'None':
            pubmed = find_pmid(person_array[i][9])
            for id in pubmed['IdList']:
                name = json.loads(fetch_summary(id))

                pmc = ''
                doi = ''

                for j in range(len(name['result'][id]['articleids'])):

                    if(name['result'][id]['articleids'][j]['idtype'] == 'doi'):
                        doi = name['result'][id]['articleids'][j]['value']

                    if(name['result'][id]['articleids'][j]['idtype'] == 'pmc'):
                        pmc = name['result'][id]['articleids'][j]['value']
                
                if doi == '':
                    doi = str(name['result'][id]['references'][0]['refsource']).split('doi: ')[1]

                if doi == person_array[i][9]:
                    person_array[i][10] = id
                    person_array[i][11] = pmc
                else:
                    pubmed = find_pmid(person_array[i][0])
                    if str(name['result'][id]['title']) == person_array[i][0]:
                        person_array[i][10] = id
                        person_array[i][11] = pmc
                break
        
        try:
            p = requests.get('https://doi.org/' + person_array[i][9], headers = headers)
        except:
            print("Error fetching doi for " + person_array[i][9] + " continuing...\n")
            continue
            
        soup = BeautifulSoup(p.text, features = "html.parser")
        text = soup.get_text()

        for item in codes:
            # if((item.replace("[Grants and Funding]",'') in text) or (item.replace("[Grants and Funding]",'').split(' ')[1] in text)):
            if(item.replace("[Grants and Funding]",'') in text):
                if(person_array[i][13] == 'None'):
                    person_array[i][13] = item.replace("[Grants and Funding]",'')
                else:
                    person_array[i][13] += ', ' + item.replace("[Grants and Funding]",'')

        # # if there was an error, put in runlog
        # except Exception as e:
        #     RunLog.append_row(['Runtime Error, continuing program...', str(e), str(date.today())])
    # append to master list the person's google scholar results
    list_of_lists.append(person_array)

driver.close()
driver.quit()
####################################################################################################################### Authors Names Sorting

def update_author_general(sheet_authors, value, row, column):
        while True:
            try:
                sheet_authors.update_cell(row - 1, column, value)
                break
            except Exception as e:
                print(e)
                pass

def new_author(sheet_authors, author_dict, name, aliases):
    while True:
        try:
            sheet_authors.append_row([name[0], author_dict['Middle Initial'], name[len(name) - 1], aliases, 1])
            break
        except:
            pass

def create_author(initial, index, pub_num, aliases):
    author = {}
    author['Middle Initial'] = initial
    author['Sheets Index'] = index
    author['Pub_Number'] = pub_num
    author['Aliases'] = aliases
    return author

f = open('sheet_name.txt','r')
sheet_name = f.readlines()[0].strip()

print("Fetching Authors..........................\n")

sheet_authors = client.open(sheet_name).worksheet('Known Authors')

known_authors = {}

index = 2

full_authors = sheet_authors.get_all_values()
time.sleep(1)

full_authors.pop(0)
# time.sleep(5)
# column1 = sheet_authors.col_values(1)
# column2 = sheet_authors.col_values(2)
# column3 = sheet_authors.col_values(3)
# column4 = sheet_authors.col_values(4)
# column5 = sheet_authors.col_values(5)

# for i in range(len(column1)):
#     if(i != 0):
#         full_authors.append([column1[i],column2[i],column3[i],column4[i],column5[i]])

for i in range(len(full_authors)):
    row_values = full_authors[i]
    if(row_values == []):
        break
    index += 1

    print(row_values)

    key = row_values[0] + ' ' + row_values[2]
    known_authors[key] = create_author(row_values[1], index, 0, [])

    if len(row_values) > 4:
        if(row_values[3] != ''):
            known_authors[key]['Aliases'] = row_values[3].split(', ')
        if(row_values[4] != ''):
            known_authors[key]['Pub_Number'] = int(row_values[4])

    print(index)

while True:
    try:
        sheet = client.open(sheet_name).worksheet('Publication Info')
        break
    except:
        pass

Display_authors = ['Display Authors', '']

tqdm_authors = tqdm(list_of_lists)

tqdm_authors.set_description("Checking Publication Author Names")

for person in tqdm_authors:
    for publication in person:
        item = publication[1]

        if item == '':
            Display_authors.append('')

        else:
            auth_array = item.split(', ')

            for person in auth_array:
                person_arr = person.split(' ')

                for i in range(len(person_arr) - 1):
                    if person_arr[i] == '':
                        person_arr.remove('')

                if len(person_arr) >= 4: 
                    first_last = person_arr[0] + ' ' + person_arr[len(person_arr) - 1] 

                    if(first_last in known_authors):

                        if(known_authors[first_last]['Middle Initial'] == '' or known_authors[first_last]['Middle Initial'] == person_arr[1][0].upper() or known_authors[first_last]['Middle Initial'] == person_arr[2][0].upper()):
                            known_authors[first_last]['Middle Initial'] = person_arr[1][0].upper()

                            update_author_general(sheet_authors, person_arr[1][0].upper(), known_authors[first_last]['Sheets Index'], 2)
                            
                            known_authors[first_last]['Pub_Number'] += 1

                            update_author_general(sheet_authors, known_authors[first_last]['Pub_Number'], known_authors[first_last]['Sheets Index'], 5)

                            time.sleep(2)
                        else:
                            known_authors[first_last] = create_author(person_arr[1][0].upper(), len(known_authors) + 3, 1, [person])
                            new_author(sheet_authors, known_authors[first_last], first_last.split(' '), person)
                            time.sleep(1)

                    # elif person_arr[0][0] + ' ' + person_arr[2] in known_authors: 
                    #     known_authors[first_last] = known_authors.pop(person_arr[0][0] + ' ' + person_arr[2])
                    #     known_authors[first_last]['Pub_Number'] += 1

                    #     known_authors[first_last]['Aliases'].append(person)

                    #     update_author_general(sheet_authors, person_arr[0], known_authors[first_last]['Sheets Index'], 1)
                    #     update_author_general(sheet_authors, known_authors[first_last]['Pub_Number'], known_authors[first_last]['Sheets Index'], 5)
                    #     update_author_general(sheet_authors, known_authors[first_last]['Aliases'], known_authors[first_last]['Sheets Index'], 4)
                    #     time.sleep(3)

                    else:
                        found = False

                        for key in known_authors:
                            split = str(key).split(' ')

                            if person_arr[len(person_arr) - 1] == split[len(split) - 1] and person_arr[0][0].upper() == split[0][0].upper():
                                array_key = key

                                if(person != key):
                                    if(person > key):
                                        array_key = person
                                        known_authors[array_key] = known_authors.pop(key)

                                        if(key not in known_authors[array_key]['Aliases']):
                                            known_authors[array_key]['Aliases'].append(key)
                                            update_author_general(sheet_authors, ', '.join(known_authors[array_key]['Aliases']), known_authors[array_key]['Sheets Index'], 4)
                                            time.sleep(1)
                                        
                                        update_author_general(sheet_authors, person_arr[0], known_authors[array_key]['Sheets Index'], 1)
                                        time.sleep(1)
                                    else:
                                        if(person not in known_authors[array_key]['Aliases']):
                                            known_authors[array_key]['Aliases'].append(person)
                                            update_author_general(sheet_authors, ', '.join(known_authors[array_key]['Aliases']), known_authors[array_key]['Sheets Index'], 4)
                                            time.sleep(1)

                                found = True
                                known_authors[array_key]['Pub_Number'] += 1
                                update_author_general(sheet_authors, known_authors[array_key]['Pub_Number'], known_authors[array_key]['Sheets Index'], 5)
                                time.sleep(1)
                                break

                        if not found:
                            known_authors[person] = create_author(person_arr[1][0].upper(), len(known_authors) + 3, 1, [person])
                            new_author(sheet_authors, known_authors[person], person_arr, person)
                            time.sleep(1)

                elif len(person_arr) == 3:
                    first_last = person_arr[0] + ' ' + person_arr[2] 

                    if(first_last in known_authors):

                        if(known_authors[first_last]['Middle Initial'] == '' or known_authors[first_last]['Middle Initial'] == person_arr[1][0].upper()):
                            known_authors[first_last]['Middle Initial'] = person_arr[1][0].upper()

                            update_author_general(sheet_authors, person_arr[1][0].upper(), known_authors[first_last]['Sheets Index'], 2)
                            
                            known_authors[first_last]['Pub_Number'] += 1

                            update_author_general(sheet_authors, known_authors[first_last]['Pub_Number'], known_authors[first_last]['Sheets Index'], 5)

                            time.sleep(2)
                        else:
                            known_authors[first_last] = create_author(person_arr[1][0].upper(), len(known_authors) + 3, 1, [person])
                            new_author(sheet_authors, known_authors[first_last], first_last.split(' '), person)
                            time.sleep(1)
                            pass # found a new author

                    # elif person_arr[0][0] + ' ' + person_arr[2] in known_authors: 
                    #     known_authors[first_last] = known_authors.pop(person_arr[0][0] + ' ' + person_arr[2])
                    #     known_authors[first_last]['Pub_Number'] += 1
                    #     update_author_general(sheet_authors, person_arr[0], known_authors[first_last]['Sheets Index'], 1)
                    #     update_author_general(sheet_authors, known_authors[first_last]['Pub_Number'], known_authors[first_last]['Sheets Index'], 5)
                    #     time.sleep(2)

                    else:
                        known_authors[first_last] = create_author(person_arr[1][0].upper(), len(known_authors) + 3, 1, [person])
                        new_author(sheet_authors, known_authors[first_last], first_last.split(' '), person)
                        time.sleep(1)

                elif len(person_arr) == 2:
                    if person_arr[0].replace('.','').isupper() or person_arr[1].replace('.','').isupper():
                        print(person + " is not full")
                        
                        if(person_arr[0].replace('.','').isupper()):
                            upper = 0
                            full = 1
                        else:
                            upper = 1
                            full = 0

                        author_found = False

                        for key in known_authors:
                            split = str(key).split(' ')

                        #     if split[1] == person_arr[full] and split[0][0] == person_arr[upper][0]:
                        #         author_found = True
                        #         known_authors[key]['Pub_Number'] += 1
                        #         update_author_general(sheet_authors, known_authors[key]['Pub_Number'], known_authors[key]['Sheets Index'], 5)
                        #         time.sleep(1)
                        #         alias_found = False

                        for alias in known_authors[key]['Aliases']:

                            if(alias == person):
                                alias_found = True
                                
                        #         if not alias_found:
                        #             print(person)
                        #             print(known_authors[key]['Aliases'])
                        #             known_authors[key]['Aliases'].append(person)
                        #             update_author_general(sheet_authors, ', '.join(known_authors[key]['Aliases']), known_authors[key]['Sheets Index'], 4)
                        #             time.sleep(1)
                        #         break

                        # if not author_found:
                        if not alias_found:
                            if(len(person_arr[upper].replace('.','')) > 1):
                                known_authors[person_arr[upper][0] + ' ' + person_arr[full]] = create_author(person_arr[upper].replace('.','')[1], len(known_authors) + 3, 1, [person])
                            else:
                                known_authors[person_arr[upper][0] + ' ' + person_arr[full]] = create_author('', len(known_authors) + 3, 1, [person])
                                
                            new_author(sheet_authors, known_authors[person_arr[upper][0] + ' ' + person_arr[full]], person_arr, person)
                            time.sleep(1)

                    else:
                        found = False

                        for key in known_authors:
                            split = str(key).split(' ')

                            if person_arr[1] == split[1] and person_arr[0][0].upper() == split[0][0].upper():
                                array_key = key

                                if(person != key):
                                    if(person > key):
                                        array_key = person
                                        known_authors[array_key] = known_authors.pop(key)

                                        if(key not in known_authors[array_key]['Aliases']):
                                            known_authors[array_key]['Aliases'].append(key)
                                            update_author_general(sheet_authors, ', '.join(known_authors[array_key]['Aliases']), known_authors[array_key]['Sheets Index'], 4)
                                            time.sleep(1)
                                        
                                        update_author_general(sheet_authors, person_arr[0], known_authors[array_key]['Sheets Index'], 1)
                                        time.sleep(1)
                                    else:
                                        if(person not in known_authors[array_key]['Aliases']):
                                            known_authors[array_key]['Aliases'].append(person)
                                            update_author_general(sheet_authors, ', '.join(known_authors[array_key]['Aliases']), known_authors[array_key]['Sheets Index'], 4)
                                            time.sleep(1)

                                found = True
                                known_authors[array_key]['Pub_Number'] += 1
                                update_author_general(sheet_authors, known_authors[array_key]['Pub_Number'], known_authors[array_key]['Sheets Index'], 5)
                                time.sleep(1)
                                break

                        if not found:
                            known_authors[person] = create_author('', len(known_authors) + 3, 1, [person])
                            new_author(sheet_authors, known_authors[person], person_arr, person)
                            time.sleep(1)

import pprint

with open('Authors_dict.txt', 'w', encoding = 'utf-8') as dict:
    pprint.pprint(known_authors, dict)

####################################################################################################################### Write to Google Sheets
final_list = []

# format list of lists to be one huge list for printing purposes
for array in list_of_lists:
    final_list += array

for publication in final_list:
    publication[2] = str('=display_authors_one_line(CONCATENATE("B", TO_TEXT(ROW(B' + str(disp_authors_start) + '))))')
    disp_authors_start += 1

# write results to the google sheet
header = ['Title', 'Raw Authors', 'Display Authors', 'Trainees', 'Journal', 'Pubtype', 'Issue', 'Volume', 'Pages', 'DOI', 'PMID', 'PMCID', 'ISSN', 'Pubdate', 'Funding', 'Abstract', 'Notes']

titles = sheet.row_values(1)
time.sleep(1)

if(str(titles) == "[]"):
    sheet.insert_row(header,1)
    time.sleep(1)

# pprint.pprint(final_list)
sheet.append_row(['Start of week ' + str(date.today() - timedelta(weeks=week)) + ' to ' + str(date.today())])
time.sleep(1)
sheet.append_rows(final_list, value_input_option='USER_ENTERED')
time.sleep(1)
sheet.append_row(['End of week ' + str(date.today() - timedelta(weeks=week)) + ' to ' + str(date.today())])
time.sleep(1)
sheet.columns_auto_resize(0, 9)
