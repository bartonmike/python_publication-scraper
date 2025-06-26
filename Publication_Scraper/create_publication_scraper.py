# Name: create_paper_scraper.py
#
# Description: Initialization program for paper_scraper.py, initializes google sheet and creates txt files
#
# See paper_scraper.py for more details


#gspread imports
import gspread # pip3 install gspread
from oauth2client.service_account import ServiceAccountCredentials # pip3 install oauth2client

scope = [
    'https://www.googleapis.com/auth/drive',
    'https://www.googleapis.com/auth/drive.file'
    ]
file_name = 'client_key.json'
creds = ServiceAccountCredentials.from_json_keyfile_name(file_name,scope)
client = gspread.authorize(creds)

sheet_name = input("What is the name of your google sheet?")

sheet = client.open(sheet_name)

pub_info = sheet.sheet1
pub_info.update_title("Publication Info")

prompts = sheet.add_worksheet("Prompts", 200, 10)

header = ['Title', 'Raw Authors', 'Display Authors', 'Trainees', 'Pubtype', 'Journal', 'Issue', 'Volume', 'Pages', 'DOI', 'PMID', 'PMCID', 'ISSN', 'Pubdate', 'Funding', 'Abstract', 'Notes']

pub_info.append_row(header)

prompts.update_cell(row = 1,col = 1, value = 'Google Scholar Profiles (Code, Name)')
prompts.update_cell(row = 3,col = 1,value = 'Trainees')
prompts.update_cell(row = 5,col = 1,value = 'Weeks Ago')
prompts.update_cell(row = 1,col = 2,value = "Funding Codes")

changelog = sheet.add_worksheet("Change Log", 200, 10)
runlog = sheet.add_worksheet("Run Log", 200, 10)
known_authors = sheet.add_worksheet("Known Authors", 200, 10)

known_authors.append_row(['First', 'Middle', 'Last', 'Aliases', '# of Publications'])

changelog.append_row(['Changes', 'Date'])
runlog.append_row(['Title', 'Message', 'Date'])

f = open('people.txt', 'x')
f.close()

f = open('trainees.txt', 'x')
f.close()

f = open('date.txt', 'x')
f.close()

f = open('sheet_name.txt', 'x')
f.close()

f = open('sheet_name.txt', 'w')
f.write(sheet_name)
f.close()
