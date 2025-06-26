# About
This Repository stores two types of python publication scrapers, the first is the "original" publication scraper, and the second is the pubmedid publication scraper

# "Original" Publication Scraper

 Description: This program takes in google scholar profile codes, along with names, and pulls data from the 
              publications within a certain amount of weeks. It utilizes the biopython and habanero to fetch the entrez
              and Crossref libraries respectively, using the PMID and name of the article to fetch the data. It also uses the three superfund funding codes
              to locate publications in the pubmed database. It writes the data it fetches to a google sheet thats name was initialized in the create_paper_scraper.py 
              program and written to a txt file named sheet_name.txt. This program develops a list of unique publications across all entered people. It also has added trainee 
              functionality, where one can provide a list of trainees and the program will automatically search for those names in the publication authors, putting any matches onto the sheet.

 Requirements: 
               Python Modules: 
                       biopython
                       fake_useragent
                       beautifulsoup4
                       habanero
                       selenium
                       gspread
                       oauth2client
                       lxml_html_clean

               files and other requirements:
                       client_key.json file in same directory, containing info on a google developer bot
                       create_paper_scraper.py ran, which creates all the needed txt files
                       bot whose details are containted in the client_key.json has edit access to given google sheet
              
 Inputs: Google scholar profile codes and names, trainee names, weeks ago one wants to search.


 Outputs: All outputs are on google sheet
           Table: Title, Authors, Journal, DOI, PMID, PMCID, Pubdate, Funding, Trainees, Notes
           Changelog: removed and added for the inputs 
           Run Log: entries found for each person/funding
