# About
This Repository stores two types of Python publication scrapers: the first is the "original" publication scraper, and the second is the PubMed ID publication scraper

# "Original" Publication Scraper

 Description: This program takes in Google Scholar profile codes, along with names, and pulls data from the 
              publications within a certain number of weeks. It utilizes biopython and habanero to fetch the Entrez
              and Crossref libraries, respectively, using the PMID and name of the article to fetch the data. It also uses the three superfund funding codes
              to locate publications in the PubMed database. It writes the data it fetches to a Google sheet whose name was initialized in the create_paper_scraper.py 
              program and is written to a txt file named sheet_name.txt. This program generates a list of unique publications for all entered individuals. It also has added trainee 
              functionality, where one can provide a list of trainees, and the program will automatically search for those names in the publication authors, putting any matches onto the sheet.

 Requirements: <br/>
               Python Modules: <br/>
                       biopython <br/>
                       fake_useragent <br/>
                       beautifulsoup4 <br/>
                       habanero <br/>
                       selenium <br/>
                       gspread <br/>
                       oauth2client <br/>
                       lxml_html_clean <br/>
               files and other requirements: <br/>
                       client_key.json file in the same directory, containing info on a Google developer bot <br/>
                       create_paper_scraper.py is run, which creates all the needed txt files <br/>
                       bot whose details are contained in the client_key.json, has edit access to the given Google sheet <br/>
              
 Inputs: Google scholar profile codes and names, trainee names, weeks ago one wants to search.<br/>


 Outputs: All outputs are on a Google sheet<br/>
           Table: Title, Authors, Journal, DOI, PMID, PMCID, Pubdate, Funding, Trainees, Notes<br/>
           Changelog: removed and added for the inputs <br/>
           Run Log: entries found for each person/funding
