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
                       lxml_html_clean <br/><br/>
               files and other requirements: <br/>
                       client_key.json file in the same directory, containing info on a Google developer bot <br/>
                       create_publication_scraper.py is run, which creates all the needed txt files <br/>
                       bot whose details are contained in the client_key.json, has edit access to the given Google sheet <br/>

All Python modules for the program are contained in the requirements.txt file <br/>
before running the program, open the Windows command line, CD into the directory of the program, and type this command followed by pressing enter: <br/>
```
pip install -r requirements.txt
```
              
 Inputs: Google scholar profile codes and names, trainee names, weeks ago one wants to search.<br/>


 Outputs: All outputs are on a Google sheet<br/>
           Table: Title, Authors, Journal, DOI, PMID, PMCID, Pubdate, Funding, Trainees, Notes<br/>
           Changelog: removed and added for the inputs <br/>
           Run Log: entries found for each person/funding

# Pubmed ID Publication Scraper

Description: Instead of utilizing google scholar profile codes and funding codes like the original, this program directly pulls information based on the pubmed ID codes that it is given. This program does not "look" for publications based off of a certain timeframe,
             it is given the publication codes and just fetches the publication data through apis such as the Pubmed Entrez api and the Crossref Habanero api. It writes the publications in the same format as the original, populating the first sheet in a google sheet with the 
             fetched data. Additionally, this program provides a "known author" sheet that lists all the unique authors present in the publications given, listed as First Name, Middle Initial, Last Name, Aliases, # of Publications, and DOI's that the author is present in.
             All outputs of the original publication scraper are present in this one. 

 Requirements: <br/>
               Python Modules: <br/>
                       biopython <br/>
                       beautifulsoup4 <br/>
                       habanero <br/>
                       gspread <br/>
                       oauth2client <br/>
                       lxml_html_clean <br/><br/>
               files and other requirements: <br/>
                       client_key.json file in the same directory, containing info on a Google developer bot <br/>
                       create_publication_scraper.py is run, which creates all the needed txt files <br/>
                       bot whose details are contained in the client_key.json, has edit access to the given Google sheet <br/>

All Python modules for the program are contained in the requirements.txt file <br/>
before running the program, open the Windows command line, CD into the directory of the program, and type this command followed by pressing enter: <br/>
```
pip install -r requirements.txt
```
              
 Inputs: Google scholar profile codes and names, trainee names, weeks ago one wants to search.<br/>


 Outputs: All outputs are on a Google sheet<br/>
           Table: Title, Authors, Journal, DOI, PMID, PMCID, Pubdate, Funding, Trainees, Notes<br/>
           Changelog: removed and added for the inputs <br/>
           Run Log: entries found for each person/funding <br/>
           Known Authors: First Name, Middle Initial, Last Name, Aliases, # of Publications, and DOI's
