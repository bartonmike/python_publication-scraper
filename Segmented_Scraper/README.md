# Segmented Publication Scraper

Overview
- This package is a modular refactor of the original `publication_scraper.py`. It splits functionality into helpers for PubMed/Crossref, Google Scholar scraping, Google Sheets access, and author processing.

Requirements
- Python 3.9 — 3.12
- Packages: biopython, fake_useragent, beautifulsoup4, habanero, selenium, gspread, oauth2client, lxml_html_clean, tqdm, requests
- Chrome + matching ChromeDriver on PATH (for Selenium headless mode)
- Google service account credentials file: `client_key.json` in `Segmented_Scraper`, with edit access to the Google Sheet

Installation

1. Create and activate a virtual environment.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install --upgrade pip
pip install biopython fake-useragent beautifulsoup4 habanero selenium gspread oauth2client lxml_html_clean tqdm requests
```

2. Put `client_key.json` in `Segmented_Scraper` and create a Google Sheet shared with that service account.

3. Run the setup script from `Segmented_Scraper` and enter the name of the existing Google Sheet:

```powershell
python create_publication_scraper.py
```

This creates the `Publication Info`, `Prompts`, `Change Log`, `Run Log`, and `Known Authors` worksheets, and creates or updates the local `sheet_name.txt`, `people.txt`, `trainees.txt`, and `date.txt` files. Existing people, trainee, and date files are preserved.

4. Ensure `ChromeDriver` is installed and available on `PATH`.

Usage

Run from the `Segmented_Scraper` folder:

```powershell
python -m segmented_scraper.main
```

Reprocess authors for existing site publications

To rebuild `Authors_dict.txt` and the author counts/aliases in the `Known Authors` worksheet from the publications already in `Publication Info`, run:

```powershell
python process_site_publication_authors.py
```

This rebuild does not scrape for new publications or modify publication rows. Existing author publication counts are recalculated from the source worksheet.

Fill missing PubMed identifiers

To search PubMed for exact title matches for `journal_article` entries without a PMID in the `Publication Info` worksheet, then fill in the PMID and any available DOI/PMCID, run from this project folder:

```powershell
python update_missing_pubmed_ids.py
```

The script reads the Google Sheet name from `sheet_name.txt` and authenticates with `client_key.json`. It searches PubMed with an explicit AND of the title's significant words in the title field, then retries with distinctive word pairs if needed. This avoids PubMed phrase-parser failures around punctuation such as hyphens. Results are accepted only when their normalized title matches the sheet title. The script fills DOI/PMCID only when those cells are empty, leaves existing identifiers unchanged, and reports the number of changed articles and their sheet row numbers. It also lists the titles and row numbers of eligible articles for which it found no PubMed match. Set the `PUBMED_EMAIL` environment variable to an email address for NCBI Entrez requests.

What the modules do
- `pubmed_utils.py`: Entrez/Crossref helpers, data parsing, duplication checks, DOI/PMCID resolution.
- `scholar_utils.py`: Google Scholar scraping (Selenium) and conversion to publication entries.
- `gsheets_utils.py`: Google Sheets authorization and worksheet handles.
- `authors_utils.py`: Known-Authors processing and Display Authors generation.
- `main.py`: Orchestrates the flow (reads prompts, searches PubMed/funding codes, scrapes Google Scholar, resolves DOIs/PMIDs, writes publication data, then batches Known Authors changes in one sheet request).

Notes & Caveats
- The segmented version keeps the same behavior as the original script but is split for maintainability. Some non-critical logging and error handling were simplified; review `segmented_scraper/main.py` for integration points.
- Running Selenium headless on CI or servers may require additional Chrome options or a display driver.
- Be cautious: the script will modify the Google Sheet referenced in `sheet_name.txt`.

Further work
- Add unit tests and CI job to verify helper functions.
- Add a `requirements.txt` or `pyproject.toml` if you want pinned dependency versions.
