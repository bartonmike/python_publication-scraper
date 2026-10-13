import json
import html
import os
import time
import unicodedata
from pathlib import Path

from Bio import Entrez

from gsheets_utils import init_gsheets


PROJECT_DIR = Path(__file__).resolve().parent
SHEET_NAME_FILE = PROJECT_DIR / "sheet_name.txt"
ENTREZ_EMAIL = os.environ.get("PUBMED_EMAIL", "")
Entrez.email = ENTREZ_EMAIL


def _normalize_title(title):
    title = unicodedata.normalize("NFKD", html.unescape(title)).casefold()
    return "".join(
        character for character in title
        if character.isalnum() and not unicodedata.combining(character)
    )


def _is_empty(value):
    return not value or value.strip().lower() in {"none", "null"}


def _find_header_columns(rows):
    required_headers = {"title", "pubtype", "doi", "pmid", "pmcid"}
    for row_number, row in enumerate(rows, start=1):
        columns = {
            value.strip().lower(): column_number
            for column_number, value in enumerate(row, start=1)
        }
        if required_headers.issubset(columns):
            return row_number, columns
    raise ValueError(
        'Could not find a header row containing Title, Pubtype, DOI, PMID, and PMCID.'
    )


def find_pubmed_match(title):
    """Return matching PMID, DOI, and PMCID, or None if PubMed has no exact title match."""
    query_title = title.replace('"', " ").replace("–", "-").strip()
    print(f"Searching PubMed for title: {query_title}")
    if not query_title:
        return None

    expected_title = _normalize_title(title)
    checked_pmids = set()


    search_handle = Entrez.esearch(
        db="pubmed",
        term=query_title + '[Title]',
        retmax="100",
        retmode="xml",
    )
    search_results = Entrez.read(search_handle)
    time.sleep(0.34)

    # print(search_results)

    print(f"Found {len(search_results.get('IdList', []))} PubMed IDs for search term: {query_title}")

    for pmid in search_results.get("IdList", []):
        if pmid in checked_pmids:
            continue
        checked_pmids.add(pmid)

        summary_handle = Entrez.esummary(db="pubmed", id=pmid, retmode="json")
        summary = json.loads(summary_handle.read())
        time.sleep(0.34)
        record = summary.get("result", {}).get(pmid, {})
        if _normalize_title(record.get("title", "")) != expected_title:
            continue

        identifiers = {"doi": "", "pmc": ""}
        for article_id in record.get("articleids", []):
            id_type = article_id.get("idtype")
            if id_type in identifiers:
                identifiers[id_type] = article_id.get("value", "")

        return pmid, identifiers["doi"], identifiers["pmc"]

    return None


def main():
    with SHEET_NAME_FILE.open("r", encoding="utf-8") as sheet_name_file:
        sheet_name = sheet_name_file.readline().strip()
    if not sheet_name:
        raise ValueError(f"{SHEET_NAME_FILE.name} does not contain a Google Sheet name.")

    client, _, _, _, _ = init_gsheets(
        client_key_path=str(PROJECT_DIR / "client_key.json"),
        sheet_name_file=str(SHEET_NAME_FILE),
    )
    worksheet = client.open(sheet_name).worksheet("Site Publications")
    rows = worksheet.get_all_values()
    header_row, columns = _find_header_columns(rows)

    updates = []
    changed_rows = []
    unmatched_articles = []
    for row_number, row in enumerate(rows[header_row:], start=header_row + 1):
        title = row[columns["title"] - 1].strip() if len(row) >= columns["title"] else ""
        pubtype = row[columns["pubtype"] - 1].strip() if len(row) >= columns["pubtype"] else ""
        pmid = row[columns["pmid"] - 1] if len(row) >= columns["pmid"] else ""
        if pubtype.lower() != "journal_article" or not _is_empty(pmid) or not title:
            continue

        match = find_pubmed_match(title)
        if match is None:
            unmatched_articles.append((row_number, title))
            continue

        matched_pmid, doi, pmcid = match
        updates.append({
            "range": f"{_column_letter(columns['pmid'])}{row_number}",
            "values": [[matched_pmid]],
        })

        for header, value in (("doi", doi), ("pmcid", pmcid)):
            current_value = row[columns[header] - 1] if len(row) >= columns[header] else ""
            if value and _is_empty(current_value):
                updates.append({
                    "range": f"{_column_letter(columns[header])}{row_number}",
                    "values": [[value]],
                })
        changed_rows.append(row_number)

    if updates:
        worksheet.batch_update(updates, value_input_option="RAW")

    print(f"Updated {len(changed_rows)} article(s).")
    if changed_rows:
        print("Changed sheet rows: " + ", ".join(str(row) for row in changed_rows))
    else:
        print("Changed sheet rows: none")

    print(f"Articles without a matching PubMed PMID: {len(unmatched_articles)}.")
    if unmatched_articles:
        print("Unmatched articles (row: title):")
        for row_number, title in unmatched_articles:
            print(f"  {row_number}: {title}")
    else:
        print("Unmatched articles: none")


def _column_letter(column_number):
    letters = ""
    while column_number:
        column_number, remainder = divmod(column_number - 1, 26)
        letters = chr(65 + remainder) + letters
    return letters


if __name__ == "__main__":
    main()
