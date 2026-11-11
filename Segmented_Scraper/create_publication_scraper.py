"""Initialize the Google Sheet and local files for the segmented scraper."""

from pathlib import Path

import gspread
from oauth2client.service_account import ServiceAccountCredentials


PROJECT_DIR = Path(__file__).resolve().parent
CLIENT_KEY_FILE = PROJECT_DIR / "client_key.json"
SHEET_NAME_FILE = PROJECT_DIR / "sheet_name.txt"


def main():
    scope = [
        "https://www.googleapis.com/auth/spreadsheets",
        "https://www.googleapis.com/auth/drive",
    ]
    credentials = ServiceAccountCredentials.from_json_keyfile_name(
        str(CLIENT_KEY_FILE), scope
    )
    client = gspread.authorize(credentials)

    sheet_name = input("What is the name of your Google Sheet? ").strip()
    if not sheet_name:
        raise ValueError("The Google Sheet name cannot be empty.")

    spreadsheet = client.open(sheet_name)

    publications = spreadsheet.sheet1
    publications.update_title("Publication Info")
    publications.append_row([
        "Title",
        "Raw Authors",
        "Display Authors",
        "Trainees",
        "Journal",
        "Pubtype",
        "Issue",
        "Volume",
        "Pages",
        "DOI",
        "PMID",
        "PMCID",
        "ISSN",
        "Pubdate",
        "Funding",
        "Abstract",
        "Notes",
    ])

    prompts = spreadsheet.add_worksheet("Prompts", rows=200, cols=10)
    for row, column, value in (
        (1, 1, "Google Scholar Profiles (Code, Name)"),
        (3, 1, "Trainees"),
        (5, 1, "Weeks Ago"),
        (1, 2, "Funding Codes"),
        (5, 2, "Custom PMIDs"),
    ):
        prompts.update_cell(row=row, col=column, value=value)

    change_log = spreadsheet.add_worksheet("Change Log", rows=200, cols=10)
    change_log.append_row(["Changes", "Date"])

    run_log = spreadsheet.add_worksheet("Run Log", rows=200, cols=10)
    run_log.append_row(["Title", "Message", "Date"])

    known_authors = spreadsheet.add_worksheet("Known Authors", rows=200, cols=10)
    known_authors.append_row([
        "First",
        "Middle",
        "Last",
        "Aliases",
        "# of Publications",
    ])

    for file_name in ("people.txt", "trainees.txt", "date.txt"):
        (PROJECT_DIR / file_name).touch(exist_ok=True)
    SHEET_NAME_FILE.write_text(sheet_name, encoding="utf-8")

    print(f"Initialized {sheet_name!r} and saved its name to {SHEET_NAME_FILE.name}.")


if __name__ == "__main__":
    main()
