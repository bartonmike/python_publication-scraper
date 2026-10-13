from authors_utils import process_publication_authors, write_known_author_updates
from gsheets_utils import init_gsheets


def _replace_author_semicolons(rows):
    updates = []
    for row_number, row in enumerate(rows, start=1):
        if not row:
            continue

        title = row[0].strip() if row else ''
        if title.lower() == 'title' or title.lower().startswith(('start of week', 'end of week')):
            continue

        raw_authors = row[1] if len(row) > 1 else ''
        display_authors = row[2] if len(row) > 2 else ''
        normalized_raw_authors = raw_authors.replace(';', ',')
        normalized_display_authors = display_authors.replace(';', ',')

        if len(row) > 1:
            row[1] = normalized_raw_authors
        if len(row) > 2:
            row[2] = normalized_display_authors

        if normalized_raw_authors != raw_authors or normalized_display_authors != display_authors:
            while len(row) < 3:
                row.append('')

            updates.append({
                'range': f'B{row_number}:C{row_number}',
                'values': [[row[1], row[2]]]
            })

    return updates


def _publication_author_lists(rows):
    publications = []
    for row in rows:
        if not row:
            continue

        title = row[0].strip() if row else ''
        if title.lower() == 'title' or title.lower().startswith(('start of week', 'end of week')):
            continue

        raw_authors = row[1].strip() if len(row) > 1 else ''
        if title or raw_authors:
            publications.append(['', raw_authors])

    if not publications:
        raise ValueError('The "Publication Info" worksheet contains no publication rows.')

    return [publications]


def main():
    client, _, _, _, _ = init_gsheets()
    with open('sheet_name.txt', 'r', encoding='utf-8') as sheet_name_file:
        sheet_name = sheet_name_file.readline().strip()

    publications_sheet = client.open(sheet_name).worksheet('Publication Info')
    rows = publications_sheet.get_all_values()
    author_column_updates = _replace_author_semicolons(rows)
    if author_column_updates:
        publications_sheet.batch_update(author_column_updates, value_input_option='RAW')

    publication_lists = _publication_author_lists(rows)
    _, known_authors_sheet, _, author_updates = process_publication_authors(
        publication_lists,
        client,
        rebuild=True
    )
    write_known_author_updates(known_authors_sheet, author_updates)
    print(f'Processed {len(publication_lists[0])} publications from "Publication Info".')


if __name__ == '__main__':
    main()
