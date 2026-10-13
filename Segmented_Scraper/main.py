from gsheets_utils import init_gsheets
from pubmed_utils import (search, fetchData, checkDupe, build_list, compare_people, search_pubmed_by_title,
                           build_pub, find_pmid, find_authors, fetch_summary, fetch_abstract, find_matches, resolve_google_entries)
from scholar_utils import scrape_google_scholar
from authors_utils import process_known_authors, process_publication_authors, write_known_author_updates
from fake_useragent import UserAgent
from datetime import date, timedelta
import time
import re
import json

def main():
    client, sheet, prompts, ChangeLog, RunLog = init_gsheets()

    # load stored files
    with open("people.txt", 'r', encoding = 'utf-8') as file:
        old_people = file.readlines()
    with open("trainees.txt", 'r', encoding = 'utf-8') as file:
        old_trainees = file.readlines()
    with open("date.txt", 'r', encoding = 'utf-8') as file:
        old_week = file.readlines()

    old_people = [p.strip() for p in old_people]
    old_trainees = [t.strip() for t in old_trainees]
    old_week = [w.strip() for w in old_week]

    people = prompts.cell(2,1).value
    if people != None:
        people = [p.strip() for p in people.split('\n')]
        people = [p.split(', ') for p in people]
    else:
        people = None

    trainees = prompts.cell(4,1).value
    if trainees != None:
        trainees = [t.strip().split(' ') for t in trainees.split('\n')]

    week = int(prompts.cell(6,1).value)

    people_diff1, people_diff2 = compare_people(people or [], old_people)
    train_diff1, train_diff2 = compare_people([' '.join(t) for t in (trainees or [])], old_trainees)

    str_diff1 = ''
    str_diff2 = ''
    week_diff = False
    if(old_week == []):
        old_week.append(None)
    if(old_week[0] != str(week)):
        week_diff = True
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

    if(str_diff1.strip() != '' and str_diff2.strip() != ''):
        ChangeLog.append_row([str_diff1 + '\n' + str_diff2, str(date.today())])
        ChangeLog.columns_auto_resize(0, 2)
    elif(str_diff1.strip() != ''):
        ChangeLog.append_row([str_diff1, str(date.today())])
        ChangeLog.columns_auto_resize(0, 2)
    elif(str_diff2.strip() != ''):
        ChangeLog.append_row([str_diff2, str(date.today())])
        ChangeLog.columns_auto_resize(0, 2)

    if(people_diff1.strip() != '' or people_diff2.strip() != ''):
        with open("people.txt", "w") as c:
            c.write('\n'.join([', '.join(p) if isinstance(p, list) else p for p in (people or [])]))

    if(train_diff2.strip() != '' or train_diff1.strip() != ''):
        with open("trainees.txt", "w") as c:
            c.write('\n'.join([' '.join(t) for t in (trainees or [])]))

    if(week_diff):
        with open("date.txt", "w") as c:
            c.write(str(week))

    # PUBMED search
    ua = UserAgent()
    headers = {"User-Agent": str(ua.chrome)}

    try:
        codes = prompts.cell(2,2).value.split('\n')
    except Exception:
        codes = []

    pubmed_array = []
    for code in codes:
        results = search(code, week)
        code_array = fetchData(results, code.replace("[Grants and Funding]",''), trainees, headers, RunLog)
        pubmed_array = checkDupe(pubmed_array, code_array)

    list_of_lists = [pubmed_array]

    # Handle custom PMIDs (copied from original logic)
    try:
        custom_pmids_raw = prompts.cell(6,2).value
    except Exception:
        custom_pmids_raw = None

    if custom_pmids_raw:
        pmid_candidates = re.split(r'[\n,;]+', str(custom_pmids_raw).strip())
        custom_entries = []
        for candidate in pmid_candidates:
            pmid = ''.join(ch for ch in candidate if ch.isdigit())
            if not pmid:
                continue
            try:
                resp = json.loads(fetch_summary(pmid))
            except Exception:
                continue

            r = resp.get('result', {}).get(pmid, {})
            if not r:
                continue

            doi = ''
            pmcid = ''
            for aid in r.get('articleids', []):
                if aid.get('idtype') == 'doi':
                    doi = aid.get('value', '')
                if aid.get('idtype') == 'pmc':
                    pmcid = aid.get('value', '')

            if doi == '':
                try:
                    doi = str(r['references'][0]['refsource']).split('doi: ')[1]
                except Exception:
                    doi = ''

            authors_fullname = ''
            if doi:
                try:
                    authors_fullname = find_authors(doi)
                except Exception:
                    authors_fullname = ''

            if authors_fullname == 'Authors Not Found' or authors_fullname == '':
                try:
                    authors_list = r.get('authors', [])
                    authors_fullname = ', '.join([a.get('name') if isinstance(a, dict) else str(a) for a in authors_list])
                except Exception:
                    authors_fullname = ''

            try:
                journal = r.get('pubtype', [''])[0]
            except Exception:
                journal = ''

            try:
                abstract = str(fetch_abstract(pmid))
            except Exception:
                abstract = ''

            entry = [
                r.get('title', ''),
                authors_fullname,
                '',
                find_matches(trainees, authors_fullname),
                journal,
                r.get('source', ''),
                r.get('volume', ''),
                r.get('issue', ''),
                r.get('pages', ''),
                doi,
                pmid,
                pmcid,
                r.get('issn', ''),
                r.get('pubdate', ''),
                'None',
                abstract,
                ''
            ]

            custom_entries.append(entry)

        if custom_entries:
            pubmed_array = checkDupe(custom_entries, pubmed_array)
            RunLog.append_row(['Number of Custom PMIDs Processed', str(len(custom_entries)), str(date.today())])
            RunLog.columns_auto_resize(0,3)

    # Google Scholar scraping
    google_lists = scrape_google_scholar(people, trainees, week, headers, codes, RunLog,
                                        lambda title: search_pubmed_by_title(title, trainees, headers),
                                        lambda doi, trainees, a, b, c: build_pub(doi, trainees, a, b, c))

    # for each person's google scholar results: dedupe against existing lists, resolve pmid/pmcid and funding codes, then append
    for gl in google_lists:
        person_array = build_list(list_of_lists, gl)
        person_array = resolve_google_entries(person_array, codes, headers)
        list_of_lists.append(person_array)

    # Author processing stages Known Authors changes for the final sheet write.
    known_authors, sheet_authors, Display_authors, author_updates = process_publication_authors(list_of_lists, client)

    # final write to sheet logic (trimmed to critical steps)
    existing_dois = set()
    existing_pmids = set()

    def _normalize_doi(s):
        if not s:
            return ''
        s = s.strip()
        if s.lower().startswith('http'):
            s = s.rstrip('/').split('/')[-1]
        return s.lower()

    try:
        existing_rows = sheet.get_all_values()
        for row in existing_rows:
            if not row:
                continue
            joined = ' '.join(row).strip()
            if joined.startswith('Start of week') or joined.startswith('End of week'):
                continue
            doi = row[9].strip() if len(row) > 9 else ''
            pmid = row[10].strip() if len(row) > 10 else ''
            norm_doi = _normalize_doi(doi)
            if norm_doi and norm_doi.lower() not in ('none',''):
                existing_dois.add(norm_doi)
            if pmid and pmid.lower() not in ('none',''):
                existing_pmids.add(pmid)
    except Exception as e:
        print('Warning reading existing Publication Info:', e)

    final_list = []
    duplicates_skipped = 0
    for array in list_of_lists:
        for pub in array:
            doi = pub[9].strip() if len(pub) > 9 else ''
            pmid = pub[10].strip() if len(pub) > 10 else ''
            if _normalize_doi(doi) in existing_dois or (pmid and pmid in existing_pmids):
                duplicates_skipped += 1
                continue
            final_list.append(pub)

    if duplicates_skipped > 0:
        try:
            RunLog.append_row(['Duplicates skipped (already on sheet)', str(duplicates_skipped), str(date.today())])
            RunLog.columns_auto_resize(0,3)
        except Exception:
            pass

    disp_authors_start = len(sheet.col_values(3)) + 2
    for publication in final_list:
        publication[2] = str('=display_authors_one_line(CONCATENATE("B", TO_TEXT(ROW(B' + str(disp_authors_start) + '))))')
        disp_authors_start += 1

    header = ['Title', 'Raw Authors', 'Display Authors', 'Trainees', 'Journal', 'Pubtype', 'Issue', 'Volume', 'Pages', 'DOI', 'PMID', 'PMCID', 'ISSN', 'Pubdate', 'Funding', 'Abstract', 'Notes']
    titles = sheet.row_values(1)
    time.sleep(1)
    if(str(titles) == "[]"):
        sheet.insert_row(header,1)
        time.sleep(1)
    sheet.append_row(['Start of week ' + str(date.today() - timedelta(weeks=week)) + ' to ' + str(date.today())])
    time.sleep(1)
    sheet.append_rows(final_list, value_input_option='USER_ENTERED')
    time.sleep(1)
    sheet.append_row(['End of week ' + str(date.today() - timedelta(weeks=week)) + ' to ' + str(date.today())])
    time.sleep(1)
    write_known_author_updates(sheet_authors, author_updates)
    sheet.columns_auto_resize(0, 9)

if __name__ == '__main__':
    main()
