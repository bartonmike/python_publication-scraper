from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException
from fake_useragent import UserAgent
from bs4 import BeautifulSoup
import requests
import time
from datetime import date, timedelta
import re
from pubmed_utils import search_pubmed_by_title, search_crossref_by_title, build_pub, find_pmid, fetch_summary

def scrape_google_scholar(people, trainees, week, headers, codes, RunLog, pubmed_search_func, build_pub_func):
    ua = UserAgent()
    headers['User-Agent'] = ua.random
    options = webdriver.ChromeOptions()
    options.add_experimental_option("excludeSwitches", ['enable-logging'])
    options.add_argument('--headless=new')
    options.add_argument('--log-level=3')
    options.set_capability("browserVersion", "117")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument(f'--user-agent={headers["User-Agent"]}')
    driver = webdriver.Chrome(options = options)

    google_master = []

    if people is None:
        driver.close()
        driver.quit()
        return google_master

    for term in people:
        google_scholar = []
        processed_links = set()
        out_of_date = False
        url = 'https://scholar.google.com/citations?view_op=list_works&hl=en&hl=en&user=' + term[0] + '&sortby=pubdate'
        driver.get(url)
        try:
            WebDriverWait(driver, 30).until(EC.presence_of_element_located((By.CLASS_NAME, "gsc_a_at")))
        except TimeoutException:
            pass

        while True:
            elements = driver.find_elements(By.CLASS_NAME, "gsc_a_at")
            found_new = False
            for element in elements:
                title_text = element.text.strip() or element.get_attribute('textContent') or ''
                if not title_text or title_text in processed_links:
                    continue
                pub_entry = None
                # Try PubMed first
                try:
                    pub_entry = pubmed_search_func(title_text)
                except Exception:
                    pub_entry = None

                # If PubMed fails, try Crossref
                if pub_entry is None:
                    try:
                        pub_entry = search_crossref_by_title(title_text, trainees)
                    except Exception:
                        pub_entry = None

                if pub_entry is not None:
                    processed_links.add(title_text)
                    min_date = date.today() - timedelta(weeks=week)
                    min_year = min_date.year
                    pubdate_str = pub_entry[13] or ''
                    year_match = re.search(r"(\d{4})", pubdate_str)
                    withindate = True
                    if year_match:
                        try:
                            pub_year = int(year_match.group(1))
                            if pub_year < min_year:
                                withindate = False
                        except Exception:
                            pass
                    if withindate:
                        google_scholar.append(pub_entry)
                    else:
                        out_of_date = True
                        break
                    continue

                link = element.get_attribute('href')
                if not link or link in processed_links:
                    continue
                found_new = True
                processed_links.add(link)
                # Rate limiting: add delay before making external requests to avoid rate-limiting
                time.sleep(2)
                p = requests.get(link, headers = headers, timeout=30)
                soup = BeautifulSoup(p.text, features = "html.parser")
                parts = soup.find_all(attrs= {'class':'gs_scl'})
                parent = soup.find(attrs={'class':'gsc_oci_title_link'})
                try:
                    title = parent.get_text()
                    url = parent['href']
                except Exception:
                    title = title_text

                Authors = ''
                Publication_date = ''
                Journal = ''
                for part in parts:
                    text = part.find(attrs={'class':'gsc_oci_field'})
                    values = part.find(attrs={'class': 'gsc_oci_value'})
                    if text and values:
                        if 'Authors' in text.get_text():
                            Authors = values.get_text()
                        if 'Publication date' in text.get_text():
                            Publication_date = values.get_text()
                        if 'Journal' in text.get_text():
                            Journal = values.get_text()

                date_today = str(date.today() - timedelta(weeks=week)).replace('-','/')
                pubcount = Publication_date.count('/') if Publication_date else 0
                withindate = False
                if Publication_date == '' or Publication_date is None:
                    withindate = False
                elif pubcount == 0:
                    if date_today.split('/')[0] <= Publication_date:
                        withindate = True
                elif pubcount == 1:
                    pub_try = Publication_date
                    if len(pub_try) < 7:
                        pub_try = pub_try.replace('/','/0') + '/00'
                    if date_today.replace('/','') <= pub_try.replace('/',''):
                        withindate = True
                elif pubcount == 2:
                    Split_pub = Publication_date.split('/')
                    if len(Split_pub) >= 3:
                        if len(Split_pub[1]) < 2:
                            Split_pub[1] = '0' + Split_pub[1]
                        if len(Split_pub[2]) < 2:
                            Split_pub[2] = '0' + Split_pub[2]
                        Split_pub = ''.join(Split_pub)
                        if Split_pub >= date_today.replace('/',''):
                            withindate = True

                if withindate:
                    link_cr = 'https://search.crossref.org/search/works?q=' + title + '&from_ui=yes'
                    time.sleep(1.5)  # Rate limit Crossref web search
                    p = requests.get(link_cr, headers = headers)
                    soup = BeautifulSoup(p.text, features = "html.parser")
                    doi_text = ''
                    doi_node = soup.find(attrs={'class':'item-links'})
                    if doi_node:
                        doi_text = doi_node.text.strip().split('\n')[0].replace('https://doi.org/','')
                    if doi_text:
                        b = build_pub_func(doi_text, trainees, Authors, Journal, Publication_date)
                        if b:
                            google_scholar.append(b)
                else:
                    out_of_date = True
                    break

            if out_of_date:
                break

            if found_new:
                try:
                    show_more = driver.find_element(By.ID, 'gsc_bpf_more')
                    if show_more.is_displayed() and show_more.is_enabled():
                        try:
                            show_more.click()
                            time.sleep(2)  # Increased delay after "Show more" clicks
                            continue
                        except Exception:
                            driver.execute_script("arguments[0].click();", show_more)
                            time.sleep(2)  # Increased delay after "Show more" clicks
                            continue
                    else:
                        break
                except Exception:
                    break
            else:
                break

        RunLog.append_row(['Number of Hits for: ' + term[1], str(len(google_scholar)), str(date.today())])
        RunLog.columns_auto_resize(0,3)
        google_master.append(google_scholar)

    driver.close()
    driver.quit()
    return google_master