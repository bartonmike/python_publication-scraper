import gspread
from oauth2client.service_account import ServiceAccountCredentials
from datetime import date

def init_gsheets(client_key_path='client_key.json', sheet_name_file='sheet_name.txt'):
    scope = [
        'https://www.googleapis.com/auth/spreadsheets',
        'https://www.googleapis.com/auth/drive'
    ]
    creds = ServiceAccountCredentials.from_json_keyfile_name(client_key_path,scope)
    client = gspread.authorize(creds)

    with open(sheet_name_file,'r', encoding='utf-8') as f:
        sheet_name = f.readlines()[0].strip()

    sheet = client.open(sheet_name).sheet1
    prompts = client.open(sheet_name).worksheet('Prompts')
    ChangeLog = client.open(sheet_name).worksheet('Change Log')
    RunLog = client.open(sheet_name).worksheet('Run Log')

    return client, sheet, prompts, ChangeLog, RunLog
