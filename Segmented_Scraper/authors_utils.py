import time
import pprint
import re

def create_author(initial, index, pub_num, aliases):
    author = {}
    author['Middle Initial'] = initial
    author['Sheets Index'] = index
    author['Pub_Number'] = pub_num
    author['Aliases'] = aliases
    return author

def _next_author_index(author_rows):
    return max(author_rows.keys(), default=1) + 2

def _normalize_first_name(value):
    if (len(value) == 2 and value.endswith('.')) or (len(value) == 4 and value.count('.') == 2) or (len(value) == 6 and value.count('.') == 3):
        return value[0]
    return value

def _parse_two_initials(value):
    match = re.fullmatch(r'([A-Z])\.?([A-Z])\.?', value)
    return match.groups() if match else None

def _set_author_value(author_rows, changed_rows, row, column, value):
    sheet_row = row - 1
    row_values = author_rows[sheet_row]
    while len(row_values) < 5:
        row_values.append('')
    if column == 1:
        value = _normalize_first_name(value)
    row_values[column - 1] = value
    changed_rows.add(sheet_row)

def _queue_new_author(author_rows, changed_rows, author_dict, name, aliases):
    sheet_row = author_dict['Sheets Index'] - 1
    alias_value = ', '.join(aliases) if isinstance(aliases, list) else aliases
    author_rows[sheet_row] = [
        _normalize_first_name(name[0]),
        author_dict['Middle Initial'],
        name[len(name) - 1],
        alias_value,
        1
    ]
    changed_rows.add(sheet_row)

def _add_author_alias(author_rows, changed_rows, author, full_name):
    aliases = author.setdefault('Aliases', [])
    if full_name not in aliases:
        aliases.append(full_name)
        _set_author_value(
            author_rows,
            changed_rows,
            author['Sheets Index'],
            4,
            ', '.join(aliases)
        )

def _update_author_if_alias_exists(author_rows, changed_rows, known_authors, full_name):
    for author in known_authors.values():
        if full_name in author.get('Aliases', []):
            _add_author_alias(author_rows, changed_rows, author, full_name)
            author['Pub_Number'] += 1
            _set_author_value(
                author_rows,
                changed_rows,
                author['Sheets Index'],
                5,
                author['Pub_Number']
            )
            return True
    return False

def write_known_author_updates(sheet_authors, author_updates):
    if author_updates:
        sheet_authors.batch_update(author_updates, value_input_option='RAW')

def process_known_authors(client, sheet_name_file='sheet_name.txt'):
    with open(sheet_name_file,'r', encoding='utf-8') as f:
        sheet_name = f.readlines()[0].strip()

    sheet_authors = client.open(sheet_name).worksheet('Known Authors')
    known_authors = {}
    index = 2
    full_authors = sheet_authors.get_all_values()
    time.sleep(1)
    full_authors.pop(0)
    for row_index, row_values in enumerate(full_authors):
        if(row_values == []):
            break
        index += 1
        key = row_values[0] + ' ' + row_values[2]
        known_authors[key] = create_author(row_values[1], index, 0, [])
        if len(row_values) > 4:
            if(row_values[3] != ''):
                known_authors[key]['Aliases'] = row_values[3].split(', ')
            if(row_values[4] != ''):
                known_authors[key]['Pub_Number'] = int(row_values[4])

    with open('Authors_dict.txt', 'w', encoding = 'utf-8') as dict:
        pprint.pprint(known_authors, dict)

    return known_authors, sheet_authors


def process_publication_authors(list_of_lists, client, sheet_name_file='sheet_name.txt', rebuild=False):
    """Ported author name sorting and Known Authors update from the original script.
    Returns author data and queued sheet updates for the caller to persist.
    """
    with open(sheet_name_file,'r', encoding='utf-8') as f:
        sheet_name = f.readlines()[0].strip()

    sheet_authors = client.open(sheet_name).worksheet('Known Authors')

    known_authors = {}
    index = 2

    full_authors = sheet_authors.get_all_values()
    time.sleep(1)
    if full_authors:
        full_authors.pop(0)

    author_rows = {
        row_index + 2: row_values[:5]
        for row_index, row_values in enumerate(full_authors)
        if row_values
    }
    changed_author_rows = set()

    for row_index, row_values in enumerate(full_authors):
        if(row_values == []):
            break
        index += 1
        key = row_values[0] + ' ' + row_values[2]
        known_authors[key] = create_author(row_values[1], index, 0, [])
        if len(row_values) > 4:
            if(row_values[3] != ''):
                known_authors[key]['Aliases'] = row_values[3].split(', ')
            if not rebuild and row_values[4] != '':
                known_authors[key]['Pub_Number'] = int(row_values[4])
        if rebuild:
            while len(author_rows[row_index + 2]) < 5:
                author_rows[row_index + 2].append('')
            author_rows[row_index + 2][4] = 0
            changed_author_rows.add(row_index + 2)

    Display_authors = ['Display Authors', '']

    from tqdm import tqdm
    tqdm_authors = tqdm(list_of_lists)
    tqdm_authors.set_description("Checking Publication Author Names")

    for person in tqdm_authors:
        for publication in person:
            item = publication[1]
            if item == '':
                Display_authors.append('')
            else:
                auth_array = item.split(', ')
                for person_name in auth_array:
                    person_arr = person_name.split(' ')
                    for i in range(len(person_arr) - 1):
                        if person_arr[i] == '':
                            person_arr.remove('')

                    if len(person_arr) >= 4:
                        first_last = person_arr[0] + ' ' + person_arr[len(person_arr) - 1]
                        if(first_last in known_authors):
                            if(known_authors[first_last]['Middle Initial'] == '' or known_authors[first_last]['Middle Initial'] == person_arr[1][0].upper() or known_authors[first_last]['Middle Initial'] == person_arr[2][0].upper()):
                                _add_author_alias(author_rows, changed_author_rows, known_authors[first_last], person_name)
                                known_authors[first_last]['Middle Initial'] = person_arr[1][0].upper()
                                _set_author_value(author_rows, changed_author_rows, known_authors[first_last]['Sheets Index'], 2, person_arr[1][0].upper())
                                known_authors[first_last]['Pub_Number'] += 1
                                _set_author_value(author_rows, changed_author_rows, known_authors[first_last]['Sheets Index'], 5, known_authors[first_last]['Pub_Number'])
                            else:
                                if not _update_author_if_alias_exists(
                                    author_rows, changed_author_rows, known_authors, person_name
                                ):
                                    known_authors[first_last] = create_author(person_arr[1][0].upper(), _next_author_index(author_rows), 1, [person_name])
                                    _queue_new_author(author_rows, changed_author_rows, known_authors[first_last], first_last.split(' '), person_name)
                        else:
                            found = False
                            for key in list(known_authors.keys()):
                                split = str(key).split(' ')
                                if person_arr[len(person_arr) - 1] == split[len(split) - 1] and person_arr[0][0].upper() == split[0][0].upper():
                                    array_key = key
                                    if(person_name != key):
                                        if(person_name > key):
                                            array_key = person_name
                                            known_authors[array_key] = known_authors.pop(key)
                                            _add_author_alias(author_rows, changed_author_rows, known_authors[array_key], key)
                                            _set_author_value(author_rows, changed_author_rows, known_authors[array_key]['Sheets Index'], 1, person_arr[0])
                                        else:
                                            _add_author_alias(author_rows, changed_author_rows, known_authors[array_key], person_name)
                                    _add_author_alias(author_rows, changed_author_rows, known_authors[array_key], person_name)
                                    found = True
                                    known_authors[array_key]['Pub_Number'] += 1
                                    _set_author_value(author_rows, changed_author_rows, known_authors[array_key]['Sheets Index'], 5, known_authors[array_key]['Pub_Number'])
                                    break

                            if not found:
                                if not _update_author_if_alias_exists(
                                    author_rows, changed_author_rows, known_authors, person_name
                                ):
                                    known_authors[person_name] = create_author(person_arr[1][0].upper(), _next_author_index(author_rows), 1, [person_name])
                                    _queue_new_author(author_rows, changed_author_rows, known_authors[person_name], person_arr, person_name)

                    elif len(person_arr) == 3:
                        first_last = person_arr[0] + ' ' + person_arr[2]
                        if(first_last in known_authors):
                            if(known_authors[first_last]['Middle Initial'] == '' or known_authors[first_last]['Middle Initial'] == person_arr[1][0].upper()):
                                _add_author_alias(author_rows, changed_author_rows, known_authors[first_last], person_name)
                                known_authors[first_last]['Middle Initial'] = person_arr[1][0].upper()
                                _set_author_value(author_rows, changed_author_rows, known_authors[first_last]['Sheets Index'], 2, person_arr[1][0].upper())
                                known_authors[first_last]['Pub_Number'] += 1
                                _set_author_value(author_rows, changed_author_rows, known_authors[first_last]['Sheets Index'], 5, known_authors[first_last]['Pub_Number'])
                            elif person_name in known_authors[first_last].get('Aliases', []):
                                _add_author_alias(author_rows, changed_author_rows, known_authors[first_last], person_name)
                                known_authors[first_last]['Pub_Number'] += 1
                                _set_author_value(
                                    author_rows,
                                    changed_author_rows,
                                    known_authors[first_last]['Sheets Index'],
                                    5,
                                    known_authors[first_last]['Pub_Number']
                                )
                            else:
                                if not _update_author_if_alias_exists(
                                    author_rows, changed_author_rows, known_authors, person_name
                                ):
                                    known_authors[first_last] = create_author(person_arr[1][0].upper(), _next_author_index(author_rows), 1, [person_name])
                                    _queue_new_author(author_rows, changed_author_rows, known_authors[first_last], first_last.split(' '), person_name)
                        else:
                            if not _update_author_if_alias_exists(
                                author_rows, changed_author_rows, known_authors, person_name
                            ):
                                known_authors[first_last] = create_author(person_arr[1][0].upper(), _next_author_index(author_rows), 1, [person_name])
                                _queue_new_author(author_rows, changed_author_rows, known_authors[first_last], first_last.split(' '), person_name)

                    elif len(person_arr) == 2:
                        compact_initials_index = None
                        compact_initials = None
                        for candidate_index in range(2):
                            parsed_initials = _parse_two_initials(person_arr[candidate_index])
                            other_token = person_arr[1 - candidate_index]
                            if parsed_initials and any(char.islower() for char in other_token):
                                compact_initials_index = candidate_index
                                compact_initials = parsed_initials
                                break

                        single_initial_index = None
                        if compact_initials_index is None:
                            for candidate_index in range(2):
                                if (
                                    re.fullmatch(r'[A-Z]\.?', person_arr[candidate_index])
                                    and any(char.islower() for char in person_arr[1 - candidate_index])
                                ):
                                    single_initial_index = candidate_index
                                    break

                        if compact_initials_index is not None or single_initial_index is not None:
                            if compact_initials_index is not None:
                                upper = compact_initials_index
                                full = 1 - upper
                                first_initial, middle_initial = compact_initials
                            elif single_initial_index == 0:
                                upper = 0
                                full = 1
                            else:
                                upper = 1
                                full = 0
                            if compact_initials_index is None:
                                initials = person_arr[upper].replace('.', '')
                                first_initial = initials[0]
                                middle_initial = initials[1:2]

                            alias_found = False
                            matching_keys = []
                            for key in list(known_authors.keys()):
                                for alias in known_authors[key].get('Aliases', []):
                                    if(alias == person_name):
                                        alias_found = True
                                        matching_keys = [key]
                                        break

                                key_parts = key.split()
                                if (
                                    len(key_parts) >= 2
                                    and key_parts[-1].casefold() == person_arr[full].casefold()
                                    and key_parts[0][0].casefold() == person_arr[upper][0].casefold()
                                ):
                                    author_middle = known_authors[key]['Middle Initial'].strip().upper()
                                    person_middle = middle_initial.upper()
                                    if not person_middle or not author_middle or person_middle == author_middle:
                                        matching_keys.append(key)

                                if alias_found:
                                    break

                            if matching_keys:
                                matching_key = matching_keys[0]
                                _add_author_alias(
                                    author_rows,
                                    changed_author_rows,
                                    known_authors[matching_key],
                                    person_name
                                )
                                author_middle = known_authors[matching_key]['Middle Initial'].strip().upper()
                                if middle_initial and not author_middle:
                                    known_authors[matching_key]['Middle Initial'] = middle_initial.upper()
                                    _set_author_value(
                                        author_rows,
                                        changed_author_rows,
                                        known_authors[matching_key]['Sheets Index'],
                                        2,
                                        middle_initial.upper()
                                    )
                                known_authors[matching_key]['Pub_Number'] += 1
                                _set_author_value(
                                    author_rows,
                                    changed_author_rows,
                                    known_authors[matching_key]['Sheets Index'],
                                    5,
                                    known_authors[matching_key]['Pub_Number']
                                )
                            elif not alias_found:
                                if not _update_author_if_alias_exists(
                                    author_rows, changed_author_rows, known_authors, person_name
                                ):
                                    author_key = first_initial + ' ' + person_arr[full]
                                    known_authors[author_key] = create_author(
                                        middle_initial.upper(),
                                        _next_author_index(author_rows),
                                        1,
                                        [person_name]
                                    )
                                    _queue_new_author(
                                        author_rows,
                                        changed_author_rows,
                                        known_authors[author_key],
                                        [first_initial, person_arr[full]],
                                        person_name
                                    )
                        else:
                            found = False
                            for key in list(known_authors.keys()):
                                split = str(key).split(' ')
                                if person_arr[1] == split[1] and person_arr[0][0].upper() == split[0][0].upper():
                                    array_key = key
                                    if(person_name != key):
                                        if(person_name > key):
                                            array_key = person_name
                                            known_authors[array_key] = known_authors.pop(key)
                                            _add_author_alias(author_rows, changed_author_rows, known_authors[array_key], key)
                                            _set_author_value(author_rows, changed_author_rows, known_authors[array_key]['Sheets Index'], 1, person_arr[0])
                                        else:
                                            _add_author_alias(author_rows, changed_author_rows, known_authors[array_key], person_name)
                                    _add_author_alias(author_rows, changed_author_rows, known_authors[array_key], person_name)
                                    found = True
                                    known_authors[array_key]['Pub_Number'] += 1
                                    _set_author_value(author_rows, changed_author_rows, known_authors[array_key]['Sheets Index'], 5, known_authors[array_key]['Pub_Number'])
                                    break

                            if not found:
                                if not _update_author_if_alias_exists(
                                    author_rows, changed_author_rows, known_authors, person_name
                                ):
                                    known_authors[person_name] = create_author('', _next_author_index(author_rows), 1, [person_name])
                                    _queue_new_author(author_rows, changed_author_rows, known_authors[person_name], person_arr, person_name)

    import pprint
    with open('Authors_dict.txt', 'w', encoding = 'utf-8') as dict:
        pprint.pprint(known_authors, dict)

    author_updates = [
        {'range': f'A{row}:E{row}', 'values': [author_rows[row]]}
        for row in sorted(changed_author_rows)
    ]
    return known_authors, sheet_authors, Display_authors, author_updates
