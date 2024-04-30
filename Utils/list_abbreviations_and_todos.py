import csv
import os

DIRECTORY = "Text"
OUT_ABBR = "List of abbreviations.csv"
OUT_TODO = "List of TODOs.txt"

def replace_multiple(string: str, repl_list, repl_by):
    for repl in repl_list:
        string = string.replace(repl, repl_by)
    return string

def scan_abbreviations_in_file(path):
    file_abbr_list = list()

    print(f"Checking {path} for abbreviations")
    with open(path, "r", encoding='utf-8') as file:
        text = file.read()
        text = replace_multiple(text, "{}()[]-", " ")
        words = text.split()

        for word in words:
            # word = "".join(char for char in word if char.isalnum())
            word = word.strip(" .,-")
            if word.startswith("\\"):  # skip commands
                pass
            elif is_abbreviation(word):
                file_abbr_list.append(word)
    
    return file_abbr_list


def is_abbreviation(word: str):
    try:
        assert len(word) >= 2  # word is at least 2 characters

        upper = 0
        for c in word:
            upper += int(c.isupper())
        
        assert upper >= 2  # at least 2 characters are uppercase
        assert upper >= len(word) - upper  # there is at least as many uppercase chars as lowercase

    except AssertionError:
        return False
    else:
        return True


def scan_todos_in_file(path):
    filename = path.split("\\")[-1]
    file_todo_list = list()

    print(f"Checking {path} for todos")
    with open(path, "r", encoding='utf-8') as file:
        lines = file.readlines()

        for linenum, line in enumerate(lines):
            if (i := line.upper().find("TODO")) >= 0:
                todo = replace_multiple(line[i:], "{}()[]-%", "").strip()
                file_todo_list.append((filename, linenum, todo))
    return file_todo_list


def find_abbr_meaning_in_files(filelist, abbr_list):
    abbr_expl_list = list()
    text = str()

    for path in filelist:
        with open(path, "r", encoding='utf-8') as file:
            text += file.read()

    text = replace_multiple(text, "{}[]-,.;\\", " ")
    words = text.split()

    idx = -1
    for abbr in abbr_list:
        all_abbr_expl = list()
        try:
            while True:
                idx = words.index(f"({abbr})", idx+1)
                potential_abbr_expl = words[max(0, idx - len(abbr) * 2) : idx]
                all_abbr_expl.append(search_for_abbr_expl_in_text(abbr, potential_abbr_expl))
        except ValueError:
            idx = -1

        abbr_expl_list.append((abbr, *all_abbr_expl))

    return abbr_expl_list
            

def search_for_abbr_expl_in_text(abbr, potential_abbr_expl):
    print(abbr, potential_abbr_expl)

    for idx in range(len(abbr), 0, -1):
        abbr_expl = potential_abbr_expl[idx:]
        if abbr_expl[0].lower().startswith(abbr.lower()[0]):
            return " ".join(abbr_expl)
    
    return None


abbr_list = list()
todo_list = list()
abbr_filelist = list()

for filename in os.listdir(DIRECTORY):
    f = os.path.join(DIRECTORY, filename)
    if os.path.isfile(f) and f.endswith(".tex"):
        if "main" not in f:
            abbr_filelist.append(f)
            abbr_list = list(set(abbr_list + scan_abbreviations_in_file(f)))
        todo_list = list(set(todo_list + scan_todos_in_file(f)))

abbr_list = find_abbr_meaning_in_files(abbr_filelist, abbr_list)

abbr_list.sort()

print(abbr_list)
print(todo_list)

# writing to csv file  
with open(DIRECTORY + "\\" + OUT_ABBR, 'w', newline='') as csvfile:  
    csvwriter = csv.writer(csvfile, dialect="excel")

    for abbr in abbr_list:
        csvwriter.writerow(abbr) 

with open(DIRECTORY + "\\" + OUT_TODO, 'w', newline='') as txtfile:
    for file, line, todo in todo_list:
        txtfile.write(f"File {file:<30} at line {line:<3}: {todo}\r\n")
