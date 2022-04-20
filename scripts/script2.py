# -*- coding:utf-8 -*-

import pandas as pd
import re
import sys
import os
import time
import pathlib
from datetime import datetime
from time import sleep
from Bio import Entrez
from Bio import SeqIO

# Function define
def timestamp():
    dt = datetime.now()
    ts = dt.strftime('%d%b%Y-%H:%M:%S')
    return ts

def elapsed_time(dt_start, dt_end):
    ds = (dt_end - dt_start).total_seconds()
    h = int(ds // 3600)
    m = int((ds % 3600) // 60)
    s = int(ds % 60)
    str_line = str(h).zfill(2) + 'h ' + str(m).zfill(2) + 'm ' + str(s).zfill(2) + 's'
    return str_line

# Open setting text file
with open('setting.txt', encoding = 'shift-jis') as st:
    txt = st.read()

# color setting
class Color:
    BLACK = '\033[30m'
    RED ='\033[31m'
    GREEN ='\033[32m'
    YELLOW ='\033[33m'
    BLUE ='\033[34m'
    PURPLE ='\033[35m'
    CYAN ='\033[36m'
    WHITE ='\033[37m'
    END ='\033[0m'
    BOLD ='\033[1m'
    UNDERLINE ='\033[4m'
    IVNISIBLE ='\033[08m'
    REVERCE ='\033[07m'

# Introduce to NCBI
Entrez.email = re.search(r'Email\s*=\s*(\S+)',txt).group(1)
Entrez.api_key = re.search(r'Api_key\s*=\s*(\S+)',txt).group(1)

print('\n' + Color.GREEN + \
    "# ========================================= #\n"\
    "# === Genbank Information download step === #\n"\
    "# ========================================= #\n"\
    + Color.END)

# Function define
def write_file(out_file, mode, info_str):
    # write_file
    with open(out_file, mode=mode) as OUT:
        OUT.write(info_str + '\n')

# ACC chunk function
def get_acc_chunks(id_list_file, max_in_chunk):
    # group specified number of accession ids in a chunk
    # return a list of chunks and total id count
    # regex to skip some lines
    comment_ptn = re.compile('^#')
    blank_ptn = re.compile(r'^\s*$')

    chunks = []
    chunk = []
    total_id_count = 0
    with open(id_list_file, mode = 'r') as ACC:
        id_count_in_a_chunk = 0
        for l in ACC :
            l = l.rstrip() # remove annoying EOL
            if(comment_ptn.match(l)) or (blank_ptn.match(l)):
                continue # skip those lines
            else:
                id_count_in_a_chunk += 1
                if id_count_in_a_chunk <= int(max_in_chunk):
                    chunk.append(l)
                else :
                    # dump a chunk
                    chunks.append(chunk)
                    total_id_count += len(chunk)
                    # start a new chunk
                    chunk = [l]
                    id_count_in_a_chunk = 1
        # append the last chunk
        chunks.append(chunk)
        total_id_count += len(chunk)

    return chunks, total_id_count

# flag function
def when_exeption_happened(flag, e, acc_chunk, display_chunk_num):
    ts = timestamp()
    if flag == 1: # error in efetch
        where = 'efetch'
    elif flag == 2 : # error in SeqIO
        where = 'SeqIO'
    # report on STDOUT
    print(' --- Error in ' + where + '. See log file. @' + ts)
    print(str(type(e)))
    print(e.args)
    if flag == 1:
        print('Will bae countinued after 1 min')
        print()
    # report to log
    log_lines = [
        'Error in ' + where + '] chunk:' + display_chunk_num + '@' + ts,
        str(type(e)),
        str(e.args), # if any
        '# accession ids in the chunk'
    ]

    # append accession ids in the chunk
    for acc in acc_chunk:
        log_lines.append(acc)
    log_lines.append('\n')
    write_file(log_file, 'a', '\n'.join(log_lines))

def main():
    print(' --- program started; ', timestamp(), '\n')
    # mark start time
    dt_start = datetime.now()

    # write log header
    log_lines = [
        '# processed by: ' + sys.argv[0],
        '# process started: ' + timestamp(),
        '# accession id file: ' + acc_id_file_path,
        '# db=' + db,
        '# rettype=' + rettype,
        '# retmode=' + retmode,
        '# id count per query:' + str(count_in_query) + '\n'
    ]

    write_file(log_file, 'w', '\n'.join(log_lines)) # create a new file

    # get accession id chunks from accession id file
    acc_chunks, total_id_count = get_acc_chunks(acc_id_file_path, count_in_query)
    print(f" --- total_id_count : {total_id_count} seqs.")
    print()

    # some buckets to monitor the processes
    chunk_num = 0
    done_id_count = 0
    error_id_count = 0
    error_chunk_count = 0

    # loop
    for acc_chunk in acc_chunks:
        chunk_num += 1

        display_chunk_num = str(chunk_num).zfill(len(str(len(acc_chunks))))
        print('      ' + acc_id_file_stem + ': chunk:' + Color.GREEN + display_chunk_num + Color.END)

        # error capture 1: efetch
        try:
            net_handle = Entrez.efetch(db=db, id=acc_chunk, rettype=rettype, retmode=retmode)
        except Exception as e:
            when_exeption_happened(1, e, acc_chunk, display_chunk_num)
            # counts
            error_chunk_count += 1
            error_id_count += len(acc_chunk)
            # wait 1 min
            time.sleep(60)
            continue
        else:
            pass # I do not want nested try-except. Pass the case to the next try-exept.

        # error capture 2: broken records?
        try:
            num = 1
            # if record is broken, it will be trapped in the next line.
            for record in SeqIO.parse(net_handle, fmt):
                # write records into a file
                out_file = './GI_Folder/' + record.id + ext
                SeqIO.write(record, out_file, fmt)
                    ##
                    #  I want to add a function to convert a .bg format to .fasta format and omit the step to output the .bg files.
                    ## 
        except Exception as e:
            when_exeption_happened(2, e, acc_chunk, display_chunk_num)
            # counts
            error_chunk_count += 1
            error_id_count += len(acc_chunk)
        else:
            # report on STDOUT
            dt_now = datetime.now()
            et = elapsed_time(dt_start, dt_now)
            print('    Done. ET:' + et)
            print()
            done_id_count += len(acc_chunk)
        finally:
            net_handle.close()
            time.sleep(1)

    # post loop works
    ts = timestamp()
    dt_now = datetime.now()
    et = elapsed_time(dt_start, dt_now)

    # log process summaray
    log_lines = [
        '# === SUMMARY ===',
        '# total accession id count: ' + str(total_id_count),
        '# chunk created: ' + str(len(acc_chunks)),
        '# error chunks: ' + str(error_chunk_count),
        '# done id count: ' + str(done_id_count),
        '# error id count: ' + str(error_id_count),
        '# process finished' + ts,
        '# total elapsed time; ' + et]
    write_file(log_file, 'a', '\n'.join(log_lines))

if __name__ == '__main__':
    # get accession id file path
    acc_id_file_path = './ID/ID.txt'

    # vars for efetch
    # 'nucleotide' will work as well
    db = re.search(r'Database2\s*=\s*(\S+)', txt).group(1)
    
    # use 'fasta', 'gbwithparts', for fasta, GenBank(full)
    rettype = re.search(r'Rettype2\s*=\s*(\S+)', txt).group(1)
    retmode = re.search(r'Retmode2\s*=\s*(\S+)', txt).group(1)
    
    # count in one query for NCBI
    count_in_query = re.search(r'Query_per_Chunk2\s*=\s*(\S+)', txt).group(1)  

    # vars for SeqIO
    # format name adjustment and file extension
    if rettype in ['gb', 'gbwithparts']:
        fmt = 'gb'
        ext = '.gb'
    elif rettype == 'fasta':
        fmt = 'fasta'
        ext = '.fa'
    else:
        print("rettype shold be 'gb', 'gbwithparts', or 'fasta'")
        sys.exit(1)

    # log file path construction
    acc_id_file_stem = os.path.splitext(os.path.basename(acc_id_file_path))[0]
    # Make file name
    runday = sys.argv[1]
    log_file = './Results/' + f'{runday}' + "_" + acc_id_file_stem + '_download_log.txt'

    # run !
    main()
