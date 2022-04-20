# -*- coding:utf-8 -*-

import pathlib
import re
import os
import shutil
import sys
import glob
from datetime import datetime
from Bio import Entrez
import numpy as np

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

# Check Results Folder
def check_result_folder():
    if not os.path.exists('./Results'):
        pathlib.Path('./Results').mkdir(exist_ok=True)
    else:
        shutil.rmtree('./Results')
        pathlib.Path('./Results').mkdir(exist_ok=True)
    # Check ID stock folder
    if not os.path.exists('./ID'):
        pathlib.Path('ID').mkdir(exist_ok=True)
    else:
        pass
    # Check Stock ID file
    if not os.path.exists('./ID/StockID.txt'):
        with open('./ID/StockID.txt', mode = 'w') as stockid:
            stockid.write('')        
    else:
        pass
    # Check ID file
    today = datetime.now().strftime('%Y%d')
    if not os.path.exists('./ID/ID.txt'):
        with open('./ID/ID.txt', mode = 'w') as i:
            i.write('')        
    else:
        pathlib.Path('./ID/ID.txt', missing_ok = True).unlink()

    # Cechk Fasta file folder
    if not os.path.exists('./Results/fasta'):
        pathlib.Path('./Results/fasta').mkdir(exist_ok=True)
    else:
        pass
    # Check Genbank ID folder
    if not os.path.exists('./GI_Folder'):
        pathlib.Path('./GI_Folder').mkdir(exist_ok=True)
    else:
        pass
    
# Get difference of two lists func
def diff_list(li1, li2):
    return list(set(li1) - set(li2))

# Omit nothing id list func.
def no_id_omit_func(acc:list):
    '''
    if downloaded acc is empty, remove 'Supplied+id+parameter+is+empty.' in list.
    '''
    if acc[0] == 'Supplied+id+parameter+is+empty.':
        return []
    else:
        return acc

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
Entrez.email = re.search(r'Email\s*=\s*(\S+)', txt).group(1)
Entrez.api_key = re.search(r'Api_key\s*=\s*(\S+)', txt).group(1)

# Check
check_result_folder()

print('\n' + Color.GREEN + \
    '# =============================== #\n' + \
    '# = Download target sequence ID = #\n' + \
    '# =============================== #\n' + \
    Color.END)
print(" --- Program started:         ", timestamp())

# Mark start time
dt_start = datetime.now()

if __name__ == '__main__':
    # vers for efetch
    db = re.search(r'Database1\s*=\s*(\S+)', txt).group(1)
    rettype = re.search(r'Rettype1\s*=\s*(\S+)', txt).group(1)
    retmode = re.search(r'Retmode1\s*=\s*(\S+)', txt).group(1)
    count_in_query = re.search(r'Count_in_query1\s*=\s*(\S+)', txt).group(1)

    before_date = datetime.strptime(re.search(r'Before_day\s*=\s*(\S+)', txt).group(1), '%Y/%m/%d').strftime('%Y/%m/%d')
    max_date = datetime.today().strftime('%Y/%m/%d')

    print(' --- Sequence download period:', before_date, '-', max_date)
    print('')

    # Teleo
    handle_teleo = Entrez.esearch(db=db, retmax=count_in_query,\
        term = f'txid32443[Organism] AND 12S[All Fields] AND (ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter] AND {before_date}[PDAT] : {max_date}[PDAT] AND 120[SLEN] : 30000[SLEN]')
    term_teleo = f'txid32443[Organism] AND 12S[All Fields]) AND (ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter] AND {before_date}[PDAT] : {max_date}[PDAT] AND 120[SLEN] : 30000[SLEN]'

    handle_comp_teleo = Entrez.esearch(db=db, retmax=count_in_query,\
        term = f'txid32443[Organism] AND complete[prop] AND (ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter] AND {before_date}[PDAT] : {max_date}[PDAT]')
    term_comp_teleo = f'txid32443[Organism] AND complete[prop] AND (ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter] AND {before_date}[PDAT] : {max_date}[PDAT]'

    print('')
    print( ' --- Teleostei nucleotide search term is')
    print(f' --- {term_teleo}')
    print(f' --- {term_comp_teleo}')

    # Chondri
    handle_chondri = Entrez.esearch(db=db, retmax=count_in_query,\
        term = f'txid7777[Organism] AND 12S[All Fields] AND (ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter] AND {before_date}[PDAT] : {max_date}[PDAT] AND 120[SLEN] : 30000[SLEN]')
    term_chondri = f'txid7777[Organism] AND 12S[All Fields]) AND (ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter] AND {before_date}[PDAT] : {max_date}[PDAT] AND 120[SLEN] : 30000[SLEN]'

    handle_comp_chondri = Entrez.esearch(db=db, retmax=count_in_query,\
        term = f'txid7777[Organism] AND complete[prop] AND (ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter] AND {before_date}[PDAT] : {max_date}[PDAT]')
    term_comp_chondri = f'txid7777[Organism] AND complete[prop] AND (ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter] AND {before_date}[PDAT] : {max_date}[PDAT]'

    print('')
    print( ' --- Chondri nucleotide search term is')
    print(f' --- {term_chondri}')
    print(f' --- {term_comp_chondri}')

    # Cyclostomata
    handle_cyclostomata = Entrez.esearch(db=db, retmax=count_in_query,\
        term = f'txid1476529[Organism] AND 12S[All Fields] AND (ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter] AND {before_date}[PDAT] : {max_date}[PDAT] AND 120[SLEN] : 30000[SLEN]')
    term_cyclostomata = f'txid1476529[Organism] AND 12S[All Fields]) AND (ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter] AND {before_date}[PDAT] : {max_date}[PDAT] AND 120[SLEN] : 30000[SLEN]'

    handle_comp_cyclostomata = Entrez.esearch(db=db, retmax=count_in_query,\
        term = f'txid1476529[Organism] AND complete[prop] AND (ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter] AND {before_date}[PDAT] : {max_date}[PDAT]')
    term_comp_cyclostomata = f'txid1476529[Organism] AND complete[prop] AND (ddbj_embl_genbank[filter] AND is_nuccore[filter] AND mitochondrion[filter] AND {before_date}[PDAT] : {max_date}[PDAT]'

    print('')
    print( ' --- Cyclostomata nucleotide search term is')
    print(f' --- {term_cyclostomata}')
    print(f' --- {term_comp_cyclostomata}')
    

    record_teleo = Entrez.read(handle_teleo)
    record_chondri = Entrez.read(handle_chondri)
    record_cyclostomata = Entrez.read(handle_cyclostomata)

    record_comp_teleo = Entrez.read(handle_comp_teleo)
    record_comp_chondri = Entrez.read(handle_comp_chondri)
    record_comp_cyclostomata = Entrez.read(handle_comp_cyclostomata)

    Partial_res1 = Entrez.efetch(db='nuccore', id=record_teleo['IdList'], rettype = "acc", retmode="text")
    Partial_res2 = Entrez.efetch(db='nuccore', id=record_chondri['IdList'], rettype = "acc", retmode="text")
    Partial_res3 = Entrez.efetch(db='nuccore', id=record_cyclostomata['IdList'], rettype = "acc", retmode="text")

    Comp_res1 = Entrez.efetch(db='nuccore', id=record_comp_teleo['IdList'], rettype = "acc", retmode="text")
    Comp_res2 = Entrez.efetch(db='nuccore', id=record_comp_chondri['IdList'], rettype = "acc", retmode="text")
    Comp_res3 = Entrez.efetch(db='nuccore', id=record_comp_cyclostomata['IdList'], rettype = "acc", retmode="text")

    seqacc  = no_id_omit_func([id.strip() for id in Partial_res1])
    seqacc += no_id_omit_func([id.strip() for id in Partial_res2])
    seqacc += no_id_omit_func([id.strip() for id in Partial_res3])

    seqacc += no_id_omit_func([id.strip() for id in Comp_res1])
    seqacc += no_id_omit_func([id.strip() for id in Comp_res2])
    seqacc += no_id_omit_func([id.strip() for id in Comp_res3])

    exists_gb = [os.path.basename(p).replace('.gb', '') for p in glob.iglob(r'./GI_Folder/*.gb')] # remove duplicate and change numpy array

    # Removing already exists acc from download acc list.
    no_exists_acc = diff_list(seqacc,exists_gb)

    # Differences of GI number bwtween StockID and ID file
    # Stock IDs
    with open('./ID/StockID.txt', mode = 'r') as sid:
        stock_seqidbf = sid.read().splitlines()

    # Dif
    DifSeqID = diff_list(seqacc,stock_seqidbf)
    print(DifSeqID)

    if len(DifSeqID) <= 0:
        print('')
        print(' --- Newly ID is nothing.')
        print(' --- Download program stop.')
        print('')
        sys.exit()
    else:
        pass

    str_ID = '\n'.join(map(str,DifSeqID))
    with open('./ID/ID.txt', mode = 'w') as f:
        f.write(str_ID + '\n')

    str_stockID = '\n'.join(map(str,DifSeqID))

    # Opne with add mode
    with open('./ID/StockID.txt', mode = 'a') as sid2:
        sid2.write(str_stockID+'\n')

    # Opne with read only mode
    with open('./ID/StockID.txt', mode = 'r') as sid3:
        stock_seqidaf = set(sid3.read().splitlines())

    print('')
    print(' --- Number of stock IDs before add is', stock_seqidbf)
    print(' --- Number of newly ID is', DifSeqID)
    print(' --- Number of add after stock ID is', stock_seqidaf)

dt_now = datetime.now()
et = elapsed_time(dt_start, dt_now)

print('')
print(' --- ID download Done.ET:' + et + '\n')
