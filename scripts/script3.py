# -*- coding:utf-8 -*-

from datetime import datetime
import re
import os
import sys
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm
import pathlib
import joblib


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

def makeSeqRecord(seq_obj, id_for_fasta):
    record = SeqRecord(seq_obj, id = id_for_fasta, description = '')
    return record

# Cechk Fasta file folder
if not os.path.exists('./Results/nodatafasta'):
    pathlib.Path('./Results/nodatafasta').mkdir(exist_ok=True)
else:
    pass

if not os.path.exists('./Results/fasta'):
    pathlib.Path('./Results/fasta').mkdir(exist_ok=True)
else:
    pass

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


if __name__ == '__main__':
    print('\n' + Color.GREEN + \
    "# ========================================================== #\n"\
    "# === Extract 12s rRNA sequence from Genbank information === #\n"\
    "# ========================================================== #\n"\
    + Color.END)

    # mark start time
    print(' --- program started; ', timestamp(), '\n')
    dt_start = datetime.now()

    gi_filepath = os.listdir('./GI_Folder')
    
    def gen_fa(gilist):
        # basket to product var
        product_var = []

        for record in SeqIO.parse(f'./GI_Folder/{gi_filepath[gilist]}', 'genbank'):

            # basket to collect rRNA slices
            rRNA_records = []
            # gat entire seq
            parental_seq = record.seq # seq obj
            # get accession
            acc = record.id
            # get organism name
            organism = record.annotations['organism']

           
            for feature in record.features:
                # pick onlay rRNA
                if feature.type == 'rRNA':
                    #  get sequence slice with location
                    seq_slice = feature.location.extract(parental_seq)
                    
                    start = feature.location.parts[0].start.position + 1
                    end = feature.location.parts[-1].end.position

                    # get product annotation
                    if 'product' in feature.qualifiers:
                        product = feature.qualifiers['product'][0]
                    else:
                        product = 'NA'

                    product_var += product.split(',')

                    rRNA_id_line = acc + '_' + organism + '_' + product
                    rRNA_id_line += '_(' + str(start) + '-' + str(end) + ')'

                    # construct SeqRecord object and put it in a busket
                    # output file name construction
                    rRNA_records.append(makeSeqRecord(seq_slice, rRNA_id_line))

                    # log DL species list
                    with open('./Results/MiFishSeq_temp_dlsplist.txt', mode = 'a') as spltemp:
                        spltemp.write(organism + '\n')

            check = 1
            for feature in record.features:
                if feature.type == 'rRNA':
                    
                    # Avoid ACC of shortgun genome sequences that do not contain sequence information.
                    key = re.compile('WGS.*')
                    if not key.search(str(record.annotations['keywords']).upper()):
                        rRNA_count = len(rRNA_records)
                        if rRNA_count == 1:
                            out_fasta = f'./Results/fasta/{acc}_all-rRNA_1.fasta'
                        elif rRNA_count >= 2:
                            out_fasta = f'./Results/fasta/{acc}_all-rRNA_{str(rRNA_count)}.fasta'
                        else:
                            out_fasta = f'./Results/nodatafasta/{acc}_No_rRNA-FOUND.fasta'
                            check = 0
                    else:
                        rRNA_count = len(rRNA_records)
                        if rRNA_count == 1:
                            out_fasta = f'./Results/nodatafasta/{acc}_all-rRNA_1_noincludeseq.fasta'
                        elif rRNA_count >= 2:
                            out_fasta = f'./Results/nodatafasta/{acc}_all-rRNA_{str(rRNA_count)}_noincludeseq.fasta'
                        else:
                            out_fasta = f'./Results/nodatafasta/{acc}_No_rRNA-FOUND_noincludeseq.fasta'
                        check = 2                                               

                    with open(out_fasta, mode = 'w', encoding = 'utf-8') as ofa:
                        if check == 1:
                            SeqIO.write(rRNA_records, ofa, 'fasta')
                        elif check == 2:
                            ofa.write('NO SEQUENCE INFORMATION FOUND')
                        else:
                            ofa.write('NO rRNA ANNOTATION FOUND')

            # make product list
            with open('./Results/product_list.tsv', mode = 'a') as prod:
                prod.write('\n'.join(map(str,set(product_var))) + '\n')

    # parallel run
    joblib.Parallel(n_jobs=-1)(joblib.delayed(gen_fa)(gilist) for gilist in range(0, len(gi_filepath)))

    # remove duplicate spcies list    
    with open('./Results/MiFishSeq_temp_dlsplist.txt', mode = 'r') as spltemp2:
        with open('./Results/MiFishSeq_dlsplist.txt', mode = 'a') as spl:
            spl.write(''.join(map(str, set(spltemp2.readlines()))))
    
    # remove temp splist file
    pathlib.Path('./Results/MiFishSeq_temp_dlsplist.txt', missing_ok=True).unlink()


dt_now = datetime.now()
et = elapsed_time(dt_start, dt_now)

print('')
print(' --- Extract 12s rRNA sequence Done.ET:' + et + '\n')
