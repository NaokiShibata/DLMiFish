# -*- coding:utf-8 -*-

from curses import COLOR_GREEN
import subprocess as sp
import pathlib
import glob
import re
import os
from datetime import datetime
import shutil
from Bio import Entrez
import sys

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

if __name__ == "__main__":
    print('\n' + Color.GREEN + \
        '# =============================== #\n' + \
        '# === Convert sequence format === #\n' + \
        '# =============================== #\n' + \
        Color.END)
    
    # mark start time
    print(' --- program started; ', timestamp(), '\n')
    dt_start_f1 = datetime.now()
    dt_start_f2 = datetime.now().strftime('%Y%m%d_%H%M%S')
    dt_start_f3 = datetime.now().strftime('%Y%m%d')

    fasta_list = glob.glob(r'./Results/fasta/*.fasta') # pattern match fasta file
    merged_fasta = './Results/fasta/rRNA.fas' # make merged fasta file

    # confirm fasta file exist
    if not os.path.isfile(merged_fasta):
        with open(merged_fasta, mode = 'w') as fme:
            pass
    else:
        renamed_file = f'./Results/fasta/rRNA{dt_start_f2}_tempfile.fas'
        os.renames(merged_fasta,renamed_file)

    # Merge fasta file by shell command and seqkit
    mcom = 'cat ./Results/fasta/*.fasta | seqkit seq -w 0 > ./Results/fasta/rRNA.fas'
    sp.run(mcom , shell = True)

    # remove each fasta sequence because already merge those
    for f in fasta_list:
        os.remove(f)

    # Make path to save for extracted 12s rRNA seq file.

    if not os.path.exists(f'./Results/fasta/{dt_start_f2}_12s'):
        pathlib.Path(f'./Results/fasta/{dt_start_f2}_12s').mkdir(exist_ok=True)
    else:
        pass

    # Extracat 12s rRNA sequence by seqkit
    ecom = f"seqkit grep -w 0 -nirp '12s|s-rrna|small|rrnS' ./Results/fasta/rRNA.fas > ./Results/fasta/{dt_start_f2}_12s/00_12srRNA_merge.fas"
    sp.run(ecom, shell = True)

    print('\n' + Color.GREEN + \
        '# ============================================ #\n' + \
        '# === Move 12s rRNA seq to MiFishDB Folder === #\n' + \
        '# ============================================ #\n' + \
        Color.END)
    print()

    # MiFish DB 
    db_dir = f'./MiFishDB/{dt_start_f2}_12s'
    pathlib.Path(f'{db_dir}').mkdir(exist_ok=True)

    # Moving 12s rRNA merged fasta file
    shutil.move(f'./Results/fasta/{dt_start_f2}_12s/00_12srRNA_merge.fas', db_dir)

    # Miving log file
    runday = sys.argv[1]
    shutil.move(f'./Results/{runday}_ID_download_log.txt',db_dir) 
    
    # product var
    shutil.move(f'./Results/product_list.tsv', db_dir)

    print('\n' + Color.GREEN + \
        '# ================================= #\n' + \
        '# === Sequence file filteration === #\n' + \
        '# ================================= #\n' + \
        Color.END)
    print() 

    # Remove unnecessary information 
    with open(f'{db_dir}/00_12srRNA_merge.fas', mode = 'r') as fa:
        fasta = fa.read()
        frp = re.sub(r'_NA_.*', '', fasta, flags = re.IGNORECASE)
        frp = re.sub(r'_rrnS_.*', '', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_mt-_.*', '', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_rnr1_.*', '', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_12S rRNA_.*', '', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_12S_.*', '', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_s-rRNA_.*', '', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_rrn12_.*', '', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_RRN_.*', '', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_12_.*', '', frp)
        frp = re.sub(r'_12SrRNA_.*', '', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_12S ribosomal RNA_.*', r'', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_12S small subunit ribosomal RNA_.*','', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_small subunit ribosomal RNA_.*','', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_12s large subunit ribosomal RNA_.*','', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_12S ribosomal RNA subunit_.*','', frp, flags = re.IGNORECASE)
        frp = re.sub(r'_12S ribosomal RNA subunit_.*','', frp, flags = re.IGNORECASE)        

        # Output
        with open(f'{db_dir}/01_12srRNA_merge2.fas', mode='w') as qe:
            qe.write(frp)

    ### option ###
    opt = re.search(r'option\s*=\s*(\S+)',txt).group(1)
    err = re.search(r'error\s*=\s*(\S+)',txt).group(1)
    minL = re.search(r'minlength\s*=\s*(\S+)',txt).group(1)

    if re.match(opt, 'yes', flags = re.IGNORECASE):
        error = 0
        cores = 24
        input = f'{db_dir}/01_12srRNA_merge2.fas'
        lastout = f'{db_dir}/02_12srRNA_merge3.fas'

        # primer set
        mfuf = 'GTCGGTAAAACTCGTGCCAGC'
        mfur = 'CATAGTGGGGTATCTAATCCCAGTTTG'
        mfev2f = 'RGTTGGTAAATCTCGTGCCAGC'
        mfev2r = 'GCATAGTGGGGTATCTAATCCTAGTTTG'
        mfu2f = 'GCCGGTAAAACTCGTGCCAGC'
        mfu2r = 'CATAGGAGGGTGTCTAATCCCCGTTTG'
        mflf = 'GCTGGTAAACCTCGTGCCAGC'
        mflr = 'CATAGCGGGGTATCTAATCCCGGTTTG'

        # primer list
        fpls = [mfuf,mfev2f,mfu2f,mflf]
        rpls = [mfur,mfev2r,mfu2r,mflr]

        pnum = 1
        # Reverse primer treatment        
        for fpl in fpls:
            if pnum == 1:
                print(f'Forward primer seq: {fpl}')

                fpcutcom = f'cutadapt -e {error} --cores {cores} -g {fpl} {input} > {db_dir}/temp{pnum}.fas'
                sp.run(fpcutcom, shell =True)
                pnum += 1
            else:
                print(f'Forward primer seq: {fpl}')

                bpnum = pnum -1
                fpcutcom = f'cutadapt -e {error} --cores {cores} -g {fpl} {db_dir}/temp{bpnum}.fas > {db_dir}/temp{pnum}.fas'
                sp.run(fpcutcom, shell =True)
                pnum += 1

        # Total primer num
        ptotal_num = len(fpls) + len(rpls)

        # Reverse primer treatment
        for rpl in rpls:

            # reverse complement
            rpseq1 = sp.run(f"echo '>seq\n{rpl}'|sed -z 's/\\n$//'|seqkit seq -t dna -prvu",shell =True, stdout=sp.PIPE, stderr=sp.STDOUT, text = True)
            # remove '>seq' and '\n'
            rp = re.sub('>seq\n', '', rpseq1.stdout.strip())

            if pnum != ptotal_num :
                bpnum = pnum -1
                
                rpcutcom = f'cutadapt -e {error} --cores {cores} -a {rp} {db_dir}/temp{bpnum}.fas > {db_dir}/temp{pnum}.fas'
                sp.run(rpcutcom, shell = True)
                
                pnum += 1
            else:
                bpnum = pnum -1

                rpcutcom =  f'cutadapt -e {error} --cores {cores} -a {rp} {db_dir}/temp{bpnum}.fas > {lastout}'
                sp.run(rpcutcom, shell = True)

        # remove sequence of less than specified length
        input1 = f'{lastout}'
        out1 = f'{db_dir}/03_12srRNA_merge_L{minL}.fas'

        sh_seq1 = f'seqkit seq -g -m {minL} -w 0 --upper-case {input1} > {out1}'
        sp.run(sh_seq1, shell = True)

        # print DB status
        sqwatchcom = f'seqkit watch -B 5 --fields ReadLen {out1}'
        sp.run(sqwatchcom, shell = True)

        dt_now = datetime.now()
        et = elapsed_time(dt_start_f1, dt_now)
        print('Option command run finished.')
        print(' --- Done.ET:' + et + '\n')

        # extract sequence information
        out2 = f'{db_dir}/03_12srRNA_merge_L{minL}.tsv'
        sqnamecom = f"seqkit fx2tab -n {out1}|sed -z 's/\_/\t/g' > {out2}"
        sp.run(sqnamecom, shell = True)

    else:
        dt_now = datetime.now()
        et = elapsed_time(dt_start_f1, dt_now)

        print('Option command passed.')
        print('')
        print(' --- Done.ET:' + et + '\n')
        pass

