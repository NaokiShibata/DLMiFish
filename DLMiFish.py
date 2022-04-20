# -*- coding:utf-8 -*-

# import module
import subprocess as sp
from datetime import datetime

dt_start_f1 = datetime.now()
today = '{0:%Y%m%d}'.format(dt_start_f1)

# com run
sp.run(["python3", "./scripts/script1.py"]) # Search target taxon Genbank IDs
sp.run(["python3", "./scripts/script2.py", f'{today}']) # Download Genbank files
sp.run(["python3", "./scripts/script3.py"]) # Extract 12s rRNA from .gb files
sp.run(["python3", "./scripts/script4.py", f'{today}']) # Filterlation seq and seq info.