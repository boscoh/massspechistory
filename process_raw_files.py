import os
import glob
import re
from datetime import datetime
import logging

from massspechistory import datafile, morpheus

"""
Finds Thermo .raw files and automatically runs a Morpheus MS/MS
search against either an E. Coli or Human database.
"""

logger = logging.getLogger('process_raw_files')


for instrument in ['qeclassic', 'qeplus']:
    raw_files_dir = os.path.join('..', instrument, "instrument_data")
    results_dir = os.path.join('..', instrument)
    web_dir = os.path.join('..', instrument, 'web')
    if not os.path.isdir(web_dir):
        os.makedirs(web_dir)
    log = os.path.join(web_dir, 'run.log')

    logging.basicConfig(
        level=logging.INFO, 
        filename=log,
        format='%(asctime)s|%(name)s|%(levelname)s|%(message)s',
        datefmt='%Y-%m-%d|%H:%M:%S')

    logger.info('Checking raw files on %s' % instrument)

    for morpheus_dir, regex, db in [
            ('hela_morpheus', r"\b(Hela|hela).*\d{12}", "db/HUMAN.fasta"),
            ('ecoli_morpheus', r"\b(E|e)coli.*\d{12}", "db/E_coli_uniprot_iRT.fasta"), 
        ]:

        out_dir = os.path.join(results_dir, morpheus_dir)
        out_dir_fn = lambda f: os.path.join(out_dir, datafile.get_base(f))

        raw_files =  datafile.glob_re(
            os.path.join(raw_files_dir, "*.raw"),
            regex)

        raw_files = [f for f in raw_files if 'corrupt' not in f.lower()]
        morpheus.batch(
            raw_files,
            out_dir_fn,
            {
                '-db': db, 
                '-ad': 'true', 
                '-mmu': 'true',
                '-precmtv': '20',
                '-precmtu': 'ppm',
                '-prodmtv': '20',
                '-prodmtu': 'ppm',
                '-pmc': 'true',
                '-minpmo': '-3',
                '-maxpmo': '+1',
            }
        )

    root = logging.getLogger()
    map(root.removeHandler, root.handlers[:])
