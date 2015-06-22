import os
from datetime import datetime
import logging

from massspechistory import datafile, morpheus

"""
Finds Thermo .raw files and automatically runs a Morpheus MS/MS
search against either an E. Coli or Human database.
"""


raw_files_dir = os.path.join('..', 'qeplus', "instrument_data")
results_dir = os.path.join('..', 'qeplus')
log = os.path.join('..', 'qeplus', 'web', 'run.log')

logging.basicConfig(filename=log,level=logging.INFO)
logger = logging.getLogger('process_raw_files')

logger.info('>>> Checking raw files on ' + str(datetime.now()))

for morpheus_dir, regex, db in [
        ('ecoli_morpheus', r'(E|e)coli*\d{12}', "E_coli_uniprot_iRT.fasta"), 
        ('hela_morpheus', r"(Hela|hela).*\d{12}", "HUMAN.fasta")]:
    out_dir = os.path.join(results_dir, morpheus_dir)
    morpheus.batch(
        datafile.glob_re(
            os.path.join(raw_files_dir, "*.raw"),
            regex),
        os.path.join("db", db),
        lambda f: os.path.join(out_dir, datafile.get_base(f)),
        dummy=False,
    )


