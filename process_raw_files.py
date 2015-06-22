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
web_dir = os.path.join('..', 'qeplus', 'web')
if not os.path.isdir(web_dir):
    os.makedirs(web_dir)
log = os.path.join(web_dir, 'run.log')

logging.basicConfig(
    level=logging.INFO, 
    filename="../qeplus/web/run.log",
    format='%(asctime)s|%(name)s|%(levelname)s|%(message)s',
    datefmt='%Y-%m-%d|%H:%M:%S')

logging.info('>>> Run Morpheus search on raw files')

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


