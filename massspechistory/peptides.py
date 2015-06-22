# -*- coding: utf-8 -*-


__doc__ = """
Identifies the 300 common peptides found in a bunch of
MS/MS search results in Morpheus summary.txt files.

For each peptide, the mean and standard deviation of the
fraction-of-ion-products and ms-intensities are saved.

Analyzes the morpheus summary files of a list
of experiments and then identifies the most commonly matched
peptides with mean and standard deviation of the fraction of
products and intensities identified and then saves this to
log/msms.json.
"""

import os
import glob
import logging

import datafile


logger = logging.getLogger('peptides')


def print_pep(pep):
    logger.debug("n_log=%d n_psm=%d, ion=%.2f(%.2f) intensity=%.2f(%.2f) seq=%s" % \
        (pep['n_log'],
         pep['n_psm'],
         pep['ion_avg'],
         pep['ion_stdv'],
         pep['intensity_avg'],
         pep['intensity_stdv'],
         pep['sequence']))


def get_peptide_by_seq(fnames):
    times = []
    peptide_by_seq = {}
    for fname in fnames:
        date = datafile.get_date_from_fname(fname)

        logger.debug("Reading peptides in %s" % date)
        seqs = []
        for entry in datafile.read_csv(fname):
            seq = entry['Peptide Sequence']
            if seq not in peptide_by_seq:
                peptide_by_seq[seq] = {
                    'sequence': seq,
                    'base_sequence': entry['Base Peptide Sequence'],
                    'intensity_fractions': [],
                    'ion_fractions': [],
                    'n_log': 0
                }
            peptide = peptide_by_seq[seq]
            peptide['intensity_fractions'].append(
                float(entry['Fraction of Intensity Matching']))
            peptide['ion_fractions'].append(
                float(entry['Ratio of Matching Products']))
            seqs.append(seq)

        for seq in set(seqs):
            peptide_by_seq[seq]['n_log'] += 1

    for peptide in peptide_by_seq.values():
        ion_avg, ion_stdv = datafile.get_avg_std(peptide['ion_fractions'])
        intensity_avg, intensity_stdv = datafile.get_avg_std(
            peptide['intensity_fractions'])
        peptide.update({
            'ion_avg': ion_avg,
            'ion_stdv': ion_stdv,
            'intensity_avg': intensity_avg,
            'intensity_stdv': intensity_stdv,
            'n_psm': len(peptide['intensity_fractions']),
        })

    return peptide_by_seq


def extract_top_peptides(peptide_by_seq, n):
    peptides = peptide_by_seq.values()
    peptides.sort(key=lambda p: p['intensity_avg'], reverse=True)
    peptides.sort(key=lambda p: p['n_log'], reverse=True)
    top_peptides = {}
    for peptide in peptides:
        if peptide['intensity_avg'] < 0.1:
            # for some entries with 0 intensities
            continue
        print_pep(peptide)
        top_peptides[peptide['sequence']] = peptide
        if len(top_peptides) >= n:
            break
    return top_peptides


def find_top_peptides(psm_fnames, out_base, n_peptide):
    psm_fnames_json = out_base + '.fnames.json'
    if os.path.isfile(psm_fnames_json):
        prev_psm_fnames = datafile.load_json(psm_fnames_json)
    else:
        datafile.write_json(psm_fnames, psm_fnames_json)
        prev_psm_fnames = []

    if psm_fnames == prev_psm_fnames:
        return

    peptide_by_seq = get_peptide_by_seq(psm_fnames)
    top_peptides = extract_top_peptides(peptide_by_seq, n_peptide)
    datafile.write_json(top_peptides, out_base + '.json')



