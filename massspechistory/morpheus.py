#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import shutil
import platform
import logging
import time

import datafile


logger = logging.getLogger('morpheus')

ROOT_DIR = os.path.dirname(__file__)
MORPHEUS_DIR = os.path.join(ROOT_DIR, 'morpheus', 'standard')
MORPHEUS_THERMO_DIR = os.path.join(ROOT_DIR, 'morpheus', 'thermo')


__doc__ = """
Morpheus proteomics search engine wrapper in python

If a 'modifications.tsv' file is found in the current directory,
will add that to the default 'modifictions.tsv' found in the
morpheus directory.

Usage: morpheus.py [options]

  -d          data (comma separated)
  -db         database

  -o          output_folder [default: ./]
  -ad         append_decoys [default: false]
  -minprecz   min_precursor_charge [default: 2]
  -maxprecz   max_precursor_charge [default: 4]
  -at         absolute_threshold [default: -1.0]
  -rt         relative_threshold [default: -1.0]
  -mp         maximum_peaks [default: 400]
  -acs        assign_charge_states [default: true]
  -di         deisotope [default: true]
  -p          protease [default: trypsin (no proline rule)" from "proteases.csv"]
  -mmc        maximum_missed_cleavages [default: 2]
  -imb        initiator_methionine_behaviour [default: Variable]; Retain or Cleave
  -fm         fixed_mods (comma separated)
  -vm         variable_mods (comma separated)
  -mvmi       max_variable_mod_isoforms_per_peptide [default: 1024]
  -precmtv    precursor_mass_tolerance_value [default: 2.1]
  -precmtu    precursor_mass_tolerance_units [default: Da]; or ppm
  -precmt     monoisotopic_precursor_mass_type [default: Monoisotopic]; or Average
  -pmc        precursor_mono_correction [default: false]
  -minpmo     min_prec_mono_offset [default: -3]
  -maxpmo     max_prec_mono_offset [default: 1]
  -prodmtv    product_mass_tolerance_value [default: 0.025]
  -prodmtu    product_mass_tolerance_units [default: Da]; or ppm
  -prodmt     product_mass_type [default: Monoisotopic]; or Average
  -fdr        maximum_fdr [default: 0.01]
  -cmu        consider_mods_unique [default: false]
  -mt         max_threads (default to number of processors)
  -mmu        minimise_memory_usage [default: false]
"""


def get_morpheus_bin(is_thermo=False):
    if not is_thermo:
        cmd = os.path.join(MORPHEUS_DIR, 'morpheus_cl.exe')
    else:
        cmd = os.path.join(MORPHEUS_THERMO_DIR, 'morpheus_tmo_cl.exe')
    if platform.system() != 'Windows':
        cmd = 'mono ' + cmd
    return cmd


def get_modifications_tsv(is_thermo=False):
    if not is_thermo:
        return os.path.join(MORPHEUS_DIR, 'modifications.tsv')
    else:
        return os.path.join(MORPHEUS_THERMO_DIR, 'modifications.tsv')


def add_modifications(modifications_tsv, extra_modifications_tsv):
    mods = list(datafile.read_csv(modifications_tsv))
    default_mods = [mod['Description'] for mod in mods]
    logger.debug("Adding %s to %s" % (extra_modifications_tsv, modifications_tsv))
    new_mod_lines = []
    with open(extra_modifications_tsv, 'Ur') as f:
        for line in f.readlines()[1:]:
            mod = line.split('\t')[0]
            if mod not in default_mods:
                new_mod_lines.append(line)
    with open(modifications_tsv, 'a') as f:
        for line in new_mod_lines:
            f.write(line)


def print_modifications(modifications_tsv):
    logger.debug('modifications.tsv: ' + os.path.relpath(modifications_tsv))
    mods = [g['Description'] for g in datafile.read_csv(modifications_tsv)]
    logger.debug('modifications: %s' % mods)


def run(options):
    is_thermo = False
    if '-d' in options:
        is_thermo = options['-d'].strip().lower().endswith('raw')

    logger.debug("options: %s" % options)

    if '-o' in options:
        d = options['-o']
        if d and not os.path.isdir(d):
            logger.debug("Creating output directory " + d)
            os.makedirs(d)

    default_modifications_tsv = get_modifications_tsv(is_thermo)

    backup_modifications_tsv = None
    if os.path.isfile('modifications.tsv'):
        backup_modifications_tsv = default_modifications_tsv + '.backup'
        shutil.copy(default_modifications_tsv, backup_modifications_tsv)
        add_modifications(default_modifications_tsv, 'modifications.tsv')

    print_modifications(default_modifications_tsv)

    option_str = ''
    for key in options:
        option_str += ' %s \'%s\'' % (key, options[key])

    cmd = get_morpheus_bin(is_thermo) + ' ' + option_str
    prompt = "morpheus: "
    if is_thermo:
        prompt = "morpheus(raw): "
    logger.info(prompt + datafile.get_base(options['-d']))
    logger.debug(cmd)
    os.system(cmd)

    if backup_modifications_tsv:
        shutil.copy(backup_modifications_tsv, default_modifications_tsv)


def get_args_from_doc(doc):
    words = [l.split()[0] for l in doc.split() if l]
    return [w for w in words if w.startswith('-')]


def get_options(words):
    options = {}
    curr_arg = None
    for word in words:
        if word.startswith('-'):
            curr_arg = word
            options[curr_arg] = ''
        elif curr_arg:
            if options[curr_arg] != '':
                options[curr_arg] += ' '
            options[curr_arg] += word
    return options


def is_good_morpheus_output(out_dir):
    fnames = os.listdir(out_dir)
    is_summary = any('summary.tsv' in f for f in fnames)
    is_psm = any('PSMs.tsv' in f for f in fnames)
    return is_summary and is_psm and len(fnames) >= 6


def batch(fnames, out_dir_fn, options={'-ad':'true','-mmu':'true'}, dummy=False):
    for fname in fnames:
        if not datafile.get_date_from_fname(fname):
            continue
        out_dir = out_dir_fn(fname)
        if os.path.isdir(out_dir):
            if is_good_morpheus_output(out_dir):
                logger.debug("Skipping " + datafile.get_base(out_dir))
                continue
        if dummy:
            continue
        try:
            start = time.time()

            if os.path.isdir(out_dir):
                shutil.rmtree(out_dir)

            params = {
              '-d': fname,
              '-o': out_dir
            }
            params.update(options)
            run(params)

            if not is_good_morpheus_output(out_dir):
                raise Exception

            end = time.time()
            c = end - start

            hours = c // 3600 % 24
            minutes = c // 60 % 60
            seconds = c % 60

            logger.info("finished in %d:%02d:%04.1f" % (hours, minutes, seconds))

        except KeyboardInterrupt:
            raise
        except:
            logger.error("failed: " + datafile.get_base(out_dir))



if __name__ == "__main__":

    if len(sys.argv) == 1:
        logger.info(__doc__)
    else:
        options = get_options(sys.argv[1:])
        args = get_args_from_doc(__doc__)
        for key in options:
            if key not in args:
                logger.info("Can't recognize option %s" % key)
                sys.exit(1)
        run(options)

