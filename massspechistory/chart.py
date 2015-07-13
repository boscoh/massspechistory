# -*- coding: utf-8 -*-


__doc__ = """
Generates chart-data from MS/MS files of Ecoli/Hela digest,
and peak properties of the iRT peptide system. As well,
takes a list of most common peptides in Ecoli/Hela digest,
and calculates the fraction of common peptides.
"""

import os
import glob
import calendar
import logging
import pytz

import datafile


logger = logging.getLogger('chart')


def parse_morpheus_summary(fname):
    result = {}
    if 'hela' in fname:
        celltype = 'Hela'
    elif 'ecoli' in fname:
        celltype = 'Ecoli'
    for entry in datafile.read_csv(fname):
        for key, val in entry.items():
            key = celltype + ' ' + key
            if not key:
                result[key] = None
            else:
                result[key] = datafile.parse_string(val)
        break
    return result


def parse_morpheus_psm(fname):
    if 'hela' in fname:
        celltype = 'Hela'
    elif 'ecoli' in fname:
        celltype = 'Ecoli'
    key = 'Precursor Mass Error (Da)'
    values = []
    for entry in datafile.read_csv(fname):
        if key not in entry:
            continue
        values.append(abs(float(entry[key])))
    if len(values) == 0:
        return {}
    else:
        param = celltype + ' ' + key
        avg, std = datafile.get_avg_std(values)
        return { param: avg }


def calculate_crt(peptides):
    dRT = None
    if peptides['pep_b']['RT'] and peptides['pep_k']['RT']:
        cRT0 = peptides['pep_b']['RT']
        cRT100 = peptides['pep_k']['RT']
        dRT = cRT100 - cRT0
    for peptide in peptides.values():
        peptide['crt'] = None
        if peptide['RT'] and dRT:
            peptide['crt'] = (peptide['RT'] - cRT0) / dRT * 100


def parse_irt_log(fname):
    format_peptide_id = lambda i: 'pep_' + i[-1]

    log = {'peptides': {}}
    params = []
    for line in open(fname):
        if line.startswith("Component Name "):
            tokens = line.split(';')
            params = line.strip().split(" ; ")[1:]
            continue

        if line.startswith("iRT-pep_"):
            tokens = line.split(" ; ")

            peptide_id = format_peptide_id(tokens[0])
            if peptide_id not in log['peptides']:
                log['peptides'][peptide_id] = {}
            peptide = log['peptides'][peptide_id]

            for param, val in zip(params, tokens[1:]):
                if 'N/A' in val:
                    val = None
                elif param == "System Suitability":
                    val = val.strip()
                else:
                    val = float(val)
                peptide[param] = val

    calculate_crt(log['peptides'])

    return log


def parse_logs(fnames, parse_fn, cache_yaml):
    if os.path.isfile(cache_yaml):
        logs = datafile.load_yaml(cache_yaml)
    else:
        logs = []

    processed_fnames = [log['fname'] for log in logs]

    for fname in fnames:
        if fname in processed_fnames:
            continue

        date = datafile.get_date_from_fname(fname)
        if date is None:
            continue

        try:
            log = {
              'fname': fname,
              'timestamp': calendar.timegm(date.timetuple()),
              'iso_date_str': date.isoformat(),
            }
            logs.append(log)
            log.update(parse_fn(fname))
        except KeyboardInterrupt:
            raise
        except Exception as E:
            logger.debug("Error parsing " + fname + ": " + str(E))
            continue

    logs.sort(key=lambda l:l['timestamp'])
    datafile.write_yaml(logs, cache_yaml)

    return logs


def get_param(log, params):
    result = log
    for param in params:
        result = result[param]
    return result


def make_chart(logs, params, title, description=''):
    result = {
        "title": title,
        "description": description,
        "chart_data": [],
    }
    for param in params:
        name, keys = param
        chart = { "key": name, "values": [] }
        for log in logs:
            try:
                try:
                    y = get_param(log, keys)
                except:
                    y = None
                x = log['timestamp']*1000  # -> milliseconds for js 
                chart['values'].append([x, y])
            except:
                pass
        result['chart_data'].append(chart)
    return result


def make_pep_id_params(pep_ids, param):
    pep_id_params = []
    for pep_id in sorted(pep_ids):
        pep_id_params.append([pep_id, ['peptides', pep_id, param]])
    return pep_id_params


