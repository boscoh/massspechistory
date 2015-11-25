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
            if not key:
                result[key] = None
            else:
                result[celltype + ' ' + key] = datafile.parse_string(val)
        break
    return result


def parse_morpheus_psm(fname):
    if 'hela' in fname:
        celltype = 'Hela'
    elif 'ecoli' in fname:
        celltype = 'Ecoli'
    values = []
    key = 'Precursor Mass Error (ppm)'
    for entry in datafile.read_csv(fname):
        if entry['Target?'].lower() != "true":
            continue
        if float(entry['Q-Value (%)']) > 1:
            continue
        value = float(entry[key])
        if abs(value) > 20:
            continue
        values.append(value)
    if len(values) == 0:
        return {}
    else:
        average_param = celltype + ' ' + key
        upper_param = celltype + ' ' + key + ' Upper'
        avg, std = datafile.get_avg_std(values)
        result = {
            average_param: avg,
            upper_param: avg + std
        }
        return result


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

def summarise_irt(log):
#    print "Summarising iRT"
    iRT_dict = {}
    for l in ['R','T','O','B','D','S','W','N','C']:
        iRT_dict[l]=0
#    print iRT_dict
#    print log['peptides'].keys()
    for pep in log['peptides'].keys():
#        print pep
        if not log['peptides'][pep]['System Suitability']:
#            print "No suitability found for: ",pep
            continue
#        print log['peptides'][pep].keys()
        if 'Failed' in log['peptides'][pep]['System Suitability']:
#            print "Found failed"
            fs = log['peptides'][pep]['System Suitability'].split("; ")[-1].split()[1:]
            for f in fs:
                iRT_dict[f]+=1
#    for k,v in iRT_dict.iteritems():   # use for removing 0 points
#        if v==0:
#            iRT_dict[k]=None
#    print iRT_dict
    return iRT_dict

def parse_irt_log(fname):
    format_peptide_id = lambda i: 'pep_' + i[-1]

    log = {'peptides': {}}
    params = []
    for line in open(fname):
        if line.startswith("Component Name "):
            tokens = line.split(';')
            params = line.strip().split(" ; ")[1:]
#            print params
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
#    print log['peptides']
    log['irt_summ'] = summarise_irt(log)
    calculate_crt(log['peptides'])
#    print log.keys()
    return log

def parse_logs(fnames, parse_fn, cache_yaml):
    if os.path.isfile(cache_yaml):
        logs = datafile.load_yaml(cache_yaml)
    else:
        logs = []
    logs = filter(lambda log: 'fname' in log, logs)
    processed_fnames = [log['fname'] for log in logs]
#    print processed_fnames

    for fname in fnames:
        if fname in processed_fnames:
            print "processed",fname
            continue
        print fname
        date = datafile.get_date_from_fname(fname)
        print date
        if date is None:
            continue

        try:
            log = {
              'fname': fname,
              'timestamp': calendar.timegm(date.timetuple()),
              'iso_date_str': date.isoformat(),
            }
            parsed_log = parse_fn(fname)
            log.update(parsed_log)
            logs.append(log)
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
        values = [pair[1] for pair in chart['values']]
        values = filter(lambda v: v is not None, values)
        if len(values) > 0:
            result['chart_data'].append(chart)
    return result


def make_pep_id_params(pep_ids, param):
    pep_id_params = []
    for pep_id in sorted(pep_ids):
        pep_id_params.append([pep_id, ['peptides', pep_id, param]])
    return pep_id_params


