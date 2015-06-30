
__doc__ = """
Utility module for various parsing, file-handling, data
manipulation routines
"""

import os
import glob
import re
import datetime
import copy
import json
import csv
import math
import shutil
from collections import defaultdict

import yaml


def write_json(logs, cache_json):
    with open(cache_json, 'w') as f:
        json.dump(logs, f)


def load_json(cache_json):
    with open(cache_json, 'Ur') as f:
        return json.load(f)


def write_jsonp(logs, cache_jsonp, callback='jsonp_callback'):
    text = '%s(\n%s\n);' % (callback, json.dumps(logs))
    with open(cache_jsonp, 'w') as f:
        f.write(text)


def load_yaml(cache_yaml):
    with open(cache_yaml, 'Ur') as f:
        return yaml.load(f)


def load_cache_yaml(cache_yaml):
    if os.path.isfile(cache_yaml):
        return load_yaml(cache_yaml)
    return {}


def to_dict(d):
    if not isinstance(d, defaultdict):
        return d
    result = {}
    for key, val in d.items():
        result[key] = to_dict(val)
    return result


def write_yaml(logs, cache_yaml):
    with open(cache_yaml, 'w') as f:
        yaml.safe_dump(
            to_dict(logs),
            f,
            encoding='utf-8',
            default_flow_style=False,
            allow_unicode=True)


def get_base(fname):
    return os.path.splitext(os.path.basename(fname))[0]


def get_date_from_fname(filename):
    strptime = datetime.datetime.strptime
    match = re.search("(\d{12})", filename)
    if match:
        return strptime(match.groups()[0], '%y%m%d%H%M%S')
    else:
        return None


float_regex_pattern = r"""
^
[-+]? # optional sign
(?:
    (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
    |
    (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
)
# followed by optional exponent part if desired
(?: [Ee] [+-]? \d+ ) ?
$
"""
float_regex = re.compile(float_regex_pattern, re.VERBOSE)


def parse_string(s):
    "Converts a string to a float or int if matches numerical pattern"
    if re.search(r'^[-+]?\d+$', s):
        return int(s)
    elif float_regex.match(s):
        return float(s)
    else:
        return s


def get_avg_std(values):
    n = len(values)
    avg = sum(values) / float(n)
    sum_sq = sum((v - avg)**2 for v in values)
    var = sum_sq / float(n)
    std = math.sqrt(var)
    return avg, std


def guess_delimiter(fname):
    if os.path.isfile(fname):
        with open(fname, 'Ur') as f:
            line = f.readline()
            if '\t' in line:
                return '\t'
            elif ',' in line:
                return ','
            else:
                raise IOError(
                    "Can't recognize delimiter in first line of '%s'" %
                    ext)
    else:
        if fname.endswith('.csv'):
            return ','
        elif fname.endswith('.txt'):
            return '\t'
        else:
            raise IOError("Can't recognize fname extension in '%s'" % ext)


def get_headers(fname):
    delimiter = guess_delimiter(fname)
    with open(fname, 'Ur') as f:
        reader = csv.reader(f, delimiter=delimiter)
        return reader.next()


def read_csv(fname):
    delimiter = guess_delimiter(fname)
    with open(fname, 'Ur') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_csv(rows, csv_fname):
    delimiter = guess_delimiter(csv_fname)
    with open(csv_fname, 'w') as f:
        writer = csv.writer(f, delimiter=delimiter)
        for row in rows:
            writer.writerow(row)


def glob_re(glob_tag, regex=None):
    fnames = glob.glob(glob_tag)
    if regex:
        fnames = [f for f in fnames if re.search(regex, f)]
    return fnames


def copy_dir(in_dir, out_dir):
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    for f in glob.glob(os.path.join(in_dir, '*')):
        if os.path.isfile(f):
            shutil.copy(f, out_dir)
        else:
            target_dir = os.path.join(out_dir, os.path.basename(f))
            copy_dir(f, target_dir)
