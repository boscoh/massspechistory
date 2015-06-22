import os
import sys
import shutil
import glob
import logging
import platform

from massspechistory import chart
from massspechistory import datafile
from massspechistory import outliers
from massspechistory import peptides



def find_top_peptides(data_dir, web_dir):
    logging.info(">>> Calculating 300 most common peptides in ecoli")
    for celltype in ['ecoli', 'hela']:
        base = os.path.join(web_dir, '%s_top_peptides' % celltype)
        search_tag = os.path.join(data_dir, '%s_morpheus/*/*PSMs.tsv' % celltype)
        peptides.find_top_peptides(datafile.glob_re(search_tag), base, 300)


def make_chart_data(data_dir, website_dir):
    logging.info(">>> Making chart data")
    if not os.path.isdir(website_dir):
        os.makedirs(website_dir)
    datafile.copy_dir('template.web', website_dir)

    charts = []

    for cell in ['Ecoli', 'Hela']:
        morpheus_dir = os.path.join(data_dir, '%s_morpheus' % cell.lower())
        morpheus_yaml = os.path.join(website_dir, '%s_msms.logs.yaml' % cell.lower())
        logs = chart.parse_logs(
            glob.glob(morpheus_dir + '/*/summary.tsv'),
            chart.parse_morpheus_summary, 
            morpheus_yaml)
        charts.append(chart.make_chart(
            logs, 
            [['MS/MS Spectra', ['MS/MS Spectra']]], 
            '%s MS/MS Spectra Count' % cell))
        charts.append(chart.make_chart(
            logs, 
            [['Target PSMs', ['Target PSMs']]], 
            '%s PSM Count' % cell))
        charts.append(chart.make_chart(
            logs, 
            [['Unique Target Peptides', ['Unique Target Peptides']]], 
            '%s Peptide Count' % cell))
        charts.append(chart.make_chart(
            logs, 
            [['Target Protein Groups', ['Target Protein Groups']]], 
            '%s Protein Count' % cell))

        top_peptides = datafile.load_json(
            os.path.join(website_dir, '%s_top_peptides.json'  % cell.lower()))
        logs = chart.parse_logs(
            glob.glob(morpheus_dir + '/*/*PSMs.tsv'),
            lambda f: chart.parse_top_peptides_of_morpheus_psm(f, top_peptides), 
            os.path.join(website_dir, '%s_top_peptides.logs.yaml' % cell.lower()))
        charts.append(chart.make_chart(
            logs, 
            [['n_top_peptide', ['n_top_peptide']]], 
            '%s Top 300 Peptide Count' % cell,
            'Peptides matching the MS/MS ions of top 300 Peptides'))

    logs = chart.parse_logs(
        datafile.glob_re(
            os.path.join(data_dir, 'instrument_data', '*.txt'), 
            r'(Hela|hela|ecoli).*iRT'),
        chart.parse_irt_log, 
        os.path.join(website_dir, 'irt_peptides.logs.yaml'))
    pep_ids = logs[0]['peptides'].keys()
    charts.append(chart.make_chart(
        logs, 
        chart.make_pep_id_params(pep_ids, 'RT'), 
        'iRT Peptides Retention Time',
        'The measured retention time of the peptides'))
    charts.append(chart.make_chart(
        logs, 
        chart.make_pep_id_params(pep_ids, 'Height'), 
        'iRT Peptides Peak Height', 
        'Height of the peak associated at the retention time'))
    charts.append(chart.make_chart(
        logs, 
        chart.make_pep_id_params(pep_ids, 'Width'), 
        'iRT Peptides Retention Width', 
        'Width of the peak associated at the retention time'))

    jsonp_fname = os.path.join(website_dir, 'load_charts.jsonp')
    datafile.write_jsonp(charts, jsonp_fname, 'load_charts')



def check_timepoints_for_outliers(website_dir, recipients=[]):
    logging.info(">>> Checking outliers")
    if platform.system() == 'Windows':
        return

    timepoints_yaml = os.path.join(website_dir, 'timepoints.yaml')
    timepoints = datafile.load_cache_yaml(timepoints_yaml)

    # Calculate limits for all variables
    limit = {}
    for cell in ['ecoli', 'hela']:
        logs = datafile.load_yaml(
            os.path.join(website_dir, '%s_msms.logs.yaml' % cell))
        params_list = [
            ['MS/MS Spectra'], 
            ['Target PSMs'], 
            ['Unique Target Peptides'],
            ['Target Protein Groups'],
        ]
        for params in params_list:
            outliers.set_lower_limit_of_param(
                limit, params, cell + "_" + params[0], logs, timepoints)

        logs = datafile.load_yaml(
            os.path.join(
                website_dir, '%s_top_peptides.logs.yaml' % cell))
        outliers.set_lower_limit_of_param(
            limit, ['n_top_peptide'], cell + "_" + 'n_top_peptide', logs, timepoints)

    logs = datafile.load_yaml(
        os.path.join(website_dir, 'irt_peptides.logs.yaml'))
    pep_ids = logs[0]['peptides'].keys()
    for pep_id in pep_ids:
        outliers.set_lower_limit_of_param(
            limit, ['peptides', pep_id, 'RT'], pep_id + '_RT', logs, timepoints)

    # Compare each data timepoint to limits
    bad_times = []
    for time in sorted(timepoints.keys()):
        time_point = timepoints[time]
        if not 'considered' in time_point:
            bad_params = [p for p in time_point if not time_point[p]]
            bad_pep_params = [p for p in bad_params if p.startswith("pep")]
            if len(bad_pep_params) <= 3:
                # if only 3 or less irt peptides are bad,
                # consider okay and remove from bad_params
                bad_params = [p for p in bad_params if not p.startswith("pep")] 
            if bad_params:
                bad_times.append(time)
            time_point['considered'] = True

    # Send reports
    if bad_times:
        message = outliers.bad_times_message(bad_times, timepoints, limit)
        outliers.report_by_email(
            message, 
            recipients)

    datafile.write_yaml(timepoints, timepoints_yaml)



# Main Loop

data_dir = "../qeplus"
web_dir = "../qeplus/web"
target_web_dir = '../web'
logging.basicConfig(
    level=logging.INFO, 
    filename="../qeplus/web/run.log",
    format='%(asctime)s|%(name)s|%(levelname)s|%(message)s',
    datefmt='%Y-%m-%d|%H:%M:%S')
logging.getLogger().addHandler(logging.StreamHandler())
recipients = [
    'apposite@gmail.com', 
    # 'oded.kleifeld@monash.edu', 
    # 'robert.goode@monash.edu', 
    # 'ralf.schittenhelm@monash.edu'
]
if not os.path.isdir(web_dir):
    os.makedirs(web_dir)
find_top_peptides(data_dir, web_dir)
make_chart_data(data_dir, web_dir)
check_timepoints_for_outliers( web_dir, recipients)
datafile.copy_dir(web_dir, target_web_dir)



