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


logger = logging.getLogger('update_website')


def make_chart_data(data_dir, website_dir, title, description):
    if not os.path.isdir(website_dir):
        os.makedirs(website_dir)
    datafile.copy_dir('massspechistory/template.web', website_dir)

    charts = []


    morpheus_yaml = os.path.join(website_dir, 'msms.logs.yaml')
    logs = chart.parse_logs(
        glob.glob(
            os.path.join(data_dir, '*_morpheus/*/summary.tsv')),
        chart.parse_morpheus_summary, 
        morpheus_yaml)

    charts.append(chart.make_chart(
        logs, 
        [['Ecoli Spectra', ['Ecoli MS/MS Spectra']],
         ['Hela Spectra', ['Hela MS/MS Spectra']]], 
        'Digest MS/MS Identification Count'))
    charts.append(chart.make_chart(
        logs, 
        [['Ecoli PSM', ['Ecoli Target PSMs']],
         ['Hela PSM', ['Hela Target PSMs']]], 
        'Digest Peptide-Spectrum-Match Count'))
    charts.append(chart.make_chart(
        logs, 
        [['Ecoli Peptides', ['Ecoli Unique Target Peptides']],
         ['Hela Peptides', ['Hela Unique Target Peptides']]], 
        'Digest Peptide Count'))
    charts.append(chart.make_chart(
        logs, 
        [['Ecoli Proteins', ['Ecoli Target Protein Groups']],
         ['Hela Proteins', ['Hela Target Protein Groups']]], 
        'Digest Protein Count'))


    morpheus_yaml = os.path.join(website_dir, 'psm.logs.yaml')
    logs = chart.parse_logs(
        glob.glob(
            os.path.join(data_dir, '*_morpheus/*/*PSMs.tsv')),
        chart.parse_morpheus_psm, 
        morpheus_yaml)
    charts.append(chart.make_chart(
        logs, 
        [['Ecoli Mass Error', ['Ecoli Precursor Mass Error (Da)']],
         ['Hela Mass Error', ['Hela Precursor Mass Error (Da)']]], 
        'Average Absolute Precursor Mass Error [Da]'))


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


    datafile.write_jsonp(
        charts, 
        os.path.join(website_dir, 'load_charts.jsonp'),
        'load_charts')

    datafile.write_jsonp(
        { 'title': title, 'description': description },
        os.path.join(website_dir, 'load_title.jsonp'),
        'load_title')


def check_timepoints_for_outliers(website_dir, recipients=[]):
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

for instrument, description in [
        ('qeplus', 'Monash Proteomics Facility. Thermo QExactive Plus'),
        ('qeclassic', 'Monash Proteomics Facility. Thermo QExactive')
        ]:
    data_dir = "../" + instrument
    web_dir = "../%s/web" % instrument
    log = "../%s/web/run.log" % instrument
    recipients = [
        'apposite@gmail.com', 
        # 'oded.kleifeld@monash.edu', 
        # 'robert.goode@monash.edu', 
        # 'ralf.schittenhelm@monash.edu'
    ]

    if not os.path.isdir(web_dir):
        os.makedirs(web_dir)

    logging.basicConfig(
        level=logging.INFO, 
        filename=log,
        format='%(asctime)s|%(name)s|%(levelname)s|%(message)s',
        datefmt='%Y-%m-%d|%H:%M:%S')

    # logging.getLogger().addHandler(logging.StreamHandler())

    logger.info("Making chart data for " + instrument)

    make_chart_data(data_dir, web_dir, instrument, description)
    # check_timepoints_for_outliers( web_dir, recipients)

    root = logging.getLogger()
    map(root.removeHandler, root.handlers[:])


