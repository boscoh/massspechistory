import os
import sys
import shutil
import glob
import logging
import platform

from massspechistory import chart
from massspechistory import datafile
from massspechistory import outliers




recipients = [
    'apposite@gmail.com', 
    'oded.kleifeld@monash.edu', 
    'robert.goode@monash.edu', 
    'david.steer@monash.edu',
]



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
        'Digest MS/MS Spectra'))
    charts.append(chart.make_chart(
        logs, 
        [['Ecoli PSM', ['Ecoli Target PSMs']],
         ['Hela PSM', ['Hela Target PSMs']]], 
        'Digest Peptide-Spectrum Matches'))
    charts.append(chart.make_chart(
        logs, 
        [['Ecoli Peptides', ['Ecoli Unique Target Peptides']],
         ['Hela Peptides', ['Hela Unique Target Peptides']]], 
        'Digest Unique Peptides'))
    charts.append(chart.make_chart(
        logs, 
        [['Ecoli Proteins', ['Ecoli Target Protein Groups']],
         ['Hela Proteins', ['Hela Target Protein Groups']]], 
        'Digest Protein Groups'))


    morpheus_yaml = os.path.join(website_dir, 'psm.logs.yaml')
    logs = chart.parse_logs(
        glob.glob(
            os.path.join(data_dir, '*_morpheus/*/*PSMs.tsv')),
        chart.parse_morpheus_psm, 
        morpheus_yaml)
    charts.append(chart.make_chart(
        logs, 
        [['Ecoli dMass avg', ['Ecoli Precursor Mass Error (ppm)']],
         ['Ecoli dMass avg+std', ['Ecoli Precursor Mass Error (ppm) Upper']],
         ['Hela dMass avg', ['Hela Precursor Mass Error (ppm)']], 
         ['Hela dMass avg+std', ['Hela Precursor Mass Error (ppm) Upper']]], 
        'Digest Precursor dMass [ppm]'))


    logs = chart.parse_logs(
        datafile.glob_re(
            os.path.join(data_dir, 'instrument_data', '*.txt'), 
            r'(Hela|hela|ecoli).*iRT'),
        chart.parse_irt_log, 
        os.path.join(website_dir, 'irt_peptides.logs.yaml'))
    if len(logs) > 0:
#        print logs[0].keys()
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
            [['Symmetry',['irt_summ','S']],
            ['Resolution',['irt_summ','R']],
            ['Peak Width',['irt_summ','W']],
            ['Tailing',['irt_summ','T']],
            ['Column Overload',['irt_summ','O']],
            ['Baseline Clipping',['irt_summ','B']],
            ['Sig to Noise',['irt_summ','N']],
            ['Concave',['irt_summ','C']],
            ['Saturation',['irt_summ','D']]],
            'System suitability flags',
            'Number of iRT peptides failing given suitability test'))
#        for log in logs[0:1]:
#            print log['peptides']

    datafile.write_jsonp(
        charts, 
        os.path.join(website_dir, 'load_charts.jsonp'),
        'load_charts')

    datafile.write_jsonp(
        { 'title': title, 'description': description },
        os.path.join(website_dir, 'load_title.jsonp'),
        'load_title')


def check_timepoints_for_outliers(website_dir, instrument, recipients=[]):
#    if platform.system() == 'Windows':
#        return

    timepoints_yaml = os.path.join(website_dir, 'timepoints.yaml')
    timepoints = datafile.load_cache_yaml(timepoints_yaml)

    # Calculate limits for all variables
    limit = {}
    logs = datafile.load_yaml(
        os.path.join(website_dir, 'msms.logs.yaml'))
    params_list = [
        ['Hela MS/MS Spectra'], 
        ['Hela Target PSMs'], 
        ['Hela Unique Target Peptides'],
        ['Hela Target Protein Groups'],
        ['Ecoli MS/MS Spectra'], 
        ['Ecoli Target PSMs'], 
        ['Ecoli Unique Target Peptides'],
        ['Ecoli Target Protein Groups'],
    ]
    for params in params_list:
        outliers.set_lower_limit_of_param(
            limit, params, params[0], logs, timepoints)

    logs = datafile.load_yaml(
        os.path.join(website_dir, 'irt_peptides.logs.yaml'))
    if len(logs) > 0:
        pep_ids = logs[0]['peptides'].keys()
        for pep_id in pep_ids:
            outliers.set_lower_limit_of_param(
                limit, ['peptides', pep_id, 'RT'], pep_id + '_RT', logs, timepoints)

    # Compare each data timepoint to limits
    bad_times = []
    for time in sorted(timepoints.keys()):
        time_point = timepoints[time]
        if 'considered' in time_point:
            continue
        bad_params = [p for p in time_point if not time_point[p]]
        bad_pep_params = [p for p in bad_params if p.startswith("pep")]
        print time,bad_params
        print time,bad_pep_params
        if len(bad_pep_params) <= 3:
            # if only 3 or less irt peptides are bad,
            # consider okay and remove from bad_params
            bad_params = [p for p in bad_params if not p.startswith("pep")] 
        if bad_params:
            bad_times.append(time)
        time_point['considered'] = True

    # Send reports
    if bad_times:
        message = outliers.bad_times_message(
            instrument, bad_times, timepoints, limit)
        print message
#        outliers.report_by_email(
#            instrument, message, recipients)

    datafile.write_yaml(timepoints, timepoints_yaml)



# Main Loop

# automatically any directories that has instrument_data for processing
for instrument in open('instruments.txt', 'Ur').read().split():

    description =  'Monash Proteomics Facility.'
    description_txt = os.path.join('..', instrument, 'description.txt')
    if os.path.isfile(description_txt):
        description = open(description_txt, 'Ur').read()

    data_dir = "../" + instrument
    web_dir = "../%s/web" % instrument
    log = "../%s/web/run.log" % instrument

    if not os.path.isdir(web_dir):
        os.makedirs(web_dir)

    logging.basicConfig(
        level=logging.INFO, 
        filename=log,
        format='%(asctime)s|%(name)s|%(levelname)s|%(message)s',
        datefmt='%Y-%m-%d|%H:%M:%S')

    # uncomment this to send ouptut to screen, otherwise silent for cron jobs
    logging.getLogger().addHandler(logging.StreamHandler())

    logger.info("Making chart data for " + instrument)

    make_chart_data(data_dir, web_dir, instrument, description)
    check_timepoints_for_outliers(web_dir, instrument, recipients)

    root = logging.getLogger()
    map(root.removeHandler, root.handlers[:])


