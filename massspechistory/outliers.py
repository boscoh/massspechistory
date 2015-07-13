import os
import platform
import collections
from pprint import pprint
from email.mime.text import MIMEText
from subprocess import Popen, PIPE
import logging

import datafile


logger = logging.getLogger('outliers')


def time_from_log(log):
    date = datafile.get_date_from_fname(log['fname'])
    return date.isoformat()


def get_param(log, params):
    result = log
    for param in params:
        result = result[param]
    return result


def set_lower_limit_of_param(limit, params, name, logs, timepoints):
    times = []
    values = []
    for log in logs:
        time = time_from_log(log)
        times.append(time)
        try:
            value = get_param(log, params)
            values.append(value)
        except:
            values.append(None)

    try:
        # first pass, includes outliers
        avg, std = datafile.get_avg_std(
            [v for v in values if v != None])

        # remove outliers
        lower_limit = avg - std
        avg, std = datafile.get_avg_std(
            [v for v in values if v != None and v > lower_limit])
    except:
        return

    limit[name] = {
        'avg': avg,
        'std': std,
        'lower_limit': avg - 2*std
    }

    for time, value in zip(times, values):
        if time not in timepoints:
            timepoints[time] = {}
        if value is None:
            timepoints[time][name] = False
        else:
            timepoints[time][name] = value > limit[name]['lower_limit']


def sendmail(to_address, from_address, subject, body):
    sendmail = "/usr/sbin/sendmail"
    if not os.path.isfile(sendmail):
        return
    msg = MIMEText(body)
    msg["To"] = to_address
    msg["From"] = from_address
    msg["Subject"] = subject
    p = Popen([sendmail, "-t", "-oi"], stdin=PIPE)
    p.communicate(msg.as_string())


def report_by_email(message, recipients=[]):
    logger.info('Sending outlier email to ' + ','.join(recipients))
    logger.debug(message)
    for recipient in recipients:
        sendmail(
            recipient, 
            'proteome@monash.edu', 
            'Warning: QC for QEPLUS', 
            message)


def bad_times_message(bad_times, timepoints, limit):
    result = "Hi, \n\n"
    result += "This is the automated QC for the QEPLUS.\n"
    result += "Full report: http://monash.edu/proteomics/qc/qeplus/index.html.\n\n"
    result += "Insufficent counts were encountered on:\n\n"

    for time in bad_times:
        result += " - %s:\n" % time.replace('T', ', ')
        for param in sorted(timepoints[time]):
            if not timepoints[time][param]:
                result += "   - %s\n" % param

    result = "\n\nThe reference counts ranges are:\n"

    for param in sorted(limit):
        result += " - %s\n" % param
        result += "    - range = %.f +/- %.f\n" % (limit[param]['avg'], limit[param]['std'])
        result += "    - lower-95%% = %.f\n" % limit[param]['lower_limit']

    return result



