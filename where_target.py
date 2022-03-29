#! /usr/bin/env python3
# https://eco.mtk.nao.ac.jp/cgi-bin/koyomi/cande/horizontal.cgi
import sys
from datetime import datetime, timedelta, timezone
import json
import argparse
import re
from pathlib import Path

import ephem
import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate

deg = 180.0/np.pi

add_targets = {}
conf_file = Path(__file__).parent/'where_target_config.json'
if conf_file.is_file():
    with open(conf_file, 'r') as f:
        for target in json.load(f):
            ttype = target['type']
            if ttype in add_targets:
                add_targets[ttype].append(target)
            else:
                add_targets[ttype] = [target]

def make_fixed_target(name, ra, dec, magnitude):
    # https://rhodesmill.org/pyephem/quick.html#catalog-format
    # http://www.clearskyinstitute.com/xephem/help/xephem.html#mozTocId468501
    line = '{},f,{},{},{},2000'.format(name, ra, dec, magnitude)
    patch = ephem.readdb(line)
    return patch

def get_targets(targets):
    res = []
    if 'all' in targets or 'sun' in targets:
        res.append(ephem.Sun())
    if 'all' in targets or 'moon' in targets:
        res.append(ephem.Moon())
    if 'all' in targets or 'planets' in targets:
        res.append(ephem.Mars())
        res.append(ephem.Jupiter())
        res.append(ephem.Saturn())
        res.append(ephem.Venus())
    for ttype in add_targets:
        if 'all' in targets or ttype in targets:
            for add_target in add_targets[ttype]:
                res.append(make_fixed_target(add_target['name'], add_target['ra'], add_target['dec'], add_target['amplitude']))
    return res

def main(args):
    obs = ephem.Observer()
    obs.lat = "-22.9579"
    obs.lon = "-67.7862"
    # obs = ephem.city('Tokyo')

    if args.start is None:
        start_time = datetime.now(tz=timezone.utc)
    elif re.match('[0-9][0-9][0-9][0-9]/[0-9][0-9]?/[0-9][0-9]?-[0-9][0-9]?:[0-9][0-9]?', args.start):
        start_time = datetime.strptime(args.start, '%Y/%m/%d-%H:%M')
    elif re.match('[0-9][0-9][0-9][0-9]/[0-9][0-9]?/[0-9][0-9]?', args.start):
        start_time = datetime.strptime(args.start, '%Y/%m/%d')
    else:
        print('invalid time format (start time)', file=sys.stderr)
        return

    if args.end is None:
        end_time = start_time + timedelta(days=1)
    elif re.match('[0-9][0-9][0-9][0-9]/[0-9][0-9]?/[0-9][0-9]?-[0-9][0-9]?:[0-9][0-9]?', args.end):
        end_time = datetime.strptime(args.end, '%Y/%m/%d-%H:%M')
    elif re.match('[0-9][0-9][0-9][0-9]/[0-9][0-9]?/[0-9][0-9]?', args.end):
        end_time = datetime.strptime(args.end, '%Y/%m/%d')
    else:
        print('invalid time format (end time)', file=sys.stderr)
        return
    assert start_time < end_time

    if args.interval is None:
        interval = 60   # 60 minutes = 1 hour
    else:
        interval = args.interval

    targets = get_targets(args.targets)
    # for plot 1
    table_label = ['time (UTC)']+[target.name for target in targets]
    table_contents = []
    # for plot 2
    targets_az = {}
    for target in targets:
        targets_az[target.name] = []

    times = np.arange(datetime.timestamp(start_time), datetime.timestamp(end_time), interval*60)
    for time in times:
        date = datetime.fromtimestamp(time, tz=start_time.tzinfo)
        obs.date = date

        line_contents = [date.strftime('%Y/%m/%d %H:%M')]
        for target in targets:
            target.compute(obs)
            line_contents.append('az: {:05.1f} deg | el: {:+05.1f} deg'.format(deg*target.az, deg*target.alt))
            targets_az[target.name].append(target.az)
        table_contents.append(line_contents)

    # time table
    table_str = tabulate(table_contents, table_label, tablefmt='orgtbl', stralign='center')
    print(table_str)
    fs = 15
    fig1 = plt.figure(figsize=(len(targets)*3, len(times)/3))
    ax11 = fig1.add_subplot(111)
    ax11.axis('off')
    table = ax11.table(cellText=table_contents, colLabels=table_label, rowLabels=None, loc='center', colLoc='center', rowLoc='center', fontsize=fs)
    table.auto_set_font_size(False)
    for pos, cell in table.get_celld().items():
        cell.set_height(1/len(times))

    fig1.tight_layout(rect=[0.05, 0.05, 0.95, 0.9])  # left, bot, right, top
    fig1.savefig('tmp/targets1.pdf')

    # azimuth changes
    fig2 = plt.figure()
    ax21 = fig2.add_subplot(111)
    for tname in targets_az:
        ax21.plot(times, np.array(targets_az[tname])*deg, '-', label=tname)
    print_times = times[[0, int(len(times)/3), int(len(times)*2/3), -1]]
    ax21.set_xticks(print_times)
    ax21.set_xticklabels([datetime.fromtimestamp(t, tz=start_time.tzinfo).strftime('%Y/%m/%d\n%H:%M') for t in print_times])
    ax21.set_ylabel('azimuth [deg]')
    fig2.legend()
    fig2.savefig('tmp/targets2.pdf')

    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-targets', help='target', default=['all'], type=str, nargs='*', choices=['all', 'sun', 'moon', 'planets']+list(add_targets.keys()))
    parser.add_argument('-start', help='start time (UTC): year/month/day or year/month/day-hour:min', type=str)
    parser.add_argument('-end', help='end time (UTC): year/month/day or year/month/day-hour:min', type=str)
    parser.add_argument('-interval', help='interval time [minutes]', type=float)
    args = parser.parse_args()
    main(args)

