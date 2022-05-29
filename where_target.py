#! /usr/bin/env python3
# https://eco.mtk.nao.ac.jp/cgi-bin/koyomi/cande/horizontal.cgi
import os
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
conf = None
conf_file = Path(__file__).parent/'where_target_config.json'
if conf_file.is_file():
    with open(conf_file, 'r') as f:
        conf = json.load(f)
        if 'targets' in conf:
            for target in conf['targets']:
                ttype = target['type']
                if ttype in add_targets:
                    add_targets[ttype].append(target)
                else:
                    add_targets[ttype] = [target]

class EphemWrapper(ephem.FixedBody):
    def __init__(self, type, raster=False):
        super().__init__()
        self.type = type
        self.raster = raster

def make_fixed_target(name, type, ra, dec, raster):
    # https://rhodesmill.org/pyephem/quick.html#catalog-format
    # http://www.clearskyinstitute.com/xephem/help/xephem.html#mozTocId468501
    # line = '{},f,{},{},{},2000'.format(name, ra, dec, magnitude)
    # star = ephem.readdb(line)
    star = EphemWrapper(type, raster)
    star.name = name
    star._ra = ephem.hours(ra)
    star._dec = ephem.degrees(dec)
    star._epoch = '2000'
    return star

def get_targets(targets):
    res = []
    if 'all' in targets or 'sun' in targets:
        _sun = ephem.Sun()
        _sun.type = 'sun'
        res.append(_sun)
    if 'all' in targets or 'moon' in targets:
        _moon = ephem.Moon()
        # _moon.type = 'moon'
        res.append(_moon)
    if 'all' in targets or 'planets' in targets:
        _mars = ephem.Mars()
        # _mars.type = 'planets'
        res.append(_mars)
        _jupiter = ephem.Jupiter()
        # _jupiter.type = 'planets'
        res.append(_jupiter)
        _sat = ephem.Saturn()
        # _sat.type = 'planets'
        res.append(_sat)
        _venus = ephem.Venus()
        # _venus.type = 'planets'
        res.append(_venus)
    for ttype in add_targets:
        if 'all' in targets or ttype in targets:
            for add_target in add_targets[ttype]:
                if 'raster' in add_target:
                    raster = add_target['raster']
                else:
                    raster = False
                res.append(make_fixed_target(add_target['name'], ttype, add_target['ra'], add_target['dec'], raster))
    return res

def main(args):
    # set observation place
    if conf is not None and 'obs' in conf:
        obs = ephem.Observer()
        obs.lat = conf['obs']['lat']
        obs.lon = conf['obs']['lon']
    else:
        obs = ephem.Observer()
        obs.lat = "-22.9579"
        obs.lon = "-67.7862"

    # set time
    if args.start is None:
        # start_time = datetime.now(tz=timezone.utc)
        start_time = datetime.utcnow()
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

    # save directory
    savedir = Path(args.output)
    if not savedir.is_dir():
        os.makedirs(savedir)
        print('make directory {}'.format(savedir))

    # calculate
    if conf is not None and 'conditions' in conf:
        cond = conf['conditions']
    else:
        cond = {}
    targets = get_targets(args.targets)
    ## for table
    table_label = ['time (UTC)']+[target.name for target in targets]
    table_contents = []
    ## for az/el plot
    targets_az = {}
    targets_el = {}
    ## observable check
    if 'el_min' in cond:
        el_min = cond['el_min']*np.pi/180.
    else:
        el_min = 30.*np.pi/180.
    if 'el_max' in cond:
        el_max = cond['el_max']*np.pi/180.
    else:
        el_max = 180.*np.pi/180.
    if 'az_min' in cond:
        az_min = cond['az_min']*np.pi/180.
    else:
        az_min = 0.*np.pi/180.
    if 'az_max' in cond:
        az_max = cond['az_max']*np.pi/180.
    else:
        az_max = 360.*np.pi/180.
    if 'sun_separation' in cond:
        sun_thd = cond['sun_separation']*np.pi/180.
    else:
        sun_thd = 5*np.pi/180.
    if 'moon_separation' in cond:
        moon_thd = cond['moon_separation']*np.pi/180.
    else:
        moon_thd = 5*np.pi/180.
    targets_obsable = {}
    ## raster scan range
    if 'raster_az_offset' in cond:
        d_az = cond['raster_az_offset']*np.pi/180.
    else:
        d_az = 2.*np.pi/180.
    if 'raster_el_offset' in cond:
        d_el = cond['raster_el_offset']*np.pi/180.
    else:
        d_el = 2.*np.pi/180.
    sun = ephem.Sun()
    moon = ephem.Moon()
    print('el_min: {:.2f} degree'.format(el_min/np.pi*180.))
    print('el_max: {:.2f} degree'.format(el_max/np.pi*180.))
    print('az_min: {:.2f} degree'.format(az_min/np.pi*180.))
    print('az_max: {:.2f} degree'.format(az_max/np.pi*180.))
    print('sun_separation: {:.2f} degree'.format(sun_thd/np.pi*180.))
    print('moon_separation: {:.2f} degree'.format(moon_thd/np.pi*180.))
    print('raster_az_offset: {:.2f} degree'.format(d_az/np.pi*180.))
    print('raster_el_offset: {:.2f} degree'.format(d_el/np.pi*180.))

    for target in targets:
        targets_az[target.name] = []
        targets_el[target.name] = []
        targets_obsable[target.name] = []

    times = np.arange(datetime.timestamp(start_time), datetime.timestamp(end_time), interval*60)
    for time in times:
        date = datetime.fromtimestamp(time, tz=start_time.tzinfo)
        obs.date = date

        line_contents = [date.strftime('%Y/%m/%d %H:%M')]
        for target in targets:
            target.compute(obs)
            line_contents.append('az: {:05.1f} deg | el: {:+05.1f} deg'.format(deg*target.az, deg*target.alt))
            targets_az[target.name].append(target.az)
            targets_el[target.name].append(target.alt)
            sun.compute(obs)
            moon.compute(obs)
            targets_obsable[target.name].append(\
                    el_min<=target.alt<=el_max and\
                    az_min<=target.az<=az_max and\
                    ephem.separation(target, sun)>sun_thd and\
                    ephem.separation(target, moon)>moon_thd)
        table_contents.append(line_contents)
    for target in targets:
        targets_az[target.name] = np.array(targets_az[target.name])
        targets_el[target.name] = np.array(targets_el[target.name])

    cmap = plt.get_cmap('tab10')
    fs = 15
    plt.rcParams['font.size'] = fs
    print_times = times[[0, int(len(times)/3), int(len(times)*2/3), -1]]
    # make table
    table_str = tabulate(table_contents, table_label, tablefmt='orgtbl', stralign='center')
    if False:
        print(table_str)
        weight = len(targets)*3+2.3
        colWidths = np.concatenate(([2.3], np.ones(len(targets))*3))/weight
        fig1 = plt.figure(figsize=(weight, len(times)/3))
        ax11 = fig1.add_subplot(111)
        ax11.axis('off')
        table = ax11.table(cellText=table_contents, colLabels=table_label, rowLabels=None, loc='center', colLoc='center', rowLoc='center', colWidths=colWidths, fontsize=fs)
        table.auto_set_font_size(False)
        for pos, cell in table.get_celld().items():
            cell.set_height(1/len(times))

        # plot table
        fig1.tight_layout(rect=[0.05, 0.05, 0.95, 0.9])  # left, bot, right, top
        fig1.savefig(savedir/'time_table.pdf')

    with open(savedir/'table.txt', 'w') as f:
        f.write(table_str)
    if args.show:
        print(table_str)

    # az plot
    fig2 = plt.figure(figsize=(16/1.5, 9/1.5))
    ax21 = fig2.add_subplot(111)
    for i,target in enumerate(targets):
        tname = target.name
        ax21.plot(times, np.array(targets_az[tname])*deg, '-', color=cmap(i), label=tname)
        if hasattr(target, 'raster') and target.raster:
            ax21.fill_between(times, (targets_az[tname]-d_az)*deg, (targets_az[tname]+d_az)*deg, alpha=0.3, color=cmap(i), label=None)
        ax21.hlines(az_min*deg, times.min(), times.max(), colors='gray', ls='--')
        ax21.hlines(az_max*deg, times.min(), times.max(), colors='gray', ls='--')
    ax21.set_xticks(print_times)
    ax21.set_xticklabels([datetime.fromtimestamp(t, tz=start_time.tzinfo).strftime('%Y/%m/%d\n%H:%M') for t in print_times])
    ax21.set_ylabel('azimuth [deg]')
    fig2.legend()
    fig2.savefig(savedir/'az_plot.pdf')

    # el plot
    fig3 = plt.figure(figsize=(16/1.5, 9/1.5))
    ax31 = fig3.add_subplot(111)
    for i,target in enumerate(targets):
        tname = target.name
        ax31.plot(times, np.array(targets_el[tname])*deg, '-', color=cmap(i), label=tname)
        if hasattr(target, 'raster') and target.raster:
            ax31.fill_between(times, (targets_el[tname]-d_el)*deg, (targets_el[tname]+d_el)*deg, alpha=0.3, color=cmap(i), label=None)
        ax31.hlines(el_min*deg, times.min(), times.max(), colors='gray', ls='--')
        ax31.hlines(el_max*deg, times.min(), times.max(), colors='gray', ls='--')
    ax31.set_xticks(print_times)
    ax31.set_xticklabels([datetime.fromtimestamp(t, tz=start_time.tzinfo).strftime('%Y/%m/%d\n%H:%M') for t in print_times])
    ax31.set_ylabel('elevation [deg]')
    fig3.legend()
    fig3.savefig(savedir/'el_plot.pdf')

    # observable flag plot
    fig4 = plt.figure(figsize=(16/1.5, 9/1.5))
    ax41 = fig4.add_subplot(111)
    for i,target in enumerate(targets):
        ax41.plot(times, np.where(targets_obsable[target.name], len(targets)-i, np.nan), '-', lw=4)
    ax41.set_xticks(print_times)
    ax41.set_xticklabels([datetime.fromtimestamp(t, tz=start_time.tzinfo).strftime('%Y/%m/%d\n%H:%M') for t in print_times])
    ax41.set_xlim([times.min(), times.max()])
    ax41.set_yticks([len(targets)-i for i in range(len(targets))])
    ax41.set_yticklabels(targets_obsable.keys())
    ax41.set_ylim([0.5, len(targets)+0.5])
    fig4.savefig(savedir/'obs_flags.pdf')

    # pointing of stars
    fig5 = plt.figure(figsize=(16/1.5, 9/1.5))
    ax51 = fig5.add_subplot(111, projection='mollweide')
    for i,target in enumerate(targets):
        lat = np.array(targets_az[target.name])
        lat = np.where(lat>np.pi, lat-2*np.pi, lat)
        lon = np.array(targets_el[target.name])
        ax51.plot(lat[1:-1], lon[1:-1], '.', color=cmap(i), label=target.name)
        # ax51.plot(lat, lon, '--', color=cmap(i), label=target.name)
        ax51.plot([lat[0]], [lon[0]], '*', color=cmap(i), ms=6)
        ax51.plot([lat[-1]], [lon[-1]], 'x', color=cmap(i), ms=6)
        if hasattr(target, 'raster') and target.raster:
            for j in range(len(lat)):
                ax51.fill_between([lat[j]-d_az, lat[j]+d_az], [lon[j]+d_el, lon[j]+d_el], [lon[j]-d_el, lon[j]-d_el], alpha=0.3, color=cmap(i))
    fig5.legend()
    ax_pos = ax51.get_position()
    fig5.text(ax_pos.x1, ax_pos.y0+0.1, 'start: *\nend: x')
    fig5.savefig(savedir/'mollweide.pdf')

    if args.show:
        plt.show()
    plt.close('all')
    print('saved files at {}'.format(savedir))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-targets', help='target', default=['all'], type=str, nargs='*', choices=['all', 'sun', 'moon', 'planets']+list(add_targets.keys()))
    parser.add_argument('-start', help='start time (UTC): year/month/day or year/month/day-hour:min', type=str)
    parser.add_argument('-end', help='end time (UTC): year/month/day or year/month/day-hour:min', type=str)
    parser.add_argument('-interval', help='interval time [minutes]', type=float)
    parser.add_argument('-show', help='show the made plots.', action='store_true')
    parser.add_argument('-o', '--output', help='output directory', type=str, default=Path(__file__).parent/'tmp')
    args = parser.parse_args()
    main(args)

