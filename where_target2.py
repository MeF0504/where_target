#! /usr/bin/env python3

import os
import sys
import json
import argparse
from datetime import datetime, timedelta, timezone
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate

# from astropy.visualization import astropy_mpl_style, quantity_support
import astropy.units as u
from astropy.coordinates import get_sun, get_moon, get_body, \
    AltAz, EarthLocation, SkyCoord, solar_system_ephemeris
from astropy.time import Time

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


class Target:
    def __init__(self, name, type, body, raster, ces):
        self.name = name
        self.type = type
        self.body = body
        self.raster = raster
        self.ces = ces


def make_fixed_target(name, type, ra, dec, raster, ces):
    body = SkyCoord(ra=ra, dec=dec, frame='icrs', unit=(u.hourangle, u.deg))
    target = Target(name, type, body, raster, ces)
    return target


def get_targets(targets, time):
    res = []
    if 'all' in targets or 'sun' in targets:
        tmp = Target('Sun', 'sun', get_sun(time), False, False)
        res.append(tmp)
    if 'all' in targets or 'moon' in targets:
        tmp = Target('Moon', 'moon', get_moon(time), False, False)
        res.append(tmp)
    if 'all' in targets or 'planets' in targets:
        # list: solar_system_ephemeris.bodies
        with solar_system_ephemeris.set('builtin'):
            mars = Target('Mars', 'planet', get_body('mars', time),
                          False, False)
            jup = Target('Jupiter', 'planet', get_body('jupiter', time),
                         False, False)
            sat = Target('Saturn', 'planet', get_body('saturn', time),
                         False, False)
            venus = Target('Venus', 'planet', get_body('venus', time),
                           False, False)
            res += [mars, jup, sat, venus]
    for ttype in add_targets:
        if 'all' in targets or ttype in targets:
            for add_target in add_targets[ttype]:
                if 'raster' in add_target:
                    raster = add_target['raster']
                else:
                    raster = False
                if 'ces' in add_target:
                    ces = add_target['ces']
                else:
                    ces = False
                tmp = make_fixed_target(add_target['name'], ttype,
                                        add_target['ra'], add_target['dec'],
                                        raster, ces)
                res.append(tmp)
    return res


def chk_deg_condition(cond, item, default):
    if item in cond:
        return cond[item]/deg
    else:
        return default/deg


def calc_CES(target, obs, times, elevation):
    L = len(times)
    ok = np.zeros(L, dtype=bool)
    az = np.zeros((2, L))
    for i, t in enumerate(times):
        altaz = target.body.transform_to(AltAz(obstime=t, location=obs))
        alt_min = altaz.alt.min()/u.deg/deg
        alt_max = altaz.alt.max()/u.deg/deg
        alts = altaz.alt/u.deg/deg
        azs = altaz.az/u.deg/deg
        sort_idx = np.argsort(alts)
        if alt_min < elevation < alt_max:
            ok[i] = True
            if elevation < alts[sort_idx[-1]] and \
               alts[sort_idx[-2]] < elevation:
                #  _/\_____el
                #  /  \
                #  \   \
                #   \  /
                #    \/
                idx1 = sort_idx[[-1, -2]]
                idx2 = sort_idx[[-1, -3]]
            elif elevation < alts[sort_idx[-2]] and \
                 alts[sort_idx[-3]] < elevation:
                #   /\
                #  /  \
                # -\   \-----el
                #   \  /
                #    \/
                idx1 =  sort_idx[[-2, -4]]
                idx2 =  sort_idx[[-1, -3]]
            elif elevation < alts[sort_idx[1]] and \
                 alts[sort_idx[0]] < elevation:
                #   /\
                #  /  \
                #  \   \
                #  _\  /_____el
                #    \/
                idx1 = sort_idx[[0, 2]]
                idx2 = sort_idx[[0, 1]]
            else:
                print('failed to set az range')
                idx1 = sort_idx[[-1, -1]]
                idx2 = sort_idx[[0, 0]]
            tmp_az1 = np.interp(elevation, [alts[idx1[0]], alts[idx1[1]]],
                                [azs[idx1[0]], azs[idx1[1]]])
            # tmp_az1 = np.arccos(np.cos(tmp_az1))
            az[0, i] = tmp_az1
            tmp_az2 = np.interp(elevation, [alts[idx2[0]], alts[idx2[1]]],
                                [azs[idx2[0]], azs[idx2[1]]])
            # tmp_az2 = np.arccos(np.cos(tmp_az2))
            az[1, i] = tmp_az2
        else:
            ok[i] = False
            az[:, i] = np.array([np.nan, np.nan])
    return ok, az


def main(args):
    # set observation place
    if conf is not None and 'obs' in conf:
        lat = float(conf['obs']['lat'])
        lon = float(conf['obs']['lon'])
        if 'height' in conf['obs']:
            h = float(conf['obs']['height'])
        else:
            h = 0.0
    else:
        #  Chile, Atacama
        lat = -22.9579
        lon = -67.7862
        h = 0.0
    obs = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=h*u.m)

    # set time
    # +0000 = UTC timezone.
    # start time
    if args.start is None:
        start_time = datetime.now(tz=timezone.utc)
    elif re.match('[0-9][0-9][0-9][0-9]/[0-9][0-9]?/[0-9][0-9]?-[0-9][0-9]?:[0-9][0-9]?', args.start):
        start_time = datetime.strptime(args.start+' +0000', '%Y/%m/%d-%H:%M %z')
    elif re.match('[0-9][0-9][0-9][0-9]/[0-9][0-9]?/[0-9][0-9]?', args.start):
        start_time = datetime.strptime(args.start+' +0000', '%Y/%m/%d %z')
    else:
        print('invalid time format (start time)', file=sys.stderr)
        return
    # end time
    if args.end is None:
        end_time = start_time + timedelta(days=1)
    elif re.match('[0-9][0-9][0-9][0-9]/[0-9][0-9]?/[0-9][0-9]?-[0-9][0-9]?:[0-9][0-9]?', args.end):
        end_time = datetime.strptime(args.end+' +0000', '%Y/%m/%d-%H:%M %z')
    elif re.match('[0-9][0-9][0-9][0-9]/[0-9][0-9]?/[0-9][0-9]?', args.end):
        end_time = datetime.strptime(args.end+' +0000', '%Y/%m/%d %z')
    else:
        print('invalid time format (end time)', file=sys.stderr)
        return
    assert start_time < end_time
    # calculation interval
    if args.interval is None:
        interval = 60   # 60 minutes = 1 hour
    else:
        interval = args.interval
    # time array
    raw_times = np.arange(datetime.timestamp(start_time),
                          datetime.timestamp(end_time),
                          interval*60)
    times = Time(raw_times, format='unix')

    # save directory
    savedir = Path(args.output)
    if not savedir.is_dir():
        os.makedirs(savedir)
        print('make directory {}'.format(savedir))

    # condition setting
    if conf is not None and 'conditions' in conf:
        cond = conf['conditions']
    else:
        cond = {}
    # observable check
    el_min = chk_deg_condition(cond, 'el_min', 30.)
    el_max = chk_deg_condition(cond, 'el_max', 90.)
    az_min = chk_deg_condition(cond, 'az_min', 0.)
    az_max = chk_deg_condition(cond, 'az_max', 360.)
    sun_thd = chk_deg_condition(cond, 'sun_separation', 5.)
    moon_thd = chk_deg_condition(cond, 'moon_separation', 5.)
    con_el = chk_deg_condition(cond, 'constant_elevation', 41.)
    targets_observable = {}
    # raster scan range
    if 'raster_az_offset' in cond:
        d_az = cond['raster_az_offset']/deg
    else:
        d_az = 2./deg
    if 'raster_el_offset' in cond:
        d_el = cond['raster_el_offset']/deg
    else:
        d_el = 2./deg
    # print the conditions
    print('el_min: {:.2f} degree'.format(el_min*deg))
    print('el_max: {:.2f} degree'.format(el_max*deg))
    print('az_min: {:.2f} degree'.format(az_min*deg))
    print('az_max: {:.2f} degree'.format(az_max*deg))
    print('sun_separation: {:.2f} degree'.format(sun_thd*deg))
    print('moon_separation: {:.2f} degree'.format(moon_thd*deg))
    print('raster_az_offset: {:.2f} degree'.format(d_az*deg))
    print('raster_el_offset: {:.2f} degree'.format(d_el*deg))
    print('constant_elevation: {:.2f} degree'.format(con_el*deg))

    # calculation preparation
    targets = get_targets(args.targets, times)
    # variables storing results.
    targets_az = {}
    targets_el = {}
    targets_observable = {}

    # calculate
    sun = get_sun(times)
    moon = get_moon(times)
    for target in targets:
        if target.ces:
            ok, az = calc_CES(target, obs, times, con_el)
            targets_az[target.name] = az
            targets_el[target.name] = np.where(ok, con_el, np.nan)
            targets_observable[target.name] = ok
            continue
        target_altaz = target.body.transform_to(AltAz(obstime=times,
                                                      location=obs))
        sun_sep = np.array(target.body.separation(sun))/deg
        moon_sep = np.array(target.body.separation(moon))/deg
        targets_az[target.name] = np.array(target_altaz.az)/deg
        targets_el[target.name] = np.array(target_altaz.alt)/deg
        targets_observable[target.name] = \
            (el_min <= targets_el[target.name]) * \
            (targets_el[target.name] <= el_max) * \
            (az_min <= targets_az[target.name]) * \
            (targets_az[target.name] <= az_max) * \
            (sun_sep > sun_thd) * \
            (moon_sep > moon_thd)

    # outputs
    # make table
    table_label = ['time (UTC)']+[target.name for target in targets]
    dates = [datetime.fromtimestamp(t, tz=start_time.tzinfo) for t in raw_times]
    table_contents = [[date.strftime('%Y/%m/%d %H:%M')] for date in dates]
    for i, line in enumerate(table_contents):
        for target in targets:
            if target.ces:
                az = targets_az[target.name][:, i]*deg
                el = targets_el[target.name][i]*deg
                if not np.isnan(az[0]):
                    ces_az_min = int(az[0])
                    ces_az_max = int(az[1])
                    line.append('az: {:03d}-{:03d} deg | el: {:+05.1f} deg'.format(ces_az_min, ces_az_max, el))
                else:
                    line.append('az: ----- deg | el: {:+05.1f} deg'.format(el))
            else:
                az = targets_az[target.name][i]*deg
                el = targets_el[target.name][i]*deg
                line.append('az: {:05.1f} deg | el: {:+05.1f} deg'.format(az, el))
    table_str = tabulate(table_contents, table_label,
                         tablefmt='orgtbl', stralign='center')
    with open(savedir/'table.txt', 'w') as f:
        f.write(table_str)
    if args.show:
        print(table_str)
    # plot figures
    cmap = plt.get_cmap('tab10')
    fs = 15
    plt.rcParams['font.size'] = fs
    print_times = raw_times[[0, int(len(raw_times)/3), int(len(raw_times)*2/3), -1]]
    xlabel = [datetime.fromtimestamp(t, tz=start_time.tzinfo).strftime('%Y/%m/%d\n%H:%M') for t in print_times]

    # az plot
    fig1 = plt.figure(figsize=(16/1.5, 9/1.5))
    ax11 = fig1.add_subplot(111)
    for i, target in enumerate(targets):
        tname = target.name
        if target.ces:
            ces_az_max = targets_az[tname][0]
            ces_az_min = targets_az[tname][1]
            ax11.fill_between(raw_times, ces_az_min*deg, ces_az_max*deg,
                              color=cmap(i), label=tname)
        else:
            ax11.plot(raw_times, np.array(targets_az[tname])*deg, '-',
                      color=cmap(i), label=tname)
            if target.raster:
                ax11.fill_between(raw_times, (targets_az[tname]-d_az)*deg,
                                  (targets_az[tname]+d_az)*deg, alpha=0.3,
                                  color=cmap(i), label=None)
    ax11.hlines(az_min*deg, raw_times.min(), raw_times.max(),
                colors='gray', ls='--')
    ax11.hlines(az_max*deg, raw_times.min(), raw_times.max(),
                colors='gray', ls='--')
    ax11.set_xticks(print_times)
    ax11.set_xticklabels(xlabel)
    ax11.set_ylabel('azimuth [deg]')
    fig1.legend()
    fig1.savefig(savedir/'az_plot.pdf')

    # el plot
    fig2 = plt.figure(figsize=(16/1.5, 9/1.5))
    ax21 = fig2.add_subplot(111)
    for i, target in enumerate(targets):
        tname = target.name
        ax21.plot(raw_times, np.array(targets_el[tname])*deg, '-',
                  color=cmap(i), label=tname)
        if target.raster:
            ax21.fill_between(raw_times, (targets_el[tname]-d_el)*deg,
                              (targets_el[tname]+d_el)*deg, alpha=0.3,
                              color=cmap(i), label=None)
    ax21.hlines(el_min*deg, raw_times.min(), raw_times.max(),
                colors='gray', ls='--')
    ax21.hlines(el_max*deg, raw_times.min(), raw_times.max(),
                colors='gray', ls='--')
    ax21.set_xticks(print_times)
    ax21.set_xticklabels(xlabel)
    ax21.set_ylabel('elevation [deg]')
    fig2.legend()
    fig2.savefig(savedir/'el_plot.pdf')

    # observable flag plot
    fig3 = plt.figure(figsize=(16/1.5, 9/1.5))
    ax31 = fig3.add_subplot(111)
    for i, target in enumerate(targets):
        ax31.plot(raw_times, np.where(targets_observable[target.name],
                                      len(targets)-i, np.nan),
                  '-', lw=4)
    ax31.set_xticks(print_times)
    ax31.set_xticklabels(xlabel)
    ax31.set_xlim([raw_times.min(), raw_times.max()])
    ax31.set_yticks([len(targets)-i for i in range(len(targets))])
    ax31.set_yticklabels(targets_observable.keys())
    ax31.set_ylim([0.5, len(targets)+0.5])
    fig3.savefig(savedir/'obs_flags.pdf')

    # pointing of stars
    fig4 = plt.figure(figsize=(16/1.5, 9/1.5))
    ax41 = fig4.add_subplot(111, projection='mollweide')
    for i, target in enumerate(targets):
        if target.ces:
            continue
            ces_label = target.name
            for j in range(len(times)):
                if np.isnan(targets_el[target.name][j]):
                    continue
                print(targets_az[target.name][:, j])
                ax41.fill_between(targets_az[target.name][:, j],
                                  [targets_el[target.name][j],
                                   targets_el[target.name][j]],
                                  [targets_el[target.name][j],
                                   targets_el[target.name][j]],
                                  alpha=0.3, color=cmap(i), label=ces_label)
                ces_label = None
        else:
            lat = np.array(targets_az[target.name])
            lat = np.where(lat > np.pi, lat-2*np.pi, lat)
            lon = np.array(targets_el[target.name])
            ax41.plot(lat[1:-1], lon[1:-1], '.',
                      color=cmap(i), label=target.name)
            # ax41.plot(lat, lon, '--', color=cmap(i), label=target.name)
            ax41.plot([lat[0]], [lon[0]], '*', color=cmap(i), ms=6)
            ax41.plot([lat[-1]], [lon[-1]], 'x', color=cmap(i), ms=6)
            if target.raster:
                for j in range(len(lat)):
                    ax41.fill_between([lat[j]-d_az, lat[j]+d_az],
                                      [lon[j]+d_el, lon[j]+d_el],
                                      [lon[j]-d_el, lon[j]-d_el],
                                      alpha=0.3, color=cmap(i))
    fig4.legend()
    ax_pos = ax41.get_position()
    fig4.text(ax_pos.x1, ax_pos.y0+0.1, 'start: *\nend: x')
    fig4.savefig(savedir/'mollweide.pdf')

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
