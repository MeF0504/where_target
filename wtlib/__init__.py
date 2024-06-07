import json
from pathlib import Path

import astropy.units as u
from astropy.coordinates import get_sun, get_body, \
    SkyCoord, solar_system_ephemeris, EarthLocation


add_targets = {}
conf = None
conf_file = Path(__file__).parent.parent/'where_target_config.json'
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


def get_obs():
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
    return obs


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
        tmp = Target('Moon', 'moon', get_bosy('moon', time), False, False)
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
