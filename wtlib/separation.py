import sys
from pathlib import Path

import numpy as np
from astropy.time import Time
from astropy.coordinates import get_sun, get_moon, AltAz

sys.path.append(str(Path(__file__).parent.parent))
from wtlib import conf, get_targets, get_obs

deg = 180.0/np.pi


def separation_x(unix_time_array, target_name, sun_moon):
    time = Time(unix_time_array, format='unix')
    targets = get_targets(target_name, time)
    obs = get_obs()

    ret = []
    for target in targets:
        if target.ces:
            print("CES is not supported")
            continue
        result = {}
        target_altaz = target.body.transform_to(AltAz(obstime=time,
                                                      location=obs))
        sep = np.array(target.body.separation(sun_moon))/deg
        result['name'] = target.name
        result['sep'] = sep
        result['az'] = np.array(target_altaz.az)/deg
        result['el'] = np.array(target_altaz.alt)/deg

        ret.append(result)

    return ret


def separation_sun(unix_time_array, target_name):
    time = Time(unix_time_array, format='unix')
    return separation_x(unix_time_array, target_name, get_sun(time))


def separation_moon(unix_time_array, target_name):
    time = Time(unix_time_array, format='unix')
    return separation_x(unix_time_array, target_name, get_moon(time))


if __name__ == '__main__':
    from datetime import datetime, timedelta

    import matplotlib.pyplot as plt
    tdy = datetime.today()
    time = np.linspace((tdy+timedelta(weeks=-2)).timestamp(),
                       (tdy+timedelta(weeks=2)).timestamp(),
                       30000)
    print(time[-1]-time[0])
    ret = separation_sun(time, ['calibration'])
    if not ret:
        print('no targets')
        exit()
    fig1 = plt.figure()
    for i, target in enumerate(ret):
        axes = fig1.add_subplot(len(ret), 1, i+1)
        # y = np.where(30/deg <= target['el'], target['sep'], np.nan)
        y = target['sep']
        axes.plot(time, y, '.', label=target['name'])
    plt.show()
