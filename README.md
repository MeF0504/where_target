# where target <!-- omit in toc -->

This is a script that calculates the location of stars, planets, and so on over the sky.
This script makes following table and plots in the output directory.
- `table.txt`
  Show the elevation and azimuth of each star and each specified time.
- `az_plot.pdf`
  Show the azimuth of each specified time.
- `el_plot.pdf`
  Show the elevation of each specified time.
- `obs_flags.pdf`
  Show each star is observable at each specified time. Here, I consider the elevation and azimuth range, and whether far away from the Sun and Moon.
- `mollweide.pdf`
  Show the pointing of each star in the sky using Mollweide method.

## Table of Contents <!-- omit in toc -->

- [Requirement](#requirement)
- [Configuration](#configuration)
  - [obs](#obs)
  - [targets](#targets)
  - [conditions](#conditions)
- [Usage](#usage)

## Requirement
- python3
- modules
    - `ephem`
    - `numpy`
    - `matplotlib`
    - `tabulate`

    You can install them from pip.
    ``` shell
    pip install ephem numpy matplotlib tabulate
    ```

## Configuration
You can set some configurations by making `where_target_config.json` at the same path as where_target.py.

Available contents are following.
Please also see the `sample_where_target_config.json`.

### obs
  Observation site configuration.
  Available options are:
  - `lat` (string)
    - latitude of the observation site.
  - `lon` (string)
    - longitude of the observation site.

### targets
  Configuration list of additional stars. Sun, Moon, and planets are already defined. These stars are considered fixed bodies.

  Available options for each item are:
  - `name` (string)
    - Name of the star. It is shown in the made plots.
  - `type` (string)
    - Type of stars. It is used in the argument to specify the displayed stars in the plot.
  - `ra` (string)
    - Right Ascension of the star (hour:min:sec).
  - `dec` (string)
    - Declination of the star (deg:arcmin:arcsec).
  - `raster` (bool)
    - The flag that the scan of this star is a raster scan or not. If true, the raster scan region is also shown in plots.

### conditions
  Observation configurations.
  Available options are:
  - `el_min` (float)
    - The minimum value of elevation (degree). It is shown in elevation plot and used in the observable flag.
  - `el_max` (float)
    - The maximum value of elevation (degree). It is shown in elevation plot and used in the observable flag.
  - `az_min` (float)
    - The minimum value of azimuth (degree). It is shown in azimuth plot and used in the observable flag.
  - `az_max` (float)
    - The maximum value of azimuth (degree). It is shown in azimuth plot and used in the observable flag.
  - `sun_separation` (float)
    - The minimum observable angle from the Sun (degree). Is is used in the observable flag.
  - `moon_separation` (float)
    - The minimum observable angle from the Moon (degree). Is is used in the observable flag.
  - `raster_el_offset` (float)
    - The elevation range of the raster scan (degree).
  - `raster_az_offset` (float)
    - The azimuth range of the raster scan (degree).

## Usage
``` shell
python3 where_target.py [options]
```
currently available options are:
- -targets
  - select targets. sun, moon, and planets are available by default.
  You can append targets by adding configurations in `where_target_config.json`.
- -start
  - start time (UTC). The format is `yyyy/MM/dd` or `yyyy/MM/dd-hh:mm`.
- -end
  - end time (UTC). The format is `yyyy/MM/dd` or `yyyy/MM/dd-hh:mm`.
- -interval
  - time interval (minutes).
- -show
  - show images after this script is executed.
- -o/--output
  - specify the output directory.
