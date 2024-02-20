from datetime import datetime, timedelta, timezone

import numpy as np
from cartopy.feature import nightshade

from demo import load_home


def hamid_equivalents() -> None:
    """
    http://www.praytimes.org/calculation/
    """
    local_now = datetime.now().astimezone()
    Lng, Lat = load_home()

    today = local_now.replace(hour=0, minute=0, second=0, microsecond=0)
    TimeZone = local_now.tzinfo.utcoffset(local_now) / timedelta(hours=1)
    jd = nightshade._julian_day(today)

    d = jd - 2451545  # jd is the given Julian date
    g = np.deg2rad(np.fmod(357.529 + 0.98560028 * d, 360))
    q = np.deg2rad(np.fmod(280.459 + 0.98564736 * d, 360))
    L = q + np.deg2rad(np.fmod(1.915 * np.sin(g) + 0.020 * np.sin(2 * g), 360))
    R = 1.00014 - 0.01671 * np.cos(g) - 0.00014 * np.cos(2 * g)
    e = np.deg2rad(np.fmod(23.439 - 0.00000036 * d, 360))
    RA = np.arctan2(
        np.cos(e) * np.sin(L),
        np.cos(L),
    ) / np.deg2rad(15)
    D = np.arcsin(np.sin(e) * np.sin(L))  # declination of the Sun
    EqT = q / np.deg2rad(15) - RA  # equation of time
    EqT = sorted(
        [
            (abs(offset), offset)
            for offset in (EqT, EqT - 24)
        ]
    )[0][1]
    Dhuhr = 12 + TimeZone - Lng / 15 - EqT

    alpha = np.deg2rad(0.83)
    T = np.arccos(
        (
            -np.sin(alpha) - np.sin(L) * np.sin(D)
        ) / np.cos(L) / np.cos(D)
    ) / np.deg2rad(15)
    sunrise = Dhuhr - T
    sunset = Dhuhr + T

    local_now = today + timedelta(hours=sunrise)
    utcnow = local_now.astimezone(timezone.utc)
