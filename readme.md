Introduction
============

This is a project to explore astro-geological calculation and visualisation for
Muslim prayer times. It's not meant to be authoritative (either scientifically
or religiously); it's only for fun and interest. It's written from a
perspective having only an introductory familiarity with the complex and
nuanced Muslim praxis and with astrodynamics.


Setup
=====

Requirements
------------

On first run the program needs access to the internet to download images from NASA.

The program needs modern Python to run, having been tested on 3.10+.

All direct dependencies are listed in `requirements.txt`.

On operating systems (typically Windows) that don't have a good built-in
timezone database, you will also 
[need to install `tzdata`](https://docs.python.org/3/library/zoneinfo.html#data-sources).

Dependencies
------------

Install third-party packages via

`pip install -r requirements.txt`


Usage
=====

Populate `.home.json` with the location of prayer, to look like

```json
{
  "lat": 40,
  "lon": -80
}
```

then run `python -m isochrones`. On first run it will download Earth surface
images.


Interpretation
==============

This describes how to interpret the program UI.

![screenshot](https://github.com/reinderien/isochrones/assets/1236420/43cde1cc-3036-4f8e-a3fe-2337dbb93588)

The plot is divided into two hemispheres. These are probably not exactly
opposite to each other and so some overlap will be apparent. The left
hemisphere is centred on "home", and the right hemisphere is centred on the
Kaaba. These points are both shown with small crosses.

From the home point, a heading from north is drawn in red to indicate the
departing angle of the _qibla_ geodesic. The geodesic is drawn in green in both
hemispheres, and the distance of this geodesic to the Kaaba is shown at the
edge of the home hemisphere.

A timestamp is shown at both hemisphere centres. At the home centre the
timezone is assumed to be the local timezone of the user's computer, and at the
Kaaba the timezone is that of Saudi Arabia. If in static plot mode, these are
the times of program start. If in animation mode, these are "now".

The globe is drawn with three shading states: day (no shading), dusk (partial
shading), and night (full shading). Night is still lighter than black to keep
the map visible.


Theory and terms
================

[Salah](https://en.wikipedia.org/wiki/Salah) (صَلَاة)
is the practice of daily Muslim prayer.

Geodesics
---------

The [Kaaba](https://en.wikipedia.org/wiki/Kaaba) (ٱلْكَعْبَة)
in [Mecca](https://en.wikipedia.org/wiki/Mecca) (مكة),
Saudi Arabia is the geographical centre of the Islam religion. Muslims face
toward it when praying. Since the Earth is an
[oblate spheroid](https://en.wikipedia.org/wiki/Equatorial_bulge), the shortest
surface path to the Kaaba is via
[geodesic](https://en.wikipedia.org/wiki/Geodesics_on_an_ellipsoid). This
geodesic's initial [heading](https://en.wikipedia.org/wiki/Heading_(navigation))
is called the [qibla](https://en.wikipedia.org/wiki/Qibla) (قِبْلَة).
This is the same path and initial heading that an airplane takes when flying to
Saudi Arabia.

The user's screen is flat, so all maps need to choose a projection to
draw the spheroidal [globe](https://en.wikipedia.org/wiki/Globe).
The [orthographic projection](https://en.wikipedia.org/wiki/Orthographic_projection)
is convenient because it is easy to understand, appears to the viewer as if
they are viewing one
[hemisphere](https://en.wikipedia.org/wiki/Hemispheres_of_Earth)
of a 3D globe, and draws all geodesics through the centre as straight lines
that  otherwise appear curved on
[non-azimuthal](https://en.wikipedia.org/wiki/Map_projection#Azimuthal_.28projections_onto_a_plane.29)
projections. Thus, this program draws two hemispheres, one centred on the
user's home location and one on the Kaaba.

Reference frames and time
-------------------------

The typical five daily salah prayers are

- fajr (صلاة الفجر, at dawn),
- zuhr (صَلَاة ٱلظُّهْر, at noon),
- asr (صلاة العصر, at afternoon),
- maghrib (صلاة المغرب, at sunset), and
- isha (صلاة العشاء, at dusk).

Each is tied to a specific [solar time](https://en.wikipedia.org/wiki/Solar_time).
Calculations relating solar time and
[coordinated time](https://en.wikipedia.org/wiki/Coordinated_Universal_Time)
require some background in terrestrial frames of reference. The Earth has at
least four different kinds of pole:

- The terrestrial poles - also called geodetic,
  [geographic](https://en.wikipedia.org/wiki/Geographic_coordinate_system) or
  [true poles](https://en.wikipedia.org/wiki/True_north) - are intersection
  points of Earth's rotational axis through its surface. We colloquially call
  these the north and south pole. They slowly
  [precess](https://en.wikipedia.org/wiki/Axial_precession) over tens of
  thousands of years.
- The [magnetic poles](https://en.wikipedia.org/wiki/Poles_of_astronomical_bodies#Magnetic_poles)
  orient a compass. When Muslims determine the _qibla_ for a given location,
  they need to account for
  [magnetic declination](https://en.wikipedia.org/wiki/Magnetic_declination) to
  translate from the magnetic frame to the terrestrial frame, in turn
  determining a heading from true north. These magnetic poles
  [wander slowly](https://en.wikipedia.org/wiki/Paleomagnetism)
  over millions of years.
- The [celestial poles](https://en.wiktionary.org/wiki/celestial_pole) sit at
  arbitrary heights above the terrestrial poles and are independent of the
  Earth's surface, instead intersecting a virtual
  [celestial sphere](https://en.wikipedia.org/wiki/Celestial_sphere).
- The [ecliptic (or orbital) poles](https://en.wikipedia.org/wiki/Orbital_pole#Ecliptic_pole)
  sit on the [normal](https://en.wikipedia.org/wiki/Normal_(geometry)) to the
  solar system [orbital plane](https://en.wikipedia.org/wiki/Ecliptic). The
  boundary of day and night will always intersect the ecliptic poles, before
  correcting for [diffraction](https://en.wikipedia.org/wiki/Dawn).
  The ecliptic moves through an annual cycle.

Isochrones
----------

At any given time, for a specific prayer, there is one curve across the Earth
where all points on that curve satisfy the solar conditions for that prayer.
This curve is an [isochrone](https://en.wiktionary.org/wiki/isochrone). We can
understand these isochrones using a series of simplifications that we
progressively have to give up.

If the celestial and ecliptic poles were aligned, then most prayer isochrones
would be [meridians](https://en.wikipedia.org/wiki/Meridian_(geography)) in the
rotational frame. However, since the Earth's rotation has an
[axial tilt](https://en.wikipedia.org/wiki/Axial_tilt) or _ecliptic obliquity_,
the poles are offset by about 23° and the isochrones are actually _ecliptic_
meridians extending to the ecliptic poles, not the terrestrial poles. This has
the effect of changing prayer times through the year.

Only the noon _zuhr_ is actually an ecliptic meridian. Since the _fajh_,
_maghrib_ and _isha_ prayer definitions need various amounts of angular
correction to account for refraction, they are not really meridians. They're
ellipsoidal sections offset (sometimes significantly, up to 20-some degrees)
from the ecliptic poles.

The afternoon _asr_ is much more complicated, and not an ellipsoidal section at
all. Since it uses the absolute value of the ecliptic latitude, it's not even
[differentiable](https://en.wikipedia.org/wiki/Differentiable_function).

Degenerate locations
--------------------

The ecliptic poles for the current day of year are
[gimbal-locked](https://en.wikipedia.org/wiki/Gimbal_lock), which implies an
undefined prayer schedule at that location: it is both always time for _zuhr_
and never time for _zuhr_. At this location, the sun circles the horizon but
never rises or sets.

There are two degenerate polar regions; whether they appear in the north or
south depends on the time of year. In one region spanning from an ecliptic pole
to the terminus of the _fajr_/_isha_ isochrones, it is always dark but never
dark enough to meet the refraction threshold for _fajr_ and _isha_, so _fajr_
and _isha_ cannot happen. In the other region bounded by the other ecliptic
pole, it is always bright, the sun never sets, and _fajr_ and _isha_ still
cannot happen.

[Fatwa 2769](https://islamqa.info/en/answers/5842/how-to-pray-and-fast-in-countries-where-the-day-or-night-is-continuous)
attempted to clear these ambiguities.


Implementation
==============

This code follows (some of the) schedule definitions from
[Hamid Zarrabi-Zadeh's prayer time guide](http://www.praytimes.org/calculation/).
However, its definition for _asr_ is incorrect, so the astronomy code also
follows
[Radhi Fadlillah's guide](https://radhifadlillah.com/blog/2020-09-06-calculating-prayer-times/).

The implementation relies heavily on [Cartopy](https://scitools.org.uk/cartopy/docs/latest)
for geographic visualisation, in turn relying heavily on
[Pyproj](https://pyproj4.github.io/pyproj/stable/) for coordinate system
transformations.

The image of the globe is from NASA's
[Blue Marble Next Generation](https://visibleearth.nasa.gov/images/73751/july-blue-marble-next-generation-w-topography-and-bathymetry/73753l)
dataset. This shows topography and bathymetry in "true colour" without
atmospheric occlusion.

Accuracy
--------

Prayer times that depend on sunlight level technically depend on surface
topology. For example, if you live in a valley and the sun needs to pass a
mountain for sunset to occur, sunset will be much earlier; but this program
does not take that into account. It assumes that the Earth is smooth.

The Earth is an oblate spheroid, but Cartopy's 
[`Orthographic` projection](https://scitools.org.uk/cartopy/docs/latest/reference/projections.html#orthographic)
doesn't support ellipsoidal coordinate systems; it can only plot with spherical
coordinate systems. For this reason, the program plots are not to scale.
However, time, angle and distance calculations still use the more accurate
[WGS84 ellipsoid globe](https://scitools.org.uk/cartopy/docs/latest/reference/generated/cartopy.crs.Globe.html#cartopy.crs.Globe).

Cartopy uses a popular generalised formula for solar times relative to
sunrise/sunset which follows the 
[sunrise equation on Wikipedia](https://en.wikipedia.org/wiki/Sunrise_equation#Generalized_equation).
This has accuracy limits not discussed here.

Cartopy calculates the [Julian date](https://en.wikipedia.org/wiki/Julian_day)
using an algorithm written in
[nightshade.py](https://github.com/SciTools/cartopy/blob/main/lib/cartopy/feature/nightshade.py#L101)
adapted from
[Vallado's Fundamentals of Astrodynamics and Applications, Algorithm 14](https://archive.org/details/FundamentalsOfAstrodynamicsAndApplications/page/n209/mode/2up).
Vallado describes this as being accurate to about 0.4 ms through the year 2100.
In the same file and from the same textbook, Cartopy adapts
[algorithm 29](https://archive.org/details/FundamentalsOfAstrodynamicsAndApplications/page/n305/mode/1up)
for the sun position vector. Vallado warns that these are not as precise as
the JPL's Sun-Earth ephemerides, but are still accurate to 0.01 degrees through
the year 2050.

The program uses a mid-resolution (5400x2700) Blue Marble image. This has about
0.2 degrees of resolution on the equator, chosen due to a modest download size
of 2.3 MB. You won't be able to resolve your swimming pool.

The program is somewhat slow, so there will be lag typically up to one second
between the depicted time and the true time. The true time depends on an
accurate operating system clock. The program does not look up the correct time
zone for 'home'; instead it assumes that the operating system's local time zone
corresponds to the correct time zone for 'home'.


Test Cases
==========

In all cases, use the Islamic Society of North America (15/15) for prayer method,
and Shafi for Asar Method.

[St. John's, Canada](https://www.salahtimes.com/canada/st-johns) during 
fajr, dhuhr, asr, maghrib, and isha:

    python -m isochrones -s -t 2024-09-08T05:02-02:30 -x -52.703471 -y 47.567119
    python -m isochrones -s -t 2024-09-08T12:59-02:30 -x -52.703471 -y 47.567119
    python -m isochrones -s -t 2024-09-08T16:34-02:30 -x -52.703471 -y 47.567119
    python -m isochrones -s -t 2024-09-08T19:27-02:30 -x -52.703471 -y 47.567119
    python -m isochrones -s -t 2024-09-08T20:55-02:30 -x -52.703471 -y 47.567119

[Trondheim, Norway during asr](https://www.salahtimes.com/norway/trondheim):

    python -m isochrones -s -t 2024-09-08T10:44-04:00 -x 10.39332 -y 63.42932

[Bloemfontein, South Africa during asr](https://www.salahtimes.com/south-africa/bloemfontein):

    python -m isochrones -s -t 2024-09-08T09:34-04:00 -x 26.21904 -y -29.12502

