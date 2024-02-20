Introduction
============

This is a project to explore astro-geological calculation and visualisation for
Muslim prayer times. It's not meant to be authoritative (either scientifically
or religiously); it's only for fun and interest. It's written from a
perspective having only an introductory familiarity with the complex and
nuanced Muslim praxis and with astrodynamics.


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

- The terrestrial poles - also called geodetic, geographic or
  [true poles](https://en.wikipedia.org/wiki/True_north) - are intersection
  points of Earth's rotational axis through its surface. We colloquially call
  these the north and south pole. They slowly
  [precess](https://en.wikipedia.org/wiki/Axial_precession) over tens of
  thousands of years.
- The [magnetic poles](https://en.wikipedia.org/wiki/Poles_of_astronomical_bodies#Magnetic_poles)
  orient a compass. When Muslims determine the qibla for a given location, they
  need to account for
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
This curve is an [isochrone](https://en.wiktionary.org/wiki/isochrone). If the
celestial and ecliptic poles were aligned, then prayer isochrones would be
[meridians](https://en.wikipedia.org/wiki/Meridian_(geography)) in the
rotational frame. However, since the Earth's rotation has an
[axial tilt](https://en.wikipedia.org/wiki/Axial_tilt) or _ecliptic obliquity_,
the poles are offset by about 23° and the isochrones are actually _ecliptic_
meridians extending to the ecliptic poles, not the terrestrial poles. This has
the effect of changing prayer times through the year.

Further, since all of the prayer definitions need various amounts of angular
correction to account for refraction, they are not really meridians. They're
geodesics offset (sometimes significantly, up to 20-some degrees) from the
ecliptic poles.

Degenerate locations
--------------------

The ecliptic poles for the current day of year are
[gimbal-locked](https://en.wikipedia.org/wiki/Gimbal_lock), which implies an
undefined prayer schedule at that location: it is both always time for zuhr and
never time for zuhr. At this location, the sun circles the horizon but never
rises or sets.

There are two degenerate polar regions; whether they appear in the north or
south depends on the time of year. In one region spanning from an ecliptic pole
to the terminus of the fajr/isha geodesics, it is always dark but never
dark enough to meet the refraction threshold for fajr and isha, so fajr and
isha cannot happen. In the other region bounded by the other ecliptic pole, it
is always bright, the sun never sets, and fajr and isha still cannot happen.

[Fatwa 2769](https://islamqa.info/en/answers/5842/how-to-pray-and-fast-in-countries-where-the-day-or-night-is-continuous)
attempted to clear these ambiguities.


Implementation
==============

This code follows (some of the) schedule definitions from
[Hamid's prayer time guide](http://www.praytimes.org/calculation/).


Setup
=====

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

then run `demo.py`.
