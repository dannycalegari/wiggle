# wiggle

A program to explore and to certify properties of the space of *wiggles*.

For a complex parameter `z` with `|z| < 1` and `|1-z| < 1`, let

```
f(x) = -zx + z,      g(x) = (z-1)x + 1
```

The attractor of the iterated function system `{f,g}` is a curve `gamma(z)`
from `0` to `1` — a *squiggle*. It is a *wiggle* when it is embedded, and

```
W = { z : gamma(z) is embedded },     R = { |z-1/2| < 1/2 } - W
```

The program draws these sets and produces machine certificates that a given
parameter lies in `W`, or in `R`, together with a rigorous radius about that
parameter on which the certificate remains valid. It is the program behind the
computer-assisted proof in the paper *Wiggle Island*, which shows that `W` is
disconnected.

## Building

Requires a C++11 compiler and X11 (for the interactive modes only).

```
make
```

On macOS the X11 headers are expected in `/opt/X11`; edit `IFLAGS` and `LFLAGS`
in the makefile if yours are elsewhere.

## Running

With no arguments the program is interactive and opens an X11 window:

```
./wiggle
```

It then asks whether to plot roots or to draw a map of the parameter plane. In
the map, click to zoom, `z`/`x` to refine or coarsen, `a`/`s` to change the
search depth, `q` to quit.

All the other modes are batch, need no display, and write their output to a
file:

| command | what it does |
| --- | --- |
| `wiggle -fig maxlen rmax N out.ppm` | roots of the differences `p_u - p_v` for all words of length at most `maxlen` lying in `\|z-1/2\| < rmax` |
| `wiggle -cert x y [depth] [cap]` | search for a stable crossing at `x+iy` and report its certified radius |
| `wiggle -w x y [depth] [cap]` | the `W` test on its own: does the list of viable pairs empty? |
| `wiggle -env x y rho k` | verify `gamma(z)` is inside the envelope for all `\|z-(x+iy)\| <= rho`, using words of length `k` |
| `wiggle -locate x y hw N [depth] [cap] [out.ppm]` | classify an `N x N` grid: `W`, `R`, or undecided |
| `wiggle -map x y hw N [depth] [cap] [threads] [out.ppm]` | as above, shaded by the length of the shortest certifying pair; multithreaded; also writes a `.dat` of the raw values |
| `wiggle -scan x y hw N [depth] [out.ppm]` | as above, shaded by the certified radius |
| `wiggle -loop x y rad [depth] [coarse] [out.txt]` | cover the circle of radius `rad` about `x+iy` by certified disks, and report whether it closes up |

Images are written as binary PPM; convert with e.g. `convert out.ppm out.png`.

To reproduce the main computation of the paper:

```
./wiggle -env  0.3409 0.43486 0.0001 18
./wiggle -loop 0.34090875 0.43484625 0.000042 300 600 loop.txt
./wiggle -w    0.34090875 0.43484625 3000 3000000
```

The first verifies the hypothesis that `gamma(z)` lies inside the envelope; the
second covers a circle around the island by 1046 certified disks with no gaps;
the third certifies that the centre of that circle lies in `W`. Together these
prove that `W` is disconnected. The file `certified_loop.txt` in this
repository is the output of the second command: one line per disk, giving the
angle, the centre, the certified radius, and the two words `u`, `v`.

## How a crossing is certified

For words `u` beginning with `f` and `v` beginning with `g`, put `A = u^{-1}v`.
Then `u gamma` meets `v gamma` exactly when `gamma` meets `A gamma`. A *stable
crossing* is certified by exhibiting four straight rays, from the four points

```
0, 1  on gamma        A(0), A(1)  on A gamma
```

whose arguments occur in alternating cyclic order, such that the rays from `0`
and `1` miss the envelope of `A gamma`, the rays from `A(0)` and `A(1)` miss
the envelope of `gamma`, and each ray of one pair is disjoint from each ray of
the other. The two resulting proper paths can then meet only along
`gamma ^ A gamma`, and a planar topology argument forces that intersection to
be non-empty.

Fixing the four directions, the certificate persists as long as eight
clearances stay positive, and each is 1-Lipschitz in the objects that move; the
derivatives are bounded above on a disk by carrying `d/dz` through the
composition in disk arithmetic. This yields an explicit radius about `z` on
which the crossing is certified.

Two caveats, stated also in the paper: the disk arithmetic absorbs floating
point error by a fixed relative inflation rather than by directed rounding, and
the `W` test is carried out in ordinary floating point. A fully rigorous
implementation would use interval arithmetic throughout.

## Files

| file | contents |
| --- | --- |
| `wiggle.cc` | main; command line dispatch |
| `graphics.cc` | X11 window handling |
| `polynomial.cc` | polynomials, evaluation, root finding |
| `roots.cc` | roots of the differences `p_u - p_v` |
| `figures.cc` | non-interactive figure generation |
| `prune.cc` | the pair enumeration, and the ball disjointness test |
| `certify.cc` | disk arithmetic, the crossing test, certified radii, maps, the loop |
| `interface.cc` | the interactive map viewer |
| `certified_loop.txt` | the certificates covering the loop around the island |

## Authorship and licence

The program was written mainly in 2022 by Danny Calegari. The additions made
for the 2026 revision of the paper — `figures.cc`, `certify.cc`, the
certified-radius and loop machinery, and the corrected crossing test — were
written by Claude (Anthropic) under the supervision of Danny Calegari.

Released under the terms of the GNU GPL.
