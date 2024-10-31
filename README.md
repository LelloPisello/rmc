# RMC
## What is it?
It is a fixed step raymarcher that marches the rays through a curved and bounded manifold.
rays are then checked for intersections against volumes.

### What's wrong with it?
I chose to make it a fixed step raymarcher because conciling SDFs with curved space seemed like a quick path to hell.
This is so far the only conscious design flaw i have made.
This is basically a speed coding project, so the code is riddled with issues both in its architecture and its legibility.

### More details
The manifold (which acts as the 'scene') is modelled via a metric tensor (and an inverse metric) which is saved onto a buffer and linearly interpolated.
Christoffel symbols are computed and again saved onto their own buffer. 
Rays are marched using the geodesic equation:

$\ddot{x}^i=-\Gamma^i_{jk}\dot{x}^j\dot{x}^k$

The thorough implementation of covectors/vectors is due to the possibility of the implementation of vector and covector fields later down the line.

