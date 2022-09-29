---
layout: post
title:  "Calculating Optimal Bounding Volumes, Pt 1: Capsules"
date:   2022-09-28 22:31:25 -0700
---

Creating some sort of bounding volume for a set of vertices is a common problem in simulations and games. Capsules, or swept spheres, are a powerful option to use as a bounding volume because of how little data is required to define them and how useful they are for collision detection. Imagine taking a sphere and sweeping it along a single line in space; the shape of this would be a capsule. We can define a capsule then as:

```
Vector3 start;
Vector3 end;
float radius;
```

These techniques can be used on any sort of point cloud, but let’s state the specifics of our case. We have a set of vertices that define some bone (right forearm, left leg, hip, etc.) and we want to find a capsule that best fits them. We’ll show how to do this generically and then how to adapt it to use Unity components.