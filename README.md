# Cyvincenty

A fast Cython implementation of the Vincenty algorithm for calculating the distance in kilometers between 2 co-ordinates.

This module is heavily inspired by [uvincenty](https://github.com/vivescere/uvincenty) - a pure C Python extension - and is just as fast, if not slightly faster despite being written in Python! (technically Cython)

## Installation

```
pip install cyvincenty
```

## Usage

```
>> from cyvincenty import vincenty

>> boston = (42.3541165, -71.0693514)
>> newyork = (40.7791472, -73.9680804)

>> cyvincenty(*boston, *newyork)
```


## Benchmarks
```
from cyvincenty import vincenty
import uvincenty
from vincenty import vincenty
import geopy.distance

boston = (42.3541165, -71.0693514)
newyork = (40.7791472, -73.9680804)

%timeit uvincenty.vincenty(*boston, *newyork)

%timeit cyvincenty(*boston, *newyork)

%timeit vincenty(boston, newyork)

%timeit geopy.distance.geodesic(boston, newyork)

```