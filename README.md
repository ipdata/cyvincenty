# Cyvincenty

A fast Cython implementation of the Vincenty algorithm for calculating the distance in kilometers between 2 co-ordinates.

This module is heavily inspired by [uvincenty](https://github.com/vivescere/uvincenty) and [vincenty](https://github.com/maurycyp/vincenty). It is just as fast if not slightly faster than uvincenty which is a pure C Python extension despite being written in Python! (technically Cython).

## Installation

```bash
pip install cyvincenty
```

## Usage

```python
>> from cyvincenty import vincenty

>> boston = (42.3541165, -71.0693514)
>> newyork = (40.7791472, -73.9680804)

>> cyvincenty(*boston, *newyork)
```


## Benchmarks

Using ipython

```python
>> import cyvincenty
>> import uvincenty
>> import vincenty
>> import geopy.distance

>> boston = (42.3541165, -71.0693514)
>> newyork = (40.7791472, -73.9680804)

>> %timeit uvincenty.vincenty(*boston, *newyork)
744 ns ± 2.58 ns per loop (mean ± std. dev. of 7 runs, 1,000,000 loops each)

>> %timeit cyvincenty.vincenty(*boston, *newyork)
736 ns ± 2.82 ns per loop (mean ± std. dev. of 7 runs, 1,000,000 loops each)

>> %timeit vincenty.vincenty(boston, newyork)
10.2 µs ± 60.1 ns per loop (mean ± std. dev. of 7 runs, 100,000 loops each)

>> %timeit geopy.distance.geodesic(boston, newyork)
191 µs ± 1.52 µs per loop (mean ± std. dev. of 7 runs, 10,000 loops each)
```