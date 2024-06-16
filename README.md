# fastDTW in Odin - Fast Dynamic Time Warping
A fast port to Odin of the fastDTW algorithm, that was original in Python.

## Description
This is a port from Python to Odin of the fastDTW library. The original library is the implementation of the fastDTW algorithm that is described in the paper "FastDTW: Toward Accurate Dynamic Time Warping in Linear Time and Space" by Stan Salvador and Philip Chan. This is a algorithm to determine the distance between 2 time series, and the point of best alignment between them in O( N ). This Odin implementation is approximately 3x and some change faster than the Python implementation. A tested it with the distance calculation and path alignment between 2 series of 2_500_000 elements each with a radius of 20. It took 15 minutes and 91 GB of memory. Normally with the normal DTW algorithm, one can only calculate the DTW distance between of 2 series of up to 3_000 elements, because the algorithm complexity is O( N * M ) where N is the len( seq_1 ), and M is the len( seq_2 ). This algorithm, fastDTW is O( N ) and the memory complexity is also O( N ).

## Original github repository in Python of this port

Github - slaypni - fastdtw <br>
Kazuaki Tanida <br>
[https://github.com/slaypni/fastdtw](https://github.com/slaypni/fastdtw) <br>
<br>
Python implementation of FastDTW, which is an approximate Dynamic Time Warping (DTW) algorithm that provides optimal or near-optimal alignments with an O(N) time and memory complexity.

## The original paper from 2007
Stan Salvador, and Philip Chan. "FastDTW: Toward accurate dynamic time warping in linear time and space." Intelligent Data Analysis 11.5 (2007): 561-580.
[https://cs.fit.edu/~pkc/papers/tdm04.pdf](https://cs.fit.edu/~pkc/papers/tdm04.pdf)

## To compile and run do

```
$ make opti

$ make run
```

## License

MIT Open Source License

## Have fun!
Best regards, <br>
João Carvalho, <br>


