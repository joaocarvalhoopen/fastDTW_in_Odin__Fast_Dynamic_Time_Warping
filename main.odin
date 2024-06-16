// Project name       : fastDTW in Odin - Fast Dynamic Time Warping
// Date               : 2024.06.16
// Author of the port : João Carvalho
// License            : MIT Open Source License
//
// Description        : This is a port from Python to Odin of the fastDTW library.
//                      The original library is the implementation of the fastDTW
//                      algorithm that is described in the paper "FastDTW: Toward
//                      Accurate Dynamic Time Warping in Linear Time and Space"
//                      by Stan Salvador and Philip Chan.
//                      This is a algorithm to determine the distance between 2 time
//                      series, and the point of best alignment between them in O( N ).
//                      This Odin implementation is approximately 3x and some change
//                      faster than the Python implementation.
//                      A tested it with the distance calculation and path alignment
//                      between 2 series of 2_500_000 elements each with a radius
//                      of 20. It took 15 minutes and 91 GB of memory.
//                      Normally with the normal DTW algorithm, one can only
//                      calculate the DTW distance between of 2 series of up to 3_000
//                      elements, because the algorithm complexity is O( N * M )
//                      where N is the len( seq_1 ), and M is the len( seq_2 ).
//                      This algorithm, fastDTW is O( N ) and the memory complexity
//                      is also O( N ).
//
// Original repository of this port:
//     Github - slaypni - fastdtw
//     Kazuaki Tanida
//
//     https://github.com/slaypni/fastdtw
//
//     Python implementation of FastDTW, which is an approximate Dynamic Time Warping (DTW)
//     algorithm that provides optimal or near-optimal alignments with an O(N) time and
//     memory complexity.
//
// The original paper from 2007:
//     Stan Salvador, and Philip Chan. "FastDTW: Toward accurate dynamic time warping in
//     linear time and space." Intelligent Data Analysis 11.5 (2007): 561-580.
//
//     https://cs.fit.edu/~pkc/papers/tdm04.pdf
//
//
// To compile do:
// 
//   $ make opti
//   $ make run
//
//

package main

import fdtw "./fastDTW"

import "core:fmt"
import "core:math"
import t "core:time"

main :: proc () {
    fmt.printfln( " Begin fastDTW - fast Dynamic Time Warping...\n" )

    a := [?]f32 { 1, 2, 3, 4, 5 }
    b := [?]f32 { 3, 4 } 

    a_slice := a[ : ]
    b_slice := b[ : ]

    // DTW
    fmt.println( "\nDTW\n" )
    dist_proc := fdtw.__difference
    dist, path := fdtw.__dtw( a_slice, b_slice, nil, dist_proc )
    fdtw.print_info( dist, path )


    // Fast DTW
    fmt.println( "\nFast DTW\n" )
    dist, path = fdtw.fastdtw( a_slice, b_slice, radius=1 ) 
    fdtw.print_info( dist, path )


    // Long Fast DTW
    num_samples := 100_000
    radius      := 20      // 100  // 1
    fmt.printfln( "\nLong Fast DTW\n\n    num_samples: %v\n" + 
                  "    radius: %v\n",
                  num_samples, radius )

    // num_samples := 1_500_000

    // 15 minuts and 45 seconds on one core of my machine,
    //    91 GB of used RAM - 7.5 GB for the OS and Open Apps.
    // num_samples := 2_500_000   // radius := 20   

  
    a_slice = fdtw.gen_slice( num_samples, 0 )
    b_slice = fdtw.gen_slice( num_samples, 20 )

    // Time START how long it took
    stopwatch: t.Stopwatch
    t.stopwatch_start( & stopwatch )

    dist, path = fdtw.fastdtw( a_slice, b_slice, radius )
  
    // Time STOP how long it took
    t.stopwatch_stop( & stopwatch )
    duration_1 := t.stopwatch_duration( stopwatch )
    t.stopwatch_reset( & stopwatch )

    fdtw.print_info( dist, path )

    fmt.printfln("\n Execution duration : %v  \n", duration_1 )

    fmt.printfln( "\n\n...end fastDTW - fast Dynamic Time Warping.\n" )

} 