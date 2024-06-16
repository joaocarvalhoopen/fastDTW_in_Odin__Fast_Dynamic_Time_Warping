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

package fastDTW

import "core:fmt"
import "core:math"
import "core:slice"
import t "core:time"


MAX_SEQ_NUM_ELEMENTS :: 10_000_000

//  function fastdtw
// 
//  return : the approximate distance between 2 time series with O(N)
//           time and memory complexity
//
//  Parameters
//  ----------
//     x : []f32 slice input 1
//
//     y : []f32 slice input 2
//
//     radius : int
//         size of neighborhood when expanding the path. A higher value will
//         increase the accuracy of the calculation but also increase time
//         and memory consumption. A radius equal to the size of x and y will
//         yield an exact DTW Dynamic Time Warping calculation.
//
//     The method for calculating the distance between x[i] and y[j] is
//     abs(x[i] - y[j]) .
//
//  Returns
//  -------
//     distance : float
//        the approximate distance between the 2 time series
//     path : list
//        list of indexes for the inputs x and y
//
//         Examples
//         --------
//         >>> import numpy as np
//         >>> import fastdtw
//         >>> x = np.array([1, 2, 3, 4, 5], dtype='float')
//         >>> y = np.array([2, 3, 4], dtype='float')
//         >>> fastdtw.fast_dtw(x, y)
//         (2.0, [(0, 0), (1, 0), (2, 1), (3, 2), (4, 2)])
//
fastdtw :: proc ( x : []f32, y : []f32, radius : int = 1 ) ->
                ( distance : f32, path : Path ) {
 
    return __fastdtw( x, y, radius )
}

__fastdtw :: proc ( x : []f32, y : []f32, radius : int ) ->
                  ( distance : f32, path : Path ) {
    
    min_time_size := radius + 2

    if len( x ) < min_time_size || len( y ) < min_time_size {
        
        distance, path = dtw( x, y, __difference )    
        return distance, path
    }
    
    x_shrinked := __reduce_by_half( x )
    y_shrinked := __reduce_by_half( y )

    distance, path = __fastdtw( x_shrinked, y_shrinked, radius=radius )
  
    window := __expand_window( path,
                               len( x ),
                               len( y ),
                               radius )
  
    delete ( x_shrinked )
    delete ( y_shrinked )
    delete ( path )

    distance, path = __dtw( x, y, window, __difference )

    delete( window )

    return distance, path
}

__reduce_by_half :: proc ( x : [ ]f32 ) -> [ ]f32 {

    res_list := make( [ dynamic ]f32, len=0, cap=len( x ) / 2 )
    if res_list == nil do panic( "ERROR: __reduce_by_half, res_list, Memory allocation failed" )
   
    len_x := len( x ) 

    for i : int = 0; i < ( len_x - len_x % 2 ); i += 2 {
        append_elem( & res_list, ( x[ i ] + x[ i + 1 ] ) / 2 )
    }

    return res_list[ : ]
}

__difference :: #force_inline proc ( a : f32, b : f32 ) -> f32 {
    return math.abs( a - b )
}


//    function dtw
//
//    return the distance between 2 time series without approximation
//
//    Parameters
//    ----------
//    x : [ ]f32 slice of array 1
//
//    y : [ ]f32 slice of array 2
//
//    The distance function is abs( x[ i ] - y[ j ] ) .
//
//    Returns
//    -------
//    distance : f32
//        the approximate distance between the 2 time series
//    path : list
//        list of indexes for the inputs x and y
//
//    Examples
//    --------
//    >>> import numpy as np
//    >>> import fastdtw
//    >>> x = np.array( [ 1, 2, 3, 4, 5 ], dtype='float' )
//    >>> y = np.array( [ 2, 3, 4 ], dtype='float' )
//    >>> fastdtw.dtw( x, y )
//    (2.0, [(0, 0), (1, 0), (2, 1), (3, 2), (4, 2)])
//

dtw :: proc ( x :[ ]f32, y : [ ]f32, dist : Dist_Proc ) ->
              ( distance : f32, path : Path ) {
    
    distance, path = __dtw( x, y, nil, dist )
    return distance, path
}

Coord :: struct {
    i : int,
    j : int,
}

Window :: [ dynamic ]Coord

c2d_to_1d :: #force_inline proc ( i : int, j : int ) -> int {
    return i * MAX_SEQ_NUM_ELEMENTS + j
}

Path :: [ dynamic ]Coord

Dist_Proc :: proc ( a : f32, b : f32 ) -> f32

D_struct :: struct {
    dist : f32,
    i    : int,
    j    : int,
}

plus2D_1 :: #force_inline proc ( coord : Coord ) -> ( i, j : int ) {
    return coord.i + 1, coord.j + 1
}

__dtw :: proc ( a_slice : []f32, b_slice : []f32, window : Window, dist_proc : Dist_Proc ) ->
              ( res_dist : f32, path : Path ) {

    window := window

    len_x := len( a_slice )
    len_y := len( b_slice )

    if window == nil {
        window = make( Window, len=0, cap=len_x * len_y )
        if window == nil do panic( "Memory allocation failed" )   
        
        for i in 0 ..< len_x {
            for j in 0 ..< len_y {
                
                append_elem( & window, Coord{ i, j } )
            }
        }
    }

    D := make( map[ int ]D_struct, len( window ) + 4 )
    if D == nil do panic( "Memory allocation failed" )
    defer delete( D )

    inf_value := max( f32 )

    d_struct := D_struct{ 0, 0, 0 }

    D[ c2d_to_1d( 0, 0 ) ] = D_struct{ 0, 0, 0 }
    for coord in window {

        i := coord.i + 1
        j := coord.j + 1

        dt := abs( a_slice[ i - 1 ] - b_slice[ j - 1 ] )

        // D[i, j] = min((D[i-1, j][0]+dt, i-1, j), (D[i, j-1][0]+dt, i, j-1),
        //               (D[i-1, j-1][0]+dt, i-1, j-1), key=lambda a: a[0])

        // ( D[ i - 1, j     ][ 0 ] + dt, i - 1, j     )
        // ( D[ i,     j - 1 ][ 0 ] + dt,     i, j - 1 )
        // ( D[ i - 1, j - 1 ][ 0 ] + dt, i - 1, j - 1 )

        value_dist_a := inf_value
        value_a, ok_a := D[ c2d_to_1d( i - 1, j ) ]
        if ok_a do value_dist_a = value_a.dist
        a_dist := value_dist_a + dt 

        value_dist_b := inf_value
        value_b, ok_b := D[ c2d_to_1d( i, j - 1 ) ]
        if ok_b do value_dist_b = value_b.dist
        b_dist := value_dist_b + dt

        value_dist_c := inf_value
        value_c, ok_c := D[ c2d_to_1d( i - 1, j - 1 ) ]
        if ok_c do value_dist_c = value_c.dist
        c_dist := value_dist_c + dt
        
        // Find the minimum distance of 3 elements.
        if a_dist < b_dist {
            if a_dist < c_dist {
                // A
 
                d_struct.dist = a_dist
                d_struct.i    = i - 1
                d_struct.j    = j
            } else {
                // C

                d_struct.dist = c_dist
                d_struct.i    = i - 1
                d_struct.j    = j - 1
            }
        } else {
            if b_dist < c_dist {
                // B

                d_struct.dist = b_dist
                d_struct.i    = i
                d_struct.j    = j - 1

            } else {
                // C

                d_struct.dist = c_dist
                d_struct.i    = i - 1
                d_struct.j    = j - 1
            }
        }

        D[ c2d_to_1d( i, j ) ] = d_struct

    }

    path = make( [ dynamic ]Coord, len=0, cap=max( len_x, len_y ) * 2 )
    if path == nil do panic( "Memory allocation failed" )

    i, j := len_x, len_y
    for  ! ( i == 0 && j == 0 ) {        
        append( &path, Coord{ i - 1, j - 1 } )
        tmp := D[ c2d_to_1d( i, j ) ]
        i, j = tmp.i, tmp.j
    }

    slice.reverse( path[ : ]  )

    res_dist = D[ c2d_to_1d(len_x, len_y ) ].dist
    return  res_dist, path
}

__expand_window :: proc ( path : Path, len_x : int , len_y : int, radius : int ) ->
                        ( windows : Window ) {

    // ====>>>>> A


    // // Time START how long it took
    // stopwatch_1 : t.Stopwatch
    // t.stopwatch_start( & stopwatch_1 )


    // When we transform into a map we loose the order of the elements.
    // path_tmp := make( map[ int ]Coord, len( path ) * 2 )

    path_tmp := make( map[ int ]Coord, len( path ) * radius )

    if path_tmp == nil do panic( "Error: __expand_window(), path_, Memory allocation failed." )
    defer delete( path_tmp )
    
    for coord in path {
        path_tmp[ c2d_to_1d( coord.i, coord.j ) ] = coord
    }


    // // Time STOP how long it took
    // t.stopwatch_stop( & stopwatch_1 )
    // duration_1 := t.stopwatch_duration( stopwatch_1 )
    // t.stopwatch_reset( & stopwatch_1 )

    // fmt.printfln("\n\n Execution duration 1  copy_path -> path_tmp : %v", duration_1 )

 
    // // Time START how long it took
    // // stopwatch_1 : t.Stopwatch
    // t.stopwatch_start( & stopwatch_1 )
    

    // for i, j in path:
    //     for a, b in ((i + a, j + b)
    //                  for a in range(-radius, radius+1)
    //                  for b in range(-radius, radius+1)):
    //         path_.add((a, b))

    for coord in path {
        i, j := coord.i, coord.j
        for a in -radius ..< radius + 1 {
            for b in -radius ..< radius + 1 {

                key := c2d_to_1d( i + a, j + b )
                if ! ( key in path_tmp ) {
                    path_tmp[ key ] = Coord{ i + a, j + b } 
                }

                // path_tmp[ c2d_to_1d( i + a, j + b ) ] = Coord{ i + a, j + b }
            }
        }
    }


    // // Time STOP how long it took
    // t.stopwatch_stop( & stopwatch_1 )
    // duration_1 = t.stopwatch_duration( stopwatch_1 )
    // t.stopwatch_reset( & stopwatch_1 )

    // fmt.printfln("Execution duration 2.a  expand radius -> path_tmp : %v", duration_1 )


    // When we transform into a map we loose the order of the elements.
    window_tmp := make( map[ int ]Coord, len( path_tmp ) * 4 )
    if window_tmp == nil do panic( "Error: __expand_window(), window_, Memory allocation failed." )
    defer delete( window_tmp )


    // // Time START how long it took
    // // stopwatch_1 : t.Stopwatch
    // t.stopwatch_start( & stopwatch_1 )


    // for i, j in path_:
    //     for a, b in ((i * 2, j * 2), (i * 2, j * 2 + 1),
    //                  (i * 2 + 1, j * 2), (i * 2 + 1, j * 2 + 1)):
    //         window_.add((a, b))

    // ( i * 2,     j * 2 ),
    // ( i * 2,     j * 2 + 1 ),
    // ( i * 2 + 1, j * 2 ),
    // ( i * 2 + 1, j * 2 + 1 )


    delta := [4]Coord { Coord{ 0, 0 },
                        Coord{ 0, 1 },
                        Coord{ 1, 0 },
                        Coord{ 1, 1 } }

    for key, coord in path_tmp {
        i, j := coord.i, coord.j

        for coord_delta in delta {
            a, b := i * 2 + coord_delta.i, j * 2 + coord_delta.j
            
            window_tmp[ c2d_to_1d( a, b ) ] = Coord{ a, b }
        }
    }


    // // Time STOP how long it took
    // t.stopwatch_stop( & stopwatch_1 )
    // duration_1 = t.stopwatch_duration( stopwatch_1 )
    // t.stopwatch_reset( & stopwatch_1 )

    // fmt.printfln("Execution duration 3  fill_window_tmp -> window_tmp : %v", duration_1 )
    
    window := make( Window, len=0, cap=len( window_tmp ) )
    if window == nil do panic( "Error: __expand_window(), window, Memory allocation failed." )

    
    // // Time START how long it took
    // // stopwatch_1 : t.Stopwatch
    // t.stopwatch_start( & stopwatch_1 )


    // for i in range(0, len_x):
    //     new_start_j = None
    //     for j in range(start_j, len_y):
    //         if (i, j) in window_:
    //             window.append((i, j))
    //             if new_start_j is None:
    //                 new_start_j = j
    //         elif new_start_j is not None:
    //             break
    //     start_j = new_start_j


    start_j : int = 0

    new_start_j_local : int = 0

    for i in 0 ..< len_x {
        new_start_j : ^int = nil
        // new_start_j : ^i32 = nil
        for j in start_j ..< len_y {
            // if (i, j) in window_:
            if c2d_to_1d( i, j ) in window_tmp {   
                append_elem( & window, Coord{ i, j } )
                if new_start_j == nil {
                    new_start_j = & new_start_j_local
                    new_start_j^ = j
                }  
            } else if new_start_j != nil {
                break
            }
        }
        start_j = new_start_j^
    }

    // // Time STOP how long it took
    // t.stopwatch_stop( & stopwatch_1 )
    // duration_1 = t.stopwatch_duration( stopwatch_1 )
    // t.stopwatch_reset( & stopwatch_1 )
    
    // fmt.printfln("Execution duration 4  fill_final_window -> window : %v\n\n\n", duration_1 )

    return window
}

print_info :: proc ( dist : f32, path : Path ) {
    fmt.printfln( "Distance: %v", dist )
    fmt.printfln( "Path:" )
    fmt.printf( "[ " )
    for coord in path {
        fmt.printf( "(%v, %v), ", coord.i, coord.j )
    }
    fmt.printfln( "]" )
    fmt.printfln( "Distance: %v", dist )
}

gen_slice :: proc ( n : int, start_value : f32 = 0 ) -> [ ]f32 {
    a := make( [ ]f32, n )
    for i in int( start_value ) ..< n {
        a[i] = f32( i )
    }
    return a
}

