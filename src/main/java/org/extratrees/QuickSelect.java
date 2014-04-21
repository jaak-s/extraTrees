package org.extratrees;

import java.util.ArrayList;

public class QuickSelect {
	/**
     * Quick selection algorithm.
     * Places the item closest to the k-th quantile in a[k*size()-1] and returns it.
     * @param a an ArrayList of Double items.
     * @param k the desired rank from 0 (smallest) to 1 (largest).
     * @return the selected value, i.e. a[k-1].
	 */
	public static double quickSelect( ArrayList<Double> a, double k) {
		int element = 1+(int)Math.round( k*a.size() );
		if (element>a.size()) {
			element = a.size();
		}
		return quickSelect(a, element );
	}
	
	
    /**
     * Quick selection algorithm.
     * Places the kth smallest item in a[k-1] and returns it.
     * @param a an ArrayList of Double items.
     * @param k the desired rank (1 is minimum) in the entire ArrayList.
     */
    public static double quickSelect( ArrayList<Double> a, int k ) {
    	if (k <= 0) {
    		return Double.NaN;
    	}
        quickSelect( a, 0, a.size() - 1, k );
        return a.get(k-1);
    }
    
    /**
     * Internal selection method that makes recursive calls.
     * Uses median-of-three partitioning and a cutoff of 10.
     * Places the kth smallest item in a[k-1].
     * @param a an ArrayList of Double items.
     * @param low the left-most index of the subarray.
     * @param high the right-most index of the subarray.
     * @param k the desired rank (1 is minimum) in the entire array.
     */
    private static void quickSelect( ArrayList<Double> a, int low, int high, int k ) {
        if( low + CUTOFF > high )
            insertionSort( a, low, high );
        else {
            // Sort low, middle, high
            int middle = ( low + high ) / 2;
            if( a.get( middle ).compareTo( a.get( low ) ) < 0 )
                swapReferences( a, low, middle );
            if( a.get( high ).compareTo( a.get( low ) ) < 0 )
                swapReferences( a, low, high );
            if( a.get( high ).compareTo( a.get( middle ) ) < 0 )
                swapReferences( a, middle, high );
            
            // Place pivot at position high - 1
            swapReferences( a, middle, high - 1 );
            Double pivot = a.get( high - 1 );
            
            // Begin partitioning
            int i, j;
            for( i = low, j = high - 1; ; ) {
                while( a.get( ++i ).compareTo( pivot ) < 0 )
                    ;
                while( pivot.compareTo( a.get( --j ) ) < 0 )
                    ;
                if( i >= j )
                    break;
                swapReferences( a, i, j );
            }
            
            // Restore pivot
            swapReferences( a, i, high - 1 );
            
            // Recurse; only this part changes
            if( k <= i )
                quickSelect( a, low, i - 1, k );
            else if( k > i + 1 )
                quickSelect( a, i + 1, high, k );
        }
    }
    
    
    /**
     * Internal insertion sort routine for subarrays
     * that is used by quicksort.
     * @param a an ArrayList of Double items.
     * @param low the left-most index of the subarray.
     * @param n the number of items to sort.
     */
    private static void insertionSort( ArrayList<Double> a, int low, int high ) {
        for( int p = low + 1; p <= high; p++ ) {
            Double tmp = a.get(p);
            int j;
            
            for( j = p; j > low && tmp.compareTo( a.get(j-1) ) < 0; j-- )
                a.set( j, a.get(j-1) );
            a.set( j, tmp );
        }
    }
    
    
    private static final int CUTOFF = 10;
    
    /**
     * Method to swap to elements in an array.
     * @param a an ArrayList of Double items.
     * @param index1 the index of the first object.
     * @param index2 the index of the second object.
     */
    public static final void swapReferences( ArrayList<Double> a, int index1, int index2 ) {
        Double tmp = a.get( index1 );
        a.set( index1, a.get( index2 ) );
        a.set( index2, tmp );
    }
}
