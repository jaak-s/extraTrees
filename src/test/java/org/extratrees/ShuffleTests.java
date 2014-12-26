package org.extratrees;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.NoSuchElementException;
import java.util.Random;

import org.junit.Test;

public class ShuffleTests {

	@Test
	public void testShuffle() {
		ArrayList<Integer> data = new ArrayList<Integer>();
		int n=10;
		for (int i=0; i<n; i++) {
			data.add(i);
		}
		ShuffledIterator<Integer> si = new ShuffledIterator<Integer>(data, new Random());
		HashSet<Integer> result = new HashSet<Integer>();
		while (si.hasNext()) {
			int x = si.next();
			assertTrue( data.contains(x) );
			result.add( x );
		}
		assertEquals("some elements are missing", n, result.size() );
		
		// no more elements should be returned:
		try {
			int x = si.next();
			fail("NoSuchElementsException should be thrown but got "+x);
		} catch(NoSuchElementException e) {
			
		}
		// resetting and iterating again:
		si.reset();
		result = new HashSet<Integer>();
		while (si.hasNext()) {
			int x = si.next();
			assertTrue( data.contains(x) );
			result.add( x );
		}
		assertEquals("some elements are missing", n, result.size() );
		
	}

}
