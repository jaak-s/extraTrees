package org.extratrees;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Random;

public class ShuffledIterator<E> implements Iterator<E> {
	ArrayList<E> data;
	int nextIndex = 0;
	Random random;
	
	/**
	 * @param c - creates iterator for Collection {@code c}, 
	 * elements are returned in random order. 
	 */
	public ShuffledIterator(Collection<E> c) {
		data = new ArrayList<E>(c);
		random = new Random();
	}

	@Override
	public boolean hasNext() {
		return nextIndex<data.size();
	}

	@Override
	public E next() {
		if (nextIndex>=data.size()) {
			throw new NoSuchElementException();
		}
		int i = random.nextInt(data.size()-nextIndex);
		// swapping elements nextIndex and nextIndex+i:
		E e1 = data.get(nextIndex+i);
		data.set(nextIndex+i, data.get(nextIndex) );
		data.set(nextIndex, e1 );
		nextIndex++;
		// returning previous nextIndex+i:
		return e1;
	}

	@Override
	public void remove() {
		throw new RuntimeException("unimplemented");
	}
	
	/**
	 * Resets the iterator so we can start from a new shuffle again. 
	 */
	public void reset() {
		nextIndex = 0;
	}

}
