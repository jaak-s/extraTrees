package org.extratrees;

/**
 * 
 * @author jaak
 *
 * @param <T> Tree type
 * @param <D>
 * @param <D> Value type
 */

public interface Aggregator<D> {
	D getResult();
	void processLeaf(D leafValue);
}
