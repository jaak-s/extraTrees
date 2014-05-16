package org.extratrees;

public interface Aggregator<T extends AbstractBinaryTree<T>, D> {
	D getResult();
	void processLeaf(T leaf);
}
