package org.extratrees;

public interface Aggregator<T extends AbstractBinaryTree> {
	void processLeaf(T leaf);
}
