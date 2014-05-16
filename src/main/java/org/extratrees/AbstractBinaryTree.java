package org.extratrees;

import java.util.Set;

/**
 * All subclasses should have their generic argument equal to itself, 
 * i.e. X extends AbstractBinaryTree<X>.
 * Otherwise getItself() will break. 
 * 
 * @author jaak
 *
 * @param <T> class that extends ABT
 * @param <D> class of value
 */
public abstract class AbstractBinaryTree <T extends AbstractBinaryTree<T, D>, D> {
	/** tree for elements below threshold.
	 * if left==null, it is a leaf node
     * if left!=null, not a leaf
	 *  */
	public T left;
	/** tree for elements equal or above threshold. */
	public T right;

	/** number of elements in the tree */
	public int    nSuccessors;
	/** feature ID used for cutting */
	public int    column=-1;
	/** threshold of cutting */
	public double threshold;
	/** value of the node (estimated by its samples).
	 *  Non-leaf nodes (may) also store value, allowing to change size of final nodes on-the-fly. */
	public D      value;

	/** tasks that are active in this thread */
	Set<Integer> tasks;
	
	public abstract D getNA();
	
	/** @return value of current node */
	public D getValue() {
		return value;
	}
	
	/**
	 * @param input the vector of input values
	 * @return the leaf node (BinaryTree) for the input
	 */
	public T getLeaf(double[] input) {
		if (left==null) {
			return getItself();
		}
		if (Double.isNaN(input[column])) {
			return null;
		}
		if (input[column]<threshold) {
			return left.getLeaf(input);
		}
		return right.getLeaf(input);
	}
	
	public D getValue(double[] input) {
		T leaf = getLeaf(input);
		if (leaf == null) return getNA();
		return leaf.getValue();
	}
	
	/**
	 * @param x
	 * @param task
	 * @return return multitask value for given input and task
	 */
	public T getLeafMT(double[] input, int task) {
		if (left==null) {
			return getItself();
		}
		if (column < 0) {
			// task cut:
			if (left.tasks.contains(task)) {
				return left.getLeafMT(input, task);
			}
			return right.getLeafMT(input, task);
		}
		if (Double.isNaN(input[column])) {
			return null;
		}
		// feature cut
		if (input[column]<threshold) {
			return left.getLeafMT(input, task);
		}
		return right.getLeafMT(input, task);

	}

	public D getValueMT(double[] input, int task) {
		T leaf = getLeafMT(input, task);
		if (leaf == null) return getNA();
		return leaf.getValue();
	}


	@SuppressWarnings("unchecked")
	public T getItself() {
		return (T)this;
	}
}
