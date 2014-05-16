package org.extratrees;

import java.util.Set;

/**
 * All subclasses should have their generic argument equal to itself, i.e X<X>.
 * Otherwise getItself() will break. 
 * 
 * @author jaak
 *
 * @param <T>
 */
public abstract class AbstractBinaryTree <T extends AbstractBinaryTree<T>> {
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

	/** tasks that are active in this thread */
	Set<Integer> tasks;
	
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


	@SuppressWarnings("unchecked")
	public T getItself() {
		return (T)this;
	}
}
