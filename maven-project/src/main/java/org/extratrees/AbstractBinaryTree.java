package org.extratrees;

import java.io.Serializable;
import java.util.Set;

import org.extratrees.data.Row;

/**
 * All subclasses should have their generic argument equal to itself, 
 * i.e. X extends AbstractBinaryTree<X>.
 * Otherwise getItself() will break. 
 * 
 * @author Jaak Simm
 *
 * @param <T> class that extends ABT
 * @param <D> class of value
 */
public abstract class AbstractBinaryTree <T extends AbstractBinaryTree<T, D>, D> implements Serializable {
	private static final long serialVersionUID = 5530548509407016040L;
	
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
	public T getLeaf(Row input) {
		if (left==null) {
			return getItself();
		}
		if (Double.isNaN(input.get(column))) {
			return null;
		}
		if (input.get(column)<threshold) {
			return left.getLeaf(input);
		}
		return right.getLeaf(input);
	}
	
	public D getValue(Row input) {
		return valueFromLeaf( getLeaf(input) );
	}
	
	/**
	 * @param x
	 * @param task
	 * @return return multitask value for given input and task
	 */
	public T getLeafMT(Row input, int task) {
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
		if (Double.isNaN( input.get(column)) ) {
			return null;
		}
		// feature cut
		if (input.get(column) < threshold) {
			return left.getLeafMT(input, task);
		}
		return right.getLeafMT(input, task);

	}

	public D getValueMT(Row input, int task) {
		return valueFromLeaf( getLeafMT(input, task) );
	}
	
	/**  uses nmin to choose the depth:
	 * @param input - input vector
	 * @param nmin  - number of elements in the final node (used for value).
	 * @return value in the tree for <b>input</b>.
	 */
	public T getLeaf(double[] input, int nmin) {
		if (this.nSuccessors <= nmin || this.left == null) {
			return getItself(); // leaf node OR below nmin
		}
		if (Double.isNaN(input[column])) {
			return null;
		}
		if (input[column] < threshold) {
			return left.getLeaf(input, nmin);
		}
		return right.getLeaf(input, nmin);
	}
	
	public D getValue(double[] input, int nmin) {
		return valueFromLeaf( getLeaf(input, nmin) );
	}

	
	/**
	 * 
	 * @param leaf
	 * @return leaf.value if leaf is not null or getNA() if leaf is null.
	 */
	private D valueFromLeaf(T leaf) {
		if (leaf == null) return getNA();
		return leaf.getValue();
	}


	@SuppressWarnings("unchecked")
	public T getItself() {
		return (T)this;
	}
}
