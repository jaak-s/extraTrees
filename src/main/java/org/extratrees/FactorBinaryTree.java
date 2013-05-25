package org.extratrees;

public class FactorBinaryTree extends AbstractBinaryTree {
	/** tree for elements below threshold.
	 * if left==null, it is a leaf node
     * if left!=null, not a leaf
	 *  */
	public FactorBinaryTree left;
	/** tree for elements equal or above threshold. */
	public FactorBinaryTree right;
	/** Value of the node (estimated by its nodes), value of the node: 
	 *  from 0 to (#numFactors-1).
	 *  Non-leaf nodes also store value, allowing to change size of final nodes on-the-fly. */
	//public double value=Double.NEGATIVE_INFINITY;
	public int value;
	
	public FactorBinaryTree() {
		
	}

	/**
	 * Returns the value from the whole tree.
	 * @param input
	 * @return
	 */
	public int getValue(double[] input) {
		if (left==null) {
			return value; // leaf node
		}
		if (input[column]<threshold) {
			return left.getValue(input);
		}
		return right.getValue(input);
	}
	

	/**
	 * @param x
	 * @param task
	 * @return return multitask value for given input and task
	 */
	public int getValueMT(double[] input, int task) {
		if (left==null) {
			return value;
		}
		if (column<0) {
			// task cut:
			if (left.tasks.contains(task)) {
				return left.getValueMT(input, task);
			}
			return right.getValueMT(input, task);
		}
		// feature cut
		if (input[column]<threshold) {
			return left.getValueMT(input, task);
		}
		return right.getValueMT(input, task);

	}


	/**
	 * Returns values for data points, each data point is a row of the matrix.
	 * @param input
	 * @return
	 */
	public int[] getValues(Matrix input) {
		int[] values = new int[input.nrows];
		double[] temp = new double[input.ncols];
		for (int row=0; row<input.nrows; row++) {
			// copying matrix row to temp:
			for (int col=0; col<input.ncols; col++) {
				temp[col] = input.get(row, col);
			}
			values[row] = getValue(temp);
		}
		return values;
	}

	/**  uses nmin to choose the depth:
	 * @param input - input vector
	 * @param nmin  - number of elements in the final node (used for value).
	 * @return value in the tree for <b>input</b>.
	 */
	public int getValue(double[] input, int nmin) {
		if (this.nSuccessors<nmin || this.left==null) {
			return value; // leaf node OR below nmin
		}
		if (input[column]<threshold) {
			return left.getValue(input, nmin);
		}
		return right.getValue(input, nmin);
	}
	


	/**
	 * uses nmin to choose the depth:
	 * @param input
	 * @param nmin
	 * @return value in the tree for all <b>inputs</b> (each row in the matrix is an input).
	 */
	public int[] getValues(Matrix input, int nmin) {
		int[] values = new int[input.nrows];
		double[] temp = new double[input.ncols];
		for (int row=0; row<input.nrows; row++) {
			// copying matrix row to temp:
			for (int col=0; col<input.ncols; col++) {
				temp[col] = input.get(row, col);
			}
			values[row] = getValue(temp, nmin);
		}
		return values;
	}
	
	/** count the number of times each column is used */
	public int[] countColumns(int ncols) {
		int[] counts = new int[ncols];
		this.countColumns( counts );
		return(counts);
	}

	/** counts the number of times each column is used recursively */
	private void countColumns(int[] counts) {
		if (left==null) {
			return;
		}
		counts[this.column]++;
		this.left.countColumns(counts);
		this.right.countColumns(counts);
	}

}
