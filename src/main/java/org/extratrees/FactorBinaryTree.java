package org.extratrees;

public class FactorBinaryTree extends AbstractBinaryTree<FactorBinaryTree, Integer> {
	public FactorBinaryTree() {
	}
	
	@Override
	public Integer getNA() {
		return -1;
	}

	/**
	 * Returns values for data points, each data point is a row of the matrix.
	 * @param input
	 * @return
	 */
	public int[] getValues(Matrix input) {
		int[] values = new int[input.nrows()];
		double[] temp = new double[input.ncols()];
		for (int row=0; row < input.nrows(); row++) {
			// copying matrix row to temp:
			for (int col=0; col < input.ncols(); col++) {
				temp[col] = input.get(row, col);
			}
			values[row] = getValue(temp);
		}
		return values;
	}

	/**
	 * uses nmin to choose the depth:
	 * @param input
	 * @param nmin
	 * @return value in the tree for all <b>inputs</b> (each row in the matrix is an input).
	 */
	public int[] getValues(Matrix input, int nmin) {
		int[] values = new int[input.nrows()];
		double[] temp = new double[input.ncols()];
		for (int row=0; row < input.nrows(); row++) {
			// copying matrix row to temp:
			for (int col=0; col<input.ncols(); col++) {
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
