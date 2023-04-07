package org.extratrees.data;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Column-oriented SparseMatrix, using HashMaps
 * 
 * @author Jaak Simm
 *
 */
public class CSparseMatrix implements Array2D {
	ArrayList<HashMap<Integer, Double>> cols;
	int nrows;
	
	/**
	 * Creates a sparse matrix with the given data.
	 * @param rows   array of row ids
	 * @param cols   array of col ids
	 * @param values array of values
	 * @param nrows  number of rows
	 * @param ncols  number of columns
	 */
	public CSparseMatrix(
			int[] rows, int[] cols, double[] values, 
			int nrows, int ncols)
	{
		if (rows.length != cols.length) {
			throw new RuntimeException("Length of rows and cols does not equal.");
		}
		if (rows.length != values.length) {
			throw new RuntimeException("Length of rows and values does not equal.");
		}
		this.cols  = new ArrayList<HashMap<Integer,Double>>( ncols );
		this.nrows = nrows;
		for (int i=0; i < ncols; i++) {
			this.cols.add( new HashMap<Integer,Double>());
		}
		/** adding data */
		for (int i=0; i < rows.length; i++) {
			if (ncols <= cols[i]) {
				throw new IllegalArgumentException(String.format("ncols (%d) is smaller or equal than cols[%d] (%d)",
						ncols, i, cols[i]));
			}
			if (nrows <= rows[i]) {
				throw new IllegalArgumentException(String.format("nrows (%d) is smaller or equal than rows[%d] (%d)",
						nrows, i, rows[i]));
			}
			this.cols.get( cols[i] ).put( rows[i], values[i] );
		}
		
	}

	@Override
	public double get(int row, int col) {
		Double d = cols.get(col).get(row);
		return d == null ? 0  : d;
	}

	@Override
	public int ncols() {
		return cols.size();
	}

	@Override
	public int nrows() {
		return nrows;
	}

	@Override
	public Row getRow(int row) {
		return new SMRow(row);
	}
	
	public class SMRow implements Row {
		int row;
		
		public SMRow(int row) { this.row = row; }
		@Override
		public double get(int col) { 
			return CSparseMatrix.this.get(row, col); 
		}
	}
	
	public int getNumNonZero() {
		int nnzero = 0;
		for (HashMap<Integer, Double> c : cols) {
			nnzero += c.size();
		}
		return nnzero;
	}

}
