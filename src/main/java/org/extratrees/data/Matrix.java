package org.extratrees.data;

public class Matrix implements Array2D {
	public double[] v;
	public int ncols, nrows;
	
	/** Matrix filled with 0s. */
	public Matrix(int nrows, int ncols) {
		this(new double[nrows*ncols], nrows, ncols);
	}
	
	/** Matrix filled with v. */
	public Matrix(double[] v, int nrows, int ncols) {
		if (v.length != nrows*ncols) {
			throw( new IllegalArgumentException(
				"Length of v ("+v.length+") is not equal to nrows*ncols ("+
				nrows+"*"+ncols+")"));
		}
		this.v = v;
		this.nrows = nrows;
		this.ncols = ncols;
	}
	
	/** row is from 0 to (nrow-1) 
	 *  col is from 0 to (ncol-1)
	 * */
	public void set(int row, int col, double value) {
		this.v[row + col*nrows] = value;
	}
	
	@Override
	public double get(int row, int col) {
		return this.v[row + col*nrows];
	}
	
	@Override
	public int ncols() {
		return ncols;
	}
	
	@Override
	public int nrows() {
		return nrows;
	}
	
	@Override
	public Row getRow(int row) {
		return new MRow(row);
	}
	
	/**
	 *  Provide access to rows of the matrix without copy-pasting
	 */
	public class MRow implements Row {
		int row;
		public MRow(int row) {
			this.row = row;
		}
		@Override
		public double get(int col) {
			return Matrix.this.get(row, col);
		}
	}
	
	/**
	 * Copies values from given row to the vector. 
	 * @param row
	 * @param vector
	 */
	public void copyRow(int row, double[] vector) {
		for (int col=0; col < this.ncols; col++) {
			vector[col] = this.get(row, col);
		}
	}

	/**
	 * Copies values from given row to the vector. 
	 * @param col
	 * @param vector
	 */
	public void copyCol(int col, double[] vector) {
		for (int row=0; row < this.nrows; row++) {
			vector[row] = this.get(row, col);
		}
	}

	/**
	 * Copies values from given row to the vector. 
	 * @param row
	 * @param vector
	 */
	public double[] getCol(int col) {
		double[] colValues = new double[ nrows ];
		copyCol( col, colValues );
		return colValues;
	}

	public void square() {
		for (int i=0; i<v.length; i++) {
			v[i] *= v[i];
		}
	}

	public String toString() {
		StringBuilder out = new StringBuilder();
		for (int row=0; row < nrows; row++) {
			for (int col=0; col < ncols; col++) {
				out.append( String.format("%1.4f", get(row,col) ) );
				out.append( " " );
			}
			out.append("\n");
		}
		return out.toString();
	}
	
	public boolean hasNaN() {
		for (int i = 0; i < v.length; i++) {
			if ( Double.isNaN(v[i]) ) {
				return true;
			}
		}
		return false;
	}
	
	public static void main(String[] args) {
		double[] v = new double[40];
		for (int i=0; i < v.length; i++) {
			v[i] = i+1;
		}
		Matrix m = new Matrix(v, 10, 4);
		System.out.println(m);
	}
}
