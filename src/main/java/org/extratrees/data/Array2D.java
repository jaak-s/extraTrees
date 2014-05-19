package org.extratrees.data;

public interface Array2D {
	public double get(int row, int col);
	public int ncols();
	public int nrows();
	public Row getRow(int row);
}
