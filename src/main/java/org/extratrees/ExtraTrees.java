package org.extratrees;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;

public class ExtraTrees extends AbstractTrees<BinaryTree> {
	Matrix input;
	double[] output;
	double[] outputSq;
	static double zero=1e-6;
	/** later shuffled and used for choosing random columns at each node */
	ArrayList<Integer> cols;
	
	// defined in AbstractTrees:
	//ArrayList<BinaryTree> trees;
	
	/** number of random cuts tried for each feature */
	int numRandomCuts = 1;
	/** whether random cuts are totally uniform or evenly uniform */
	boolean evenCuts = false;

	public ExtraTrees(Matrix input, double[] output) {
		if (input.nrows!=output.length) {
			throw(new IllegalArgumentException("Input and output do not have same length."));
		}
		this.input = input;
		this.output = output;
		this.outputSq = new double[output.length];
		for (int i=0; i<output.length; i++) {
			this.outputSq[i] = this.output[i]*this.output[i]; 
		}
		// making cols list for later use:
		this.cols = new ArrayList<Integer>(input.ncols);
		for (int i=0; i<input.ncols; i++) {
			cols.add(i);
		}
	}
	
	public boolean isEvenCuts() {
		return evenCuts;
	}
	
	/**
	 * @param evenCuts - whether the random cuts (if more than 1) are
	 * sampled from fixed even intervals (true) 
	 * or just sampled ordinary uniform way (false)
	 */
	public void setEvenCuts(boolean evenCuts) {
		this.evenCuts = evenCuts;
	}
	
	public int getNumRandomCuts() {
		return numRandomCuts;
	}
	
	public void setNumRandomCuts(int numRandomCuts) {
		this.numRandomCuts = numRandomCuts;
	}
	
	/** Builds trees with ids */
	public ArrayList<BinaryTree> buildTrees(int nmin, int K, int nTrees, int[] ids) {
		ArrayList<BinaryTree> trees = new ArrayList<BinaryTree>(nTrees);
		ShuffledIterator<Integer> cols = new ShuffledIterator<Integer>(this.cols);
		for (int t=0; t<nTrees; t++) {
			trees.add( this.buildTree(nmin, K, ids, cols) );
		}
		return trees;
	}
	
	/** Average of several trees: */
	public static double getValue(ArrayList<BinaryTree> trees, double[] input) {
		double output = 0;
		for(BinaryTree t : trees) {
			output += t.getValue(input);
		}
		return output/trees.size();
	}

	/** Average of several trees, using nmin as depth */
	public static double getValue(ArrayList<BinaryTree> trees, double[] input, int nmin) {
		double output = 0;
		for(BinaryTree t : trees) {
			output += t.getValue(input, nmin);
		}
		return output/trees.size();
	}
	
	/**
	 * @param input
	 * @return matrix of predictions where
	 * output[i, j] gives prediction made for i-th row of input by j-th tree. 
	 */
	public Matrix getAllValues(Matrix input) {
		Matrix out = new Matrix( input.nrows, trees.size() );
		// temporary vector:
		double[] temp = new double[input.ncols];
		for (int row=0; row<input.nrows; row++) {
			input.copyRow(row, temp);
			for (int j=0; j<trees.size(); j++) {
				out.set( row, j, trees.get(j).getValue(temp) );
			}
		}
		return out;
	}
	

	/**
	 * Object method, using the trees stored by learnTrees(...) method.
	 * @param input
	 * @return
	 */
	public double[] getValues(Matrix input) {
		return getValues(this.trees, input);
	}

	/** Average of several trees for many samples */
	public static double[] getValues(ArrayList<BinaryTree> trees, Matrix input) {
		double[] values = new double[input.nrows];
		double[] temp = new double[input.ncols];
		for (int row=0; row<input.nrows; row++) {
			// copying matrix row to temp:
			input.copyRow(row, temp);
			values[row] = getValue(trees, temp);
		}
		return values;
	}
	

	/**
	 * @param nmin - number of elements in leaf node
	 * @param K    - number of choices
	 */
	@Override
	public BinaryTree buildTree(int nmin, int K) {
		// generating full list of ids:
		int[]    ids = new int[output.length];
		for (int i=0; i<ids.length; i++) {
			ids[i] = i;
		}
		ShuffledIterator<Integer> cols = new ShuffledIterator<Integer>(this.cols);
		return buildTree(nmin, K, ids, cols);
	}
	
	/**
	 * 
	 * @param nmin
	 * @param K
	 * @param ids
	 * @param randomCols - passed to save memory (maybe not needed)
	 * @return
	 */
	public BinaryTree buildTree(int nmin, int K, int[] ids, ShuffledIterator<Integer> randomCols) {
		if (ids.length<nmin) {
			return makeLeaf(ids);
		}
		// doing a shuffle of cols:
		randomCols.reset();
		
		// trying K trees or the number of non-constant columns,
		// whichever is smaller:
		int k = 0, col_best=-1;
		double score_best = Double.NEGATIVE_INFINITY;
		boolean leftConst = false, rightConst = false;
		int countLeftBest = 0, countRightBest = 0;
		double t_best=Double.NaN;
		//for (int i=0; i<randomCols.size(); i++) {
		while( randomCols.hasNext() ) {
			int col = randomCols.next();
			// calculating columns min and max:
			double col_min = Double.POSITIVE_INFINITY;
			double col_max = Double.NEGATIVE_INFINITY;
			for (int n=0; n<ids.length; n++) {
				double v = input.get(ids[n], col);
				if ( v<col_min ) { col_min = v; }
				if ( v>col_max ) { col_max = v; }
			}
			if (col_max-col_min < zero) {
				// skipping, because column is constant
				continue;
			}
			// picking random test point numRepeatTries:
			double diff = (col_max-col_min);
			for (int repeat=0; repeat<this.numRandomCuts; repeat++) {
				double t;
				if (evenCuts) {
					double iStart = col_min + repeat*diff/numRandomCuts;
					double iStop  = col_min + (repeat+1)*diff/numRandomCuts;
					t = Math.random()*(iStop-iStart) + iStart;
				} else {
					t = Math.random()*diff + col_min;
				}
				
				// calculating score:
				int countLeft=0, countRight=0;
				double sumLeft=0, sumRight=0;
				double sumSqLeft=0, sumSqRight=0;
				for (int n=0; n<ids.length; n++) {
					if (input.get(ids[n], col) < t) {
						countLeft++;
						sumLeft   += output[  ids[n]];
						sumSqLeft += outputSq[ids[n]];
					} else {
						countRight++;
						sumRight   += output[  ids[n]];
						sumSqRight += outputSq[ids[n]];
					}
				}
				// calculating score:
				double varLeft  = sumSqLeft/countLeft  - (sumLeft/countLeft)*(sumLeft/countLeft);
				double varRight = sumSqRight/countRight- (sumRight/countRight)*(sumRight/countRight);
				double var = (sumSqLeft+sumSqRight)/ids.length - Math.pow((sumLeft+sumRight)/ids.length, 2.0);
				double score = 1 - (countLeft*varLeft + countRight*varRight) / ids.length / var;
				
				// if variance is 0
				if (var<zero*zero) {
					return makeLeaf(ids);
				}
				
				if (score>score_best) {
					score_best = score;
					col_best   = col;
					t_best     = t;
					leftConst  = (varLeft<zero*zero);
					rightConst = (varRight<zero*zero);
					countLeftBest  = countLeft;
					countRightBest = countRight;
				}
			}

			k++;
			if (k>=K) {
				// checked enough columns, stopping:
				break;
			}
		}
		// no score has been found, all inputs are constant:
		if (col_best<0) {
			return makeLeaf(ids);
		}
		
		// outputting the tree using the best score cut:
		int[] idsLeft  = new int[countLeftBest];
		int[] idsRight = new int[countRightBest];
		int nLeft=0, nRight=0;
		for (int n=0; n<ids.length; n++) {
			if (input.get(ids[n], col_best) < t_best) {
				// element goes to the left tree:
				idsLeft[nLeft] = ids[n];
				nLeft++;
			} else {
				// element goes to the right tree:
				idsRight[nRight] = ids[n];
				nRight++;
			}
		}
		BinaryTree bt = new BinaryTree();
		bt.column    = col_best;
		bt.threshold = t_best;
		bt.nSuccessors = ids.length;
		if (leftConst) { 
			bt.left = makeLeaf(idsLeft); // left child's output is constant 
		} else {  
			bt.left  = this.buildTree(nmin, K, idsLeft, randomCols); 
		}
		if (rightConst) {
			bt.right = makeLeaf(idsRight); // right child's output is constant
		} else {
			bt.right = this.buildTree(nmin, K, idsRight, randomCols);
		}
		// this value is used only for CV:
		bt.value  = bt.left.value*bt.left.nSuccessors + bt.right.value*bt.right.nSuccessors;
		bt.value /= bt.nSuccessors;
		return bt;
	}
	
	/**
	 * @param ids
	 * @return builds a leaf node and returns it with the given ids.
	 */
	public BinaryTree makeLeaf(int[] ids) {
		// terminal node:
		BinaryTree bt = new BinaryTree();
		bt.value = 0;
		bt.nSuccessors = ids.length;
		for (int n=0; n<ids.length; n++) {
			bt.value += output[ids[n]];
		}
		bt.value /= ids.length;
		return(bt);
	}
	
	
	
	/**
	 * @param ndata
	 * @param ndim
	 * @return
	 */
	
	
	public static ExtraTrees getSampleData(int ndata, int ndim) {
		double[] output = new double[ndata];
		double[] v = new double[ndata*ndim];
		for (int i=0; i<v.length; i++) {
			v[i] = Math.random();
		}
		Matrix m = new Matrix(v, ndata, ndim);
		// generate values for all outputs
		for (int row=0; row<output.length; row++) {
			m.set(row, 2, 0.5);
			output[row] = m.get(row, 1)+0.2*m.get(row, 3);
		}
		ExtraTrees et = new ExtraTrees(m, output);
		return et;
	}
	
	/**
	 * 
	 * @param trees
	 * @param testInput
	 * @param testOutput
	 * @return mean squared error on the test input and output.
	 */
	public static double getMeanSqError(ArrayList<BinaryTree> trees, Matrix testInput, double[] testOutput) {
		double error = 0;
		//double[] output_hat = getValues(trees, testInput);
		double[] temp = new double[testInput.ncols];
		for (int row=0; row<testInput.nrows; row++) {
			// copying matrix row to temp:
			for (int col=0; col<testInput.ncols; col++) {
				temp[col] = testInput.get(row, col);
			}
			error += Math.pow(getValue(trees, temp) - testOutput[row], 2);
		}

		return error/testOutput.length;
	}

	/** 
	 * mean squared error for CV
	 * @return mean squared error calculated only 
	 *         using data rows in testIds as test input-output.
	 **/
	public static double getMeanSqError(ArrayList<BinaryTree> trees, 
			Matrix testInput, double[] testOutput, 
			int nmin, int[] testIds)
	{
		double error = 0;
		//double[] output_hat = getValues(trees, testInput);
		double[] temp = new double[testInput.ncols];
		for (int n=0; n<testIds.length; n++) {
			int row = testIds[n];
			// copying matrix row to temp:
			for (int col=0; col<testInput.ncols; col++) {
				temp[col] = testInput.get(row, col);
			}
			error += Math.pow(getValue(trees, temp, nmin) - testOutput[row], 2);
		}

		return error/testIds.length;
	}

	public static double getMeanAbsError(ArrayList<BinaryTree> trees, Matrix testInput, double[] testOutput) {
		double error = 0;
		double[] output_hat = getValues(trees, testInput);
		for (int n=0; n<testOutput.length; n++) {
			error += Math.abs(output_hat[n] - testOutput[n]);
		}
		return error/testOutput.length;
	}
	
	public ArrayList<BinaryTree> buildTreeCV(int K, int nTrees) {
		int[] nmins = {2, 3, 5, 9, 14};
		int trainSize = (int)(2d/3d*output.length);
		
		// randomizing data:
		Integer[] ids = new Integer[this.output.length];
		for (int n=0; n<ids.length; n++) { ids[n]=n; }
		Collections.shuffle( Arrays.asList(ids) );
		
		// choosing trainSize of them as training data, others as test:
		int[] idsTrain = new int[trainSize];
		int[] idsTest  = new int[output.length - trainSize];
		for (int n=0; n<idsTrain.length; n++) {
			idsTrain[n] = ids[n];
		}
		for (int n=0; n<idsTest.length; n++) {
			idsTest[n]  = ids[n+idsTrain.length];
		}
		
		// building model:
		ArrayList<BinaryTree> t = this.buildTrees(2, K, nTrees, idsTrain);
		double[] errors = new double[nmins.length];
		double error_best = Double.POSITIVE_INFINITY;
		int nmin_best = nmins[0];
		for (int i=0; i<nmins.length; i++) {
			errors[i] = getMeanSqError(t, this.input, this.output, nmins[i], idsTest);
			if (errors[i]<error_best) {
				nmin_best = nmins[i];
				error_best = errors[i];
			}
		}
		//System.out.println( Arrays.toString(errors) );
		//System.out.println( String.format("Using nmin=%d based on CV.", nmin_best) );
		// building full model:
		ArrayList<BinaryTree> trees = this.buildTrees(nmin_best, K, nTrees);
		return trees;
	}
	

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		int ndata  = 10000, ndim=7;
		int nTrees = 15;
		
		ExtraTrees et = getSampleData(ndata, ndim);
		Date t3 = new Date();
//		ArrayList<BinaryTree> m = et.buildTreeCV(ndim, nTrees);
		Date t4 = new Date();
		System.out.println( "Took: " + (t4.getTime()-t3.getTime())/1000.0 + "s");
		//int x=1;
		//if (x==1) return;
		Date t1 = new Date();
		et.learnTrees(2, 6, nTrees);
		ArrayList<BinaryTree> trees = et.trees;
		Date t2 = new Date();
		
		ExtraTrees et2 = getSampleData(1000, ndim);
		double[] output_hat = getValues(trees, et2.input);
		for (int row=0; row<et2.output.length; row++) {
			System.out.print( String.format("%d\t%1.3f %1.3f", row, et2.output[row], output_hat[row]) );
			System.out.println();
		}
		System.out.println( "Took: " + (t2.getTime()-t1.getTime())/1000.0 + "s");

		int[] ids = new int[et2.output.length];
		for (int n=0; n<ids.length; n++) { ids[n]=n; }
		double e1 = getMeanSqError(trees, et2.input, et2.output);
		double e2 = getMeanSqError(trees, et2.input, et2.output, 5, ids );
		System.out.println( "Error: " + e1 );
		System.out.println( "Error: " + e2 );
		

//		System.out.println( Arrays.toString( bt.countColumns(m.ncols)) );
	}

}
