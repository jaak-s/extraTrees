package org.extratrees;
import java.util.ArrayList;


public class FactorExtraTrees extends AbstractTrees<FactorBinaryTree> {
	Matrix input;
	int[] output;
	/** number of factors: */
	int nFactors;
	//String[] factorNames;
	static double zero=1e-7;
	
	/** later shuffled and used for choosing random columns at each node: */
	ArrayList<Integer> cols;
	
	//ArrayList<FactorBinaryTree> trees;

	/** number of random cuts tried for each feature */
	int numRandomCuts = 1;
	/** whether random cuts are totally uniform or evenly uniform */
	boolean evenCuts = false;
	
	/**
	 * @param input    - matrix of inputs, each row is an input vector
	 * @param output   - array of ints from 0 to nFactors-1 (class label)
	 * @param nFactors - number of factors (class labels)
	 */
	public FactorExtraTrees(Matrix input, int[] output) {
		if (input.nrows!=output.length) {
			throw(new IllegalArgumentException("Input and output do not have same length."));
		}
		this.input = input;
		this.output = output;
		//this.outputSq = new double[output.length];
		//for (int i=0; i<output.length; i++) {
		//	this.outputSq[i] = this.output[i]*this.output[i]; 
		//}
		
		this.nFactors = 1;
		for (int i=0; i<output.length; i++) {
			if (output[i]<0) { 
				throw new RuntimeException("Bug: negative output (factor) values.");
			}
			if (nFactors<=output[i]) {
				nFactors = output[i]+1;
			}
		}
		
		// making cols list for later use:
		this.cols = new ArrayList<Integer>(input.ncols);
		for (int i=0; i<input.ncols; i++) {
			cols.add(i);
		}
	}
	
	public int getnFactors() {
		return nFactors;
	}
	
	public void setnFactors(int nFactors) {
		this.nFactors = nFactors;
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
	
	/**
	 * @param selection
	 * @return new ExtraTrees object with the same input and output data with
	 * only the selected trees specified by {@code selection}.
	 */
	public FactorExtraTrees selectTrees(boolean[] selection) {
		FactorExtraTrees newET = new FactorExtraTrees(input, output);
		newET.trees = new ArrayList<FactorBinaryTree>();
		for (int i=0; i<selection.length; i++) {
			if (!selection[i]) {
				continue;
			}
			newET.trees.add(this.trees.get(i));
		}
		return newET;
	}

	
	/** Builds trees with ids */
	/*
	public ArrayList<FactorBinaryTree> buildTrees(int nmin, int K, int nTrees, int[] ids) {
		ArrayList<FactorBinaryTree> trees = new ArrayList<FactorBinaryTree>(nTrees);
		ShuffledIterator<Integer> cols = new ShuffledIterator<Integer>(this.cols);
		for (int t=0; t<nTrees; t++) {
			trees.add( this.buildTree(nmin, K, ids, null, cols) );
		}
		return trees;		
	}*/

	/** Average of several trees: */
	public static int getValue(ArrayList<FactorBinaryTree> trees, double[] input, int nFactors) {
		int[] counts = new int[nFactors];
		for(FactorBinaryTree t : trees) {
			counts[ t.getValue(input) ]++;
		}
		return getMaxIndex(counts);
	}
	
	/** Average of several trees, using nmin as depth */
	/*
	public static double getValue(ArrayList<FactorBinaryTree> trees, double[] input, int nmin, int nFactors) {
		int[] counts = new int[nFactors];
		for(FactorBinaryTree t : trees) {
			counts[ t.getValue(input, nmin) ]++;
		}
		return getMaxIndex(counts);
	}*/

	/**
	 * @param values
	 * @return index of the max value (first one if there are many)
	 */
	public static int getMaxIndex(int[] values) {
		int maxIndex = -1;
		int maxValue = Integer.MIN_VALUE;
		// adding counts:
		for (int i=0; i<values.length; i++) {
			if (values[i]>maxValue) {
				maxValue = values[i];
				maxIndex = i;
			}
		}
		return maxIndex;
	}
	
	/**
	 * @param input
	 * @return matrix of predictions where
	 * output[i, j] gives prediction made for i-th row of input by j-th tree.
	 * All values are integers. 
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
	public int[] getValues(Matrix input) {
		return getValues(this.trees, input, this.nFactors);
	}

	/** Average of several trees for many samples */
	public static int[] getValues(ArrayList<FactorBinaryTree> trees, Matrix input, int nFactors) {
		int[]  values = new int[input.nrows];
		double[] temp = new double[input.ncols];
		for (int row=0; row<input.nrows; row++) {
			// copying matrix row to temp:
			for (int col=0; col<input.ncols; col++) {
				temp[col] = input.get(row, col);
			}
			values[row] = getValue(trees, temp, nFactors);
		}
		return values;
	}

	/**
	 * @param nmin - number of elements in leaf node
	 * @param K    - number of choices
	 */
	public FactorBinaryTree buildTree(int nmin, int K) {
		// generating full list of ids:
		int[]    ids = new int[output.length];
		for (int i=0; i<ids.length; i++) {
			ids[i] = i;
		}
		int[] unlabeledIds = null;
		if (unlabeled!=null) {
			unlabeledIds = new int[unlabeled.nrows];
			for (int i=0; i<unlabeledIds.length; i++) {
				unlabeledIds[i] = i;
			}
		}
		ShuffledIterator<Integer> cols = new ShuffledIterator<Integer>(this.cols);
		return buildTree(nmin, K, ids, unlabeledIds, cols);
	}
	
	/**
	 * @param counts
	 * @return Gini index ( 1 - sum(f_i) ), where f_i is 
	 *         the proportion of label i.<br>
	 *         Value 0 implies pure, high values mean noisy.  
	 */
	public static double getGiniIndex(int[] counts) {
		int sum = 0;
		int total = 0;
		for (int i=0; i<counts.length; i++) {
			sum   += counts[i]*counts[i];
			total += counts[i];
		}
		return 1 - sum / (double)( total*total );
	}
	
	public FactorBinaryTree buildTree(int nmin, int K, int[] ids, int[] unlabeledIds,
			ShuffledIterator<Integer> randomCols) 
	{
		if (ids.length<nmin) {
			return makeLeaf(ids);
		}
		// doing a shuffle of cols:
		randomCols.reset();
		
		// trying K trees or the number of non-constant columns,
		// whichever is smaller:
		int k = 0, col_best=-1;
		double score_best = Double.POSITIVE_INFINITY;
		boolean leftConst = false, rightConst = false;
		int countLeftBest = 0, countRightBest = 0;
		double t_best=Double.NaN;
		//for (int i=0; i<cols.size(); i++) {
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
			double diff = (col_max-col_min);
			for (int repeat=0; repeat<this.numRandomCuts; repeat++) {
				// picking random test point:
				double t;
				if (evenCuts) {
					double iStart = col_min + repeat*diff/numRandomCuts;
					double iStop  = col_min + (repeat+1)*diff/numRandomCuts;
					t = Math.random()*(iStop-iStart) + iStart;
				} else {
					t = Math.random()*diff + col_min;
				}
			
				// calculating QINI impurity index (0 - pure, 0 - noisy):
				int countLeft=0, countRight=0;
				int[] factorCountLeft  = new int[nFactors];
				int[] factorCountRight = new int[nFactors];
				//double sumLeft=0, sumRight=0;
				//double sumSqLeft=0, sumSqRight=0;
				for (int n=0; n<ids.length; n++) {
					if (input.get(ids[n], col) < t) {
						countLeft++;
						factorCountLeft[ output[ids[n]] ]++;
					} else {
						countRight++;
						factorCountRight[ output[ids[n]] ]++;
					}
				}
				// calculating score:
				double giniLeft  = getGiniIndex(factorCountLeft);
				double giniRight = getGiniIndex(factorCountRight);
				
				double score = (giniLeft*countLeft + giniRight*countRight) / ids.length;
				double scoreUSL = this.calculateUSL(unlabeled, unlabeledIds, col, t);
				score += scoreUSL;
				// pure node:
	//			if (score<zero*zero) {
	//				return makeLeaf(ids);
	//			}
				
				if (score < score_best) {
					score_best = score;
					col_best   = col;
					t_best     = t;
					leftConst  = (giniLeft<zero*zero);
					rightConst = (giniRight<zero*zero);
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
		FactorBinaryTree bt = new FactorBinaryTree();
		bt.column    = col_best;
		bt.threshold = t_best;
		bt.nSuccessors = ids.length;
		int[] uIdsLeft = null;
		int[] uIdsRight= null;
		// splitting unlabeled IDs:
		if (unlabeled!=null) {
			// finding ids:
			int[][] temp = AbstractTrees.splitIds(unlabeled, ids, bt.column, bt.threshold);
			uIdsLeft  = temp[0];
			uIdsRight = temp[1];
		}
		if (leftConst) { 
			bt.left = makeLeaf(idsLeft); // left child's output is constant 
		} else {  
			bt.left  = this.buildTree(nmin, K, idsLeft, uIdsLeft, randomCols); 
		}
		if (rightConst) {
			bt.right = makeLeaf(idsRight); // right child's output is constant
		} else {
			bt.right = this.buildTree(nmin, K, idsRight, uIdsRight, randomCols);
		}
		// this value is used only for CV:
		// TODO: add code for calculating value for intermediate nodes:
		//bt.value  = bt.left.value*bt.left.nSuccessors + bt.right.value*bt.right.nSuccessors;
		//bt.value /= bt.nSuccessors;
		return bt;
	}
	
	/**
	 * @param ids
	 * @return builds a leaf node and returns it with the given ids.
	 */
	public FactorBinaryTree makeLeaf(int[] ids) {
		// terminal node:
		FactorBinaryTree bt = new FactorBinaryTree();
		bt.value = 0;
		bt.nSuccessors = ids.length;
		// counting the factors:
		int[] counts = new int[nFactors];
		for (int n=0; n<ids.length; n++) {
			counts[ output[ids[n]] ]++;
		}
		bt.value = getMaxIndex(counts);
		return(bt);
	}
}
