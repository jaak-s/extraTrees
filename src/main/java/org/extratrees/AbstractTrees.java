package org.extratrees;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public abstract class AbstractTrees<E extends AbstractBinaryTree> {
	Matrix input;
	Random random = new Random();
	double[] weights;
	boolean useWeights;
	boolean hasNaN = false;
	int[] bagSizes = null;
	int[][] bagElems = null;
	protected final static double zero=1e-7;
	/** value to be returned when there is no answer, i.e., not available */
	protected final static double NA = Double.NaN;


	/** for multi-task learning, stores task indeces (null if not present) */
	int[] tasks;
	
	/** number of tasks: (tasks are indexed from 0 to (nTasks-1) */
	int nTasks;

	ArrayList<E> trees;

	/** number of threads */
	int numThreads = 1;
	
	/** number of times task cut is tried */
	int numRandomTaskCuts = 1;
	
	/** probability of trying task cutting (lambda) */
	double probOfTaskCuts = 1.0;

	/** number of random cuts tried for each feature */
	int numRandomCuts = 1;
	/** whether random cuts are totally uniform or evenly uniform */
	boolean evenCuts = false;

	/** later shuffled and used for choosing random columns at each node: */
	ArrayList<Integer> cols;
	
	public boolean isHasNaN() {
		return hasNaN;
	}
	
	public void setHasNaN(boolean hasNaN) {
		this.hasNaN = hasNaN;
	}

	public int getNumThreads() {
		return numThreads;
	}
	
	public void setNumThreads(int numThreads) {
		this.numThreads = numThreads;
	}
	
	public int getNumRandomTaskCuts() {
		return numRandomTaskCuts;
	}
	/** sets the number of times random cut for tasks is tried (assuming multitask setting).  
	 *  default is 1.
	 * */
	public void setNumRandomTaskCuts(int numRandomTaskCuts) {
		this.numRandomTaskCuts = numRandomTaskCuts;
	}
	
	public double getProbOfTaskCuts() {
		return probOfTaskCuts;
	}
	/**
	 * sets the probability of trying task cutting.
	 * @param probOfTaskCuts
	 */
	public void setProbOfTaskCuts(double probOfTaskCuts) {
		this.probOfTaskCuts = probOfTaskCuts;
	}
	
	/**
	 * @return number of trees used in ExtraTrees
	 */
	public int getNumTrees() {
		return trees.size();
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
	
	public void setInput(Matrix input) {
		this.input = input;
		// making cols list for later use:
		this.cols = new ArrayList<Integer>(input.ncols);
		for (int i=0; i<input.ncols; i++) {
			cols.add(i);
		}
	}
	
	public void setTasks(int[] tasks) {
		// making a list of tasks:
		this.tasks  = tasks; 
		this.nTasks = 1;
		if (this.tasks!=null) {
			for (int i=0; i<tasks.length; i++) {
				//taskNames.add(tasks[i]);
				if (nTasks < tasks[i] + 1) {
					nTasks = tasks[i] + 1;
				}
			}
		}
	}
	
	public void setWeights(double[] weights) {
		if (weights != null && input.nrows != weights.length) {
			throw(new IllegalArgumentException("Input and weights do not have the same number of data points."));
		}
		this.weights = weights;
		this.useWeights = (weights!=null);
	}
	
	/**
	 * Sets bag size to bagSize, so each tree is only built with bagSize (randomly selected) samples.
	 * @param bagSize
	 */
	public void setBagging(int bagSize) {
		if (bagSize > input.nrows) {
			throw( new IllegalArgumentException("Supplied bagSize exceeds the number of samples.") );
		}
		this.bagSizes = new int[]{ bagSize };
	}

	/**
	 * Sets bag sizes for each bag label group. 
	 * @param bagSizes   int[] size of the bag for each label
	 * @param bagLabels  int[] bag label for each sample, all from 0 to (Nbags - 1)
	 */
	public void setBagging(int[] bagSizes, int[] bagLabels) {
		if (bagLabels.length != input.nrows) {
			throw( new IllegalArgumentException("size of bagLabels has to equal to the number of input rows.") );
		}
		this.bagSizes = bagSizes;
		int[] counts = new int[bagSizes.length];
		for (int i=0; i<input.nrows; i++) {
			counts[ bagLabels[i] ]++;
		}
		// making sure all bags have enough elements:
		this.bagElems = new int[counts.length][];
		for (int bag=0; bag < counts.length; bag++) {
			if (counts[bag] < bagSizes[bag]) {
				throw( new IllegalArgumentException( 
						String.format("Bag %d has less elements (%d) than requested by bag size (%d).",
								bag, counts[bag], bagSizes[bag]
						)));
			}
			this.bagElems[bag] = new int[ counts[bag] ];
		}
		// adding elements to the appropriate bags:
		for (int i=0; i < input.nrows; i++) {
			int bag = bagLabels[i];
			int j = this.bagElems[bag].length - counts[bag];
			counts[bag]--;
			this.bagElems[bag][j] = i;
		}

	}

	/**
	 * @param col_min
	 * @param diff
	 * @param repeat  only used when evenCuts==true.
	 * @return random cut from col_min to col_min+diff.
	 */
	protected double getRandomCut(double col_min, double diff, int repeat) {
		double t;
		if (evenCuts) {
			double iStart = col_min + repeat*diff/numRandomCuts;
			double iStop  = col_min + (repeat+1)*diff/numRandomCuts;
			t = getRandom()*(iStop-iStart) + iStart;
		} else {
			t = getRandom()*diff + col_min;
		}
		return t;
	}
	
	/**
	 * @param ids
	 * @param col
	 * @param input
	 * @return array of min and max values.
	 */
	protected static double[] getRange(int[] ids, int col, Matrix input) {
		double[] range = new double[2];
		range[0] = Double.POSITIVE_INFINITY;
		range[1] = Double.NEGATIVE_INFINITY;
		for (int n=0; n<ids.length; n++) {
			double v = input.get(ids[n], col);
			if ( v<range[0] ) { range[0] = v; }
			if ( v>range[1] ) { range[1] = v; }
		}
		return range;
	}

	/**
	 * @param input
	 * @return array of size two: min and max values.
	 */
	protected static double[] getRange(double[] input) {
		double[] range = new double[2];
		range[0] = Double.POSITIVE_INFINITY;
		range[1] = Double.NEGATIVE_INFINITY;
		for (int n=0; n<input.length; n++) {
			if ( input[n]<range[0] ) { range[0] = input[n]; }
			if ( input[n]>range[1] ) { range[1] = input[n]; }
		}
		return range;
	}
	
	static protected class CutResult {
		double score;
		boolean leftConst;
		boolean rightConst;
		int countLeft;
		int countRight;
		double nanWeigth;
		
		public CutResult() {}
		
		public int getCountLeft() {
			return countLeft;
		}
		public void setCountLeft(int countLeft) {
			this.countLeft = countLeft;
		}
		public int getCountRight() {
			return countRight;
		}
		public void setCountRight(int countRight) {
			this.countRight = countRight;
		}
		
		public boolean isLeftConstant() {
			return leftConst;
		}
		public void setLeftConstant(boolean leftConstant) {
			this.leftConst = leftConstant;
		}
		
		public boolean isRightConstant() {
			return rightConst;
		}
		public void setRightConstant(boolean rightConstant) {
			this.rightConst = rightConstant;
		}
		
		public double getScore() {
			return score;
		}
		public void setScore(double score) {
			this.score = score;
		}
		
	}
	
	abstract public E makeLeaf(int[] ids, Set<Integer> leftTaskSet);
	
	/**
	 * Same as buildTrees() except computes in parallel.
	 * @param nmin
	 * @param K
	 * @param nTrees
	 * @return
	 */
	public ArrayList<E> buildTreesParallel(int nmin, int K, int nTrees) {
		// creating a thread pool and using it to compute nTrees:
		ExecutorService executor = Executors.newFixedThreadPool(numThreads);
		List<TreeCallable> callables = new ArrayList<TreeCallable>();
		
		// adding tasks: using the same task for each tree
		TreeCallable task = new TreeCallable(nmin, K);
		for (int i=0; i<nTrees; i++) {
			callables.add( task );
		}
		// computing and fetching results:
		List<Future<E>> results;
		try {
			results = executor.invokeAll(callables);
		} catch (InterruptedException e) {
			// not solving this error here:
			throw new RuntimeException(e);
		}
		// fetching all BinaryTrees and storing them:
		ArrayList<E> trees = new ArrayList<E>(nTrees);
		for (Future<E> f : results) {
			try {
				trees.add( f.get() );
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			} catch (ExecutionException e) {
				throw new RuntimeException(e);
			}
		}
		executor.shutdown();
		return trees;
	}
	
	
	/** Nested class for making BinaryTrees */
	public class TreeCallable implements Callable<E> {
		int nmin, K;
		public TreeCallable(int nmin, int K) {
			this.nmin = nmin;
			this.K    = K;
		}
		@Override
		public E call() throws Exception {
			return AbstractTrees.this.buildTree(nmin, K);
		}
	}

	/**
 	 * good values:
	 * n_min = 2 (size of tree element)
	 * K = 5     (# of random choices)
	 * M = 50    (# of trees)
	 * if n_min is chosen by CV, then we have pruned version

	 * @param nmin   - size of tree element
	 * @param K      - # of random choices
	 * @param nTrees - # of trees
	 * Single threaded computation.
	 * @return learned trees
	 */
	public ArrayList<E> buildTrees(int nmin, int K, int nTrees) {
		ArrayList<E> trees = new ArrayList<E>(nTrees);
		// single-threading:
		for (int t=0; t<nTrees; t++) {
			trees.add( this.buildTree(nmin, K) );
		}
		return trees;
	}

	/**
	 * stores trees with the AbstractTrees object. 
	 * Uses multiple threads if set.
	 * @param nmin
	 * @param K
	 * @param nTrees
	 */
	public void learnTrees(int nmin, int K, int nTrees) {
		if (numThreads<=1) {	//if (numThreads<=1) {
			this.trees = buildTrees(nmin, K, nTrees);
		} else {
			this.trees = buildTreesParallel(nmin, K, nTrees);
		}
	}
	
	/**
	 * @param n
	 * @return int array of [0, 1, ..., n-1].
	 */
	public static int[] seq(int n) {
		int[] seq = new int[n];
		for (int i=0; i < n; i++) {
			seq[i] = i;
		}
		return seq;
	}

	/**
	 * @param nStart
	 * @param nEnd
	 * @return int array of [nStart, nStart+1, ..., nEnd-1].
	 */
	public static int[] seq(int nStart, int nEnd) {
		int[] seq = new int[nEnd - nStart];
		for (int i=nStart; i < nEnd; i++) {
			seq[i-nStart] = i;
		}
		return seq;
	}

	/**
	 * Main method that performs tree training.
	 * @param nmin
	 * @param K
	 * @return
	 */
	public E buildTree(int nmin, int K) {
		int[] ids = getInitialSamples();
		ShuffledIterator<Integer> cols = new ShuffledIterator<Integer>(this.cols);
		
		// finding task set:
		HashSet<Integer> taskSet = getSequenceSet(nTasks);
		return buildTree(nmin, K, ids, cols, taskSet );
	}
	
	protected static ArrayList<Integer> arrayToList(int[] array) {
		ArrayList<Integer> list = new ArrayList<Integer>(array.length);
		for (int value : array) {
			list.add(value);
		}
		return list;
	}

	protected static int[] listToArray(ArrayList<Integer> list) {
		int[] a = new int[ list.size() ];
		for (int i=0; i < a.length; i++) {
			a[i] = list.get(i);
		}
		return a;
	}

	protected int[] getInitialSamples() {
		if (bagSizes == null) {
			return seq( input.nrows );
		}
		if (bagSizes.length == 1) {
			ArrayList<Integer> allIds = arrayToList( seq( input.nrows ) );
			ShuffledIterator<Integer> shuffle = new ShuffledIterator<Integer>(allIds);
			
			int[] bag = new int[ bagSizes[0] ];
			for (int i=0; i < bag.length; i++) {
				bag[i] = shuffle.next();
			}
			return bag;
		}
		// selecting random samples from each bag:
		int[] bag = new int[ sum(bagSizes) ];
		int i = 0;
		for (int b=0; b < bagSizes.length; b++) {
			ArrayList<Integer> ids = arrayToList( bagElems[b] );
			ShuffledIterator<Integer> shuffle = new ShuffledIterator<Integer>(ids);

			// filling with elements from bag[b]:
			for (int n = i + bagSizes[b]; i < n; i++) {
				bag[i] = shuffle.next();
			}
			
		}
		
		return bag; 
	}

	/** @return penalty for NaN values for the given ids */
	protected abstract double get1NaNScore(int[] ids);
	
	protected abstract void calculateCutScore(int[] ids, int col, double t, CutResult result);
	
	/**
	 * @param ids
	 * @param bestScore
	 * @return array of booleans specifying the cut (true for left tree, false for left).
	 *         Returns null if no cut better (smaller score) than bestScore was found.
	 */
	protected abstract TaskCutResult getTaskCut(int[] ids, 
			Set<Integer> tasks, double bestScore);


	/**
	 * @param m    data matrix
	 * @param ids  row numbers to be split
	 * @param dim  which dimension to use for splitting
	 * @param cut  splitting value
	 * @return filters ids <b>ids</b> into two arrays, one whose values in data values
	 * {@code m[id, dim]} are below cut and others whose values are higher. 
	 */
	public static int[][] splitIds(Matrix m, int[] ids, int dim, double cut) {
		int[][] out = new int[2][];
		int lenLower = 0;
		for (int i=0; i<ids.length; i++) {
			if (m.get(ids[i], dim)<cut) {
				lenLower++;
			}
		}
		// two vectors: lower and higher
		out[0] = new int[lenLower];
		out[1] = new int[ids.length - lenLower];
		
		int i0 = 0;
		int i1 = 0;
		for (int i=0; i<ids.length; i++) {
			if (m.get(ids[i], dim)<cut) {
				out[0][i0] = ids[i];
				i0++;
			} else {
				out[1][i1] = ids[i];
				i1++;
			}
		}
		
		return out;
	}
	
	/**
	 * 
	 * @param ids
	 * @param leftTasks
	 * @return filters <b>ids</b> into two arrays, one whose tasks[i] is is 
	 * in leftTasks and others whose values are not.
	 */
	public int[][] splitIdsByTask(int[] ids, Set<Integer> leftTasks) {
		int[][] out = new int[2][];
		int lenLower = 0;
		for (int i=0; i<ids.length; i++) {
			if ( leftTasks.contains(tasks[ids[i]]) ) {
				lenLower++;
			}
		}
		// two vectors: lower and higher
		out[0] = new int[lenLower];
		out[1] = new int[ids.length - lenLower];
		
		int i0 = 0;
		int i1 = 0;
		
		for (int i=0; i<ids.length; i++) {
			if ( leftTasks.contains(tasks[ids[i]]) ) {
				out[0][i0] = ids[i];
				i0++;
			} else {
				out[1][i1] = ids[i];
				i1++;
			}
		}
		return out;
	}
	
	public static int sum(int[] array) {
		int s = 0;
		for (int i=0; i<array.length; i++) {
			s += array[i];
		}
		return s;
	}

	public static double sum(double[] array) {
		double s = 0;
		for (int i=0; i<array.length; i++) {
			s += array[i];
		}
		return s;
	}

	public static HashSet<Integer> getSequenceSet(int n) {
		HashSet<Integer> taskSet;
		taskSet = new HashSet<Integer>(n);
		for (int i=0; i<n; i++) {
			taskSet.add(i);
		}
		return taskSet;
	}

	/**
	 * @return random between 0.0 and 1.0.
	 */
	protected double getRandom() {
		return random.nextDouble();
	}
	
	public void setSeed(long seed) {
		this.random = new Random(seed);
	}
	
	/**
	 * @param xmin
	 * @param xmax
	 * @return uniformly random value between xmin and xmax.
	 */
	protected double getRandom(double xmin, double xmax) {
		return xmin + getRandom()*(xmax - xmin);
	}

	abstract protected E makeFilledTree(E leftTree, E rightTree, 
			int col_best, double t_best,
			int nSuccessors);
	
	/**
	 * 
	 * @param nmin
	 * @param K
	 * @param ids
	 * @param randomCols - passed to save memory (maybe not needed)
	 * @return
	 */
	public E buildTree(int nmin, int K, int[] ids, 
			ShuffledIterator<Integer> randomCols, Set<Integer> taskSet) {
		if (ids.length<nmin) {
			return makeLeaf(ids, taskSet);
		}
		// doing a shuffle of cols:
		randomCols.reset();
		
		// trying K trees or the number of non-constant columns,
		// whichever is smaller:
		int k = 0, col_best=-1;
		double t_best=Double.NaN;
		double nanPenalty = hasNaN ?get1NaNScore(ids) :0;
		CutResult bestResult = new CutResult();
		bestResult.score = Double.POSITIVE_INFINITY;
		
		loopThroughColumns:
		while( randomCols.hasNext() ) {
			int col = randomCols.next();
			// calculating columns min and max:
			double[] range = getRange(ids, col, input);
			if (range[1]-range[0] < zero) {
				// skipping, because column is constant
				continue;
			}
			// picking random test point numRepeatTries:
			double diff = (range[1]-range[0]);
			for (int repeat=0; repeat<this.numRandomCuts; repeat++) {
				double t;
				t = getRandomCut(range[0], diff, repeat);
				
				CutResult result = new CutResult();
				calculateCutScore(ids, col, t, result);
				if (hasNaN) {
					result.score += result.nanWeigth * nanPenalty;
				}
				
				if (result.score < bestResult.score) {
					// found a better scoring cut
					col_best   = col;
					t_best     = t;
					
					bestResult.score = result.score;
					bestResult.leftConst  = result.leftConst;
					bestResult.rightConst = result.rightConst;
					bestResult.countLeft  = result.countLeft;
					bestResult.countRight = result.countRight;
					if (bestResult.leftConst && bestResult.rightConst && ! hasNaN) {
						// the result cannot be improved:
						break loopThroughColumns;
					}
				}
				
			}

			k++;
			if (k>=K) {
				// checked enough columns, stopping:
				break;
			}
		}
		
		// multi-task learning:
		TaskCutResult taskCutResult = null;
		if (taskSet.size() > 1) {
			// checking whether to perform task splitting
			if ( probOfTaskCuts > getRandom() ) {
				// calculating task order:
				taskCutResult = getTaskCut(ids, taskSet, bestResult.score);
			}
		}
		
		// outputting the tree using the best score cut:
		int[] idsLeft  = new int[bestResult.countLeft];
		int[] idsRight = new int[bestResult.countRight];
		
		if (col_best < 0 && taskCutResult==null) {
			// no feature or task split found
			return makeLeaf(ids, taskSet);
		}
		
		int[][] split;
		Set<Integer> leftTaskSet, rightTaskSet;
		if (taskCutResult!=null) {
			// task cut:
			split = splitIdsByTask(ids, taskCutResult.leftTasks);
			leftTaskSet  = taskCutResult.leftTasks;
			rightTaskSet = taskCutResult.rightTasks;
			col_best = -1;
			t_best = Double.NaN;
		} else {
			// splitting according to feature:
			split = splitIds(input, ids, col_best, t_best);
			leftTaskSet  = taskSet;
			rightTaskSet = taskSet;
		}
		idsLeft  = split[0];
		idsRight = split[1];
		
		E leftTree, rightTree;
		if (bestResult.leftConst) {
			 // left child's output is constant
			leftTree = makeLeaf(idsLeft, leftTaskSet); 
		} else {
			leftTree  = this.buildTree(nmin, K, idsLeft, randomCols, leftTaskSet); 
		}
		if (bestResult.rightConst) {
			// right child's output is constant
			rightTree = makeLeaf(idsRight, rightTaskSet);
		} else {
			rightTree = this.buildTree(nmin, K, idsRight, randomCols, rightTaskSet);
		}
		
		E bt = makeFilledTree(leftTree, rightTree, 
				col_best, t_best, ids.length);
		bt.tasks = taskSet;
		return bt;
	}
	


}
