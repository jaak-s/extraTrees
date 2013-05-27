package org.extratrees;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import org.extratrees.AbstractTrees.CutResult;

public abstract class AbstractTrees<E extends AbstractBinaryTree> {
	Matrix input;
	protected final static double zero=1e-7;

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


	// for semi-supervised learning:
	/** unsupervised points, used for semi-supervised learning */
	//Matrix unlabeled;
	/** weight of unsupervised learning (used to make semi-supervised signal) */
	//double weightOfUSL = 1.0;
	
//	HardMargin margin = new HardMargin();
	
	/*
	public void setUnlabeled(Matrix unlabeled) {
		this.unlabeled = unlabeled;
	}
	
	public Matrix getUnlabeled() {
		return unlabeled;
	}*/
	
//	public void setWeightOfUSL(double weightOfUSL) {
//		this.weightOfUSL = weightOfUSL;
//	}
//	
//	public double getWeightOfUSL() {
//		return weightOfUSL;
//	}

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
	
	public void setTasks(int[] tasks) {
		// making a list of tasks:
		this.tasks  = tasks; 
		this.nTasks = 1;
		if (this.tasks!=null) {
			for (int i=0; i<tasks.length; i++) {
				//taskNames.add(tasks[i]);
				nTasks = tasks[i] + 1;
			}
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
	protected double[] getRange(int[] ids, int col, Matrix input) {
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
	protected double[] getRange(double[] input) {
		double[] range = new double[2];
		range[0] = Double.POSITIVE_INFINITY;
		range[1] = Double.NEGATIVE_INFINITY;
		for (int n=0; n<input.length; n++) {
			if ( input[n]<range[0] ) { range[0] = input[n]; }
			if ( input[n]>range[1] ) { range[1] = input[n]; }
		}
		return range;
	}

	
	/**
	 * @param inputs
	 * @param ids
	 * @param dim
	 * @param cut
	 * @return score of unsupervised error (the bigger the worse), i.e., 
	 * how far the margin is from the best margin (0.5).
	 * Gives 0 if unlabeled is null. 
	 * Uses <b>weightOfSSL</b> to adjust the result.
	 */
	/*
	public double calculateUSL(Matrix inputs, int[] ids, int dim, double cut) {
		if (unlabeled==null) {
			return 0;
		}
		// best margin is 0.5
		return weightOfUSL * ( 0.5 - margin.getCriteria(inputs, ids, dim, cut) );
	}
	*/
	
	static protected class CutResult {
		double score;
		boolean leftConst;
		boolean rightConst;
		int countLeft;
		int countRight;
		Object value;
		
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
	

	public E buildTree(int nmin, int K) {
		int[]    ids = new int[input.nrows];
		for (int i=0; i<ids.length; i++) {
			ids[i] = i;
		}
		ShuffledIterator<Integer> cols = new ShuffledIterator<Integer>(this.cols);
		
		// finding task set:
		HashSet<Integer> taskSet = getSequenceSet(nTasks);
		return buildTree(nmin, K, ids, cols, taskSet );
	}

	protected abstract void calculateCutScore(int[] ids, int col, double t, CutResult result);
	//protected abstract double[] getTaskScores(int[] ids);
	/**
	 * @param ids
	 * @param bestScore
	 * @return array of booleans specifying the cut (true for left tree, false for left).
	 *         Returns null if no cut better (smaller score) than bestScore was found.
	 */
	protected abstract TaskCutResult getTaskCut(int[] ids, 
			HashSet<Integer> tasks, double bestScore);


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
		return Math.random();
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
		CutResult bestResult = new CutResult();
		bestResult.score = Double.POSITIVE_INFINITY;
		
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
				
				if (result.score < bestResult.score) {
					// found a better scoring cut
					col_best   = col;
					t_best     = t;
					
					bestResult.score = result.score;
					bestResult.leftConst  = result.leftConst;
					bestResult.rightConst = result.rightConst;
					bestResult.countLeft  = result.countLeft;
					bestResult.countRight = result.countRight;
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
				taskCutResult = getTaskCut(ids, null, bestResult.score);
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
