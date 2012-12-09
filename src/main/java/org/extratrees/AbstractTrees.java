package org.extratrees;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

public abstract class AbstractTrees<E> {
	ArrayList<E> trees;

	/** number of threads */
	int numThreads = 1;

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

	public abstract E buildTree(int nmin, int k);


}
