package org.extratrees;

import java.util.Date;

public class Timer {
	static long t1;
	static long t2;
	
	public static void tic() {
		t1 = new Date().getTime();
	}

	public static long toc(String function) {
		t2 = new Date().getTime();
		System.out.println(String.format("%s took %dms.", function, t2-t1));
		return t2 - t1;
	}

}
