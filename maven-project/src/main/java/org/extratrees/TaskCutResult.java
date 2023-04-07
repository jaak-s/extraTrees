package org.extratrees;

import java.util.Set;

import org.extratrees.AbstractTrees.CutResult;

public class TaskCutResult extends CutResult {
	Set<Integer> leftTasks;
	Set<Integer> rightTasks;
}
