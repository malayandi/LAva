package matrix_calculator;

import java.util.ArrayList;

public class Vector {
    public Vector(int... args) {
        values = args;
        numRows = args.length;
        normalized = new ArrayList<>();
        int counter = 0;
        for (int value : values) {
            counter += value * value;
        }
        magnitude = counter;
    }

    public ArrayList<Double> normalize() {
        if (!normalized.isEmpty()) {
            return normalized;
        }
        for (int value : values) {
            normalized.add((double) value / magnitude);
        }
        return normalized;
    }

    public int numRows() {
        return numRows;
    }

    public int magnitude() {
        return magnitude;
    }

    private int[] values;
    private int numRows;
    private int magnitude;
    private ArrayList<Double> normalized;
}
