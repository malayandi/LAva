package matrix_calculator;

import java.util.ArrayList;

/** A class representing a set of Vectors.
 * 
 * @author Andrew Pau */
public class VectorSet {
    public VectorSet(Vector... args) {
        _vectors = new ArrayList<>();
        for (Vector vector : args) {
            _vectors.add(vector);
        }
    }

    /** Performs the Gram-Schmidt process on this vector set. */
    public VectorSet gramSchmidt(int count, VectorSet initial) {
        Vector u = _vectors.get(count);
        for (int counter = 0; counter < count; counter++) {
            Vector next = _vectors.get(counter);
            double coefficient = 0;
            for (int index = 0; index < u.values().length; index++) {
                coefficient += u.values()[index] * next.values()[index];
            }
            coefficient = coefficient / next.magnitude();
            for (int index = 0; index < next.numRows(); index++) {
                u.values()[index] -= coefficient * next.values()[index];
            }
        }
        if (count == 0) {
            initial = new VectorSet();
        }
        initial.add(u);
        if (initial.size() != _vectors.size()) {
            return gramSchmidt(count + 1, initial);
        } else {
            return initial;
        }
    }

    /** Returns an ArrayList containing the vectors within this set. */
    public int size() {
        return _vectors.size();
    }

    /** Adds a vector to this vector set. */
    public void add(Vector vector) {
        _vectors.add(vector);
    }

    /** Gets a vector from this vector set. */
    public Vector get(int index) {
        return _vectors.get(index);
    }

    /** An ArrayList containing all vectors within this set. */
    private ArrayList<Vector> _vectors;
}
