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
        double[] entries = _vectors.get(count).values();
        double[] entries2 = new double[entries.length];
        System.arraycopy(entries, 0, entries2, 0, entries.length);
        Vector u = new Vector(entries2);
        for (int counter = 0; counter < count; counter++) {
            Vector next = initial.get(counter);
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
    
    /** Returns a Matrix containing the vectors of this set as its columns. 
     * @throws MatrixException */
    public Matrix toMatrix() throws MatrixException {
        int h = _vectors.get(0).numRows();
        int w = _vectors.size();
        Matrix matrix = new Matrix(h, w, new double[h][w]);
        for (int c = 0; c < w; c++) {
            Vector v = get(c);
            for (int r = 0; r < h; r++) {
                matrix.set(r + 1, c + 1, v.get(r));
            }
        }
        return matrix;
    }
    
    /** Returns a Square Matrix containing the vectors of this set as its columns. 
     * @throws MatrixException */
    public SquareMatrix toSquareMatrix() throws MatrixException {
        int h = _vectors.get(0).numRows();
        int w = _vectors.size();
        SquareMatrix matrix = new SquareMatrix(h, new double[h][w]);
        for (int c = 0; c < w; c++) {
            Vector v = get(c);
            for (int r = 0; r < h; r++) {
                matrix.set(r + 1, c + 1, v.get(r));
            }
        }
        return matrix;
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
