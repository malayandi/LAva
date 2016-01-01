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
    

    /** Converts this vector set into a matrix.
     * 
     * @throws MatrixException */
    public Matrix matricize() throws MatrixException {
        double[][] contents = new double[_vectors.get(0).numRows()][size()];
        int count = 0;
        for (Vector vector : _vectors) {
            for (int row = 0; row < vector.numRows(); row++) {
                contents[row][count] = vector.values()[row];
            }
            count++;
        }
        Matrix matrix;
        if (_vectors.get(0).numRows() == size()) {
            matrix = new SquareMatrix(size(), contents);
        } else {
            matrix = new Matrix(_vectors.get(0).numRows(), size(), contents);
        }
        return matrix;
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
            coefficient = coefficient / (next.magnitude() * next.magnitude());
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
    
    /** Returns a VectorSet containing orthonormalised vectors from this set. */
    public VectorSet orthonormalise() {
        Vector[] result = new Vector[size()];
        for (int i = 0; i < size(); i++) {
            Vector r = new Vector(_vectors.get(i).normalize());
            result[i] = r;
        }
        return new VectorSet(result);
    }
    
    /** Returns a Matrix containing the vectors of this set as its columns,
     * sending all zero vectors to the end of this matrix. 
     * @throws MatrixException */
    public Matrix toMatrix() throws MatrixException {
        ArrayList<Vector> vectors = new ArrayList<>();
        for (Vector v : _vectors) {
            vectors.add(v);
        }
        for (Vector v : _vectors) {
            if (v.isZero()) {
                vectors.remove(v);
                vectors.add(v);
            }
        }
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
    
    /** Returns a Square Matrix containing the vectors of this set as its columns,
     * sending all zero vectors to the end of this matrix. 
     * @throws MatrixException */
    public SquareMatrix toSquareMatrix() throws MatrixException {
        ArrayList<Vector> vectors = new ArrayList<>();
        for (Vector v : _vectors) {
            vectors.add(v);
        }
        for (Vector v : _vectors) {
            if (v.isZero()) {
                vectors.remove(v);
                vectors.add(v);
            }
        }
        int h = _vectors.get(0).numRows();
        SquareMatrix matrix = new SquareMatrix(h, new double[h][h]);
        for (int c = 0; c < h; c++) {
            Vector v = get(c);
            for (int r = 0; r < h; r++) {
                matrix.set(r + 1, c + 1, v.get(r));
            }
        }
        return matrix;
    }
    
    /** Returns the coordinates of x given coordinate-vector VECTOR and this
     * vector set. */
    public Vector coordinates(Vector coordinateVector) {
        // assert that this vector set is a basis
        double[] values = new double[coordinateVector.numRows()];
        int count = 0;
        for (Vector b : _vectors) {
            for (int counter = 0; counter < coordinateVector.numRows(); counter++) {
                values[counter] += coordinateVector.values()[count] * b.values()[counter];
                System.out.println(values[counter]);
            }
            count++;
        }
        Vector coordinates = new Vector(values);
        return coordinates;
    }

    /** Returns the coordinate-vector of Vector VECTOR with respect to this
     * vector set.
     * 
     * @throws MatrixException */
    public Matrix coordinateVector(Vector vector) throws MatrixException {
        SquareMatrix matrix = (SquareMatrix) matricize();
        SquareMatrix inverse = matrix.getInverse();
        Matrix coordinateVector = vector.matricize();
        coordinateVector = Operations.matrixMult(inverse, coordinateVector);
        return coordinateVector;
    }

    /** If this vector set is a basis w/ respect to c, given coordinateVector
     * [x]b, return [x]c. */
    public Vector cCoordinate(Vector coordinateVector) {
        return coordinates(coordinateVector);
    }

    /** Returns the change of basis matrix from this VectorSet to another
     * VectorSet SET. */
    public Matrix changeOfBasis(VectorSet set) throws MatrixException {
        Matrix B = matricize();
        Matrix C = set.matricize();
        Matrix changeOfBasis = B.changeOfBasis(C);
        return changeOfBasis;
    }

    /** Prints this VectorSet in the standard output. */
    public void print() {
        if (size() == 0) {
            System.out.println("This set is empty");
        } else {
            for (int index = 0; index < _vectors.get(0).numRows(); index++) {
                for (Vector vector : _vectors) {
                    System.out.print("[ " + Matrix.df.format(vector.get(index)) + "]");
                }
                System.out.println("");
            }
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
