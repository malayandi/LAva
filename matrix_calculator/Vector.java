package matrix_calculator;

import java.util.ArrayList;

/** A class representing the general Vector object.
 * 
 * @author Andrew Pau */
public class Vector {
    public Vector(double... args) {
        _values = args;
        _numRows = args.length;
        _normalization = new double[_numRows];
        _normalized = false;
    }
    
    /** Converts this vector into an nx1 Matrix.
     * 
     * @throws MatrixException */
    public Matrix matricize() throws MatrixException {
        double[][] contents = new double[_numRows][1];
        for (int row = 0; row < _numRows; row++) {
            contents[row][0] = _values[row];
        }
        Matrix matrix = new Matrix(_numRows, 1, contents);
        return matrix;
    }
    
    /** Returns an array representing the normalized vector. */
    public double[] normalize() {
        if (!_normalized) {
            int count = 0;
            for (double value : _values) {
                _normalization[count] = ((double) value / _magnitude);
                count++;
            }
            _normalized = true;
        }
        return _normalization;
    }
    
    /** Prints out this Matrix on the standard output. */
    public void print() {
        for (int row = 0; row < _numRows; row++) {
            System.out.print("[");
            System.out.print(_values[row]);
            System.out.print("]");
            System.out.println("");
        }
    }

    /** Returns the dot product of this vector and vector V. */
    public double dotProduct(Vector v) {
        double result = 0;
        for (int i = 0; i < numRows(); i++) {
            result += (get(i) * v.get(i));
        }
        return result;
    }

    /** Returns the number of entries in this vector. */
    public int numRows() {
        return _numRows;
    }

    /** Returns the magnitude of this vector. */
    public double magnitude() {
       _magnitude = 0;
       for (double value : _values) {
           _magnitude += value * value;
            }
        return _magnitude;
    }

    /** Returns the entries of this vector. */
    public double[] values() {
        return _values;
    }
    
    /** Returns the value at the specified index INDEX. */
    public double get(int index) {
        return _values[index];
    }

    /** Boolean indicating whether the vector has been normalized. */
    private boolean _normalized;
    /** The values of this vector. */
    private double[] _values;
    /** The number of entries in this vector. */
    private int _numRows;
    /** The magnitude of this vector. */
    private double _magnitude;
    /** The normalized entries of this vector. */
    private double[] _normalization;
}
