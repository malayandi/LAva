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
            if (_magnitude == null) {
                magnitude();
            }
            for (double value : _values) {
                if (Math.abs(value) < Matrix.epsilon) {
                    _normalization[count] = 0;
                } else {
                    _normalization[count] = ((double) value / _magnitude);
                }
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
            System.out.print(Matrix.df.format(_values[row]));
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
        _magnitude = (double) 0;
        for (double value : _values) {
            if (Math.abs(value) < Matrix.epsilon) {
                continue;
            }
            _magnitude += value * value;
        }
        _magnitude = Math.sqrt(_magnitude);
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

    /** Returns true if all the entries of this vector are 0. */
    public boolean isZero() {
        for (double v : _values) {
            if (v >= Matrix.epsilon) {
                return false;
            }
        }
        return true;
    }
    
    /** Returns true if this vector is equal to vector V. */
    public boolean equals(Vector v) {
        for (int i = 0; i < v._numRows; i++) {
            if (Math.abs(get(i) - v.get(i)) >= Matrix.epsilon) {
                return false;
            }
        }
        return true;
    }
    
    /** Returns true if this vector is equal to vector V, within a scale factor. */
    public boolean equalScale(Vector v) {
        if(v.values().length != values().length){
            System.out.println("The vector sizes do not match.");
            return false;
        }
        double factor = v.get(0) / get(0);
        for (int i = 1; i < numRows(); i++) {
            if (Math.abs(factor - (v.get(i) / get(i))) >= Matrix.epsilon) {
                return false;
            }
        }
        return true;
    }
    
    /** Scales every entry in this vector by K. */
    public void scale(double k) {
        for (int i = 0; i < _values.length; i++) {
            _values[i] = _values[i] * k;
        }
    }
    
    /** Scales this vector to contain only whole numbers, if possible.
     * If not possible, does nothing. */
    public void scaleWholeNum() {
        double[] copy = new double[_values.length];
        System.arraycopy(_values, 0, copy, 0, _values.length);
        Double factor = Math.abs(_values[0]);
        for (double v : _values) {
            if (Math.abs(v) < factor) {
                if (Math.abs((int) v - v) < Matrix.epsilon) {
                    continue;
                }
                factor = Math.abs(v);
            }
        }
        scale(1/factor);
        for (double v : _values) {
            if (Math.abs(v - (int) v) >= Matrix.epsilon) {
                _values = copy;
                break;
            }
        }
    }

    /** Boolean indicating whether the vector has been normalized. */
    private boolean _normalized;
    /** The values of this vector. */
    private double[] _values;
    /** The number of entries in this vector. */
    private int _numRows;
    /** The magnitude of this vector. */
    private Double _magnitude;
    /** The normalized entries of this vector. */
    private double[] _normalization;
}
