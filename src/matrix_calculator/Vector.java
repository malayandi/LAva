package matrix_calculator;

/** A class representing the general Vector object.
 * 
 * @author Andrew Pau */
public class Vector {
    public Vector(double... args) {
        _values = args;
        _numRows = args.length;
        _normalization = new double[_numRows];
        double counter = 0;
        for (double value : _values) {
            counter += value * value;
        }
        _magnitude = counter;
        _normalized = false;
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

    /** Returns the number of entries in this vector. */
    public int numRows() {
        return _numRows;
    }

    /** Returns the magnitude of this vector. */
    public double magnitude() {
        return _magnitude;
    }

    /** Returns the entries of this vector. */
    public double[] values() {
        return _values;
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
