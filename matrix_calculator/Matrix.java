package matrix_calculator;

import java.util.ArrayList;

/** A class representing the general Matrix object.
 * 
 * @author AndyPalan
 *
 */

public class Matrix {

    public Matrix(int row, int col, Double[][] contents) {
        _dim.add(row);
        _dim.add(col);
        _contents = contents;
    }

    /** Returns the double at row R and col C. */
    public double get(int r, int c) {
        _contents[r - 1][c - 1];
    }
    
    /** Sets the entry at row R and col C to be the double K. */
    public void set(int r, int c, double k) {
        _contents[r - 1][c - 1] = k;
    }
    
    /** Returns an ArrayList containing the dimension of the matrix, with
     * the height and the 0th index and the width at the 1st. */
    public ArrayList<Integer> getDimension() {
        return _dim;
    }

    /** Returns the integer height (number of rows) of the Matrix. */
    public Integer getHeight() {
        return _dim.get(0);
    }

    /** Returns the integer width (number of columns) of the Matrix. */
    public Integer getWidth() {
        return _dim.get(1);
    }

    /** The contents of this Matrix. */
    private Double[][] _contents;
    
    /** The dimension of this Matrix. */
    private ArrayList<Integer> _dim;

    /** A boolean that returns true if the columns of this Matrix are linearly
     * independent. */
    private boolean _linInd;
    
    /** The row reduced form of this Matrix. */
    private Matrix _rowRed;
    
    /** The row reduced echelon form of this Matrix. */
    private Matrix _rowRedEF;
    
    /** The rank of this Matrix. */
    private Integer _rank;
    
    /** The dimension of the null space of this Matrix. */
    private Integer _nullity;
    
    /** A basis for the range of this Matrix. */
    private VectorSet _basisNull;
    
    /** A basis for the null space of this Matrix. */
    private VectorSet _basisRange;
}
