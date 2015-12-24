package matrix_calculator;

import java.util.ArrayList;

/** A class representing the general Matrix object.
 *
 * @author AndyPalan */

public class Matrix {

    public Matrix(int row, int col, double[][] contents) {
        try {
            if (!(contents.length == row && contents[1].length == col)) {
                throw new MatrixException("Incorrect dimension.");
            }
            _dim = new ArrayList<Integer>();
            _dim.add(row);
            _dim.add(col);
            _contents = contents;
        } catch (MatrixException excp) {
            System.out.println(excp.getMessage());
        }
    }

    /** Returns the double at row R and col C. */
    public double get(int r, int c) {
        return _contents[r - 1][c - 1];
    }

    /** Sets the entry at row R and col C to be the double K. */
    public void set(int r, int c, double k) {
        _contents[r - 1][c - 1] = k;
    }

    /** Returns an ArrayList containing the dimension of the matrix, with the
     * height and the 0th index and the width at the 1st. */
    public ArrayList<Integer> getDimension() {
        return _dim;
    }

    /** Returns the integer height (number of rows) of the Matrix. */
    public int getHeight() {
        return _dim.get(0);
    }

    /** Returns the integer width (number of columns) of the Matrix. */
    public int getWidth() {
        return _dim.get(1);
    }

    /** Returns true if this Matrix is equal to Matrix A. */
    public Boolean equals(Matrix A) {
        for (int r = 1; r <= A.getHeight(); r++) {
            for (int c = 1; c <= A.getWidth(); c++) {
                if (get(r, c) != A.get(r, c)) {
                    return false;
                }
            }
        }
        return true;
    }

    /** Prints out this Matrix on the standard output. */
    public void print() {
        for (int r = 1; r <= getHeight(); r++) {
            System.out.print("[ ");
            for (int c = 1; c <= getWidth(); c++) {
                System.out.print(get(r, c) + " ");
            }
            System.out.print("]");
            System.out.println("");
        }
    }

    /** Sets _rowRed to be this Matrix in row reduced form if EF is false and
     * sets _rowRedEF to be this Matrix in row reduced echelon form if EF is
     * true.
     *
     * @throws MatrixException */
    public void rowReduction(Boolean EF) throws MatrixException { // Consider
        // used
        // _rowRed to
        // calculate
        // _rowRedEF if takes too long to
        // calculate _rowRedEF on its own
        Matrix B = Operations.matrixCopy(this);
        int pivot = 1;
        for (int c = 1; c <= B.getWidth(); c++) {
            if (B.count(c, pivot) == 0) {
                continue; // Non-pivot column
            } else if (B.get(c, pivot) == 0) {
                int k = c + 1;
                while (B.get(k, c) == 0) {
                    k++;
                }
                B.switchRow(c, k);
            }
            if (EF == true) {
                for (int r = 0; r <= B.getHeight(); r++) {
                    if (r == pivot) {
                        continue;
                    }
                    B.add(pivot, r, -1 * B.get(r, c) / B.get(pivot, c));
                }
            } else {
                for (int r = pivot + 1; r <= B.getHeight(); r++) {
                    B.add(pivot, r, -1 * B.get(r, c) / B.get(pivot, c));
                }
            }
            pivot++;
        }
        if (EF == true) {
            _rowRedEF = B;
        } else {
            _rowRed = B;
        }
        if (_rank == null) {
            _rank = pivot;
            _nullity = getWidth() - _rank;
            _linInd = (_rank == getWidth());
            _surjective = (_rank == getHeight());
            _injective = (_nullity == 0);
        }
    }

    /** Returns the row reduced form of this Matrix.
     *
     * @throws MatrixException */
    public Matrix getRowRed() throws MatrixException {
        if (_rowRed == null) {
            rowReduction(false);
        }
        return _rowRed;
    }

    /** Returns the row reduced echelon form of this Matrix.
     *
     * @throws MatrixException */
    public Matrix getRowRedEF() throws MatrixException {
        if (_rowRedEF == null) {
            rowReduction(true);
        }
        return _rowRedEF;
    }

    /** Returns the rank of this Matrix.
     *
     * @throws MatrixException */
    public int getRank() throws MatrixException {
        if (_rank == null) {
            getRowRed();
        }
        return _rank;
    }

    /** Returns the dimension of the Null Space of this Matrix.
     *
     * @throws MatrixException */
    public int getNullity() throws MatrixException {
        if (_nullity == null) {
            getRowRed();
        }
        return _nullity;
    }

    /** Returns true if the columns of this Matrix are linearly independent.
     *
     * @throws MatrixException */
    public Boolean isLinInd() throws MatrixException {
        if (_linInd == null) {
            getRowRed();
        }
        return _linInd;
    }

    /** Returns true if this matrix is surjective i.e. if columns of this Matrix
     * span R^n, where n is the dimension of the domain.
     *
     * @throws MatrixException */
    public Boolean isSurjective() throws MatrixException {
        if (_surjective == null) {
            getRowRed();
        }
        return _surjective;
    }

    /** Returns true if this matrix is injective i.e. if columns of this Matrix
     * span R^n, where n is the dimension of the domain of the Matrix.
     *
     * @throws MatrixException */
    public Boolean isInjective() throws MatrixException {
        if (_injective == null) {
            getRowRed();
        }
        return _injective;
    }

    /** Scalar multiplies row R of this Matrix by a constant K.
     *
     * @throws MatrixException */
    public void scalarMultRow(int r, double k) throws MatrixException {
        if (r < 1 || r > getHeight()) {
            throw new MatrixException("Row " + r + " is not a valid row.");
        }
        for (int c = 1; c <= getWidth(); c++) {
            set(r, c, get(r, c) * k);
        }
    }

    /** Scalar multiplies the entire matrix by a constant K. */
    public void scalarMult(double k) throws MatrixException {
        for (int r = 1; r <= getHeight(); r++) {
            scalarMultRow(r, k);
        }
    }

    /** Switches rows R1 and row R2 in this Matrix.
     *
     * @throws MatrixException */
    public void switchRow(int R1, int R2) throws MatrixException {
        if (R1 < 1 || R1 > getHeight()) {
            throw new MatrixException("Row " + R1 + " is not a valid row.");
        } else if (R2 < 1 || R2 > getHeight()) {
            throw new MatrixException("Row " + R2 + " is now a valid row.");
        }
        for (int c = 1; c <= getWidth(); c++) {
            double storeR1 = get(R1, c);
            set(R1, c, get(R2, c));
            set(R2, c, storeR1);
        }
    }

    /** Adds K * row R1 to row R2 of this Matrix.
     *
     * @throws MatrixException */
    public void add(int R1, int R2, double k) throws MatrixException {
        if (R1 < 1 || R1 > getHeight()) {
            throw new MatrixException("Row " + R1 + " is not a valid row.");
        } else if (R2 < 1 || R2 > getHeight()) {
            throw new MatrixException("Row " + R2 + " is now a valid row.");
        }
        for (int c = 1; c <= getWidth(); c++) {
            set(R2, c, get(R2, c) + get(R1, c) * k);
        }
    }

    /** Returns the number of non-zero entries in column C of this Matrix.
     *
     * @throws MatrixException */
    public int count(int c) throws MatrixException {
        return count(c, 1);
    }

    /** Returns the number of non-zero entries in column C of this Matrix,
     * starting at row R.
     *
     * @throws MatrixException */
    public int count(int c, int r) throws MatrixException {
        if (c < 1 || c > getWidth()) {
            throw new MatrixException("Column " + c + "is not a valid row.");
        } else if (r < 1 || r > getHeight()) {
            throw new MatrixException("Column " + c + "is not a valid row.");
        }
        int num = 0;
        for (int i = r; i <= getHeight(); i++) {
            if (get(i, c) != 0) {
                num++;
            }
        }
        return num;
    }

    /** The contents of this Matrix. */
    private double[][] _contents;

    /** The dimension of this Matrix. */
    private ArrayList<Integer> _dim;

    /** A Boolean that is true if the columns of this Matrix are linearly
     * independent. */
    private Boolean _linInd;

    /** A Boolean that is true if this Matrix is surjective, i.e. if columns of
     * this Matrix span R^n, where n is the dimension of the domain. */
    private Boolean _surjective;

    /** A Boolean that is true if this Matrix is injective. */
    private Boolean _injective;

    /** The row reduced form of this Matrix. */
    private Matrix _rowRed;

    /** The row reduced echelon form of this Matrix. */
    private Matrix _rowRedEF;

    /** The rank of this Matrix. */
    private Integer _rank;

    /** The dimension of the null space of this Matrix. */
    private Integer _nullity;
}
