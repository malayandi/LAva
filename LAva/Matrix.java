package matrix_calculator;

import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;

/** A class representing the general Matrix object.
 *
 * @author AndyPalan */

public class Matrix {

    /** Creates a new ROW x COL square Matrix with contents CONTENTS. */
    public Matrix(int row, int col, double[][] contents) throws MatrixException {
        if (!(contents.length == row && contents[1].length == col)) {
            throw new MatrixException("Incorrect dimension.");
        }
        _dim = new ArrayList<Integer>();
        _dim.add(row);
        _dim.add(col);
        _contents = contents;
        df.setRoundingMode(RoundingMode.HALF_UP);
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
                if (Math.abs(get(r, c) - A.get(r, c)) >= epsilon) {
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
                double entry = get (r, c) + 0.0;
                System.out.print(df.format(entry) + " ");
            }
            System.out.print("]");
            System.out.println("");
        }
    }

    /** Sets _transpose to be the transpose of this Matrix. */
    public void transpose() throws MatrixException {
        Matrix T = new Matrix(getWidth(), getHeight(), new double[getWidth()][getHeight()]);
        for (int r = 1; r <= getHeight(); r++) {
            for (int c = 1; c <= getWidth(); c++) {
                T.set(c, r, get(r, c));
            }
        }
        _transpose = T;
    }

    /** Sets _rowRed to be this Matrix in row reduced form if EF is false and
     * sets _rowRedEF to be this Matrix in row reduced echelon form if EF is
     * true. Simultaneously computes rank/nullity while checking for linear
     * independence of columns and injectivity/surjectivity of this Matrix.
     *
     * @throws MatrixException */
    public void rowReduction(Boolean EF) throws MatrixException {
        Matrix B = Operations.matrixCopy(this);
        _pivotCols = new ArrayList<Integer>();
        _pivotRows = new ArrayList<Integer>();
        int pivot = 1;
        for (int c = 1; c <= B.getWidth(); c++) {
            if (pivot > B.getHeight()) {
                break;
            }
            if (B.count(c, pivot) == 0) {
                continue;
            } else if (B.get(pivot, c) == 0) {
                int k = pivot + 1;
                while (B.get(k, c) == 0) {
                    k++;
                }
                B.switchRow(pivot, k);
            }
            if (EF == true) {
                for (int r = 1; r <= B.getHeight(); r++) {
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
            B.scalarMultRow(pivot, 1 / B.get(pivot, c));
            _pivotCols.add(c);
            _pivotRows.add(pivot);
            pivot++;
        }
        if (EF == true) {
            _rowRedEF = B;
        } else {
            _rowRed = B;
        }
        if (_rank == null) {
            _rank = pivot - 1;
            _nullity = getWidth() - _rank;
            _linInd = (_rank == getWidth());
            _surjective = (_rank == getHeight());
            _injective = (_nullity == 0);
        }
    }

    /** Returns the change of basis Matrix given another Matrix MATRIX.
     * 
     * @throws MatrixException */
    public Matrix changeOfBasis(Matrix matrix) throws MatrixException {
        if (matrix.getHeight() != getHeight() || matrix.getWidth() != getWidth()) {
            throw new MatrixException("The numbers of vectors do not match.");
        }
        Matrix B = Operations.matrixCopy(this);
        int pivot = 1;
        for (int c = 1; c <= B.getWidth(); c++) {
            if (matrix.get(pivot, c) == 0) {
                int k = pivot + 1;
                while (matrix.get(k, c) == 0) {
                    k++;
                }
                B.switchRow(pivot, k);
                matrix.switchRow(pivot, k);
            }
            for (int r = pivot + 1; r <= B.getHeight(); r++) {
                B.add(pivot, r, -1 * matrix.get(r, c) / matrix.get(pivot, c));
                matrix.add(pivot, r, -1 * matrix.get(r, c) / matrix.get(pivot, c));
            }
            B.scalarMultRow(pivot, 1 / matrix.get(pivot, c));
            matrix.scalarMultRow(pivot, 1 / matrix.get(pivot, c));
            
            pivot++;
        }
        for (int r = 2; r <= getHeight(); r++) {
            for (int row = r - 1; row > 0; row--) {
                B.add(r, row, -1 * matrix.get(row, r));
                matrix.add(r, row, -1 * matrix.get(row, r));
            }
        }
        return B;
    }

    /** Returns the basis for the column space of this matrix.
     * 
     * @throws MatrixException */
    public VectorSet columnSpace() throws MatrixException {
        ArrayList<Integer> columns = new ArrayList<>();
        VectorSet basis = new VectorSet();
        int count;
        for (double[] column : getRowRed()._contents) {
            count = 0;
            for (double value : column) {
                if (value != 0) {
                    columns.add(count);
                    break;
                }
                count++;
            }
        }
        for (int column : columns) {
            double[] vector = new double[_contents.length];
            for (int index = 0; index < _contents.length; index++) {
                vector[index] = get(index + 1, column + 1);
            }
            basis.add(new Vector(vector));
        }
        return basis;
    }
    
    /** Returns the basis for the null space of this Matrix as a VectorSet. */
    public VectorSet nullSpace() throws MatrixException {
        if (isLinInd()) {
            double[] zeroes = new double[getHeight()];
            Vector zero = new Vector(zeroes);
            return new VectorSet(zero);
        }
        // rowReduction(true);
        // the code only adds stuff into _pivots if I run the above line but
        // I shouldn't have to since the line just below this should automatically
        // call the above line
        Matrix RREF = getRowRedEF();
        Vector[] result = new Vector[getNullity()];
        int count = 0;
        for (int c = 1; c <= getWidth(); c++) {
            double[] current = new double[getWidth()];
            if (_pivotCols.contains(c)) {
                continue;
            }
            current[c - 1] = 1;
//          CAN POSSIBLY BE REMOVED;
//          THIS EDGE CASE SHOULD BE COVERED BY REGULAR ALGORITHM
//            if (RREF.isZero(c)) {
//                result[count] = new Vector(current);
//                count++;
//                continue;
//            }
            int r = 1;
            for (int p : _pivotCols) {
                while (RREF.get(r, p) < epsilon) {
                    r++;
                }
                current[p - 1] = -RREF.get(r, c);
            }
            result[count] = new Vector(current);
            count++;
        }
        return new VectorSet(result);
    }
    
    
    /** Returns a vector which is one possible solution to the system
     * of linear equations, Ax=b, where A is this matrix and B is
     * a vector. Stores the general solution to this system in SOLSET. */
    public Vector solve(Vector b, String[] solset) throws MatrixException {
        double[][] contents = new double[getHeight()][];
        for(int c = 0; c < getHeight(); c++){
            double[] matrix = _contents[c];
            contents[c] = new double[getWidth() + 1];
            System.arraycopy(matrix, 0, contents[c], 0, getWidth());
            contents[c][getWidth()] = b.values()[c];
        }
        Matrix augmented = new Matrix(getHeight(), getWidth() + 1, contents);
        Matrix RREF = augmented.getRowRedEF();
        double[] result = new double[getWidth()];
        if (solset == null) {
            solset = new String[getWidth()];
        }
        for (int r = 1; r <= getHeight(); r++) {
            if (RREF.inconsistent(r)) {
                System.out.println("No solution exists. This system is inconsistent.");
                return null;
            } 
        }
        for (int col = 1; col <= getWidth(); col++) {
            if (RREF.getPivotCols().contains(col)) {
                int r = RREF.getPivotRows().get(RREF.getPivotCols().indexOf(col));
                String sol = df.format(RREF.get(col, RREF.getWidth())) + "";
                for (int c = 1; c <= getWidth(); c++) {
                    if (!RREF.getPivotCols().contains(c)) {
                        String toAdd;
                        double value = -1 * RREF.get(r, c);
                        String entry = df.format(Math.abs(value));
                        if (df.format(Math.abs(value)).equals("1")) {
                            entry = "";
                        }
                        if (Math.abs(RREF.get(r, c)) < epsilon) {
                            toAdd = "";
                            continue;
                        } else if (value > 0) {
                            toAdd = " + ";
                        } else {
                            toAdd = " - ";
                        }
                        toAdd += entry + "x" + c;
                        sol += toAdd;
                    }
                }
                solset[col - 1] = sol;
                result[col - 1] = RREF.get(r, RREF.getWidth());
            } else {
                solset[col - 1] = "Free";
                result[col - 1] = 0;
            }
        }
        Vector x = new Vector(result);
        return x;
    }

    /** Returns a vector which is one possible solution to the system
     * of linear equations, Ax=b, where A is this matrix and B is
     * a vector. */
    public Vector solve(Vector b) throws MatrixException {
        return solve(b, null);
    }
    
    /** Returns a string array containing a general solution set to the
     * system of linear equations, Ax=b, where A is this matrix and B is
     * a vector. Prints out this general solution.
     * @throws MatrixException */
    public String[] generalSolution(Vector b) throws MatrixException {
        String[] solset = new String[getWidth()];
        Vector test = solve(b, solset);
        if (test == null) {
            return null;
        }
        for (int i = 1; i <= getWidth(); i++) {
            String entry = "[ x" + i + " ]   ";
            if (i == getWidth()/2) {
                entry += "=  [ " + solset[i - 1] + " ]";
            } else {
                entry += "   [ " + solset[i - 1] + " ]";
            }
            System.out.println(entry);
        }
        return solset;
    }
    
    /** Returns an ArrayList containing the index of the pivot columns of this
     * Matrix (index starting at 1). 
     * @throws MatrixException */
    public ArrayList<Integer> getPivotCols() throws MatrixException {
        if (_pivotCols == null) {
            getRowRed();
        }
        return _pivotCols;
    }
    
    /** Returns an ArrayList containing the index of the pivot row of this
     * Matrix (index starting at 1). 
     * @throws MatrixException */
    public ArrayList<Integer> getPivotRows() throws MatrixException {
        if (_pivotRows == null) {
            getRowRed();
        }
        return _pivotRows;
    }


    /** Returns true if row R is inconsistent, i.e. consists of all zeroes,
     * except in the last column. */
    public boolean inconsistent(int r) {
        for (int c = 1; c < getWidth(); c++) {
            if (get(r, c) >= epsilon) {
                return false;
            }
        }
        if (get(r, getWidth()) >= epsilon) {
            return true;
        } else {
            return false;
        }
    }
    
    /** Returns true if a column contains only zeroes. */
    public boolean isZero(int c) {
        for (int r = 1; r <= getHeight(); r++) {
            if (get(r, c) >= epsilon) {
                return false;
            }
        }
        return true;
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

    /** Returns the transpose of this Matrix. */
    public Matrix getTranspose() throws MatrixException {
        if (_transpose == null) {
            transpose();
        }
        return _transpose;
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

    /** Returns the number of non-zero entries in column C of this Matrix
     * (values smaller than epsilon are treated as 0).
     * @throws MatrixException */
    public int count(int c) throws MatrixException {
        return count(c, 1);
    }

    /** Returns the number of non-zero entries in column C of this Matrix,
     * starting at row R (values smaller than epsilon are treated as 0).
     *
     * @throws MatrixException */
    public int count(int c, int r) throws MatrixException {
        if (c < 1 || c > getWidth()) {
            throw new MatrixException("Column " + c + " is not a valid row.");
        } else if (r < 1 || r > getHeight()) {
            throw new MatrixException("Row " + r + " is not a valid row.");
        }
        int num = 0;
        for (int i = r; i <= getHeight(); i++) {
            if (Math.abs(get(i, c)) >= epsilon) {
                num++;
            }
        }
        return num;
    }
    
    /** Returns a vector set containing the columns of this matrix as vectors. */
    public VectorSet vectorSet() {
        Vector[] vectors = new Vector[getWidth()];
        for (int c = 0; c < getWidth(); c++) {
            double[] values = new double[getHeight()];
            for (int r = 0; r < getHeight(); r++) {
                values[r] = get(r + 1, c + 1);
            }
            vectors[c] = new Vector(values);
        }
        return new VectorSet(vectors);
    }

    /** The contents of this Matrix. */
    private double[][] _contents;

    /** The dimension of this Matrix. */
    private ArrayList<Integer> _dim;
    
    /** The pivot columns of this Matrix. */
    protected ArrayList<Integer> _pivotCols;
    
    /** The pivot rows of this Matrix. */
    protected ArrayList<Integer> _pivotRows;

    /** The rank of this Matrix. */
    protected Integer _rank;

    /** The dimension of the null space of this Matrix. */
    protected Integer _nullity;

    /** A Boolean that is true if the columns of this Matrix are linearly
     * independent. */
    protected Boolean _linInd;

    /** A Boolean that is true if this Matrix is surjective, i.e. if columns of
     * this Matrix span R^m, where m is the dimension of the range. */
    protected Boolean _surjective;

    /** A Boolean that is true if this Matrix is injective. */
    protected Boolean _injective;

    /** The row reduced form of this Matrix. */
    protected Matrix _rowRed;

    /** The row reduced echelon form of this Matrix. */
    protected Matrix _rowRedEF;

    /** The transpose of this Matrix. */
    protected Matrix _transpose;
    
    /** Two doubles are considered equal if they are within this margin. */
    protected static final double epsilon = 0.000001;
        
    /** The format of output for entries in the matrix. */
    protected static final DecimalFormat df = new DecimalFormat("#");

}
