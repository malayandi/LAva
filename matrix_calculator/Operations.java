package matrix_calculator;

/** Executes matrix operations.
 *
 * @author AndyPalan */

public class Operations {

    /** Returns the Matrix that is the result of adding together A and B (A +
     * B).
     *
     * @throws MatrixException */
    public static Matrix add(Matrix A, Matrix B) throws MatrixException {
        if (!A.getDimension().equals(B.getDimension())) { // Check if ArrayList
                                                          // equality works
            throw new MatrixException("These two matrices do not have"
                            + " the appropriate dimension to be added/subtracted together.");
        }
        double[][] contents = new double[A.getHeight()][A.getWidth()];
        for (int r = 1; r <= A.getHeight(); r++) {
            for (int c = 1; c <= A.getWidth(); c++) {
                contents[r - 1][c - 1] = A.get(r, c) + B.get(r, c);
            }
        }
        return new Matrix(A.getHeight(), A.getWidth(), contents);
    }

    /** Returns the Matrix that is the result of subtracting B from A (A - B).
     * 
     * @throws MatrixException */
    public static Matrix subtract(Matrix A, Matrix B) throws MatrixException {
        Matrix C = matrixCopy(B);
        C.scalarMult(-1);
        return add(A, C);
    }

    /** Returns the Matrix that has been scalar multiplied by K. 
     * @throws MatrixException */
    public static Matrix scalarMult(Matrix A, double k) throws MatrixException {
        double[][] contents = new double[A.getHeight()][A.getWidth()];
        for (int r = 1; r <= A.getHeight(); r++) {
            for (int c = 1; c <= A.getWidth(); c++) {
                contents[r - 1][c - 1] = k * A.get(r, c);
            }
        }
        return new Matrix(A.getHeight(), A.getWidth(), contents);
    }

    /** Returns the Matrix that is the result of Matrix multiplying A and B in
     * the form A x B.
     *
     * @throws MatrixException */
    public static Matrix matrixMult(Matrix A, Matrix B) throws MatrixException {
        if (A.getWidth() != B.getHeight()) {
            throw new MatrixException("These two matrices do not have"
                            + " the appropriate dimension to be multiplied together.");
        }
        double[][] contents = new double[A.getHeight()][B.getWidth()];
        for (int i = 0, j = 0; i < A.getHeight(); j++) {
            double entry = 0;
            for (int c = 1; c <= A.getWidth(); c++) {
                entry += ((A.get(i + 1, c) * B.get(c, j + 1)));
            }
            contents[i][j] = entry;
            if (j == (B.getWidth() - 1)) {
                j = -1;
                i++;
            }
        }
        return new Matrix(A.getHeight(), B.getWidth(), contents);
    }
    
    /** Returns the Square Matrix that is the result of Matrix multiplying A
     * and B in the form A x B.
     *
     * @throws MatrixException */
    public static SquareMatrix matrixMult(SquareMatrix A, SquareMatrix B) throws MatrixException {
        if (A.getWidth() != B.getHeight()) {
            throw new MatrixException("These two matrices do not have"
                            + " the appropriate dimension to be multiplied together.");
        }
        double[][] contents = new double[A.getHeight()][A.getWidth()];
        for (int i = 0, j = 0; i < A.getHeight(); j++) {
            double entry = 0;
            for (int c = 1; c <= A.getWidth(); c++) {
                entry += ((A.get(i + 1, c) * B.get(c, j + 1)));
            }
            contents[i][j] = entry;
            if (j == (B.getWidth() - 1)) {
                j = -1;
                i++;
            }
        }
        return new SquareMatrix(A.getHeight(), contents);
    }


    /** Copies the contents of Matrix A into a new Matrix and returns that
     * Matrix. 
     * @throws MatrixException */
    public static Matrix matrixCopy(Matrix A) throws MatrixException {
        Matrix B = new Matrix(A.getHeight(), A.getWidth(),
                        new double[A.getHeight()][A.getWidth()]);
        for (int r = 1; r <= A.getHeight(); r++) {
            for (int c = 1; c <= A.getWidth(); c++) {
                B.set(r, c, A.get(r, c));
            }
        }
        return B;
    }
    
    /** Copies the contents of Matrix A into a new Matrix and returns that
     * Matrix. 
     * @throws MatrixException */
    public static SquareMatrix matrixCopy(SquareMatrix A) throws MatrixException {
        SquareMatrix B = new SquareMatrix(A.getHeight(),
                        new double[A.getHeight()][A.getWidth()]);
        for (int r = 1; r <= A.getHeight(); r++) {
            for (int c = 1; c <= A.getWidth(); c++) {
                B.set(r, c, A.get(r, c));
            }
        }
        return B;
    }


}
