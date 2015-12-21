package matrix_calculator;

/** Executes matrix operations.
 * 
 * @author AndyPalan
 *
 */

public static class Artihmetic {

    /** Returns the Matrix that is the result of adding together A and B (A + B). */
    public static Matrix add(Matrix A, Matrix B) {
        if (!A.getDimension().equals(B.getDimension())) {         // Check if ArrayList equality works
            throw new MatrixException("These two matrices do not have"
                            + " the appropriate dimension to be added/subtracted together.");
        }
        Double[][] contents = new Double[A.getHeight()][A.getWidth()];
        for (int r = 1; r <= A.getHeight(); r++) {
            for (int c = 1; c <= A.getWidth(); c++) {
                contents[r-1][c-1] = A.get(r, c) + B.get(r, c);
            }
        }
        return new Matrix(A.getHeight(), A.getWidth(), contents);
    }

    /** Returns the Matrix that is the result of subtracting B from A (A - B). */
    public static Matrix subtract(Matrix A, Matrix B) {
        return add(A, scalarMult(B, -1));
    }

    /** Returns the Matrix that has been scalar multiplied by K. */
    public static Matrix scalarMult(Matrix A, Double k) {
        Double[][] contents = new Double[A.getHeight()][A.getWidth()];
        for (int r = 1; r <= A.getHeight(); r++) {
            for (int c = 1; c <= A.getWidth(); c++) {
                contents[r-1][c-1] = k * A.get(r, c);
            }
        }
        return new Matrix(A.getHeight(), A.getWidth(), contents);
    }

    /** Returns the Matrix that is the result of Matrix multiplying A and
     * B in the form A x B. */
    public static Matrix matrixMult(Matrix A, Matrix B) {
        if (A.getWidth() != B.getHeight()) {
            throw new MatrixException("These two matrices do not have"
                            + " the appropriate dimension to be multiplied together.");
        }
        Double[][] contents = new Double[A.getHeight()][A.getWidth()];
        for (int i = 0, j = 0; i < A.getHeight(); i++, j++) {
            entry = 0;
            for (int r = 1; r <= B.getHeight(); r++) {
                for (int c = 1; c <= A.getWidth(); c++) {
                    entry += (A.get(r, c) * B.get(c, r));
                }
            }
            contents[i][j] = entry;
            if (j == (B.getWidth - 1)) {
                j = -1;
            }
        }
        return new Matrix(A.getHeight(), B.getWidth(), contents);
    }

    /** Copies the contents of Matrix A into a new Matrix and returns
     * that Matrix. */
    public static Matrix matrixCopy(Matrix A) {
        Matrix B = new Matrix(A.getHeight(), A.getWidth(), new Double[][]);
        for (int r = 1; r <= A.getHeight(); r++) {
            for (int c = 1; c <= A.getWidth(); c++) {
                B.set(r, c, A.get(r, c));
            }
        }
        return B;
    }

}
