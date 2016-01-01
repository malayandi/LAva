package matrix_calculator;

import java.util.ArrayList;

/** A class representing a square Matrix object (i.e. an n x n Matrix).
 * 
 * @author AndyPalan */
public class SquareMatrix extends Matrix {

    /** Creates a new N x N square Matrix with contents CONTENTS. */
    public SquareMatrix(int n, double[][] contents) throws MatrixException {
        super(n, n, contents);
    }

    /** Sets _rowRed to be this Matrix in row reduced form if EF is false and
     * sets _rowRedEF to be this Matrix in row reduced echelon form if EF is
     * true. Simultaneously computes rank/nullity/determinant and inverse
     * while checking for linear independence of columns and injectivity/
     * surjectivity of this Matrix.
     *
     * @throws MatrixException */
    public void squareRowReduction(Boolean EF) throws MatrixException {
        SquareMatrix B = Operations.matrixCopy(this);
        SquareMatrix I = new SquareMatrix(getHeight(), new double[getHeight()][getHeight()]);
        for (int r = 1; r <= I.getHeight(); r++) {
            I.set(r, r, 1);
        }
        int pivot = 1;
        double det = 1;
        for (int c = 1; c <= B.getWidth(); c++) {
            if (B.count(c, pivot) == 0) {
                continue;
            } else if (B.get(pivot, c) == 0) {
                int k = pivot + 1;
                while (B.get(k, c) == 0) {
                    k++;
                }
                B.switchRow(pivot, k);
                I.switchRow(pivot, k);
                det *= -1;
            }
            if (EF == true) {
                for (int r = 1; r <= B.getHeight(); r++) {
                    if (r == pivot) {
                        continue;
                    }
                    double factor = -1 * B.get(r, c) / B.get(pivot, c);
                    B.add(pivot, r, factor);
                    I.add(pivot, r, factor);
                }
            } else {
                for (int r = pivot + 1; r <= B.getHeight(); r++) {
                    B.add(pivot, r, -1 * B.get(r, c) / B.get(pivot, c));
                }
            }
            double factor = 1 / B.get(pivot, c);
            B.scalarMultRow(pivot, factor);
            I.scalarMultRow(pivot, factor);
            det *= (1 / factor);
            pivot++;
        }
        if (EF == true) {
            _rowRedEF = B;
            _inverse = I;
        } else {
            _rowRed = B;
        }
        _det = det;
        if (_rank == null) {
            _rank = pivot - 1;
            _nullity = getWidth() - _rank;
            _linInd = (_rank == getWidth());
            _surjective = (_rank == getHeight());
            _injective = (_nullity == 0);
        }
    }

    /** Returns the row reduced form of this Matrix.
     *
     * @throws MatrixException */
    public SquareMatrix getRowRed() throws MatrixException {
        if (_rowRed == null) {
            squareRowReduction(false);
        }
        return _rowRed;
    }

    /** Returns the row reduced echelon form of this Matrix.
     *
     * @throws MatrixException */
    public SquareMatrix getRowRedEF() throws MatrixException {
        if (_rowRedEF == null) {
            squareRowReduction(true);
        }
        return _rowRedEF;
    }

    /** Returns the determinant of this Matrix. */
    public double getDet() throws MatrixException {
        if (_det == null) {
            getRowRed();
        }
        return _det;
    }

    /** Returns the trace of this Matrix. */
    public double getTrace() {
        if (_trace == null) {
            double trace = 0;
            for (int i = 1; i <= getHeight(); i++) {
                trace += get(i, i);
            }
            _trace = trace;
        }
        return _trace;
    }

    /** Returns the inverse of this Matrix. 
     * 
     * @throws MatrixException */
    public SquareMatrix getInverse() throws MatrixException {
        if (_inverse == null) {
            getRowRedEF();
        }
        return _inverse;
    }

    /** Returns the transpose of this Matrix. */
    public SquareMatrix getTranspose() throws MatrixException {
        if (_transpose == null) {
            transpose();
        }
        return _transpose;
    }

    /** Sets _transpose to be the transpose of this Matrix. */
    public void transpose() throws MatrixException {
        SquareMatrix T = new SquareMatrix(getHeight(), new double[getHeight()][getHeight()]);
        for (int r = 1; r <= getHeight(); r++) {
            for (int c = 1; c <= getWidth(); c++) {
                T.set(c, r, get(r, c));
            }
        }
        _transpose = T;
    }

    /** Sets _Q and _R to be the QR factorised form of this Matrix. 
     * @throws MatrixException */
    public void QR() throws MatrixException {
        VectorSet columns = vectorSet();
        VectorSet orthogonalised = columns.gramSchmidt(0, columns);
        VectorSet orthonormalised = orthogonalised.orthonormalise();
        SquareMatrix Q = orthonormalised.toSquareMatrix();
        SquareMatrix R = new SquareMatrix(Q.getHeight(), new double[getHeight()][getHeight()]);
        for (int c = 1; c <= getHeight(); c++) {
            for (int r = 1; r <= c; r++) {
                R.set(r, c, columns.get(c - 1).dotProduct(orthonormalised.get(r - 1)));
            }
        }
        _Q = Q;
        _R = R;
    }

    /** Returns an array list containing the QR factorised form of this Matrix. 
     * @throws MatrixException */
    public ArrayList<SquareMatrix> getQR() throws MatrixException {
        if (_Q == null) {
            QR();
        }
        ArrayList<SquareMatrix> result = new ArrayList<SquareMatrix>();
        result.add(_Q);
        result.add(_R);
        return result;
    }

    /** Returns Q from the QR factorised form of this Matrix. 
     * @throws MatrixException */
    public SquareMatrix getQ() throws MatrixException {
        if (_Q == null) {
            QR();
        }
        return _Q;
    }

    /** Returns R from the QR factorised form of this Matrix. 
     * @throws MatrixException */
    public SquareMatrix getR() throws MatrixException {
        if (_R == null) {
            QR();
        }
        return _R;
    }


    /** Returns true if this Matrix is diagonal (values smaller than epsilon are treated as 0). */
    public boolean isDiagonal() {
        for (int r = 1; r <= getHeight(); r++) {
            for (int c = 1; c <= getWidth(); c++) {
                if (r == c) {
                    continue;
                } else if (Math.abs(get(r, c)) >= epsilon) {
                    return false;
                }
            }
        }
        return true;
    }

    /** Returns true if this Matrix is upper triangular. 
     * 
     * @throws MatrixException */
    public boolean isUpperTriangular() throws MatrixException {
        for (int i = 1; i < getHeight(); i++) {
            if (count(i, i + 1) != 0) {
                return false;
            }
        }
        return true;
    }

    /** Returns true if this Matrix is lower triangular. */
    public boolean isLowerTriangular() throws MatrixException {
        return getTranspose().isUpperTriangular();
    }

    /** Returns true if this Matrix is triangular. */
    public boolean isTriangular() throws MatrixException {
        return isUpperTriangular() || isLowerTriangular();
    }

    /** Sets _eigenvalues to contain the eigenvalues of this Matrix. 
     * 
     * @throws MatrixException */
    public void eigenvalues() throws MatrixException {
        _eigenvalues = new ArrayList<Double>();
        SquareMatrix A = this;
        if (!(isTriangular() || isDiagonal())) {
            SquareMatrix Q = getQ();
            SquareMatrix R = getR();
            A = Operations.matrixMult(R, Q);
            while (!A.isTriangular()) {
                Q = A.getQ();
                R = A.getR();
                A = Operations.matrixMult(R, Q);
            }
        }
        for (int i = 1; i <= getHeight(); i++) {
            _eigenvalues.add(A.get(i, i));
        }
    }
    
    public ArrayList<Double> getEigenvalues() throws MatrixException {
        if (_eigenvalues == null) {
            eigenvalues();
        }
        return _eigenvalues;
    }

    /** The row reduced form of this Matrix. */
    private SquareMatrix _rowRed;

    /** The row reduced echelon form of this Matrix. */
    private SquareMatrix _rowRedEF;

    /** The determinant of this Matrix. */
    private Double _det;

    /** The trace of this Matrix. */ 
    private Double _trace;

    /** The inverse of this Matrix. */
    private SquareMatrix _inverse;

    /** The transpose of this Matrix. */
    protected SquareMatrix _transpose;

    /** The eigenvalues of this Matrix. */
    private ArrayList<Double> _eigenvalues;

    /** The eigvenvectors of this Matrix. */
    private ArrayList<Vector> _eigenvectors;

    /** True if the Matrix is diagonalisable. */
    private Boolean _diagonalisable;

    /** The diagonalised form of this Matrix. */
    private SquareMatrix _diagonalised;

    /** Q from the QR-Factorised form of this Matrix. */
    private SquareMatrix _Q;

    /** R from the QR-Factorised form of this Matrix. */
    private SquareMatrix _R;
}
