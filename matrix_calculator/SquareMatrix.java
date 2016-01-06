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
        _pivotCols = new ArrayList<Integer>();
        _pivotRows = new ArrayList<Integer>();
        for (int r = 1; r <= I.getHeight(); r++) {
            I.set(r, r, 1);
        }
        int pivot = 1;
        double det = 1;
        for (int c = 1; c <= B.getWidth(); c++) {
            // Check for free column
            if (B.count(c, pivot) == 0) {
                continue;
            }
            
            // Swap pivot with largest absolute value
            int max = pivot;
            for (int i = pivot + 1; i <= B.getHeight(); i++) {
                if (Math.abs(B.get(i, c)) > Math.abs(B.get(max, c))) {
                    max = i;
                }
            }
            if (max != pivot) {
                B.switchRow(pivot, max);
                I.switchRow(pivot, max);
                det *= -1;
//                B.print();
//                System.out.println("");
            }

            // Scale the row for the pivot to have value 1
            double scalefactor = 1 / B.get(pivot, c);
            B.scalarMultRow(pivot, scalefactor);
            I.scalarMultRow(pivot, scalefactor);
            det *= (1 / scalefactor);
//            B.print();
//            System.out.println("");
            

            // Elimination
            if (EF == true) {
                for (int r = 1; r <= B.getHeight(); r++) {
                    if (r == pivot) {
                        continue;
                    }
                    double elimfactor = -1 * B.get(r, c) / B.get(pivot, c);
                    B.add(pivot, r, elimfactor);
                    I.add(pivot, r, elimfactor);
//                    B.print();
//                    System.out.println("");
                }
            } else {
                for (int r = pivot + 1; r <= B.getHeight(); r++) {
                    double elimfactor = -1 * B.get(r, c) / B.get(pivot, c);
                    B.add(pivot, r, elimfactor);
//                    B.print();
//                    System.out.println("");
                }
            }
            _pivotCols.add(c);
            _pivotRows.add(pivot);
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

    
    /** Returns true if this matrix is similar to matrix MATRIX. 
     * @throws MatrixException */
    public boolean similar(SquareMatrix matrix) throws MatrixException{
        if (_eigenvalues == null) {
            eigenvalues();
        }
        ArrayList<Double> eigenvalues = new ArrayList<>();
        for(double eigenvalue : matrix.getEigenvalues()){
            eigenvalues.add(eigenvalue);
        }
        ArrayList<Double> checked = new ArrayList<>();
        for (double eigenvalue : _eigenvalues){
            for (double e : matrix.getEigenvalues()) {
                if (Math.abs(eigenvalue - e) < Matrix.epsilon * 100 && !checked.contains(eigenvalue)) {
                    System.out.println(eigenvalue + " " + e);
                    checked.add(eigenvalue);
                    eigenvalues.remove(e);
                    continue;
                }
            }
        }
        if(eigenvalues.size() > 0){
            return false;
        }
        return true;
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
//                Q.print();
//                System.out.println("");
                R = A.getR();
//                R.print();
//                System.out.println("");
                A = Operations.matrixMult(R, Q);
//                A.print();
//                System.out.println("");
            }
        }
        for (int i = 1; i <= getHeight(); i++) {
            _eigenvalues.add(A.get(i, i));
        }
    }
    
    /** Returns an ArrayList containing the eigenvalues of this Matrix. */
    public ArrayList<Double> getEigenvalues() throws MatrixException {
        if (_eigenvalues == null) {
            eigenvalues();
        }
        return _eigenvalues;
    }
    
    /** Sets _eigenvectors to contain the eigenvectors of this Matrix. */
    public void eigenvectors() throws MatrixException {
        ArrayList<Double> eigenvalues = getEigenvalues();
        ArrayList<Double> computed = new ArrayList<>();
        ArrayList<Vector> eigenvectors = new ArrayList<>();
        for (double value : eigenvalues) {
            if (computed.contains(value)) {
                continue;
            }
            SquareMatrix copy = Operations.matrixCopy(this);
            for (int i = 1; i <= getHeight(); i++) {
                copy.set(i, i, get(i, i) - value);
            }
            VectorSet vectors = copy.nullSpace();
            for (int i = 0; i < vectors.size(); i++) {
                eigenvectors.add(vectors.get(i));
            }
            computed.add(value);
        }
        for (Vector v : eigenvectors) {
            v.scaleWholeNum();
        }
        _eigenvectors = eigenvectors;
    }
    
    /** Returns an ArrayList of vectors containing the eigenvectors of this Matrix. */
    public ArrayList<Vector> getEigenvectors() throws MatrixException {
        if (_eigenvectors == null) {
            eigenvectors();
        }
        return _eigenvectors;
    }
    
    /** Returns true if this Matrix is diagonalisable. 
     * 
     * @throws MatrixException */
    public boolean isDiagonalisable() throws MatrixException {
        return (getEigenvectors().size() == getWidth());
    }
    
    /** Sets _diagonalised to contain P, D, P-1 (in that order) where
     * this matrix, A is expressed as A = PDP-1 and D is a diagonal matrix. */
    public void diagonalise() throws MatrixException {
        Vector[] vectors = new Vector[getWidth()];
        SquareMatrix D = new SquareMatrix(getWidth(), new double[getWidth()][getWidth()]);
        ArrayList<SquareMatrix> result = new ArrayList<>();

        for (int i = 0; i < getWidth(); i++) {
            vectors[i] = getEigenvectors().get(i);
            D.set(i + 1, i + 1, getEigenvalues().get(i));
        }
        
        VectorSet eigenvectors = new VectorSet(vectors);
        SquareMatrix P = eigenvectors.toSquareMatrix(); 
        
        result.add(P);
        result.add(D);
        result.add(P.getInverse());
        
        _diagonalised = result;
    }
    
    /** Returns an ArrayList containing P, D, P-1 (in that order) where
     * this matrix, A is expressed as A = PDP-1 and D is a diagonal matrix. */
    public ArrayList<SquareMatrix> getDiagonalised() throws MatrixException {
        if (!isDiagonalisable()) {
            System.out.println("This matrix is not diagonalisable.");
            return null;
        }
        if (_diagonalised == null) {
            diagonalise();
        }
        return _diagonalised;
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

    /** The eigenvectors of this Matrix. */
    private ArrayList<Vector> _eigenvectors;

    /** The matrices, P, D and P-1 where this matrix, A = PDP-1 and D is diagonal. */
    private ArrayList<SquareMatrix> _diagonalised;

    /** Q from the QR-Factorised form of this Matrix. */
    private SquareMatrix _Q;

    /** R from the QR-Factorised form of this Matrix. */
    private SquareMatrix _R;
}
