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
    
    /** The determinant of this Matrix. */
    private double _det;
    
    /** The trace of this Matrix. */ 
    private double _trace;
    
    /** The inverse of this Matrix. */
    private Matrix _inverse;
    
    /** The eigenvalues of this Matrix. */
    private ArrayList<Double> _eigenvalues;
    
    /** The eigvenvectors of this Matrix. */
    private ArrayList<Vector> _eigenvectors;
    
    
    
    

}
