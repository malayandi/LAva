package matrix_calculator;

/** A class representing a square Matrix object (i.e. an n x n Matrix).
 * 
 * @author AndyPalan */
public class SquareMatrix extends Matrix {
    
    /** Creates a new N x N square Matrix with contents CONTENTS. */
    public SquareMatrix(int n, double[][] contents) throws MatrixException {
        super(n, n, contents);
    }
    
    

}
