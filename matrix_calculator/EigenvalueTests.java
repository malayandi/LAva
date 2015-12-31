package matrix_calculator;

import static org.junit.Assert.*;

import org.junit.Test;

public class EigenvalueTests {

    // to test: QR
    
    @Test
    public void conversion() throws MatrixException {
        double[][] contentsA = {
                        { 5, 4, 6, 7 },
                        { 1, 0, 3, 8 },
                        { 0, 0, 7, 9 },
                        { 3, 0, 0, 2 }

        };
        SquareMatrix A = new SquareMatrix(4, contentsA);
        
        VectorSet vectorSetA = A.vectorSet();
        assertEquals(vectorSetA.get(0).get(1), 1, 0);
        assertEquals(vectorSetA.get(0).dotProduct(vectorSetA.get(3)), 49, 0);
        
        Matrix convertedMatrixA = vectorSetA.toMatrix();
        Matrix convertedSquareMatrixA = vectorSetA.toSquareMatrix();
        
        assertTrue(A.equals(convertedMatrixA));
        assertTrue(A.equals(convertedSquareMatrixA));
    }
    
    @Test
    public void triangular() throws MatrixException {
        double[][] contentsA = {
                        { 1, 2 },
                        { 0.000001, 4 }
        };
        SquareMatrix A = new SquareMatrix(2, contentsA);

        assertTrue(A.isUpperTriangular());
        assertTrue(A.isTriangular());
        
        double[][] contentsB = {
                        { 1, 0.000001 },
                        { 2, 4 }
        };
        SquareMatrix B = new SquareMatrix(2, contentsB);
        assertTrue(B.isLowerTriangular());
        assertTrue(B.isTriangular());
        
        double[][] contentsC = {
                        { 1, 0.000001 },
                        { 0.000001, 4 }
        };
        SquareMatrix C = new SquareMatrix(2, contentsC);
        assertTrue(C.isDiagonal());
    }
    
    @Test
    public void QR() throws MatrixException {
        double[][] contentsA = {
                        { 12, -51, 4 },
                        { 6, 167, -68 },
                        { -4, 24, -41 }
        };
        SquareMatrix A = new SquareMatrix(3, contentsA);
        
        double[][] contentsQ = {
                        { (double) 6/7, (double) -69/175, (double) -58/175 },
                        { (double) 3/7, (double) 158/175, (double) 6/175 },
                        { (double) -2/7, (double) 6/35, (double) -33/35 }
        };
        SquareMatrix Q = new SquareMatrix(3, contentsQ);
        
        double[][] contentsR = {
                        { 14, 21, -14 },
                        { 0, 175, -70 },
                        { 0, 0, 35 }
        };
        SquareMatrix R = new SquareMatrix(3, contentsR);

        A.getQ().print();
        System.out.println("hello");
        A.getR().print();
    }

}
