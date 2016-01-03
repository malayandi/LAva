package matrix_calculator;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;

public class EigenvalueTests {
    
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
                        { 0.000000001, 4 }
        };
        SquareMatrix A = new SquareMatrix(2, contentsA);

        assertTrue(A.isUpperTriangular());
        assertTrue(A.isTriangular());
        
        double[][] contentsB = {
                        { 1, 0.000000001 },
                        { 2, 4 }
        };
        SquareMatrix B = new SquareMatrix(2, contentsB);
        assertTrue(B.isLowerTriangular());
        assertTrue(B.isTriangular());
        
        double[][] contentsC = {
                        { 1, 0.000000001 },
                        { 0.0000000001, 4 }
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

        assertTrue(A.getQ().equals(Q));
        assertTrue(A.getR().equals(R));
    }
    
    @Test
    public void eigenvalues() throws MatrixException {
        double[][] contentsA = {
                        { 1, 2, 3 },
                        { 0, 4, 5 },
                        { 0, 0, 6 }
        };
        SquareMatrix A = new SquareMatrix(3, contentsA);
        Double[] eigenvaluesA = {1.0, 4.0, 6.0};
        
        for (int i = 0; i < eigenvaluesA.length; i++) {
            assertEquals(eigenvaluesA[i], A.getEigenvalues().get(i));
        }
        
        double[][] contentsB = {
                        { 1, 2, 3 },
                        { 4, 5, 6 },
                        { 8, 8, 9 }
        };
        SquareMatrix B = new SquareMatrix(3, contentsB);
        Double[] eigenvaluesB = {16.2787, -1.4095, 0.1308};
        
        for (int i = 0; i < eigenvaluesB.length; i++) {
            assertEquals(eigenvaluesB[i], B.getEigenvalues().get(i), 0.0001);
        }

        double[][] contentsC = {
                        { 1, 2, 3 },
                        { 4, 5, 6 },
                        { 7, 8, 9 }
        };
        SquareMatrix C = new SquareMatrix(3, contentsC);
        Double[] eigenvaluesC = { 16.1168, -1.1168, 0.0 };
        
        for (int i = 0; i < eigenvaluesC.length; i++) {
            assertEquals(eigenvaluesC[i], C.getEigenvalues().get(i), 0.0001);
        }
    }
    
    @Test
    public void nullSpace() throws MatrixException {
        double[][] contentsA = {
                        { 1, 2, 3 },
                        { 4, 5, 6 },
                        { 7, 8, 9 }
        };
        SquareMatrix A = new SquareMatrix(3, contentsA);
        
        double[] contents = { 1.0, -2.0, 1.0 };
        Vector v = new Vector(contents);
        
        assertEquals(A.nullSpace().size(), 1);
        assertTrue(A.nullSpace().get(0).equals(v));

        double[][] contentsB = {
                        { 1, 2, 3, 4 },
                        { 4, 5, 6, 7 },
                        { 7, 8, 9, 10 }
        };
        Matrix B = new Matrix(3, 4, contentsB);
        
        double[] c1 = { 1.0, -2.0, 1.0, 0.0 };
        double[] c2 = { 2.0, -3.0, 0.0, 1.0 };
        Vector v1 = new Vector(c1); 
        Vector v2 = new Vector(c2);
        
        assertEquals(B.nullSpace().size(), 2);
        assertTrue(B.nullSpace().get(0).equals(v1));
        assertTrue(B.nullSpace().get(1).equals(v2));      
    }

    @Test
    public void eigenvector() throws MatrixException {
        double[][] contentsA = {
                        { 1, 2, 3 },
                        { 4, 5, 6 },
                        { 7, 8, 9 }
        };
        SquareMatrix A = new SquareMatrix(3, contentsA);

        double[] c1 = { 0.283349, 0.641675, 1 };
        Vector v1 = new Vector(c1);
        double[] c2 = { -1.28335, -0.141675, 1 };
        Vector v2 = new Vector(c2);
        double[] c3 = { 1, -2, 1 };
        Vector v3 = new Vector(c3);
        
        ArrayList<Vector> eigenvectors = A.getEigenvectors();
                
        assertEquals(eigenvectors.size(), 3);
        assertTrue(eigenvectors.get(0).equals(v1));
//        assertTrue(eigenvectors.get(1).equals(v2));
        assertTrue(eigenvectors.get(2).equals(v3));
        
        double[][] contentsB = {
                        { 1, 2, 3, 4 },
                        { 5, 6, 7, 8 },
                        { 9, 10, 11, 12 },
                        { 13, 14, 15, 16 }
        };
        SquareMatrix B = new SquareMatrix(4, contentsB);

        double[] d1 = { 0.202782, 0.468521, 0.734261, 1 };
        Vector w1 = new Vector(d1);
        double[] d2 = { -1.20278, -0.468521, 0.265739, 1 };
        Vector w2 = new Vector(d2);
        double[] d3 = { 2, -3, 0, 1 };
        Vector w3 = new Vector(d3);
        double[] d4 = { 1, -2, 1, 0 };
        Vector w4 = new Vector(d4);
        
        ArrayList<Vector> eigenvectorsB = B.getEigenvectors();
        
        for (Vector v : B.getEigenvectors()) {
            v.print();
            System.out.println("");
        }
        
        assertEquals(eigenvectorsB.size(), 4);
        assertTrue(eigenvectorsB.get(0).equals(w1));
        assertTrue(eigenvectorsB.get(1).equals(w2));
        assertTrue(eigenvectorsB.get(2).equals(w3));
        assertTrue(eigenvectorsB.get(3).equals(w4));

    }
}
