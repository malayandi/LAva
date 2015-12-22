package matrix_calculator;

import static org.junit.Assert.*;
import org.junit.Test;

public class MatrixTest {

    @Test
    public void basicFunctionality() throws MatrixException {
        double[][] contents = {
                        { 5, 4, 6 },
                        { 6, 3, 2 }
        };
        Matrix A = new Matrix(2, 3, contents);
        
        assertEquals(A.getHeight(), 2);
        assertEquals(A.getWidth(), 3);
        assertEquals(A.get(1, 2), 4, 0);
        
        Matrix B = new Matrix(2, 3, new double[2][3]);
        B.set(1, 1, 5);
        B.set(1, 2, 4);
        B.set(1, 3, 6);
        B.set(2, 1, 6);
        B.set(2, 2, 3);
        B.set(2, 3, 2);
        assertTrue(A.equals(B));

        A.set(2, 3, 15.9);
        assertEquals(A.get(2, 3), 15.9, 0);
        A.print();        
    }
    
    @Test
    public void elementaryRowOperations() throws MatrixException {
        double[][] contents = {
                        { 5, 4, 6 },
                        { 6, 3, 2 }
        };
        Matrix A = new Matrix(2, 3, contents);

        System.out.println("");
        A.scalarMultRow(1, 3);
        A.print();
        
        System.out.println("");
        A.scalarMult(2);
        A.print();
        
        System.out.println("");
        A.switchRow(1, 2);
        A.print();
        
        System.out.println("");
        A.add(1, 2, 2);
        A.print();
    }
    
    @Test
    public void count() throws MatrixException {
        double[][] contents = {
                        { 5, 4, 6 },
                        { 1, 0, 3 },
                        { 0, 0, 7 },
                        { 3, 0, 0 }

        };
        Matrix A = new Matrix(4, 3, contents);
        
        assertEquals(A.count(1), 3);
        assertEquals(A.count(2), 1);
        assertEquals(A.count(3), 3);
        
        assertEquals(A.count(1, 3), 1);
        assertEquals(A.count(2, 2), 0);
        assertEquals(A.count(2, 3), 0);
    }

}
