package matrix_calculator;

import static org.junit.Assert.*;

import org.junit.Test;

public class DiagonalisationDemo {

    @Test
    public void diagonalised() throws MatrixException {
        // Defining A
        double[][] contents = {
                        { 1, 2, 0 },
                        { 0, 3, 0 },
                        { 2, 4, -2 }
        };
        System.out.println("Non-diagonalised:");
        long startTime = System.nanoTime();
        SquareMatrix A = new SquareMatrix(3, contents);
        
        // Finding A^10
        SquareMatrix Ak = Operations.exp(A, 10);
        
        System.out.println("");
        System.out.println("A^10 = ");
        Ak.print();
        long endTime = System.nanoTime();
        System.out.println("");
        System.out.println("It took " + (endTime - startTime)/1000 + " miliseconds"
                        + " to find the 10th exponent of A");
        
        
        System.out.println("");
        System.out.println("Diagonalised:");
        startTime = System.nanoTime();
        
        // Finding D - diagonalised form of A
        SquareMatrix D = A.getDiagonalised().get(1);
        
        System.out.println("");
        System.out.println("D = ");
        D.print();
        System.out.println("");
        System.out.println("P = ");
        
        // Finding P
        A.getDiagonalised().get(0).print();
        
        // Finding D^10
        SquareMatrix Dk = Operations.exp(D, 10);
        
        // Finding P * D^10 * P-1
        Ak = Operations.matrixMult(A.getDiagonalised().get(0), Dk);
        Ak = Operations.matrixMult(Ak, A.getDiagonalised().get(2));
        
        System.out.println("");
        System.out.println("A^10 = ");
        Ak.print();
        endTime = System.nanoTime();
        System.out.println("");
        System.out.println("It took " + (endTime - startTime)/1000 + " miliseconds"
                        + " to find the 10th exponent of A");
    }

}
