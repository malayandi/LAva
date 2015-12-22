package matrix_calculator;

/** A new type of exception for the Matrix Calculator
 *
 * @author AndyPalan */
@SuppressWarnings("serial")
public class MatrixException extends Exception {

    public MatrixException() {
    }

    public MatrixException(String msg) {
        super(msg);
    }

}
