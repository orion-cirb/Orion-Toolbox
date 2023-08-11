
package Orion_Toolbox;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;


/**
 * Print stream that prints nothing
 * @author Sylvain Hallé
 */
public class NullPrintStream extends PrintStream {
    
    /**
     * Creates a new null print stream
     */
    public NullPrintStream() {
      super(new NullByteArrayOutputStream());
    }
    
    /**
     * Output stream that prints nothing
     */
    private static class NullByteArrayOutputStream extends ByteArrayOutputStream {

        @Override
        public void write(int b) {
            // do nothing
        }

        @Override
        public void write(byte[] b, int off, int len) {
            // do nothing
        }

        @Override
        public void writeTo(OutputStream out) throws IOException {
            // do nothing
        }
    }
}


