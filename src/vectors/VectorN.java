/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vectors;

/**
 *
 * @author bowen
 */
public class VectorN implements Vector {

    private final double[] content;
    
    public VectorN(int n) {
        content = new double[n];
    }
    
    private VectorN(double[] content) {
        this.content = new double[content.length];
        System.arraycopy(content, 0, this.content, 0, content.length);
    }
    
    
    @Override
    public Vector copy() {
        return new VectorN(content);
    }

    @Override
    public double[] content() {
        return content;
    }

    @Override
    public Vector set(int i, double d) {
        content[i] = d;
        return this;
    }
    
}
