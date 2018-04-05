/*
 * The MIT License
 *
 * Copyright 2017 Bowen Peng.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package vectors;

import java.util.function.BiFunction;
import java.util.function.BinaryOperator;

/**
 *
 * @author bowen
 */
public abstract class Vectors {
    
    /**
     * Finds the minimum elements between the two vectors
     * @param <T> Type of vector
     * @param vector First vector
     * @param vector2 Second vector
     * @return A new vector with its elements equal to min(vector, vector2)
     */
    public static <T extends Vector<T>> T min(T vector, T vector2) {
        if (vector.size() == vector2.size()) {
            T newVector = vector.shell();
            for (int i=0; i<vector.size(); i++) {
                newVector.set(i, Math.min(vector.get(i), vector2.get(i)));
            }
            return newVector;
        }
        throw new IllegalArgumentException("Different vector sizes.");
    }
    /**
     * Finds the maximum elements between the two vectors
     * @param <T> Type of vector
     * @param vector First vector
     * @param vector2 Second vector
     * @return A new vector with its elements equal to max(vector, vector2)
     */
    public static <T extends Vector<T>> T max(T vector, T vector2) {
        if (vector.size() == vector2.size()) {
            T newVector = vector.shell();
            for (int i=0; i<vector.size(); i++) {
                newVector.set(i, Math.min(vector.get(i), vector2.get(i)));
            }
            return newVector;
        }
        throw new IllegalArgumentException("Different vector sizes.");
    }
    /**
     * Performs element-wise addition of two vectors of same size
     * @param <T> Type of vector
     * @param vector First vector
     * @param vector2 Second vector
     * @return A new vector with its elements equal to vector + vector2
     */
    public static <T extends Vector<T>> T add(T vector, T vector2) {
        if (vector.size() == vector2.size()) {
            return vector.copy().add(vector2);
        }
        throw new IllegalArgumentException("Different vector sizes.");
    }
    /**
     * Performs element-wise addition of a vector with a value
     * @param <T> Type of vector
     * @param vector Vector
     * @param d Value
     * @return A new vector with its elements equal to vector + d
     */
    public static <T extends Vector<T>> T add(T vector, double d) {
        return vector.copy().add(d);
    }
    
    /**
     * Performs element-wise subtraction of two vectors of same size
     * @param <T> Type of vector
     * @param vector First vector
     * @param vector2 Second vector
     * @return A new vector with its elements equal to vector - vector2
     */
    public static <T extends Vector<T>> T sub(T vector, T vector2) {
        if (vector.size() == vector2.size()) {
            return vector.copy().sub(vector2);
        }
        throw new IllegalArgumentException("Different vector sizes.");
    }
    public static <T extends Vector<T>> T sub(T vector, double d) {
        return vector.copy().sub(d);
    }
    
    /**
     * Performs element-wise multiplication of two vectors of same size
     * @param <T> Type of vector
     * @param vector First vector
     * @param vector2 Second vector
     * @return A new vector with its elements equal to vector * vector2
     */
    public static <T extends Vector<T>> T mulElem(T vector, T vector2) {
        if (vector.size() == vector2.size()) {
            return vector.copy().mulElem(vector2);
        }
        throw new IllegalArgumentException("Different vector sizes.");
    }
    public static <T extends Vector<T>> T mulElem(T vector, double d) {
        return vector.copy().mulElem(d);
    }
    
    /**
     * Performs element-wise division of two vectors of same size
     * @param <T> Type of vector
     * @param vector First vector
     * @param vector2 Second vector
     * @return A new vector with its elements equal to vector / vector2
     */
    public static <T extends Vector<T>> T div(T vector, T vector2) {
        if (vector.size() == vector2.size()) {
            return vector.copy().div(vector2);
        }
        throw new IllegalArgumentException("Different vector sizes.");
    }
    public static <T extends Vector<T>> T div(T vector, double d) {
        return vector.copy().div(d);
    }
    
    /**
     * Performs the dot-product of two vectors of same size
     * @param <T> Type of vector
     * @param vector First vector
     * @param vector2 Second vector
     * @return A new vector with its elements equal to vector (dot) vector2
     */
    public static <T extends Vector<T>> double dot(T vector, T vector2) {
        if (vector.size() == vector2.size()) {
            double ans = 0;
            for (int i=0; i<vector.size(); i++) {
                ans += vector.get(i) * vector2.get(i);
            }
            return ans;
        }
        throw new IllegalArgumentException("Different vector sizes.");
    }
    
    /**
     * Performs the scalar-product of two vectors of same size
     * @param <T> Type of vector
     * @param vector First vector
     * @param vector2 Second vector
     * @return A new vector with its elements equal to ( vector (dot) vector2 ) / norm of vector2
     */
    public static <T extends Vector<T>> double scalar(T vector, T vector2) {
        return dot(vector, vector2)/vector2.norm();
    }
    
    /**
     * Performs the projection of a vector into another vector of same size
     * @param <T> Type of vector
     * @param vector First vector
     * @param vector2 Second vector
     * @return A new vector with its elements equal to vector (proj) vector2
     */
    public static <T extends Vector<T>> T proj(T vector, T vector2) {
        double denom = Vectors.dot(vector2, vector2);
        if (denom != 0) {
            double scalar = Vectors.dot(vector, vector2)/denom;
            return Vectors.mulElem(vector, scalar);
        }
        return vector;
    }
    
    /**
     * Performs the rejection of a vector into another vector of same size
     * @param <T> Type of vector
     * @param vector First vector
     * @param vector2 Second vector
     * @return A new vector with its elements equal to vector (rej) vector2
     */
    public static <T extends Vector<T>> T rej(T vector, T vector2) {
        double denom = Vectors.dot(vector2, vector2);
        if (denom != 0) {
            double scalar = Vectors.dot(vector, vector2)/denom;
            return Vectors.sub(vector, Vectors.mulElem(vector, scalar));
        }
        return vector;
    }
    
    /**
     * Performs the cross product of two 3-vectors
     * @param vector First vector
     * @param vector2 Second vector
     * @return A new 3-vector with its elements equal to vector (cross) vector2
     */
    public static Vector3 cross(Vector3 vector, Vector3 vector2) {
        double s0 = vector.get(1)*vector2.get(2) - vector.get(2)*vector2.get(1);
        double s1 = vector.get(2)*vector2.get(0) - vector.get(0)*vector2.get(2);
        double s2 = vector.get(0)*vector2.get(1) - vector.get(1)*vector2.get(0);
        return new Vector3(s0, s1, s2);
    }
    
    /**
     * Performs a double cross product of two vectors
     * The vectors must be of size 2 or 3
     * @param a First vector
     * @param b Second vector
     * @return A new vector with its elements equal to a (cross) b (cross) b
     */
    public static Vector cross_AxBxB(Vector a, Vector b) {
        if (a instanceof Vector2 && b instanceof Vector2) {
            return cross_AxBxB((Vector2) a, (Vector2) b);
        } else if (a instanceof Vector3 && b instanceof Vector3) {
            return cross_AxBxB((Vector3) a, (Vector3) b);
        }
        throw new IllegalArgumentException("Cannot perform the cross product on a vector of this size");
    }
    
    /**
     * Performs a double cross product of two 2-vectors
     * @param a First vector
     * @param b Second vector
     * @return A new 2-vector with its elements equal to a (cross) b (cross) b
     */
    public static Vector2 cross_AxBxB(Vector2 a, Vector2 b) {
        double AxB_k = (a.get(0) * b.get(1) - a.get(1) * b.get(0));
        double AxBxB_i = -b.get(1) * AxB_k;
        double AxBxB_j = b.get(0) * AxB_k;
        return new Vector2(AxBxB_i, AxBxB_j);
    }
    
    /**
     * Performs a double cross product of two 3-vectors
     * @param a First vector
     * @param b Second vector
     * @return A new 3-vector with its elements equal to a (cross) b (cross) b
     */
    public static Vector3 cross_AxBxB(Vector3 a, Vector3 b) {
        return cross(cross(a, b), b);
    }
    
    /**
     * Performs a double cross product of two vectors
     * The vectors must be of size 2 or 3
     * @param a First vector
     * @param b Second vector
     * @return A new vector with its elements equal to a (cross) b (cross) a
     */
    public static Vector cross_AxBxA(Vector a, Vector b) {
        if (a instanceof Vector2 && b instanceof Vector2) {
            return cross_AxBxA((Vector2) a, (Vector2) b);
        } else if (a instanceof Vector3 && b instanceof Vector3) {
            return cross_AxBxA((Vector3) a, (Vector3) b);
        }
        throw new IllegalArgumentException("Cannot perform the cross product on a vector of this size");
    }
    
    /**
     * Performs a double cross product of two 2-vectors
     * @param a First vector
     * @param b Second vector
     * @return A new 2-vector with its elements equal to a (cross) b (cross) a
     */
    public static Vector2 cross_AxBxA(Vector2 a, Vector2 b) {
        double AxB_k = (a.get(0) * b.get(1) - a.get(1) * b.get(0));
        double AxBxA_i = -a.get(1) * AxB_k;
        double AxBxA_j = a.get(0) * AxB_k;
        return new Vector2(AxBxA_i, AxBxA_j);
    }
    
    /**
     * Performs a double cross product of two 3-vectors
     * @param a First vector
     * @param b Second vector
     * @return A new 3-vector with its elements equal to a (cross) b (cross) a
     */
    public static Vector3 cross_AxBxA(Vector3 a, Vector3 b) {
        return cross(cross(a, b), a);
    }
    
    
    /**
     * Computes the arithmetic mean
     * @param <T>
     * @param vectors
     * @return Mean of vectors, also called Center of mass or Centroid
     */
    public static <T extends Vector<T>> T mean(T... vectors) {
        return meanGen(1, vectors);
    }
    
    /**
     * Computes the generalized mean
     * @param <T>
     * @param p
     * @param vectors
     * @return Generalized mean of vectors
     */
    public static <T extends Vector<T>> T meanGen(double p, T... vectors) {
        if (vectors.length < 1) {
            return null;
        }
        T sumVector = vectors[0].shell();
        int n = vectors.length;
        if (p == 1) { //Arithmetic mean
            for (T v : vectors) {
                sumVector.add(v);
            }
            sumVector.div(n);
            return sumVector;
        } else if (p == 0) { //Geometric mean
            for (T v : vectors) {
                sumVector.mulElem(v);
            }
            sumVector.pow(1d/n);
            return sumVector;
        } else if (p == -1) { //Harmonic mean
            for (T v : vectors) {
                sumVector.add(v.copy().inv());
            }
            return vectors[0].shell().fill(n).div(sumVector);
        } else if (p == Double.POSITIVE_INFINITY) {//TODO: UNOPTIMIZED, Useless comparisons with first vector, repeat comparisons and useless vector instantiation
            sumVector = vectors[0].copy();
            for (T v : vectors) {
                sumVector = max(sumVector, v);
            }
            return sumVector;
        } else if (p == Double.NEGATIVE_INFINITY) {
            sumVector = vectors[0].copy();
            for (T v : vectors) {
                sumVector = min(sumVector, v);
            }
            return sumVector;
        } else {
            for (T v : vectors) {
                sumVector.add(v.copy().pow(p));
            }
            sumVector.div(n).pow(1d/p);
            return sumVector;
        }
    }
    
    /**
     * Finds the vector that has the shortest distance to the target
     * @param <T>
     * @param target
     * @param choices
     * @return Closest vector in Euclidian space
     */
    public static <T extends Vector<T>> T findMin(T target, T... choices) {
        return findMinGen(DistanceMetric.EUCLIDIAN, target, choices);
    }
    
    /**
     * Finds the vector that has the shortest distance to the target in general spaces
     * @param <T>
     * @param distanceMetric
     * @param target
     * @param choices
     * @return Closest vector based on distanceMetric
     */
    public static <T extends Vector<T>> T findMinGen(DistanceMetric distanceMetric, T target, T... choices) {
        double minScore = Double.POSITIVE_INFINITY;
        T bestVector = null;
        
        for (T v : choices) {
            double score = distanceMetric.apply(target, v);
            if (score < minScore) {
                minScore = score;
                bestVector = v;
            }
        }
        
        return bestVector;
    }
    
    /**
     * Computes the euclidian geometric median, this function is slow and can be approximated using findMin and arithmetic mean
     * @param <T>
     * @param vectors
     * @return Geometric median of vectors
     */
    public static <T extends Vector<T>> T median(T... vectors) {
        return medianGen(DistanceMetric.EUCLIDIAN, vectors);
    }
    /**
     * Computes the median in Lp space, this function is slow and can be approximated using findMinGen and mean
     * @param <T>
     * @param p
     * @param vectors
     * @return Geometric median of vectors in Lp space
     */
    public static <T extends Vector<T>> T medianLp(double p, T... vectors) {
        return medianGen(DistanceMetric.getLpSpaceMetric(p));
    }
    /**
     * Computes the generalized median, this function is slow and can be approximated using findMinGen and mean
     * @param <T>
     * @param distanceMetric
     * @param vectors
     * @return Generalized median of vectors
     */
    public static <T extends Vector<T>> T medianGen(DistanceMetric distanceMetric, T... vectors) {
        double minScore = Double.POSITIVE_INFINITY;
        T bestVector = null;
        
        for (T t : vectors) {
            double scoreSum = 0;
            for (Vector v : vectors) {
                scoreSum += distanceMetric.apply(t, v);
            }
            if (scoreSum < minScore) {
                minScore = scoreSum;
                bestVector = t;
            }
        }
        
        return bestVector;
    }
    
    /**
     * Computes L2 norm of a vector
     * @param vector
     * @return Norm of vector
     */
    public static double norm(Vector<?> vector) {
        return normLp(vector, 2);
    }
    /**
     * Computes Lp norm of a vector
     * @param vector
     * @return Lp Norm of vector
     */
    public static double normLp(Vector<?> vector, double p) {
        if (p == 0) { //"zero" norm
            int count = 0;
            for (double value : vector.content()) {
                if (value != 0) {
                    count++;
                }
            }
            return count;
        } else if (p == 1) {
            double sum = 0;
            for (double value : vector.content()) {
                sum += Math.abs(value);
            }
            return sum;
        } else if (p == 2) {
            double sumSqr = 0;
            for (double value : vector.content()) {
                sumSqr += value*value;
            }
            return Math.sqrt(sumSqr);
        } else if (Double.isInfinite(p)) {
            double max = 0;
            
            for (double value : vector.content()) {
                value = Math.abs(value);
                if (value > max) {
                    max = value;
                }
            }
            return max;
        } else {
            double sumP = 0;
            for (double value : vector.content()) {
                sumP += Math.pow(Math.abs(value), p);
            }
            return Math.pow(sumP, 1d/p);
        }
    }
    
    /**
     * Computes sum of squares of a vector
     * @param vector
     * @return Norm of vector
     */
    public static double distanceSum(Vector<?> vector) {
        return distanceSumLp(vector, 2);
    }
    public static double distanceSumLp(Vector<?> vector, double p) {
        if (p == 0) { //"zero" distance
            int count = 0;
            for (double value : vector.content()) {
                if (value != 0) {
                    count++;
                }
            }
            return count;
        } else if (p == 1) {
            double sum = 0;
            for (double value : vector.content()) {
                sum += Math.abs(value);
            }
            return sum;
        } else if (p == 2) {
            double sumSqr = 0;
            for (double value : vector.content()) {
                sumSqr += value*value;
            }
            return sumSqr;
        } else if (Double.isInfinite(p)) {
            for (double value : vector.content()) {
                if (Math.abs(value) >= 1) {
                    return Double.POSITIVE_INFINITY;
                }
            }
            return 0;
        } else {
            double sumP = 0;
            for (double value : vector.content()) {
                sumP += Math.pow(Math.abs(value), p);
            }
            return sumP;
        }
    }
}
