package fishmodel.enkf;

import java.util.Random;

public class Util {
    public static double[] reshapeArray(double[][][] arr) {
        // Calculate the total number of elements in the 3D array
        int totalElements = arr.length * arr[0].length * arr[0][0].length;

        // Create a new 1D array with the same number of elements
        double[] reshapedArray = new double[totalElements];

        // Reshape the 3D array into a 1D array
        int index = 0;
        for (int i = 0; i < arr.length; i++) {
            for (int j = 0; j < arr[i].length; j++) {
                for (int k = 0; k < arr[i][j].length; k++) {
                    reshapedArray[index] = arr[i][j][k];
                    index++;
                }
            }
        }

        return reshapedArray;
    }

    /**
     * Populates the given 3D array with the elements of a 1D array.
     * @param elements
     * @param arr
     * @param dims
     */
    public static void reshapeInto3d(double[] elements, double[][][] arr, int[] dims) {

        for (int i = 0; i < dims[0]; i++) {
            for (int j = 0; j < dims[1]; j++) {
                for (int k = 0; k < dims[2]; k++) {
                    arr[i][j][k] = elements[getStateIndex(i, j, k, dims)];
                    /*if (Double.isNaN(arr[i][j][k])) {
                        System.out.println("NaN: i="+i+", j="+j+", k="+k);
                        System.exit(0);
                    }*/
                }
            }
        }
    }

    public static int getStateIndex(int i, int j, int k, int[] cageDims) {
        return i*cageDims[1]*cageDims[2] + j*cageDims[2] + k;
    }

    public static double[] getStateCoords(int index, int[] dims) {
        double[] res = new double[3];
        res[0] = index / (dims[1]*dims[2]);
        index -= res[0]*dims[1]*dims[2];
        res[1] = index / dims[2];
        index -= res[1]*dims[2];
        res[2] = index;
        return res;
    }


    public static double updateGaussMarkov(double prev, double beta, double sigma, double dt, Random rnd) {
        double f = Math.exp(-beta * dt);
        return f * prev + Math.sqrt(1. - f * f) * rnd.nextGaussian() * sigma;
    }

    public static double getGaussValue(double sigma, Random rnd) {
        return sigma*rnd.nextGaussian();
    }

    public static void main(String[] args) {
        int[] dims = new int[] {5, 5, 3};
        int index = getStateIndex(0,2,2, dims);
        System.out.println("Index: "+index);
        double[] coord = getStateCoords(index, dims);
        System.out.println("Coord:");
        for (int i = 0; i < coord.length; i++) {
            double i1 = coord[i];
            System.out.println(i1);
        }
    }
}
